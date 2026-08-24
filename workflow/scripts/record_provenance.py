"""Verify the acquired GEO artifacts and write the provenance records.

Owner task: P0-T2. Driven by rule p0t2_record_provenance.

Re-hashes every artifact from disk rather than trusting the sidecars written at
download time — the sidecar says what arrived, this says what is there now, and
those are different claims.

Verification is two-tier (ADR 0005):

  fixed     author-uploaded supplementary files, immutable since 2022.
            A digest mismatch is a HARD FAILURE. These bytes are the data.

  volatile  GEO-generated metadata (the SOFT family file), regenerated
            server-side on request. Its digest is recorded as an observation.
            Drift WARNS and re-records; failing here would break the build for
            reasons that have nothing to do with the data. The content that
            matters is asserted on by P0-T3, which is the right place for it.

Outputs:
  resources/checksums.sha256  standard `<digest>  <path>` — deliberately in the
                              format `shasum -a 256 -c` understands, so the
                              claim is checkable without this pipeline.
  resources/provenance.tsv    url, access date, size, digest, tier, pin status.
"""

import hashlib
import json
from pathlib import Path

CHUNK = 1 << 20

artifacts = snakemake.params.artifacts
# Resolve each artifact by the path Snakemake actually handed us, rather than
# rebuilding it from a directory prefix — that keeps the record honest if the
# paths ever move, and keeps `snakemake --lint` clean.
path_by_dest = {Path(p).name: p for p in snakemake.input.artifacts}
checksums_path = Path(snakemake.output.checksums)
provenance_path = Path(snakemake.output.provenance)
log_path = Path(snakemake.log[0])
log_path.parent.mkdir(parents=True, exist_ok=True)


def sha256_of(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(CHUNK), b""):
            digest.update(chunk)
    return digest.hexdigest()


with open(log_path, "w", encoding="utf-8") as log:

    def emit(msg):
        print(msg, file=log, flush=True)

    # Sidecars are keyed by dest filename; index them so the config manifest
    # drives the ordering (deterministic output, independent of glob order).
    sidecars = {}
    for meta_file in snakemake.input.meta:
        record = json.loads(Path(meta_file).read_text(encoding="utf-8"))
        sidecars[record["artifact"]] = record

    rows = []
    failures = []
    unpinned = []

    for key in sorted(artifacts):
        spec = artifacts[key]
        dest = spec["dest"]
        path = path_by_dest[dest]
        record = sidecars[dest]

        observed = sha256_of(path)
        pinned = spec["sha256"]
        tier = spec["stability"]

        if pinned is None:
            status = "UNPINNED"
            unpinned.append((key, tier, observed))
        elif pinned == observed:
            status = "verified"
        elif tier == "fixed":
            status = "MISMATCH"
            failures.append((key, pinned, observed))
        else:
            status = "drifted"

        emit(f"{key:<18} {tier:<8} {status:<9} {observed}")

        if status == "drifted":
            emit(
                f"  WARNING: {key} is `volatile` and its digest changed "
                f"(pinned {pinned}). Re-recorded, not failed. Update the pin in "
                f"config/config.yaml if you want the record to match."
            )

        rows.append(
            {
                "artifact": key,
                "file": dest,
                "path": path,
                "url": record["url"],
                "accessed_utc": record["accessed_utc"],
                "last_modified": record["last_modified"] or "",
                "bytes": record["bytes"],
                "sha256": observed,
                "stability": tier,
                "pin_status": status,
            }
        )

    if failures:
        detail = "\n".join(
            f"  {key}\n    pinned   {pinned}\n    observed {observed}"
            for key, pinned, observed in failures
        )
        raise RuntimeError(
            "SHA-256 mismatch on `fixed` artifact(s) — these bytes are the data, "
            "so this is a hard failure, not a warning:\n"
            f"{detail}\n"
            "Either the download is corrupt (delete it and re-fetch) or the "
            "remote changed (investigate before touching the pin). Do not "
            "'fix' this by pasting the observed digest into config.yaml."
        )

    columns = [
        "artifact", "file", "url", "accessed_utc", "last_modified",
        "bytes", "sha256", "stability", "pin_status",
    ]
    with open(provenance_path, "w", encoding="utf-8") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in rows:
            handle.write("\t".join(str(row[c]) for c in columns) + "\n")

    # Paths relative to the repo root, so `shasum -a 256 -c resources/checksums.sha256`
    # works when run from the root.
    with open(checksums_path, "w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(f"{row['sha256']}  {row['path']}\n")

    emit(f"\nwrote {provenance_path}")
    emit(f"wrote {checksums_path}")

    if unpinned:
        banner = "\n".join(
            f'      sha256: "{observed}"   # {key} ({tier})'
            for key, tier, observed in unpinned
        )
        message = (
            "\n"
            + "=" * 72
            + "\nUNPINNED ARTIFACTS — provenance is NOT yet enforced for these.\n"
            + "Paste each digest into config/config.yaml under acquire.artifacts,\n"
            + "then re-run to confirm the verification path passes. A `fixed`\n"
            + "artifact left null is a provenance control that looks enforced and\n"
            + "is not; do not commit in this state.\n\n"
            + banner
            + "\n"
            + "=" * 72
        )
        emit(message)
        print(message)

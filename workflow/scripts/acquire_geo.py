"""Fetch one GEO artifact into resources/raw/, with a provenance sidecar.

Owner task: P0-T2. Driven by rule p0t2_fetch_geo, one job per artifact in
config["acquire"]["artifacts"].

This script downloads and records; it does **not** verify against the pinned
digest. Verification belongs to p0t2_record_provenance, and the split is
deliberate: the artifact outputs here are protected(), so if the pinned SHA-256
were a param of this rule then adding a pin would mark the job out of date and
Snakemake would fail trying to overwrite a write-protected file. Keeping the
pins out of this rule means adding one re-runs only the cheap verification.

Integrity at download time is still checked, against the server's own
Content-Length — a truncated transfer fails here rather than surfacing as a
confusing digest mismatch downstream.

Stdlib only (urllib): a stdlib HTTP GET does not justify a new conda env and a
regenerated lock, and workflow/envs/* are provisional until P6-T1 anyway.
"""

import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import Request, urlopen

CHUNK = 1 << 20  # 1 MiB

url = snakemake.params.url
dest = Path(snakemake.output.artifact)
meta_path = Path(snakemake.output.meta)
log_path = Path(snakemake.log[0])

log_path.parent.mkdir(parents=True, exist_ok=True)
dest.parent.mkdir(parents=True, exist_ok=True)
meta_path.parent.mkdir(parents=True, exist_ok=True)

with open(log_path, "w", encoding="utf-8") as log:

    def emit(msg):
        print(msg, file=log, flush=True)

    accessed = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    emit(f"artifact  : {snakemake.wildcards.artifact}")
    emit(f"url       : {url}")
    emit(f"accessed  : {accessed}")

    # Identify ourselves. NCBI asks automated clients to, and an anonymous
    # urllib default is the kind of thing that gets rate-limited without
    # explanation.
    request = Request(url, headers={"User-Agent": "nsclc-brainmet-time/P0-T2 (Snakemake)"})

    tmp = dest.with_name(dest.name + ".part")
    digest = hashlib.sha256()
    written = 0

    try:
        with urlopen(request) as response:
            declared = response.headers.get("Content-Length")
            last_modified = response.headers.get("Last-Modified")
            declared = int(declared) if declared is not None else None
            emit(f"declared  : {declared} bytes")
            emit(f"modified  : {last_modified}")

            with open(tmp, "wb") as handle:
                while True:
                    chunk = response.read(CHUNK)
                    if not chunk:
                        break
                    handle.write(chunk)
                    digest.update(chunk)
                    written += len(chunk)

        if declared is not None and written != declared:
            raise RuntimeError(
                f"truncated download: got {written} bytes, server declared {declared}. "
                f"Nothing was written to {dest}."
            )

        os.replace(tmp, dest)
    finally:
        # A failed transfer must not leave a .part behind masquerading as
        # progress, and must never leave a partial file at the real path.
        if tmp.exists():
            tmp.unlink()

    sha256 = digest.hexdigest()
    emit(f"bytes     : {written}")
    emit(f"sha256    : {sha256}")

    meta_path.write_text(
        json.dumps(
            {
                "artifact": snakemake.wildcards.artifact,
                "key": snakemake.params.key,
                "url": url,
                "accessed_utc": accessed,
                "bytes": written,
                "sha256": sha256,
                "last_modified": last_modified,
                "stability": snakemake.params.stability,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    emit(f"sidecar   : {meta_path}")

---
description: Scaffold a Snakemake rule + script + log/benchmark wiring
argument-hint: <phase number> <rule name>
allowed-tools: Bash(ls:*), Bash(snakemake -n:*), Bash(snakemake --lint:*), Read, Write, Edit
---

Scaffold a Snakemake rule: **$ARGUMENTS**

1. **Pick the right `.smk`** in `workflow/rules/` from the phase number
   (`00_acquire`, `01_annotate`, `02_qc`, `03_landscape`, `04_contexture`,
   `05_checkpoints`, `06_crosstalk`, `07_survival`, `99_fair`).
2. **Write the rule.** Every rule declares all of these — no exceptions (§4.3):

```python
rule <name>:
    input:
        ...,
    output:
        ...,
    params:
        seed=config["seed"],      # only if the step is stochastic
    log:
        "results/logs/<name>.log",
    benchmark:
        "results/benchmarks/<name>.tsv",
    conda:
        "../envs/<env>.yaml",
    threads: 1
    script:
        "../scripts/<name>.<py|R>"
```

3. **Create the script stub** in `workflow/scripts/`, reading its inputs from
   `snakemake.input` / `snakemake@input` rather than hardcoded paths.
4. **Add the output to that phase's `TARGETS_<PHASE>` list**, and turn the
   phase's switch on in `config/config.yaml` in the same change. A rule that
   exists but is not reachable from a target is not wired in.
5. **Verify:** `snakemake --lint` and `snakemake -n` must both be clean.

Choose the env by language boundary, never by convenience: `r-geomx` for
GeoMx/`standR` work, `r-stats` for mixed models, `py-analysis` for AnnData,
plotting, enrichment, survival and packaging. Do not reimplement `standR` or
`GeomxTools` in Python (§4.1).

If the step is stochastic, pass `config["seed"]` explicitly. No implicit RNG.

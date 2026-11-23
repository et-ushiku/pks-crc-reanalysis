# code/

## Runs from this repository, with no other input

| Script | What it does | Needs |
|---|---|---|
| `analyse.py` | Reproduces the primary result: recomputes breadth from the depth histograms, rarefies to a common read count, applies the pre-specified gate, prints per-cohort tables and the Mantel-Haenszel meta-analysis | `data/` only. Python 3, standard library |

```
python3 code/analyse.py            # default target 20,000,000 reads
python3 code/analyse.py 10000000   # try a different common depth
```

It asserts that the calls recomputed from `data/depth_hist.tsv` match the stored
calls in `data/per_sample.tsv`, so a failure means the shipped data is inconsistent.

## How data/ was generated

These read the original working tree (BAM depth tables under
`/mnt/c/crc_strain_pilot/pks/`), which is not part of this repository. They are
included so the derivation of `data/` is auditable, not because they will run here.

| Script | Output |
|---|---|
| `make_table.py` | `data/depth_hist.tsv`, `data/per_sample.tsv` |
| `make_dup_table.py` | `data/missing_split.tsv` |

## Archival: analyses reported in RESULTS.md

Also written against the original working tree. Kept for reference; the numbers they
produced are quoted in [../RESULTS.md](../RESULTS.md).

| Script | What it produced |
|---|---|
| `verify_claims.py` | The checks that refuted items 1 and 2 of [../CORRECTIONS.md](../CORRECTIONS.md) |
| `sens_tataka.py` | The calling-rule sensitivity analysis (five alternative rules) |
| `de_carrier.py` | The marker-density carrier screen on the DE-CRC cohort |

Their comments and printed output are in Japanese.

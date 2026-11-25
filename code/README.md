# code/

## Runs from this repository, with no other input

Python 3, standard library only. No arguments needed unless shown.

| Script | What it does |
|---|---|
| `analyse.py` | Reproduces the primary result: recomputes breadth from the depth histograms, rarefies to a common read count, applies the pre-specified gate, prints per-cohort tables and the Mantel-Haenszel meta-analysis |
| `verify_extra.py` | Six checks that each used to be impossible from `data/` alone: the DE-CRC depth difference, the density-ratio offset, the depth-adjusted regression, the carrier table, the assembly junction, and the strain typing |
| `verify_rarefy.py` | Checks the analytic rarefying formula against actual subsampling, and that the disagreement is smaller than the seed-to-seed spread |

```
python3 code/analyse.py            # default target 20,000,000 reads
python3 code/analyse.py 10000000   # try a different common depth
python3 code/verify_extra.py
python3 code/verify_rarefy.py
```

`analyse.py` asserts that the calls recomputed from `data/depth_hist.tsv` match the
stored calls in `data/per_sample.tsv`, so a failure means the shipped data is
inconsistent. `verify_extra.py` likewise asserts that the `kp_plus` contig bin is
exactly `kp_bin` plus every island-bearing contig.

## Needs network access

| Script | What it does |
|---|---|
| `verify_repeat.py` | Fetches AM229678.1 (~55 kb) from the ENA browser API and recomputes the duplication: its extent, its identity, the spacing between the copies, and the absence of any other repeat at comparable spacing |

```
python3 code/verify_repeat.py            # fetches the sequence
python3 code/verify_repeat.py local.fa   # or use a local copy
```

This is the only script here that reaches the network, and the only thing it sends is
the accession. The sequence itself is not redistributed in this repository.

## How data/ was generated

These read the original working tree (BAM depth tables, assemblies and alignments
under `/mnt/c/crc_strain_pilot/`), which is **not** part of this repository. They are
included so the derivation of `data/` is auditable, not because they will run here.
The newer ones take their source directory as an argument, so the paths they assume
are visible rather than buried in the file.

| Script | Output |
|---|---|
| `make_table.py` | `data/depth_hist.tsv`, `data/per_sample.tsv` |
| `make_dup_table.py` | `data/missing_split.tsv` |
| `make_carrier_table.py` | `data/carrier_six_genomes.tsv` |
| `make_junction_table.py` | `data/assembly_alignment.tsv`, `data/assembly_junction_depth.tsv` |
| `make_kleborate_table.py` | `data/kleborate/*` - the typing output of both runs, plus the contig assignment behind it |
| `make_clb_in_genomes.py` | `data/clb_in_six_genomes.tsv` - which reference genomes carry a *clb* island of their own |

`make_clb_in_genomes.py` is the exception to "will not run here": it needs only an
island FASTA, `data/clb_genes.bed` and the six-genome FASTA, and does exact 40-mer
counting with no external tool.

## Archival: analyses reported in RESULTS.md

Also written against the original working tree. Kept for reference; the numbers they
produced are quoted in [../RESULTS.md](../RESULTS.md).

| Script | What it produced |
|---|---|
| `verify_claims.py` | The checks that refuted items 1 and 2 of [../CORRECTIONS.md](../CORRECTIONS.md) |
| `sens_tataka.py` | The calling-rule sensitivity analysis (five alternative rules) |
| `de_carrier.py` | The marker-density carrier screen on the DE-CRC cohort |

Their comments and printed output are in Japanese, as are those of the `make_*`
generators.

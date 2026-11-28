# code/

## Runs from this repository, with no other input

Python 3, standard library only. No arguments needed unless shown.

| Script | What it does |
|---|---|
| `analyse.py` | Reproduces the primary result: recomputes breadth from the depth histograms, rarefies to a common read count, applies the pre-specified gate, prints per-cohort tables and the Mantel-Haenszel meta-analysis |
| `verify_extra.py` | Seven checks that each used to be impossible from `data/` alone: the DE-CRC depth difference, the density-ratio offset, the depth-adjusted regression, the carrier table, the assembly junction, the strain typing, and whether the island contigs belong to that genome |
| `verify_rarefy.py` | Checks the analytic rarefying formula against actual subsampling, and that the disagreement is smaller than the seed-to-seed spread |
| `verify_numbers.py` | Recomputes 147 of the numbers written in the prose - breadth ceiling, duplication cost, background, carrier ratios, the DE-CRC depth confound, the *K. pneumoniae* abundance column, both sensitivity tables and the island contig tiling - and exits non-zero if any disagrees |
| `verify_marker_scaling.py` | Tests the one assumption behind the abundance columns - that marker reads scale with marker length - against the six-genome index, which does not involve marker length. Asserts the interval contains 1 and excludes the unscaled alternative |

```
python3 code/analyse.py            # default target 20,000,000 reads
python3 code/analyse.py 10000000   # try a different common depth
python3 code/verify_extra.py
python3 code/verify_rarefy.py
python3 code/verify_numbers.py
python3 code/verify_marker_scaling.py
```

`analyse.py` asserts that the calls recomputed from `data/depth_hist.tsv` match the
stored calls in `data/per_sample.tsv`, so a failure means the shipped data is
inconsistent. `verify_extra.py` likewise asserts that the `kp_plus` contig bin is
exactly `kp_bin` plus every island-bearing contig.

`verify_numbers.py` is the one that checks the **documents** rather than the data: it
holds each published figure as a literal next to the computation that should produce
it. It exists because nothing did that until 2026-11-26, which is how the stale
figures in [../CORRECTIONS.md](../CORRECTIONS.md) items 13 and 17 survived. It does
not cover numbers quoted from papers, nor numbers derived from inputs we do not
ship - both limits are listed in its docstring.

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
| `make_backbone_table.py` | `data/kleborate/island_contig_backbone.tsv`, `data/kleborate/contig_depth.tsv` - each island contig's non-island portion and its coverage depth |
| `make_marker_table.py` | `data/marker_counts.tsv` - raw marker read counts per sample per SGB, which the abundance columns are computed from and which were not shipped before |
| `make_locus_table.py` | `data/kleborate/island_tiling.tsv`, `data/kleborate/island_locus.tsv` - how the island contigs overlap each other on the island reference, and where the element they belong to sits in the *Klebsiella* chromosome |
| `make_asm282_table.py` | `data/assembly_282_island.tsv` - where the island landed in `cancer_ERR1018282`'s assembly, with every hit split into *clb* gene bases and non-island flanking bases before any species is named |
| `make_mapq_table.py` | `data/mapq_split.tsv` - reads per gene with and without the MAPQ 20 filter, from `pipeline/mapq_split.sh`; the measurement that replaced item 10's reasoning |
| `make_ratio_scale.py` | `data/island_ratio_scale.tsv` - the island/species depth ratio for all 14 carriers and the samples that calibrate its scale, selected without looking at the ratio |
| `make_link_table.py` | `data/kleborate/contig_read_depth.tsv`, `data/kleborate/contig_links.tsv` - per-contig read depth (not megahit `multi`) and read pairs spanning two contigs, from `pipeline/asm_linkage.sh` |

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

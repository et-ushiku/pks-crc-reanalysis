# pks island (colibactin) carriage and colorectal cancer: a reanalysis of 397 public stool metagenomes

**Carriage of the *pks* island is not significantly associated with colorectal cancer
across four public cohorts.** The direction is consistent (more carriers in CRC in
three of four cohorts) but the confidence interval includes 1.

| | |
|---|---|
| Samples | **397** stool metagenomes, 4 cohorts, 4 countries |
| Judgeable under pre-specified criteria | **167** |
| Pre-specified primary result | **Mantel-Haenszel OR 1.91, 95% CI 0.76-4.75, P = 0.166** |
| Heterogeneity | I2 = 17% |

This agrees in both direction and significance with the largest published analysis,
Piccinno et al. 2025 (*Nat Med*, 18 cohorts, 3,741 samples), which reports:

> Although not significant, the carriage of colibactin-producing genes by
> *E. coli* and *Klebsiella* spp. was increased in CRC

A result that does not reach significance at 3,741 samples will not reach it at 397.
**This is not evidence that there is no effect; it is evidence that any effect is small.**

**This is not a diagnostic test and must not be used for clinical decisions.**

---

## What is actually new here

The negative result is a second, independent confirmation of something already
published, so its novelty is modest. The part that may save other people time is
[PITFALLS.md](PITFALLS.md) - four ways that *pks* detection from shotgun stool
metagenomes goes wrong, none of which we found documented anywhere:

1. **The reference contains a 1,248 bp duplication, so breadth of coverage can never
   reach 1.0.** In AM229678.1 the same 1,248 bp appears at 16,066-17,313 (inside
   *clbK*) and 21,719-22,966 (inside *clbJ*), 99.84% identical. Reads there are
   multi-mapping and are discarded at MAPQ 20. The observed ceiling is **~0.95**
   (maximum seen in 397 samples: 0.9607). A sample at 0.95 is not a partial island.
2. **Raw read counts against the island are not a measure of island abundance.**
   The 55,140 bp reference contains 5,685 bp that is not *clb* gene sequence.
   Island-negative samples still accumulate reads there (median background
   pks/*E. coli* ratio 0.012 across 82 such samples), with breadth near 0.004.
3. **A marker-gene density ratio carries a large systematic offset.** Ours was
   ~1.9x, decomposed as 1.115 (island length convention) x 1.14 (marker fragment
   edge loss) x ~1.50 (residual marker-vs-island mappability). Ratios are useful for
   exclusion ("this species cannot account for the signal"), not for attribution.
4. **Breadth is depth-limited below about 10x island coverage.** 34 of 44 positive
   samples sat below that. Group differences in sequencing depth therefore masquerade
   as group differences in carriage - in one of our cohorts the depth difference
   between groups is enormous (P = 1.5e-14). Equalise depth before calling.

## Reproducing the result without the raw data

The four cohorts total roughly 2.5 TB of FASTQ. You do not need it.

Every call in this analysis is a function of the per-position depth over the 18 *clb*
genes, and the expected breadth after rarefying to a fraction `f` is

```
E[breadth] = mean_i [ 1 - (1-f)^d_i ]
```

which depends only on the *multiset* of depths, not their order. So this repository
ships the depth histograms (753 KB) instead of the reads:

```
python3 code/analyse.py          # pure stdlib, no scipy, under a second
python3 code/verify_rarefy.py    # checks the analytic rarefaction against
                                 # actually rarefied reads (5 runs, 2 samples)
```

This recomputes breadth from `data/depth_hist.tsv`, checks it against the stored
per-sample values, rarefies every sample to 20M reads, applies the pre-specified
abundance gate, and prints the per-cohort tables and the meta-analysis. It reproduces
the headline numbers above exactly.

The formula is exact by linearity of expectation, but a derivation is not a check. For
two samples we also rarefied the FASTQ for real with `seqtk sample`, remapped, and
recomputed breadth; `verify_rarefy.py` puts the two side by side. The largest gap over
five runs is **0.0211**, which is smaller than the **0.0182** spread between three
different seeds at the same depth, and **all five calls agree**. See
[METHODS.md](METHODS.md) for the table and for what this check does not settle.

## Contents

| Path | What it is |
|---|---|
| [RESULTS.md](RESULTS.md) | All result tables, per cohort and pooled, plus the carrier analysis |
| [METHODS.md](METHODS.md) | Pipeline, reference, pre-specified thresholds, rarefaction derivation |
| [PITFALLS.md](PITFALLS.md) | The four failure modes above, with the data that shows each |
| [CORRECTIONS.md](CORRECTIONS.md) | Every claim we made and later had to withdraw, and why |
| `data/depth_hist.tsv` | Per-sample, per-gene depth histograms - **the reproducible dataset** |
| `data/per_sample.tsv` | 397 rows: breadth, genes called, call, read count, *E. coli* %, *K. pneumoniae* % (**NA where that species was never measured** - see CORRECTIONS.md item 7) |
| `data/missing_split.tsv` | Uncovered bases split inside/outside the duplication |
| `data/clb_genes.bed` | The 18 gene intervals on AM229678.1 |
| `code/analyse.py` | Reproduces the primary result from `data/` alone |
| `code/verify_rarefy.py` | Compares the analytic rarefaction against actually rarefied reads; largest gap 0.0211, all 5 calls agree |
| `data/rarefy_check_hist.tsv` | Depth histograms from 5 runs where the FASTQ really was subsampled and remapped |
| `code/make_table.py`, `code/make_dup_table.py` | How `data/` was generated from BAM depth |
| `pipeline/` | The mapping and download scripts, for anyone starting from FASTQ |

## Data sources

All sequence data is public, from the European Nucleotide Archive. Nothing in this
repository redistributes reads.

| Cohort | Study | ENA |
|---|---|---|
| FR-CRC (40) | Zeller et al. 2014, *Mol Syst Biol*, French cohort | PRJEB6070 |
| AT-CRC (109) | Feng et al. 2015, *Nat Commun* | PRJEB7774 |
| CN-CRC (128) | Yu et al. 2017, *Gut* | PRJEB10878 |
| DE-CRC (120) | Wirbel et al. 2019, *Nat Med* (82) + Zeller et al. 2014 German cohort (38) | PRJEB27928, PRJEB6070 |

Reference: *E. coli* IHE3034 *pks* island, GenBank **AM229678.1** (55,140 bp).

## Limitations

- **No covariates.** We did not adjust for age, sex or study, which Piccinno 2025
  does. Our cohorts are compared as strata, not as a covariate-adjusted model.
- **Direction is unstable.** With 6-8 positives per cohort the pooled direction
  changes when the analysis is varied. In CN-CRC the odds ratio moves from 0.85
  (unrarefied) to 2.13 (rarefied). Neither is significant; the honest statement is
  "no significant difference", not "fewer in cancer" or "more in cancer".
- **Carrier attribution fails outright in 4 of 14 samples.** All 14 CN-CRC positives
  were mapped to six complete Enterobacteriaceae genomes in one competitive index.
  Four could not be assigned to any species, because *E. coli* and *K. pneumoniae*
  were both abundant and both ratios were consistent with carriage. Three samples are
  best explained by *K. pneumoniae* rather than *E. coli*, but only one of those was
  confirmed by assembly, where the island and the *Klebsiella* chromosome co-occur on
  one contig; the other two rest on coverage-ratio arguments, which are suggestive,
  not proof. **Two of the three are controls.**
- **Abundance is not carriage.** Two samples carry an abundant *K. pneumoniae*
  (16.7x and **109x** genome depth, breadth >0.88) where the island nonetheless tracks
  *E. coli*. We initially read the 109x sample as *Klebsiella* on abundance alone and
  were wrong (CORRECTIONS.md item 8).
- **No wet-lab work.** No isolates, no PCR, no colibactin measurement.
- One operator, one analysis, no independent replication of the code.

## License

Code: MIT (see [LICENSE](LICENSE)).
Documentation and derived data tables: CC BY 4.0.

Underlying sequence data belongs to the original studies and remains under their terms.

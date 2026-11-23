# Methods

## Pre-specified decision rules

These were fixed before the cohorts were analysed and **were not changed after seeing
any result**. They are hard-coded at the top of `code/analyse.py`.

| Parameter | Value | Meaning |
|---|---|---|
| `MIN_DEPTH` | 1 | a position counts as covered at depth >= 1 |
| `GENE_COVER` | 0.5 | a *clb* gene counts as present if >= 50% of its bases are covered |
| `CALL_BREADTH` | 0.5 | island breadth >= 0.5 ... |
| `CALL_GENES` | 14 | ... **and** >= 14 of 18 genes present -> **present** |
| `SUSPECT_BREADTH` | 0.1 | breadth in [0.1, 0.5) -> **suspect** (not counted as positive) |
| `ECOLI_MIN` | 0.1 | samples with *E. coli* < 0.1% are **not judgeable** and are excluded |
| target depth | 20,000,000 reads | samples are rarefied to this before calling |
| MAPQ | 20 | applied at every depth computation |

The abundance gate exists because below ~0.1% *E. coli* there is not enough
enterobacterial sequence for "island absent" to mean anything. Removing the gate
changes which samples are counted, so the ungated tables in [RESULTS.md](RESULTS.md)
are reported as reference only, never as the primary analysis.

## Reference and mapping

- Island reference: *E. coli* IHE3034 *pks* island, GenBank **AM229678.1**, 55,140 bp.
- Gene intervals: `data/clb_genes.bed`, 18 genes totalling **49,455 bp**. The
  remaining 5,685 bp is intergenic/backbone and is excluded from breadth.
- Aligner: `bowtie2 --very-sensitive --no-unal`, single-end (both mates supplied
  independently as `-U`, so counts are reads, not fragments).
- Depth: `samtools depth -a -Q 20`.
- Species abundance: MetaPhlAn marker sets mapped as a combined index.
  `t__SGB10068` = *E. coli* (96,300 bp / 122 fragments),
  `t__SGB10115` = *K. pneumoniae* (6,450 bp / 11 fragments).
  Conversion: **313 reads per 1% per million reads** (`PER_PCT`), calibrated once and
  then held fixed.
- Carrier confirmation: six complete Enterobacteriaceae genomes in **one combined
  index**, so reads compete and cross-species bleed is removed by MAPQ 20
  (`pipeline/asm_genomes.sh`). Species: *E. coli* IHE3034, *K. pneumoniae* HS11286,
  *C. freundii* CFNIH1, *C. koseri* BAA-895, *E. cloacae* ATCC13047,
  *K. aerogenes* KCTC2190.

## Rarefying without the reads

Sequencing depth differs between groups in some cohorts (see
[PITFALLS.md](PITFALLS.md) item 4), so calls must be made at a common depth. The raw
FASTQ is ~2.5 TB and was deleted after mapping, so we rarefy analytically.

Let position *i* of the island have depth *d_i* at full depth, and let *f* be the
fraction of reads retained. Each read is kept independently with probability *f*, so
the probability that position *i* retains at least one read is

```
P(covered | f) = 1 - (1 - f)^d_i
```

and by linearity of expectation the expected breadth is

```
E[breadth | f] = (1 / L) * sum_i [ 1 - (1 - f)^d_i ]
```

This is **exact**, not an approximation, and it does not depend on how reads are
distributed along the island, because expectation is linear regardless of correlation
between positions. Read-length and fragment structure affect the *variance* of
breadth, not its expectation.

Because the formula depends only on the multiset of depths, the per-gene depth
**histogram** is sufficient. That is what `data/depth_hist.tsv` contains, and it is
why this repository can reproduce every number without distributing reads.

### What the method does not give you

The call is a threshold on breadth, so for a sample sitting near 0.5 the expectation
alone does not determine the call - the realised breadth could fall either side.
`code/analyse.py` therefore recomputes calls, not just breadths. Samples with fewer
than the target reads cannot be rarefied up and are excluded, which is reported in the
output.

### Checking the analytic rarefaction against actually rarefied reads

A derivation is not a check. For two samples we also did it the slow way: subsampled
the FASTQ with `seqtk sample`, remapped, and recomputed breadth from the result. Run
`python3 code/verify_rarefy.py` to reproduce this from `data/` alone:

| Sample | Target | Seed | f | Full | Analytic | Actual | Call agrees |
|---|---|---|---|---|---|---|---|
| cancer_CCIS77252613ST | 20M | 1 | 0.2235 | 0.6347 | 0.2332 | 0.2455 | yes |
| cancer_CCIS77252613ST | 40M | 1 | 0.4470 | 0.6347 | 0.3969 | 0.4064 | yes |
| cancer_CCIS77252613ST | 40M | 2 | 0.4470 | 0.6347 | 0.3969 | 0.4074 | yes |
| cancer_CCIS77252613ST | 40M | 3 | 0.4470 | 0.6347 | 0.3969 | 0.3892 | yes |
| cancer_CCIS91228662ST | 40M | 1 | 0.4252 | 0.7924 | 0.5495 | 0.5706 | yes |

- **Largest gap 0.0211** over the five runs, and **all five calls agree**.
- The three runs at 40M on the same sample exist to give the gap a scale. The actual
  breadth varies between seeds from 0.3892 to 0.4074, a **spread of 0.0182**, and the
  analytic value 0.3969 falls inside that range. An expectation cannot be asked to do
  more than land inside the range that actual subsampling produces.
- **The analytic value is below the actual in 4 of the 5 runs.** On five runs that is
  not significant (two-sided binomial P = 0.38), so we are not claiming a bias, but we
  are not claiming its absence either: five runs on two samples cannot settle it.
  Anyone extending this should add runs rather than trust the sign.
- **Gene counts diverge more than breadth does** (4 analytic vs 7 actual in one run),
  for the reason given above: the per-gene 0.5 coverage threshold is a step function,
  and a gene whose expected coverage sits just under it realises above it about half
  the time. This is why the breadth threshold and the gene threshold are applied
  together rather than the gene count alone.

## Statistics

- Per cohort: Fisher's exact test, two-sided (own implementation, `fisher()` in
  `code/analyse.py`; no SciPy dependency).
- Pooled: **Mantel-Haenszel** odds ratio with cohorts as strata, variance by
  **Robins-Breslow-Greenland**. Cohorts are not pooled by simple addition, because
  control-group carriage differs between them (0% to 27%).
- Heterogeneity: Cochran's Q and I2 on the fixed-effect log odds ratios with a
  continuity correction of 0.5.
- Multiplicity: the smallest per-cohort P is reported with Bonferroni (x4) and Sidak
  corrections, since four cohorts were tested.
- Odds ratios printed per cohort use a 0.5 continuity correction so that zero cells
  remain finite; the Mantel-Haenszel estimate does not use one.

## What was not done

- No adjustment for age, sex, BMI or study. Piccinno et al. 2025 adjust for age and
  sex; we do not have those covariates assembled for all four cohorts.
- No per-strain analysis, no assembly for most samples, no wet-lab confirmation.
- No correction for human-read contamination beyond an upper bound (the cohorts
  differ in host-removal preprocessing; one cohort's exact preprocessing is
  unpublished).

# Results

All tables below are produced by `python3 code/analyse.py` from `data/` alone.
The raw output is appended at the end of this file.

---

## Primary analysis (pre-specified: rarefied to 20M reads, *E. coli* >= 0.1%)

| Cohort | Study | CRC | control | OR | Fisher P |
|---|---|---|---|---|---|
| FR-CRC | Zeller 2014 (France) | 1/11 | 3/12 | 0.39 | 0.590 |
| AT-CRC | Feng 2015 (Austria) | 4/21 | 3/21 | 1.36 | 1.000 |
| CN-CRC | Yu 2017 (China) | 6/44 | 2/33 | 2.13 | 0.454 |
| DE-CRC | Wirbel 2019 + Zeller 2014 (Germany) | 6/16 | 0/9 | 11.76 | 0.057 |
| **Pooled (Mantel-Haenszel)** | | | | **1.91** | **0.166** |

95% CI 0.76-4.75. Heterogeneity Q = 3.60, **I2 = 17%** (low).
Smallest cohort P = 0.0571; Bonferroni x4 = **0.228**.

**167 of 397 samples were judgeable.** The 230 excluded are: below 20M reads (cannot
be rarefied to the common depth), or *E. coli* below 0.1% (not enough enterobacterial
sequence for "absent" to mean anything).

### Reading this

The direction is the same in three of four cohorts, and the pooled point estimate is
about 1.9, but the interval runs from 0.76 to 4.75. **There is no significant
association.** With 42 positive calls in total, the pooled direction is not stable
either - see the CN-CRC note below.

---

## The DE-CRC cohort is the only one near significance, and we do not rest on it

| View | CRC | control | P |
|---|---|---|---|
| **Pre-specified** (rarefied + gated) | 6/16 | 0/9 | **0.0571** |
| Ungated (**not** the pre-specified analysis) | 6/60 | 0/55 | 0.0281 |
| Bonferroni across the four cohorts | | | **0.228** |

Three reasons not to treat this as a finding:

1. **Only 9 judgeable controls.** If the true control carriage were 15%, the
   probability of observing 0 of 9 is 0.23. The zero cell is unremarkable.
2. **Four cohorts were tested.** The corrected P is 0.228.
3. **The ungated P = 0.0281 is not the pre-specified analysis.** The gate was fixed in
   advance; dropping it after seeing the result would be choosing the analysis to fit
   the answer. Without the gate, "island absent" and "no enterobacteria sequenced"
   are indistinguishable.

### The depth confound in DE-CRC is extreme, but abundance is not

```
reads       CRC median 69.5M (22.2-120.4M)   control 34.9M (16.3-53.7M)  P = 1.5e-14
E. coli %   CRC median 0.0146%               control 0.0066%             P = 0.172
K. pneu. %                                                               P = 0.502
gate passed CRC 16/60                        control 10/60
```

The depth difference is handled by rarefying. The *E. coli* abundance difference is
**not significant**, so "more *E. coli* in the cancer group" cannot explain the
signal - though a 2.2x median ratio means it cannot be fully excluded either.

---

## CN-CRC: the direction flips when depth is equalised

| | CRC | control | OR | P |
|---|---|---|---|---|
| Not rarefied | 7/44 | 6/33 | **0.85** | 1.000 |
| Rarefied to 20M (primary) | 6/44 | 2/33 | **2.13** | 0.454 |

Controls in this cohort were sequenced slightly deeper, so rarefying removes more of
their coverage and pushes borderline samples from *present* to *suspect*. Neither
result is significant, so "no significant difference" holds either way, but the
**direction is not determined** at this number of positives.

---

## Sensitivity to the calling rule

Tataka et al. 2026 report that a gene-count rule performs worse than single-gene or
summed-coverage rules. We kept our pre-registered rule and added the alternatives as a
post-hoc sensitivity analysis.

### CN-CRC

| Rule | CRC | control | P |
|---|---|---|---|
| Pre-registered (breadth >= 0.5 and >= 14/18 genes) | 7/44 | 6/33 | 1.000 |
| *clbB* coverage >= 0.5 | 8/44 | 7/33 | 0.777 |
| *clbB* coverage >= 0.9 | 6/44 | 5/33 | 1.000 |
| >= 8 of 18 genes | 8/44 | 6/33 | 1.000 |
| Read rate >= 41.7 per million | 7/44 | 6/33 | 1.000 |

### AT-CRC

All five rules gave the identical table, 4/20 vs 2/20, P = 0.661.

**No rule produces a significant difference.** The conclusion does not depend on the
calling rule.

---

## Who carries the island

For each positive sample we asked which species could account for the island's depth,
first by marker-gene density ratio, then - for a subset - by mapping to six complete
Enterobacteriaceae genomes in one competitive index. The readout is
**island depth / that species' genome depth**: near 1 means "most of that species
carries one copy".

### DE-CRC, all 9 positives

Eight are consistent with *E. coli* (R 1.37-2.01, which after dividing out the ~1.9x
offset is 0.72-1.06, i.e. one island per genome). One is not:

| Sample | Group | Breadth | Genes | Best fit | R |
|---|---|---|---|---|---|
| **n_CCMD39882286ST** | **control** | 0.778 | 16/18 | ***K. pneumoniae*** | 0.79 (*E. coli*: **723**) |

An *E. coli* R of 723 means *E. coli* cannot account for the island by a factor of
700; the ~1.9x calibration offset is irrelevant at that scale.

### CN-CRC, all 14 positives measured against six complete genomes

Each sample was mapped to six complete Enterobacteriaceae genomes in a single
competitive index at MAPQ 20. "ratio" is island depth divided by that species' genome
depth, from the same read set: near 1 means most of that species carries one copy, and
a ratio well above 1 means that species cannot account for the island at all. Genome
breadth is given because a species with breadth below ~0.1 is not really present.

| Sample | Group | Island depth | *E. coli* depth (breadth) | ratio | *K. pneumoniae* depth (breadth) | ratio | Reading |
|---|---|---|---|---|---|---|---|
| cancer_ERR1018307 | CRC | 391.2 | 407.2 (0.92) | **0.96** | 0.31 (0.02) | 1256 | *E. coli* |
| n_ERR1018246 | control | 14.1 | 13.8 (0.86) | **1.03** | 0.04 (0.01) | 340 | *E. coli* |
| cancer_ERR1018229 | CRC | 4.0 | 4.2 (0.91) | **0.95** | 1.44 (0.67) | 2.75 | *E. coli* |
| cancer_ERR1018213 | CRC | 54.5 | 63.9 (0.95) | **0.85** | 0.14 (0.13) | 387 | *E. coli* |
| cancer_ERR1018198 | CRC | 4.7 | 5.2 (0.86) | **0.91** | 16.7 (0.94) | 0.28 | *E. coli* - **despite abundant Kp** |
| cancer_ERR1018282 | CRC | 3.0 | 3.5 (0.84) | **0.85** | 109.4 (0.89) | 0.03 | *E. coli* - **despite 109x Kp** |
| cancer_ERR1018224 | CRC | 5.7 | 20.3 (0.89) | 0.28 | 0.005 (0.004) | - | *E. coli* by exclusion (nothing else present) |
| n_ERR1018245 | control | 3.9 | 27.6 (0.95) | 0.14 | 1.09 (0.59) | 3.57 | *E. coli* by exclusion (Kp too shallow) |
| **cancer_ERR1018186** | **CRC** | 13.5 | 0.41 (0.31) | **33** | 10.3 (0.88) | 1.30 | ***K. pneumoniae*, assembly-confirmed** ST23/K1/ICEKp10/*clb* |
| **n_ERR1018237** | **control** | 1.7 | **64.7 (0.81)** | **0.03** | 2.04 (0.76) | **0.84** | ***E. coli* cannot explain it** |
| n_ERR1018240 | control | 2.7 | 3.5 (0.72) | 0.76 | 2.74 (0.83) | 0.97 | **undetermined - both fit** |
| n_ERR1018238 | control | 2.5 | 8.3 (0.79) | 0.30 | 8.34 (0.91) | 0.30 | undetermined |
| cancer_ERR1018219 | CRC | 1.6 | 11.5 (0.91) | 0.14 | 7.82 (0.90) | 0.20 | undetermined |
| n_ERR1018236 | control | 2.5 | 11.8 (0.78) | 0.22 | 2.06 (0.76) | 1.23 | undetermined |

**Four of the 14 cannot be assigned a carrier species at all**, even with six complete
genomes competing in one index. This is the part of the carrier analysis that we did
not anticipate and that we think is worth reporting on its own.

`n_ERR1018240` is the clearest case: *E. coli* at ratio 0.76 and *K. pneumoniae* at
ratio 0.97, both with genome breadth above 0.7. **Both species are abundant enough to
carry the island, and the ratio fits both.** There is no coverage-based argument that
separates them. `cancer_ERR1018219` is the same situation in reverse - *E. coli* 11.5x
and Kp 7.8x, both at breadth ~0.9, with the island at only 1.6x, so either species
could carry it in a minority of cells.

**`n_ERR1018237` is the one case where *E. coli* is excluded by coverage.** *E. coli*
sits at 64.7x genome depth while the island is at 1.7x, a ratio of 0.03: if *E. coli*
carried it, only 1 in 38 *E. coli* cells would. *K. pneumoniae* is genuinely present
(2.04x, genome breadth 0.76) at a ratio of 0.84. The marker genes agree independently:
Kp markers total 6,450 bp, so 2.04x predicts ~140 reads and 127 were observed;
*E. coli* markers total 96,300 bp, so 64.7x predicts ~66,000 and 64,656 were observed.

**`n_ERR1018236` and `n_ERR1018238` are undetermined, not *Klebsiella*.** Their
*E. coli* ratios are 0.22 and 0.30, which only requires that 1 in 4 or 1 in 3 *E. coli*
cells carry the island - entirely ordinary for a gut carrying several *E. coli*
strains. Their Kp ratios (1.23 and 0.30) do not fit better. We call these undetermined
rather than counting them toward either species. `n_ERR1018236` also has **no
*K. pneumoniae* marker measurement at all** (NA, see [CORRECTIONS.md](CORRECTIONS.md)
item 7), so there is no independent cross-check for it.

**`cancer_ERR1018198` and `cancer_ERR1018282` are the counterexamples that keep this
honest.** In `cancer_ERR1018198`, *K. pneumoniae* is present at 16.7x with 0.94 genome
breadth - abundant and unambiguous - and is still not the carrier; *E. coli* is, at
ratio 0.91. `cancer_ERR1018282` is more extreme: Kp at **109x** against *E. coli* at
3.5x, and the island still tracks *E. coli* (ratio 0.85) rather than Kp (ratio 0.03,
i.e. 1 Kp cell in 33). **Abundance of a species is not carriage by that species**, and
any method that attributes the island to whichever *Klebsiella* happens to be abundant
will get both of these samples wrong. We made exactly that mistake with
`cancer_ERR1018282` before its genome measurement existed; see
[CORRECTIONS.md](CORRECTIONS.md) item 8.

### Summary of carriers

Carriers were examined in **23 positive samples across two cohorts**: all 14 CN-CRC
positives, against six complete genomes in one competitive index, and all 9 DE-CRC
positives, by marker-gene density ratio only. **The two methods are not equivalent.**
The ratio can exclude a species ("R >> 1, so it cannot account for this") but cannot
establish one, and it cannot produce the "undetermined" verdict at all - so the
4-of-14 figure below belongs to the six-genome subset and has no counterpart among
the 9.

| Attribution | n | Samples |
|---|---|---|
| *E. coli*, ratio 0.85-1.03 with no competing species | 4 | `cancer_ERR1018213`, `cancer_ERR1018229`, `cancer_ERR1018307`, `n_ERR1018246` |
| *E. coli*, despite an abundant *Klebsiella* | 2 | `cancer_ERR1018198`, `cancer_ERR1018282` |
| *E. coli* by exclusion, minority of cells | 2 | `cancer_ERR1018224`, `n_ERR1018245` |
| ***K. pneumoniae***, assembly-confirmed | 1 | `cancer_ERR1018186` (CRC) |
| ***K. pneumoniae***, *E. coli* excluded by coverage | 2 | `n_ERR1018237` (CN control), `n_CCMD39882286ST` (DE control) |
| **Undetermined** (six-genome subset only) | 4 | `cancer_ERR1018219`, `n_ERR1018236`, `n_ERR1018238`, `n_ERR1018240` |
| *E. coli* by marker ratio, DE-CRC (R 1.37-2.01, i.e. 0.72-1.06 after the ~1.9x offset) | 8 | the other 8 DE-CRC positives |

- **Only one sample is confirmed by assembly**: `cancer_ERR1018186`, where the island
  and the *Klebsiella* chromosome appear on a single contig. The other two *Klebsiella*
  attributions rest on coverage ratios, which are consistent with carriage but do not
  prove it. Those FASTQ files were deleted after measurement, so the linkage check has
  not been run on them.
- **Two of the three *Klebsiella* attributions are controls**, which is the opposite of
  what a disease association would predict. Three of 23 positives is 13%.
- **Four of 14 are undetermined.** Species attribution of the island from stool
  metagenomes fails outright in roughly a quarter to a third of positive samples, and
  it fails precisely where two Enterobacteriaceae are co-abundant. Any paper that
  reports pks carriage broken down by species from stool reads should state what
  fraction of its positives it could not assign; we have not seen one that does.

This is concordant with Putze et al. 2009, which found the island in *E. coli*,
*K. pneumoniae*, *E. aerogenes* and *C. koseri* and described it as "highly conserved
among these species". The island is not an *E. coli*-only feature, and attributing it
to *E. coli* by default is a mistake.

---

## Raw output of `code/analyse.py`

```
pks island carriage and colorectal cancer -- reanalysis of 4 public cohorts
samples: 397   target depth: 20M reads

check: breadth recomputed from depth_hist.tsv matches per_sample.tsv
       largest difference over 397 samples = 5.00e-07

judgeable after rarefying and gating: 167 of 397

======================================================================
PRIMARY  rarefied to 20M reads, E. coli >= 0.1% (both pre-specified)
======================================================================
cohort              CRC      control        OR        P
FR-CRC         1/11         3/12          0.39   0.5901
AT-CRC         4/21         3/21          1.36   1.0000
CN-CRC         6/44         2/33          2.13   0.4544
DE-CRC         6/16         0/9          11.76   0.0571

Mantel-Haenszel  OR 1.91  95% CI 0.76-4.75  P = 0.1663
heterogeneity    Q = 3.60  I2 = 17%
smallest cohort P = 0.0571  ->  Bonferroni x4 = 0.2283

======================================================================
reference  rarefied, NO abundance gate (not the pre-specified analysis)
======================================================================
cohort              CRC      control        OR        P
FR-CRC         1/20         3/20          0.38   0.6050
AT-CRC         4/46         3/63          1.83   0.4524
CN-CRC         7/73         2/54          2.37   0.2995
DE-CRC         6/60         0/55         13.24   0.0281

Mantel-Haenszel  OR 2.31  95% CI 0.97-5.51  P = 0.0590
heterogeneity    Q = 4.20  I2 = 29%
smallest cohort P = 0.0281  ->  Bonferroni x4 = 0.1124

======================================================================
reference  gated but NOT rarefied (not the pre-specified analysis)
======================================================================
cohort              CRC      control        OR        P
FR-CRC         3/11         3/12          1.12   1.0000
AT-CRC         5/21         3/21          1.76   0.6965
CN-CRC         7/44         6/33          0.85   1.0000
DE-CRC         6/16         0/10         13.00   0.0532

Mantel-Haenszel  OR 1.59  95% CI 0.73-3.49  P = 0.2440
heterogeneity    Q = 2.97  I2 = 0%
smallest cohort P = 0.0532  ->  Bonferroni x4 = 0.2128
```

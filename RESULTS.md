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

Tataka et al. 2026 (*BMC Microbiol* 26:930,
[doi:10.1186/s12866-026-05460-3](https://doi.org/10.1186/s12866-026-05460-3)) report
that a gene-count rule performs worse than single-gene or summed-coverage rules
(*clbB* alone AUC 0.826, the 19-gene read sum 0.838, ">= 8 genes" 0.734-0.773). We kept our pre-registered rule and added the alternatives as a
post-hoc sensitivity analysis.

### CN-CRC

| Rule | CRC | control | P |
|---|---|---|---|
| Pre-registered (breadth >= 0.5 and >= 14/18 genes) | 7/44 | 6/33 | 1.000 |
| *clbB* coverage >= 0.5 | 8/44 | 7/33 | 0.777 |
| *clbB* coverage >= 0.9 | 6/44 | 5/33 | 1.000 |
| >= 8 of 18 genes | 8/44 | 6/33 | 1.000 |
| Read rate >= 56.7 per million | 7/44 | 6/33 | 1.000 |

The read-rate rule has **no external threshold**: no published cut-off is
transferable to our read counts, so the threshold is set at the point that yields the
same number of positives as the pre-registered rule. It therefore tests whether the
*ranking* reshuffles cases between the groups, not whether a published cut-off
replicates. An earlier version of this table gave the threshold as 41.7 per million,
taken from the mapping log's read count rather than the `island_reads` column shipped
here; the two counts differ, and from the shipped column the threshold is 56.7. The
counts and P are unaffected, because the threshold is defined by the positive count.

### AT-CRC

| Rule | CRC | control | P |
|---|---|---|---|
| Pre-registered (breadth >= 0.5 and >= 14/18 genes) | 5/21 | 3/21 | 0.697 |
| *clbB* coverage >= 0.5 | 5/21 | 3/21 | 0.697 |
| *clbB* coverage >= 0.9 | 4/21 | 3/21 | 1.000 |
| >= 8 of 18 genes | 5/21 | 3/21 | 0.697 |
| Read rate >= 56.3 per million | 5/21 | 3/21 | 0.697 |

Four of the five rules give the same table; *clbB* >= 0.9 moves one CRC sample out.
An earlier version of this section said "all five rules gave the identical table,
4/20 vs 2/20, P = 0.661", which was stale on every number - the denominators are 21
and 21, not 20 and 20 - see [CORRECTIONS.md](CORRECTIONS.md) item 17.

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

Eight are consistent with *E. coli* (R 1.43-2.12, median 1.91, which after dividing
out the ~1.9x
offset is 0.72-1.06, i.e. one island per genome). One is not:

| Sample | Group | Breadth | Genes | Best fit | R |
|---|---|---|---|---|---|
| **n_CCMD39882286ST** | **control** | 0.778 | 16/18 | ***K. pneumoniae*** | 0.79 (*E. coli*: **764**) |

An *E. coli* R in the hundreds means *E. coli* cannot account for the island by that
factor; the ~1.9x calibration offset is irrelevant at that scale. Recomputed from the
shipped columns it is **764**; the figure originally published here was 723, from
MAPQ-filtered island reads that were not shipped. Both are reproduced by
`python3 code/verify_extra.py`, which also recomputes every other R in this table.

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
| cancer_ERR1018282 | CRC | 3.0 | 3.5 (0.84) | **0.85** | 109.4 (0.89) | 0.0273 | *E. coli* on parsimony - Kp not excluded |
| cancer_ERR1018224 | CRC | 5.7 | 20.3 (0.89) | 0.28 | 0.005 (0.004) | - | *E. coli* by exclusion (nothing else present) |
| n_ERR1018245 | control | 3.9 | 27.6 (0.95) | 0.14 | 1.09 (0.59) | 3.57 | *E. coli* by exclusion (Kp too shallow) |
| **cancer_ERR1018186** | **CRC** | 13.5 | 0.41 (0.31) | **33** | 10.3 (0.88) | 1.30 | ***K. pneumoniae*, assembly-confirmed** ST23/K1 (automatic bin); ICEKp10 and a complete *clb* 2 only with the island contigs added, which are placed in this genome by the island tiling and by the element's two flanks rather than by the typing - see [CORRECTIONS](CORRECTIONS.md) item 16 |
| **n_ERR1018237** | **control** | 1.7 | **64.7 (0.81)** | **0.03** | 2.04 (0.76) | **0.84** | ***E. coli* cannot explain it** |
| n_ERR1018240 | control | 2.7 | 3.5 (0.72) | 0.76 | 2.74 (0.83) | 0.97 | **undetermined - both fit** |
| n_ERR1018238 | control | 2.5 | 8.3 (0.79) | 0.30 | 8.34 (0.91) | 0.30 | undetermined |
| cancer_ERR1018219 | CRC | 1.6 | 11.5 (0.91) | 0.14 | 7.82 (0.90) | 0.20 | undetermined |
| n_ERR1018236 | control | 2.5 | 11.8 (0.78) | 0.22 | 2.06 (0.76) | 1.23 | undetermined |

**The other four genomes were checked and none of them is ever a better fit.** The
table shows only *E. coli* and *K. pneumoniae* because those two are the only
candidates that ever come close, but all six were measured in the same index. In no
sample does any of the other four have a genome breadth above 0.3 **and** a ratio
closer to 1 than the species reported here; the closest any of them comes is
*C. freundii* CFNIH1 in `cancer_ERR1018198`, breadth 0.496 at ratio 1.37 against
*E. coli*'s 0.91. *C. koseri* BAA-895 needs a separate note: it **carries a *clb*
island of its own** (94.88% of this reference's 40-mers, all 18 genes), and its genome
breadth never exceeds 0.028 in any sample - about the 1.1% of its genome that the
island accounts for. Its apparent presence is the island, not the organism.
`python3 code/verify_extra.py` prints all six per sample, and
`data/clb_in_six_genomes.tsv` holds the island content of each reference genome.

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

**`cancer_ERR1018198` and `cancer_ERR1018282` are the samples where abundance pulls
hardest in the wrong direction.** In `cancer_ERR1018198`, *K. pneumoniae* is present at
16.7x with 0.94 genome breadth - abundant and unambiguous - while the island sits at
4.7x and *E. coli* at 5.2x, ratio **0.911**. In `cancer_ERR1018282` the pull is extreme:
Kp at **109x** against *E. coli* at 3.5x, and the island at 3.0x, ratio **0.850** on
*E. coli*. **Abundance of a species is not carriage by that species**, and a method that
attributes the island to whichever *Klebsiella* happens to be abundant will read both of
these samples wrong. We made exactly that mistake with `cancer_ERR1018282`; see
[CORRECTIONS.md](CORRECTIONS.md) item 8.

**But the depth ratio does not exclude *Klebsiella* in either sample, and this paragraph
used to say it did.** The ratio on Kp is 0.0273 in `cancer_ERR1018282` - 1 *Klebsiella*
cell in 37 - and 0.283 in `cancer_ERR1018198` - 1 cell in 3.5. Those are ordinary
frequencies for a mobile genetic element, not impossibilities, and our own item 16 shows
an island on a *Klebsiella* chromosome in a sample where *Klebsiella* is abundant.
`cancer_ERR1018198` is worse still: *C. freundii* is at 3.4x, ratio **1.37**, so three
species fit its depths. (The text here read "1 Kp cell in 33", which comes from
inverting the rounded 0.03 rather than the measured 0.0273.)

**What survives, and what does not.** What survives is the negative half: ranking
candidates by abundance is the wrong operation, and these two samples demonstrate it.
What does not survive is the positive half - the attribution of these two samples to
*E. coli* rests on parsimony, because *E. coli* is the only candidate that needs no
minority subpopulation. Parsimony is not evidence. The scale of this ratio has now been
measured (`data/island_ratio_scale.tsv`, [CORRECTIONS.md](CORRECTIONS.md) item 8
amendment): across the four samples where five of the six species are excluded at
>= 100x it runs 0.278-1.026, so a value near 1 is where a fully-carrying population
sits - but the band spans 3.7x on n = 4, and **there is no calibrator for *Klebsiella* at
all**. Settling these two needs what settled `cancer_ERR1018186`: an assembly that puts
the island on a chromosome. See [CLAIMS.md](CLAIMS.md) row 8.

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
| *E. coli* on parsimony, with an abundant *Klebsiella* **not excluded** | 2 | `cancer_ERR1018198`, `cancer_ERR1018282` - the depth ratio fits a *Klebsiella* minority (1 cell in 3.5 and 1 in 37) as well as it fits *E. coli*, and in `cancer_ERR1018198` *C. freundii* fits too (ratio 1.37). Graded **C** in [CLAIMS.md](CLAIMS.md) |
| *E. coli* by exclusion, minority of cells | 2 | `cancer_ERR1018224`, `n_ERR1018245` |
| ***K. pneumoniae***, assembly-confirmed | 1 | `cancer_ERR1018186` (CRC) - contig `k119_45171`, 386,079 bp: *Klebsiella* chromosome to ~354,276 then the island to the end, with read depth never falling below 12 across the junction. The other three island contigs tile the same island with 119 bp overlaps - exactly the assembly's k - so they are one island with this contig, and the whole 137.7 kb insertion is flanked on both sides by contigs the automatic assignment placed in *Klebsiella*. See `data/assembly_alignment.tsv`, `data/kleborate/island_tiling.tsv`, `data/kleborate/island_locus.tsv` and `python3 code/verify_extra.py` |
| ***K. pneumoniae***, *E. coli* excluded by coverage | 2 | `n_ERR1018237` (CN control), `n_CCMD39882286ST` (DE control) |
| **Undetermined** (six-genome subset only) | 4 | `cancer_ERR1018219`, `n_ERR1018236`, `n_ERR1018238`, `n_ERR1018240` |
| *E. coli* by marker ratio, DE-CRC (R 1.43-2.12, i.e. ~0.75-1.11 after the ~1.9x offset) | 8 | the other 8 DE-CRC positives |

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

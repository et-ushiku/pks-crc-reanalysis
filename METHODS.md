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
  Conversion: **313 reads per 1% per million reads** (`PER_PCT`), calibrated once on
  the *E. coli* marker set and then held fixed. Reads landing on a marker set scale
  with **abundance x marker length**, so the constant is scaled per species by
  `marker_bp / 96,300`: for *K. pneumoniae*'s 6,450 bp that is 20.96 reads per 1% per
  million, not 313. Using 313 for both understates *K. pneumoniae* by
  96,300 / 6,450 = **14.93x**, which is what the published column did until
  2026-11-26 - see [CORRECTIONS.md](CORRECTIONS.md) item 17.
- Whether that length scaling is right was **measured, not assumed**
  (`code/verify_marker_scaling.py`). The six-genome competitive index measures the same
  two species in the same samples by a route that does not involve marker length, so
  each species' marker read density per effective marker base can be compared with its
  own genome's read density per base. The length scaling requires that this multiplier
  be **the same for both species**. Over the 5 samples where both genomes are covered
  at breadth >= 0.7 it is: *K. pneumoniae* / *E. coli* = **0.929, 95% CI 0.662-1.304**,
  against **14.93** if the scaling should not have been applied. The interval contains
  1 and excludes 14.93 by a factor of 11. The multipliers themselves are near 1 for
  both species once the end-to-end edge loss is accounted for (*E. coli* 0.92-1.23,
  *K. pneumoniae* 0.83-1.27), which is the physical reason the scaling holds: the
  marker sets recruit reads at about the rate their own genome does. **The interval
  spans a factor of 2.0**, so the scaling is confirmed only to that coarseness - this
  column must not be used to argue a two-fold difference - and it is confirmed only
  for these two marker sets, in CN-CRC, in 5 samples.
- Marker edge loss: `bowtie2` is end-to-end, so a read overhanging a marker fragment
  loses its alignment; the effective marker length is `marker_bp - fragments x
  (read_len - 1)`. The two sets differ here (122 fragments over 96,300 bp against 11
  over 6,450 bp), and including the term moves the multiplier ratio from 0.886 to
  0.929. The conclusion does not depend on the read length: at 84, 94 and 100 bp the
  ratio is 0.924, 0.929, 0.933.
- Group comparisons of reads and abundance: two-sided Mann-Whitney U.
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

## How the calling criteria compare with Piccinno et al. 2025

Piccinno et al. 2025 is the comparator in the README, so it matters how its calls are
made. We looked, because we wanted to state the comparison honestly. **The two methods
cannot be compared on sensitivity, because theirs does not have a number.**

Their paper says only:

> Anpan first filters samples to remove any without enough species-specific gene
> coverage to accurately assess gene-level effects.

No threshold is given there or anywhere else in the text. The rule is in the tool.
From Anpan v.0.3.0 (`R/filtering.R`, `filter_with_kmeans`; defaults
`filtering_method = "kmeans"`, `discard_poorly_covered_samples = TRUE`):

- For each species and sample it computes `n_nz`, the number of gene families with
  nonzero abundance, and `q50`, the median log abundance.
- It then runs **2-means clustering** on those two standardised statistics and keeps
  the samples that land in the "present" cluster. **The decision boundary is fitted to
  whatever samples are in hand.**
- If a reference-genomes file is supplied there is an additional hard floor on `n_nz`:
  `mean(gene count) - 2 x SD(gene count)` when there are at least 5 genomes, otherwise
  `(2/3) x mean(gene count)`. Samples below it are forced to "absent".
- With more than 50,000 gene families, `n_nz` is `log1p`-transformed first.

This is a reasonable design - it adapts to each species rather than imposing one
number on all of them - but it has a consequence worth stating: **the same metagenome
can be called differently depending on which other samples are in the study**, because
the k-means centres move. A threshold like ours is cruder and almost certainly less
sensitive, but it is a property of the sample. Theirs is a property of the dataset.

For the same reason, "they called 3,741 samples and we called 397" is not a
like-for-like statement, and we do not make it. What we claim is agreement of
conclusions, not of methods.

Three further things we could not resolve from their text, recorded so nobody repeats
the search:

- **Which genes.** The word "colibactin" appears exactly twice in the full text. The
  specific UniRef90 gene families used for the colibactin result are not named; they
  would be in Supplementary Table 9, which we have not opened.
- **The significance threshold is stated three ways.** Methods, for the Anpan gene
  model: "a significance threshold of an absolute coefficient of 2 and q < 0.05".
  Results and the Fig. 5a legend: "FDR global q < 0.05 and absolute coefficient
  estimate >1". The colibactin sentence itself: "q > 0.01 and abs(coefficient
  estimate) < 1". The first two are criteria for calling something significant; the
  third is the bound reported for a null.
- **How close the colibactin result was.** Their null does hold under their own rule,
  because `abs(coefficient) < 1` fails a criterion that requires it to exceed 1 (or 2,
  per Methods). But `q > 0.01` corresponds to no threshold used anywhere else in the
  paper, and **the actual q value and coefficient are never reported**, so the distance
  from significance cannot be read off. We rely on their conclusion, not on a number.

## Carrier attribution, assembly and typing

This is the part of the work RESULTS.md leads with, and until recently its methods
were not written down here at all (CORRECTIONS.md items 15 and 16).

**Step 1, density ratio.** For every positive sample, R = (island reads / 49,455 bp)
divided by (marker reads / marker bp) for the candidate species. R near 1 means that
species could carry one island copy per genome; R in the hundreds means it cannot
account for the island at all. R carries a systematic offset of about 1.9x, decomposed
in `code/verify_extra.py`; the offset is irrelevant at the scale that matters for
exclusion.

**Step 2, six competing genomes.** For the 14 CN-CRC positives, the reads were mapped
to a **single** index holding six complete Enterobacteriaceae genomes - *E. coli*
IHE3034 (CP001969.1), *K. pneumoniae* HS11286 (CP003200.1), *C. freundii* CFNIH1
(CP007557.1), *C. koseri* BAA-895 (CP000822.1), *E. cloacae* ATCC13047 (CP001918.1),
*K. aerogenes* KCTC2190 (CP002824.1) - at MAPQ 20, so the species compete and
cross-species bleed is filtered rather than counted twice. The readout is island depth
divided by that species' genome depth, from the same read set. `pipeline/asm_genomes.sh`.

**Two properties of the reference set a reader needs.** Only *E. coli* IHE3034 (100%)
and *C. koseri* BAA-895 (94.88%) carry a *clb* island; the other four share under
0.15% of its 40-mers and none of its 18 genes (`data/clb_in_six_genomes.tsv`). So
(a) *C. koseri*'s apparent genome breadth in these samples - never above 0.028 - is
mostly the 1.1% of its genome that the island accounts for, not the organism, and
(b) a contig made of *clb* can never be assigned to *Klebsiella* by this index. Point
(b) is what forced the second Kleborate run described below.

**Step 3, assembly.** `pipeline/asm_assemble.sh` assembles one sample with megahit,
aligns the island to the contigs to find those carrying it, and aligns those contigs
to the six genomes with `minimap2 -x asm10`. Per contig, matched bases are summed per
species and the species with the most wins, with the winner required to cover at least
half the contig. Summing before dividing matters: one contig aligns in several
fragments, so a per-row coverage test rejects everything.

**Step 4, is the join real.** A contig that carries chromosome and island contiguously
could still be an assembler chimera. `pipeline/asm_junction.sh` remaps the reads to
that contig and counts pairs and single reads spanning the junction, and the depth
across it. A chimera gives zero spanning evidence and a depth collapse; this one has a
minimum depth of 12 and **no base at zero** (`data/assembly_junction_depth.tsv`).

**Step 5, typing.** Kleborate 3.2.4, `-p kpsc`, on the species-assigned bin
(`pipeline/asm_kleborate.sh`), and then on that bin plus the island-bearing contigs
(`pipeline/asm_kleborate_plus.sh`). Both outputs are in `data/kleborate/`. **They do
not agree, and the difference matters**: species, ST and K type come from the
automatic bin, while ICEKp10 and a complete *clb* 2 appear only in the second run.
A typing from a hand-assembled bin is not independent evidence of what the added
contigs belong to - that evidence is step 4, and two further measurements added on
2026-11-28 (`code/make_locus_table.py`, checked by `verify_extra.py` sections 8 and 9):

- **Island tiling.** The four island contigs are laid out on the island reference
  AM229678.1. Junctions between different contigs are reported as overlaps; an overlap
  of exactly k (=119, this assembly's final k-mer length) means the two contigs are
  consecutive nodes of one path in the assembly graph. All four such junctions are
  exactly 119 bp, and both breaks fall on the two copies of the duplication of
  PITFALLS item 1. Thresholds: alignments at >= 99.0% identity over >= 1,000 bp.
- **The element's flanks.** The element containing the island is tiled on the only
  reference that carries the whole of it (*C. koseri* BAA-895, used as a coordinate
  ruler and not as a statement of origin). Its extent is grown from the island-bearing
  contigs by overlap, allowing gaps of up to 200 bp, rather than from a window chosen
  by hand. The contigs reaching its two ends are then checked for membership of the
  automatic bin, and their *K. pneumoniae* chromosome intervals for adjacency.

- **Read depth and read-pair linkage.** The reads are mapped back to the whole
  assembly - all 21,759 contigs in one bowtie2 index, so a read goes to its best place
  anywhere in the assembly rather than being forced onto the bin - and filtered at
  MAPQ 20 (`pipeline/asm_linkage.sh`). Per-contig read depth replaces megahit's
  `multi`, which is a k=119 graph multiplicity against 100 bp reads and is not a read
  depth. Pairs whose mates land on different contigs are counted, which requires not
  filtering on proper-pair: a pair spanning two contigs is by definition not
  concordant. Partners are restricted, before counting, to the 160 contigs the
  six-genome step assigned at breadth >= 0.5. Each added contig is then compared with
  the bin's own contigs **in the same length band**, because below 10 kb the bin's own
  contigs span the whole range and the test has no power.

See CORRECTIONS.md item 16.

**Step 6, do the added contigs belong to that genome.** Alignments are split at the
island boundary, where the boundary is set by exact 40-mer matches grouped into
contiguous blocks - not first-hit-to-last-hit, which inflates a scatter of hits into a
span. Taken alone this step is **negative**: `k119_70945`'s non-island portion matches
*C. koseri* better than *K. pneumoniae* (0.928 against 0.818 of matched bases over
aligned length - a PAF ratio that includes gaps, so usable for comparing species on
one contig but not an identity), so similarity does not place it. What does place it is megahit's per-contig *k*-mer coverage depth compared
with the bin's own: the bin averages 18.08, the island contigs sit at 18.01, 18.99 and
19.83, the collapsed duplication at 36.00, and *C. koseri* is absent from the sample.
That is coverage plus the absence of the alternative, not similarity, and it is
weaker than step 4. `code/make_backbone_table.py`, checked by section 7 of
`code/verify_extra.py`.

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
- No per-strain analysis, no assembly for most samples, no wet-lab confirmation. One
  sample was assembled and typed; its ICEKp10 call depends on a bin we assembled by
  hand, and nothing here is confirmed in culture.
- No correction for human-read contamination beyond an upper bound (the cohorts
  differ in host-removal preprocessing; one cohort's exact preprocessing is
  unpublished).

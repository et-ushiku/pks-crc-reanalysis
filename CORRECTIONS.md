# Claims we made and then had to withdraw

This file exists because the interesting part of a negative result is whether the
person reporting it was willing to be wrong in public. Each item is a statement we
wrote down, believed, and later refuted with our own data. Nothing here was deleted
from the record; the superseded statements are kept so the reasoning can be audited.

---

## 1. "The breadth ceiling is 0.9589"

**What we said.** That 0.9589 was the maximum attainable breadth, because the
duplicated region is always lost.

**Why it was wrong.** A sample came in at 0.9598. Then others above that, up to
0.9607.

**Root cause.** We took the highest value we had seen in a handful of deep samples and
promoted it to a derived constant. It was an observation, not a derivation.

**What is actually derivable.** If both copies of the duplication are lost entirely,
breadth is `(49455 - 2574) / 49455 = 0.9480`. (This paragraph originally said 1,248 bp
and 0.9495; the repeat was later measured at 1,287 bp - see item 9. The denominator
49,455 is itself 24 bp too large, because four adjacent gene pairs overlap - see item
18; against the 49,431 distinct bases the figure is 0.9479.) Reads that happen to span a
position where the two copies differ (they are 99.840% identical, so 0.16% of
positions) can still be placed, which lifts breadth above that figure by a
sample-dependent amount. **The ceiling is not a fixed number.** The defensible
statements are: 0.9480 is what you get if the duplication is lost completely, and no
sample in 397 exceeded 0.99.

(That sentence said **0.9495** until 2026-11-28, four lines below the paragraph that
had already been corrected to 0.9480 and three sections above item 9 saying in so many
words "not 0.9495 as stated in item 1 of this file". The retraction was written and the
sentence it retracted was left standing. Nothing turned on it - no sample in 397 falls
in the 1.5-thousandth gap between the two values, so the counts in item 2 are the same
either way, which is checked - but the file contradicted itself for three commits.)

---

## 2. "In a carrier, the uncovered bases lie only inside the duplication"

**What we said.** That for any sample with breadth >= 0.5, every uncovered base would
fall inside the duplicated region.

**Why it was wrong.** We checked every positive sample. **Only 11 satisfied it.**
Twenty-six had breadth below the derived ceiling while still passing the >= 14/18 gene
rule, which is impossible if the duplication is the only thing missing. (The counts
here were originally given as 9 and 25 out of 44 positives; the audit in item 13 found
the correct figures to be 11 and 26 out of 42.)

**Root cause.** Same as item 1: we generalised from the three or four deepest samples
in hand, which happened to be the ones where it was true.

**What replaced it.** The duplication loss is near-constant and depth-independent
(1,920-2,140 bp in the 13 samples at >= 10x island depth), while everything outside
the duplication is depth-limited. At >= 10x, 11 of 13 samples have zero missing bases
outside the duplication; below 2x the median missing runs to thousands of bases.
This is now [PITFALLS.md](PITFALLS.md) item 4, and it is the reason the rarefied
analysis is primary rather than a sensitivity check.

---

## 3. The density ratio R was treated as if it were a population fraction

**What we said.** We computed `R = (island reads / 49,455) / (marker reads / marker
bp)` and read values near 1 as "this species carries one island per genome", without
having established the scale.

**Why it was wrong.** R clustered near 1.9 for samples where *E. coli* was plainly the
carrier. A ratio with a 1.9x offset cannot be read as a fraction.

**Root cause.** We validated the ratio's ordering and never its scale.

**What we did.** Decomposed the offset: 1.115 (reads counted over the 55,140 bp
reference but divided by the 49,455 bp gene total) x 1.14 (marker fragments lose
coverage at their ends) x ~1.50 (residual mappability difference), product ~1.9,
stable across read lengths of 84, 94 and 100 bp. R is now used only for the one-sided
inference "R >> 1 means that species cannot account for the island", which a 1.9x
offset cannot overturn. R never entered the pre-specified calling criteria, so no
result changed.

---

## 4. "There was no difference" in CN-CRC, stated with the unrarefied numbers

**What we said**, at the head of the CN-CRC write-up: no difference; CRC 15.9% (7/44)
vs control 18.2% (6/33), OR 0.85, Fisher P = 1.00.

**Why it was wrong.** Those are the **unrarefied** numbers, and our own pre-specified
primary analysis rarefies to a common depth first. Rarefied, the same cohort gives
6/44 vs 2/33, OR **2.13**, P = 0.454. **The direction reverses.**

**Root cause.** The headline was written before the rarefied analysis existed and was
never revisited when it did. The correction was recorded further down the document
while the opening paragraph still asserted the old number, which is the version
anyone would actually read.

**What is correct.** Neither result is significant, so "no significant difference"
stands. But "there was no difference" overstates it, and "fewer carriers in cancer" is
not supported. At 6-8 positives the direction is not determined.

---

## 5. We looked at a group comparison before the cohort had finished

**What happened.** We had committed in writing to completing all 120 DE-CRC samples
before looking at any group difference. While testing that the cohort was wired into
the analysis script correctly, we ran the script, which prints per-cohort P values, at
77 of 120 samples. So we saw the group split early.

**What we did about it.** Recorded it in the project log with a timestamp. No threshold
was changed, no sample was added or dropped, and the cohort was run to completion. But
it is a real breach of the procedure we had set, and a reader is entitled to weigh the
DE-CRC result accordingly: it is the cohort nearest to significance.

---

## 6. A download bug that would have produced a half-depth sample

Not a claim, but a near-miss worth recording. Our fetch script expanded a
semicolon-separated URL list by piping `printf` output through `tr`. With no trailing
newline, `wc -l` returns 1 instead of 2 and a `while read` loop discards the final
line. The script would have downloaded only the first mate file, counted it as
complete, written its completion marker, and handed a half-depth sample to the
pipeline as if it were whole.

It was caught before any sample was measured, because the file count was checked
against the ENA manifest rather than against the script's own loop. Fix: expand URLs
into a file, one per line, and count with `grep -c .`. The lesson we drew is to verify
completeness against the external source of truth, never against the code that just
did the work.

---

---

## 7. "Not measured" was published as "zero" (found while preparing this release)

**What was wrong.** `data/per_sample.tsv` reported `kpneumoniae_pct = 0.000000` for 27
CN-CRC samples that have **no marker table at all**. Those samples were run through a
dedicated *E. coli* index before the combined Enterobacteriaceae index existed, so the
*K. pneumoniae* marker counts were never computed. The generator read a missing file
as zero hits and divided by the read count, producing a confident 0.000000.

**Root cause.** `marker_pct()` initialised `hit = 0` and only overwrote it if the
counts file existed. A missing input and a true zero took the same path.

**Why it matters.** "No *Klebsiella* detected" and "*Klebsiella* never looked for" are
different statements, and the second one was being published as the first. Anyone
computing *Klebsiella* prevalence from that column would have reported 27 false
negatives.

**What changed.** `marker_pct()` now returns NA when the counts file is absent. The
column is NA for 176 samples (149 that were never run through the combined index at
all, plus these 27) and a real value for 96.

**Effect on the result: none.** The abundance gate uses *E. coli*, which is measured
for all 397 samples from a dedicated path (0 NA). The primary analysis is unchanged:
167 judgeable, Mantel-Haenszel OR 1.91, 95% CI 0.76-4.75, P = 0.1663, I2 = 17%,
verified by re-running `code/analyse.py` after the fix.

**How it surfaced.** Two control samples, `n_ERR1018236` and `n_ERR1018237`, had almost
identical *K. pneumoniae* genome coverage (depth 2.063 / breadth 0.7643 and
2.036 / 0.7638) while their marker columns read 0.000000 and 0.006483. (That second
figure is now 0.096798: the column was rescaled in item 17. The zero is now `NA`.) A genome signal
that strong with no marker reads should be impossible, so one of the two numbers had to
be wrong. It was the zero.

## 8. We read `cancer_ERR1018282` as *Klebsiella* because *Klebsiella* was abundant

**What we said.** That `cancer_ERR1018282` leaned toward *K. pneumoniae*, on the
grounds that its Kp relative abundance was 6.83% against *E. coli* at 0.178% - a factor
of 38.

**Why it was wrong.** The six-genome measurement puts Kp at 109x genome depth and
*E. coli* at 3.5x, while the island sits at 3.0x. The island therefore matches
*E. coli* almost exactly (ratio 0.85, i.e. ~85% of *E. coli* cells carrying one copy)
and matches Kp not at all (ratio 0.03, i.e. 1 Kp cell in 33). **The more abundant
species is the one the island does not track.**

**Root cause.** We ranked candidates by how much of each species was present, when the
question is what fraction of each species carries the island. Those are different
quantities, and for a mobile element they can point in opposite directions. Abundance
was the only evidence available for this sample at the time, which is a reason to have
marked it undetermined rather than leaning.

**What replaced it.** `cancer_ERR1018282` is now attributed to *E. coli*, and it is
reported in [RESULTS.md](RESULTS.md) alongside `cancer_ERR1018198` as a counterexample
to abundance-based attribution - the case where the pull of abundance is strongest
(109x) and still wrong. Our *Klebsiella* count for CN-CRC is unchanged at 2, because
the sample we removed from the *Klebsiella* side was never counted there in the final
tally; what changed is that it is now positively assigned rather than left leaning.

### Amended 2026-11-28: a depth ratio cannot say "85% of cells", and it does not exclude *Klebsiella*

Two statements above claim more than a depth ratio can carry, and the first of them is
item 3 of this file repeated.

**"Ratio 0.85, i.e. ~85% of *E. coli* cells carrying one copy."** Item 3 is the
correction where we read an uncalibrated ratio as a population fraction, found it was
offset by 1.9x, and wrote that the ratio would from then on be used only for the
one-sided inference. That was the marker-based ratio R. The ratio in item 8 is a
different quantity - island read depth over that species' genome depth, both from the
six-genome competitive index - and we read *it* as a population fraction five sections
later. Same mistake, different numerator.

So we measured its scale, the way item 3 measured R's offset
(`code/make_ratio_scale.py`, `data/island_ratio_scale.tsv`). The rule was fixed before
looking: a species is excluded as the carrier when island depth is >= **100x** its
genome depth, and a **calibrator** is a sample where five of the six species are
excluded, so the sixth is the carrier beyond doubt. **No condition is placed on the
ratio itself** - conditioning on it would select values near 1 and make the calibration
circular.

Four samples calibrate *E. coli*:

| sample | island depth | *E. coli* depth | island / *E. coli* | other five excluded |
|---|---|---|---|---|
| `cancer_ERR1018224` | 5.646 | 20.278 | **0.278** | yes |
| `cancer_ERR1018213` | 54.509 | 63.859 | 0.854 | yes |
| `cancer_ERR1018307` | 391.191 | 407.178 | 0.961 | yes |
| `n_ERR1018246` | 14.129 | 13.775 | **1.026** | yes |

**The scale is near 1, not 1.9.** The ratio never runs past 1.03 in a sample where
*E. coli* must be the carrier, so unlike R this quantity has no large multiplicative
offset. One offset in it is measurable and small: the island's depth is averaged over
all 49,431 gene bases, including the ~2,000 in the repeat interior that MAPQ 20 can
never cover, which biases it **down** by 1.041-1.046x in these four. Correcting only
that, the ceiling for a population in which every cell carries one copy sits at about
**1.07**.

**But the band is 0.278-1.026 - a spread of 3.7x on n = 4.** 0.85 is 17% below the
ceiling, and reading that 17% as "15% of cells do not carry it" needs a precision this
measurement does not have. The 0.278 matters in the other direction: it shows the ratio
does carry real information below 1 and is not all offset. **"~85% of cells" was never
measurable, and is withdrawn.**

**There is no calibrator for *Klebsiella*. Zero samples.** Not one of the 14 has
*E. coli* and the other four excluded, so this ratio's scale is measured for *E. coli*
only. Any reading of a *Klebsiella* ratio as copies per cell - including the 1.30 in
the 2026-11-28 amendment to item 16 - is an assumption that two genomes behave alike
under this measurement, not a measurement. Item 16's argument does not rest on it: its
depth evidence was replaced by overlap, position and linkage precisely because depth
could not carry the claim.

**The second statement is the worse one. "Matches Kp not at all (ratio 0.03, i.e. 1 Kp
cell in 33)" is not an exclusion of anything.** One *Klebsiella* cell in 37 carrying a
mobile genetic element is an entirely ordinary frequency for a mobile genetic element.
It is low, not impossible. And item 16 of this file shows an island sitting on a
*Klebsiella* chromosome contig in a sample where *Klebsiella* is abundant, so a minority
*Klebsiella* subpopulation is not a hypothesis we are entitled to wave away. (The
figure is also slightly wrong: the ratio is 0.0273, so 1 in 36.6, not 1 in 33. "1 in
33" is what you get by inverting the rounded 0.03 - item 17's class of error in
miniature.)

**What this means for `cancer_ERR1018282`.** The depth measurement fits both
candidates: 85-100% of a 3.5x *E. coli* population carrying one copy, or 2.7% of a 109x
*Klebsiella* population carrying one copy. **The ratio excludes neither.** The
attribution to *E. coli* rests on parsimony - one hypothesis needs no subpopulation -
and parsimony is not evidence. Item 8 is still right about what was wrong: abundance
was the wrong quantity and leaning on it was wrong. What is weaker than the text claims
is the replacement verdict.

**So we assembled it, and the assembly does not settle it either.** The reads for
`cancer_ERR1018282` are still on disk, so we ran the same pipeline that settled item 16
(`pipeline/asm_assemble.sh`, 63 minutes): 44,116 contigs, longest 416,583 bp.
*Klebsiella* assembled beautifully, as 109x genome depth should - a 185,848 bp contig at
91% *K. pneumoniae* and a 39,109 bp contig at 98%, at megahit multiplicity 91.3 and
83.0.

**Eight contigs hit AM229678.1, and four of them carry no *clb* gene base at all.** The
reference is 55,140 bp of which the 18 genes are 49,431; the other 5,685 bp is
non-island flanking sequence (PITFALLS item 2). Those four contigs' hits are entirely
inside it - island coordinates 136-1,672 and 53,807-54,032 - and **two of the four are
the big *Klebsiella* contigs**:

| contig | bp | multi | island hit | of which *clb* gene | top species | verdict |
|---|---|---|---|---|---|---|
| `k119_42126` | 185,848 | 91.3 | 1,537 | **0** | Kp, 91% | *Klebsiella* chromosome, **backbone hit only** |
| `k119_49689` | 39,109 | 83.0 | 226 | **0** | Kp, 98% | *Klebsiella* chromosome, **backbone hit only** |
| `k119_2268` | 2,056 | 5.0 | 3,343 | 3,333 | *E. coli* / *C. koseri* | island fragment (the repeat) |
| `k119_77516` | 1,577 | 4.0 | 1,577 | 1,577 | *E. coli* / *C. koseri* | island fragment |
| `k119_24080` | 1,483 | 3.0 | 1,515 | 1,098 | *E. coli* / *C. koseri* | island fragment |
| `k119_58287` | 1,278 | 3.0 | 1,156 | 1,156 | *C. koseri* / *E. coli* | island fragment |
| `k119_115023` | 1,764 | 3.0 | 1,039 | **0** | *C. koseri* | backbone hit only |
| `k119_38310` | 1,373 | 3.0 | 78 | **0** | *E. coli* | backbone hit only |

**Counting the first two would have produced "the island is on *Klebsiella*" - in the
sample where we most wanted an answer - out of nothing but the reference's flanking
sequence.** That is PITFALLS item 2 exactly, and the only thing that stopped it was
splitting each hit into gene bases and non-gene bases before looking at the species
column. It is worth being explicit that the temptation was there and the check is what
caught it.

**The four contigs that do carry genes cover 7,164 of 49,431 gene bases - 14.5% - in
fragments of 1,278-2,056 bp with no chromosome attached to any of them.** There is no
island-chromosome junction anywhere in the assembly, so nothing can be attributed.

**And the multiplicities say why, which is the useful part.** The island's contigs sit
at **3-5** against the *Klebsiella* chromosome's **83-91**. At 3x, a 55 kb element
cannot be assembled through and a junction carried by a minority of molecules cannot be
bridged. The same pipeline on `cancer_ERR1018186`, where the island was at 13.5x,
produced a 386 kb *Klebsiella* contig carrying the island and a 96.5% tiling. **The
difference between the two samples is the island's own depth, not the carrier's
abundance** - and a reader should take that as the general rule for when this method
can be expected to work.

**What the assembly does establish: the island is not in the bulk *Klebsiella*
population.** If every *Klebsiella* cell carried it, its contigs would sit near 90 like
the chromosome. They sit at 3-5. So a *Klebsiella* explanation requires a minority of
roughly 3% - which is what the depth ratio said in the first place. The assembly
confirms the depth ratio and does not break the tie.

**Where we are leaving it, and the number for a reader who disagrees.** The attribution
stays *E. coli*, on parsimony, in a category of its own that says so - because the ratio
does positively put the island at the *E. coli* single-copy ceiling (0.850 against a
measured ceiling near 1.0) while the *Klebsiella* explanation requires a subpopulation
nobody has seen. That is weaker than exclusion and we are not calling it more. By the
stricter standard this item's own root cause sets - "a reason to have marked it
undetermined rather than leaning" - these two samples belong in **undetermined**, which
would move RESULTS.md's "4 unassignable of 14" to **6 of 14**. The figure is there so a
reader can apply that standard without recomputing anything.

**Settling it needs what we do not have**: long reads over the insertion, or a library
deep enough to put the island above about 10x in this sample. Neither exists, and we
are not going to get them from a public archive.

**What does not change.** Our *Klebsiella* count for CN-CRC is still 2 - this sample was
never counted on the *Klebsiella* side, and an undetermined sample is not a *Klebsiella*
sample. The primary analysis does not use carrier attribution at all: 167 judgeable, MH
OR 1.91, 95% CI 0.76-4.75, P = 0.1663.

**Shipped with this**: `data/island_ratio_scale.tsv` and
`data/assembly_282_island.tsv`, generated by `code/make_ratio_scale.py` and
`code/make_asm282_table.py`, and recomputed with their assertions by
`python3 code/verify_extra.py` (sections 11 and 12).

---

## 9. "None of these four pitfalls is documented anywhere"

**What we said**, in the README and at the head of PITFALLS.md: that the four failure
modes were ones we could not find documented in the literature. **We published that
claim before checking it properly.**

**Why it was wrong.** A targeted literature check done after publication found
published antecedents for three of the four:

- **Item 1, the duplication.** Auvray et al. 2021 (*Microb Genom* 7:000579) report
  twelve strains with a 5,651 bp deletion fusing *clbJ* and *clbK*, and state it
  "presumably resulted from recombination between two copies of a 1480 bp homologous
  sequence located in *clbJ* and *clbK*". That is the repeat. Our own measurement of
  the spacing between the two copies is 5,653 bp, within 2 bp of their deletion, so it
  is unmistakably the same feature. They say nothing about read mapping, and their one
  caution is about assembly, so the mapping consequence survives as ours - but the
  existence of the repeat was never ours to claim.
- **Item 2, read counts are not abundance.** Breadth-based presence calling is
  standard in metagenomics methods work, and Blanco-Miguez et al. 2023 (*mSystems*
  8:e00079-23) state of this very reference that "there are genes not belonging to the
  island distributed in the inner part of the cluster".
- **Item 4, depth-limited breadth.** That shallow sequencing cannot support gene
  presence/absence calls, and that sequencing effort confounded with treatment group
  produces false positives, are both established, with rarefaction the standard remedy.

Only item 3, the ~1.9x density-ratio offset, turned up nothing - and even there the
general class of bias (mappability, average genome size) is well known.

**Root cause.** We searched the *pks*-and-colorectal-cancer literature, found the
pitfalls absent from it, and generalised that to "absent from the literature". The
antecedents are in the methods and bacterial-genomics literature, which we had not
searched. A claim of novelty is a claim about the whole literature and needs a search
shaped like one.

**What replaced it.** PITFALLS.md now opens with a table of what was already known,
quotes both papers, and states the narrower claim: these are places where known
methodological facts are not applied to *pks* detection. The check that produced this
correction was itself targeted rather than systematic - six searches and two full-text
reads - so a further item may yet prove to be published.

**A second error found by the same check.** Measuring the repeat properly to compare
with Auvray et al. showed our own figure was short. We had **1,248 bp**
(16,066-17,313); the repeat extends to **1,287 bp** (16,066-17,352) by exact 30-mer
matching, at 99.845% identity, and to 1,292 bp at >=98% windowed identity. An
exhaustive 30-mer scan confirms no further indel-shifted continuation: of 1,206
sequences occurring more than once in the island, 1,198 lie at a spacing of exactly
5,653 bp and no other spacing between 5,000 and 6,500 bp exists. The derived
"everything lost" breadth therefore becomes `(49455 - 2574) / 49455 = 0.9480`, not
0.9495 as stated in item 1 of this file. `data/missing_split.tsv` was regenerated with
the corrected interval. **The measured loss is unchanged** at 1,920-2,140 bp with 11
of 13 deep samples having nothing missing outside the repeat, because the bases that
actually go uncovered lie in the interior of the repeat, bounded by read length - so
widening the interval by 39 bp at each end moves nothing.

---

## 10. "Omitting a MAPQ filter inflates RPKM for clbJ and clbK"

**What we said**, in PITFALLS.md and in the commit message that introduced the worked
example: that a pipeline without a mapping-quality filter retains the duplicated
reads and so inflates the per-gene RPKM of *clbJ* and *clbK*.

**Why it was wrong.** An aligner in its default mode reports one primary alignment per
read. A read from the repeat is counted **once**, placed arbitrarily at one of the two
copies. The total over *clbJ* + *clbK* is therefore about preserved, and nothing is
inflated. What is unreliable is only the split between the two genes. (This paragraph
said "BWA-MEM reports one primary alignment per read" until 2026-11-28. We have never
run BWA-MEM - this pipeline is bowtie2 throughout - and we had not measured either
aligner. See the amendment below, where it is measured.)

**Root cause.** We reasoned from our own filter's behaviour - we drop those reads, so
they must be being double-counted by someone who does not - without checking what an
unfiltered aligner actually emits.

**What this means, and it cuts against us.** A method that averages read counts over
the *clb* genes is **barely affected by the duplication at all**. A breadth method
loses ~19% of *clbK* and ~16% of *clbJ* permanently. So pitfall 1 is not a defect we
found in other people's pipelines; it is a property of the reference that anyone
calling presence from coverage - **us included** - has to know, and we met it because
of a design choice we made. PITFALLS.md now says this, with the per-gene breadths at
>= 10x island depth as the evidence: 0.8065 for *clbK*, 0.8450 for *clbJ*, 0.9982 mean
for the other 16 genes.

**A second, related overreach.** The preamble to PITFALLS.md framed all four items as
"known methodological facts not being applied to *pks* detection". Reading Nooij et al.
2021 in full showed that framing is too broad: that paper subsamples every metagenome
to a common 10 million read pairs, states outright that "metagenomics-based assays
require high sequencing depths to reliably detect single genes", and **excludes marker
genes after detecting homology-driven uneven coverage in them** - which is the
reasoning behind our items 1 and 3, applied before we did it. Items 3 and 4 are, in
substance, anticipated there. Only item 1 is absent from it. The preamble and the
item 4 row were rewritten accordingly.

---

### Amended 2026-11-28: we corrected a reasoning error with more reasoning, and have now measured it

The correction above is right, and the way it was arrived at was the same mistake it
was correcting. We had claimed other people's read counts were inflated by reasoning
from our own filter's behaviour. We then withdrew that by reasoning about how an
aligner behaves - naming **BWA-MEM**, which this project has never run, and measuring
nothing. Items 15 and 16 of this file are the same shape: a real statement with no
measurement under it.

So we measured it (`pipeline/mapq_split.sh`, `data/mapq_split.tsv`). One carrier,
`cancer_ERR1018186` - the only one whose reads are still on disk, island breadth 0.9578
- mapped to the island reference with **bowtie2**, the aligner this pipeline actually
uses, called the same way `asm_measure.sh` calls it.

**Is a read ever counted twice? No.**

| | |
|---|---|
| primary alignments | **16,442** |
| distinct read names | **16,442** |
| secondary alignments | 0 |
| supplementary alignments | 0 |

Equal to the last digit, with nothing secondary. That is the whole of the "nothing is
inflated" claim, and it holds.

**What does the filter do to the two genes?** The opposite of inflating them:

| gene | bp | unfiltered | MAPQ 20 | lost |
|---|---|---|---|---|
| *clbK* | 6,465 | 1,904 | 1,560 | **18.1%** |
| *clbJ* | 6,501 | 2,042 | 1,665 | **18.5%** |
| the other 16 | - | - | - | 0.2-2.1%, median 1.2% (except *clbR*, 213 bp, 11.5% on 87 reads) |

So the MAPQ filter - ours - takes about 18% off each of those two genes, and leaves
every other gene alone. **The deficit belongs to the filtered pipeline, not to the
unfiltered one.**

**And the split really is arbitrary.** The two copies of the repeat take 384 and 405
unfiltered reads, differing by 5.5%, of which only 64 and 51 survive MAPQ 20 -
115 of 789, or 14.6%. Two sequences 99.840% identical, each taking about half the reads
that could belong to either: that is what arbitrary apportionment looks like, and it is
why the per-gene split cannot be trusted while the sum can.

**Where the correction overstated, by about 5%.** "The total over *clbJ* + *clbK* is
therefore preserved" is exactly true of the counting mechanism and not quite true of
the number. Unfiltered, those two genes carry 3,946 reads over 12,966 bp = **0.3043
reads/bp**, against **0.3218** for the three genes over 4 kb that have no repeat in
them (*clbN*, *clbH*, *clbB*) - so **94.6%**, a 5.4% residual deficit with no filter at
all. At MAPQ 20 it is 77.3%. "Preserved" holds to about 5%, not exactly, and the word
was doing more work than any measurement supported. The conclusion is unaffected: 5.4%
is nothing beside the 19% a breadth method loses permanently, so a read-count method is
still barely affected by the duplication.

**What is still not measured.** BWA-MEM itself, which is what the withdrawn sentence
named and what the published pipelines we were reasoning about actually use. It is not
installed here and we did not install it to find out. One primary alignment per read is
its documented default and the mechanism above does not depend on which aligner has it,
but the sentence should not have named a program we never ran, and no longer does.

**Shipped with this**: `data/mapq_split.tsv`, measured by `pipeline/mapq_split.sh`,
tabulated by `code/make_mapq_table.py`, checked in `code/verify_numbers.py` section 15.

---

## 11. We predicted the duplication would be visible in someone else's correlation matrix

**What we said**, in the paragraph added with item 10: that omitting a mapping-quality
filter "does make their per-gene RPKM values for those two genes arbitrary, which
matters for the gene-to-gene correlation matrix they use to argue the island travels
intact". A testable prediction: *clbJ* and *clbK* should correlate worse than other
*clb* pairs in Nooij et al.'s Supplementary Figure 1.

**Why it was wrong.** We extracted that figure from the PDF and read it. ***clbJ* vs
*clbK* = 0.92**, which is in the middle of the range for *clb* pairs in the same matrix
(roughly 0.79 to 0.99 excluding *clbS*). There is no detectable degradation.

**Root cause.** We reasoned about the ambiguous region without weighing it against the
gene. The repeat is 1,287 bp inside genes of 6,465 and 6,501 bp - about a fifth. The
other four fifths are unique and dominate the per-gene RPKM, and arbitrary
apportionment of the ambiguous reads averages out across a library. The effect we
predicted is real in principle and negligible in size.

**What this means.** Taken with item 10, the duplication leaves **no trace at all** in
a read-count method: not in the summed abundance, and not in the gene-to-gene
correlations. It is a problem exclusively for coverage-based calling. Our pitfall 1
is a warning to people doing what we did, and we have now been wrong about its
consequences for other people's methods twice, in the same direction - assuming our
own method's failure mode must be everyone's.

---

## 12. The carrier analysis was reported against the wrong denominator

**What we said**, in the summary table of RESULTS.md and then on the front page of the
README: that carriers were examined "of 14 CN-CRC positives and 1 DE-CRC positive",
and that 3 samples of **15** were better explained by *K. pneumoniae* than by
*E. coli*.

**Why it was wrong.** DE-CRC has **9** positives, not 1, and all 9 were examined - the
DE-CRC section of the same file says so in its own heading ("DE-CRC, all 9 positives")
and reports that 8 of them fit *E. coli*. Only the one exception was carried into the
summary. The correct denominator is **23**.

**Effect.** Three *Klebsiella* attributions out of 15 is 20%; out of 23 it is **13%**.
We published the larger figure. No individual sample's attribution changes.

**Root cause.** The summary table was written while finishing the CN-CRC work and
counted the DE-CRC cohort by the one sample that was interesting rather than by the
samples that were examined. The error then propagated to the README when the carrier
analysis was promoted to the front page.

**A second thing this exposed.** The two cohorts were not examined by the same method:
the 14 CN-CRC samples went through the six-genome competitive index, while the 9
DE-CRC samples were assessed only by marker-gene density ratio, which can exclude a
species but cannot establish one and **cannot return the "undetermined" verdict at
all**. Pooling them into one denominator without saying so implied a uniformity that
does not exist. Both files now state the two methods separately, and the
"4 unassignable" figure is reported against 14 rather than 23, because the other 9
could not have produced it.

---

## 13. An audit of every checkable number found six more stale figures

After item 12, we recomputed **every** quantity in these documents from `data/`
rather than fixing the one that had been noticed. Six more were wrong, all of them
stale values from earlier runs that were never recomputed when the tables were
regenerated:

| Claim | Published | Correct |
|---|---|---|
| Positive samples | 44 | **42** |
| Positives below 10x island depth | 34 of 44 | **29 of 42** |
| Positives whose uncovered bases lie only in the duplication (item 2 above) | 9 | **11** |
| Positives below the derived ceiling despite >= 14 genes (item 2 above) | 25 | **26** |
| Mean bases missing outside the duplication at 2-10x | 2,597 | **2,585** |
| Island-negative samples with >= 50 reads on the reference | 82 | **106** |
| Their median breadth | 0.0038 | **0.0018** |

None of these changes a conclusion. The duplication still caps breadth, the background
is if anything **larger** than we reported, and depth still limits most positive calls.
But a repository that asks readers to check its numbers has to have numbers that
survive being checked.

**Root cause.** The data tables were regenerated three times - when the
*K. pneumoniae* NA bug was fixed (item 7), when the duplication interval was corrected
(item 9), and when `island_reads` was added here - and the prose was not recomputed
against them. Each fix was verified not to change the primary result, which was true,
and that was mistaken for verifying that it changed nothing.

**A worse finding than any of the individual numbers.** Pitfall 2's figures were
**not checkable from the repository at all.** They depend on how many reads mapped to
the 55,140 bp reference, and that column was not in any shipped table - breadth and
depth are computed over the gene intervals only, which is precisely the distinction
the pitfall is about. A reader could not have caught this, and neither could we until
we tried. `data/per_sample.tsv` now carries an `island_reads` column, and
`code/make_table.py` documents what it is.

**What we did not do at the time.** Three figures were left unverifiable because
their inputs were not here: the ~1.9x ratio decomposition, the DE-CRC between-group
depth difference (P = 1.5e-14), and the earlier regression whose odds ratio moved
after adjusting for depth. All three turned out to be computable once `island_reads`
was shipped, and `code/verify_extra.py` now does so. See item 14.

---

## 14. The three "unverifiable" figures were verifiable, and one of them was wrong

Item 13 ended by saying three figures could not be checked from the shipped data and
would be flagged in place. That was the wrong call: with `island_reads` shipped, all
three are computable, and `code/verify_extra.py` computes them. Leaving a number in a
document with a note saying nobody can check it is not much better than leaving it
wrong.

**Two reproduced exactly.**

- The DE-CRC depth difference: medians 69.5M and 34.9M reads, Mann-Whitney
  **P = 1.51e-14** against the published 1.5e-14.
- The ratio decomposition: 55140/49455 = **1.1150** for the length convention, and
  **1.1434** for marker fragment edge loss from a mean fragment of 789.3 bp at a
  100 bp read. The third factor, ~1.50, is **not an independent measurement** - it is
  the residual after dividing the observed ratio by the other two, and the script now
  says so where the document implied three measured components.

**One was wrong.** PITFALLS.md said that adjusting for sequencing depth moved the
AT-CRC odds ratio "from 3.11 to 1.31 and P from 0.080 to 0.71". Recomputed, the
unadjusted model gives **OR 3.11, P = 0.0798** - exactly as published - but adding
log10(total reads) gives **OR 1.78, P = 0.4033**, not 1.31 and 0.71. Depth itself is
significant at P = 0.0103. The 1.31 figure came from a third model with a further
covariate whose specification was never recorded, so it is dropped rather than
reconstructed. The point the paragraph makes is unaffected and arguably clearer: one
covariate halves the odds ratio and takes P from 0.08 to 0.40.

**The R values shifted slightly**, because `island_reads` is the pre-MAPQ count from
the bowtie2 log while the published figures used MAPQ-filtered reads. DE-CRC's eight
*E. coli*-consistent positives are **R 1.43-2.12, median 1.91** rather than 1.37-2.01,
and the *Klebsiella* outlier is **764** rather than 723. Both versions are stated
where they appear, and the conclusions - a tight cluster near 1.9 and one sample two
and a half orders of magnitude away - do not move.

---

## 15. The front-page finding was not reproducible from this repository either

Asked whether anything else was unverifiable, we swept the documents by what data each
claim depends on rather than by reading them. **The carrier analysis - the finding we
had just promoted to the front page - was not in the shipped data at all.** Neither
was the assembly confirmation that is the strongest single observation in the project.
Both are now shipped; both are tiny.

| Claim | Depended on | Now in |
|---|---|---|
| The 14-sample carrier table | `*.g6.depth.tsv`, 14 KB, not shipped | `data/carrier_six_genomes.tsv` (6 KB, 98 rows) |
| "The island and the *Klebsiella* chromosome on a single contig" | a 386 kb assembly contig, not shipped | `data/assembly_alignment.tsv` (4 rows) |
| That the contig is not a chimera | read depth across the junction, not shipped | `data/assembly_junction_depth.tsv` (5,001 rows) |
| The 1,287 bp repeat at 99.845% | AM229678.1, which we do not redistribute | `code/verify_repeat.py`, which fetches it |

The contig is 386,079 bp. Positions 25-354,276 align to the *K. pneumoniae* HS11286
chromosome; positions 354,259-386,079 match AM229678.1 exactly over 31,389 50-mers.
Read depth across the junction has a minimum of 12 and **no base at zero**, which is
what distinguishes a real join from an assembler's mistake. None of that was checkable
by a reader until now.

**Root cause, and it is the same one as item 13.** We shipped what the *primary*
analysis needed - the depth histograms - and assumed that covered the repository's
obligations. Every subsequent finding was written up from files that stayed on one
laptop. The check that would have caught it is not "are the numbers right" but "for
each claim, which shipped file would a reader open", and we ran it only when asked.

**What is still not verifiable from this repository**, listed so it is not discovered
a fourth time:

- ~~**The *Klebsiella* strain typing** of `cancer_ERR1018186` (ST23/K1/ICEKp10) comes
  from Kleborate output that is not shipped.~~ **Now shipped** under
  `data/kleborate/`, which turned out to matter more than we expected - see item 16.
- **The marker-index constants** (96,300 bp over 122 fragments for *E. coli*;
  `PER_PCT = 313`) are stated in `pipeline/entero_ref.py` but reproducing them means
  rebuilding the index.
- **The correlations read off Nooij et al.'s Supplementary Figure 1** were read from a
  published image, not from their data, and PITFALLS.md says so where they appear.
- **Everything quoted from a paper** is verifiable only by reading that paper; the
  DOIs and PMC identifiers are given.

---

## 16. The strain typing was one line made of two runs

Asked to ship the Kleborate output as well, we did, and the output does not say what
the one-line typing in RESULTS.md said. **"ST23/K1/ICEKp10/*clb*" is two runs put
together**, and only the first half comes from the bin that was assembled
automatically.

| | `kp_bin` (62 contigs the assignment step gave to *K. pneumoniae*) | `kp_plus` (that bin plus the 4 island-bearing contigs) |
|---|---|---|
| species | *Klebsiella pneumoniae*, strong | same |
| MLST | **ST23** | same |
| K type | **K1** | same |
| yersiniabactin | **none detected** | **ybt 1; ICEKp10** |
| CbST | 26-6LV | 29 |
| colibactin | **clb 2 (incomplete)**, *clbA*-*clbI* only | **clb 2**, complete |

So the species, the ST and the K type are robust: they come from 5,082,743 bp of
chromosome that the six-genome assignment placed in *K. pneumoniae* on its own.
**ICEKp10 and the complete *clb* 2 do not.** They appear only after four contigs are
added by hand, and the assignment step had given all four to *C. koseri* BAA-895.

### Why those contigs were added, and why it still limits the claim

Not because the sample needed it. **The reference set makes it impossible to assign a
*clb* contig to *Klebsiella*:**

| Reference genome | Share of this island's 40-mers | *clb* genes present |
|---|---|---|
| *E. coli* IHE3034 | 100.00% | 18/18 |
| *C. koseri* BAA-895 | **94.88%** | **18/18** |
| *K. aerogenes* KCTC2190 | 0.13% | 0/18 |
| ***K. pneumoniae* HS11286** | **0.11%** | **0/18** |
| *C. freundii* CFNIH1 | 0.08% | 0/18 |
| *E. cloacae* ATCC13047 | 0.08% | 0/18 |

HS11286 carries no *clb* island at all, so any contig made of *clb* can only ever be
won by *E. coli* or *C. koseri*. Adding the island contigs back was forced by the
reference set, not chosen after seeing the answer, and the rule is checkable:
`kp_plus` is **exactly** `kp_bin` plus every contig carrying any part of the island,
which `python3 code/verify_extra.py` asserts.

But the consequence stands either way: **a typing produced from a bin we assembled by
hand is not independent evidence that those contigs belong to this *Klebsiella*
genome.** We presented ICEKp10 as though it were.

The independent evidence is the one from item 15, and it does not use the added
contigs: `k119_45171` carries 354 kb of *K. pneumoniae* chromosome and 32 kb of the
island contiguously, with read depth never dropping to zero across the join - and that
contig is in the automatic bin. Consistently with this, the automatic bin alone types
as *clbA*-*clbI*, which are the clb genes that contig physically carries: reference
positions 17,195-54,756 span *clbA* through *clbJ* and part of *clbK* (the island's
gene order runs *clbQ* first in this reference), and Kleborate types neither *clbJ*
nor *clbK*.

**What changed in the text:** RESULTS.md now says which run each part of the typing
comes from, instead of running them together.

### Amended one commit later: the similarity argument is worse, and a different one works

The paragraph above said a typing from a hand-assembled bin "is not independent
evidence that those contigs belong to this *Klebsiella* genome". That stands, and
once we split the alignments at the island boundary it is **worse** than stated.

Our own working notes from the time claimed the opposite - that with the island region
excluded, `k119_70945`'s backbone "matches only *K. pneumoniae*, *C. koseri* zero".
**That does not reproduce.** There are only two alignments for that contig, and
separating island from backbone gives:

| region of `k119_70945` | *C. koseri* BAA-895 | *K. pneumoniae* HS11286 |
|---|---|---|
| backbone (the 49.6 kb that is not island) | 46,008 / 49,568 = **0.928** | 37,315 / 45,623 = **0.818**, and only over 45.6 kb |

Those are matched bases over aligned length from minimap2's PAF. **They include gaps
and are not identities** - a contig at 99% identity can read as 0.81 this way, a
mistake we made earlier in this project and recorded in `pipeline/asm_kleborate.sh`.
They are used here only to compare the two species against each other on one contig.

By similarity the backbone of that contig favours *C. koseri*, not *Klebsiella*. It
behaves like an ICE-like element with homologues in both genera, which is exactly the
kind of sequence whose species cannot be read off similarity. We had written down a
number that said otherwise and did not recheck it before quoting the conclusion.

**There is independent evidence, of a different kind, and we had not used it.**
megahit reports a k-mer coverage depth per contig, and these are all from one sample:

| contig | bp | depth | island content |
|---|---|---|---|
| the bin, length-weighted over 62 contigs | 5,082,743 | **18.08** | - |
| `k119_45171` | 386,079 | **18.01** | 31,821 bp, *clbA*-*clbJ* end |
| `k119_70945` | 63,927 | **18.99** | 14,343 bp, the *clbQ*-*clbL* end |
| `k119_70596` | 4,643 | **19.83** | all island |
| `k119_32464` | 1,248 | **36.00** | all island - **almost exactly twice**, which is what a collapsed two-copy repeat must be |
| `k119_25170` | 4,341 | **8.87** | 137 bp - **does not fit, and is barely island at all** |

The island-bearing contigs sit at the *Klebsiella* genome's own coverage depth; the
collapsed duplication sits at double it, as a two-copy repeat has to; and *C. koseri*
is not in this sample to be their source - its genome breadth is 0.014, which is
itself mostly the *clb* island it carries. **The argument is coverage depth plus the
absence of the alternative, not sequence similarity.** It is weaker than the single
contig of item 15 and we are not presenting it as equivalent, but it is evidence, and
it is now shipped: `data/kleborate/island_contig_backbone.tsv` and
`data/kleborate/contig_depth.tsv`, recomputed by `python3 code/verify_extra.py`.

Two things fall out of it that were not visible before:

- **`k119_25170` should probably never have been in the island set.** Depth 8.87,
  about half the genome's, and only 137 bp of exact island *k*-mers.
  (The depth half of this is wrong - see the 2026-11-28 amendment below. In read depth
  that contig is at 0.874 of the bin, not half; 8.87 is a megahit `multi` value and
  `multi` misreports this one contig by a factor of 1.9. The 137 bp stands and is the
  reason.)
- **`k119_25144` is not an island contig either**, though it was in the island FASTA
  and in the automatic bin: 1,024 bp of island in four scattered blocks across 214 kb.

**Root cause of this one.** We checked a published claim against a note we had written
rather than against the data. The note's figure came from one hand-typed command with
an ad hoc coordinate cut and no record of how it was computed. The first version of
the table above had its own version of the same bug - taking the island's extent as
first-hit-to-last-hit, which turns `k119_25144`'s four scattered blocks into "39,448
bp of island" - and that is recorded in `code/make_backbone_table.py`.

### Amended 2026-11-28: the depth argument was replaceable, and we replaced it

The section above ended by calling its own evidence "coverage depth plus the absence
of the alternative", "weaker than the single contig of item 15". Asked to deal with
that, we found two things wrong with the argument and one thing we had never looked
at.

**The depths in that table are not read depths.** Every number in it is megahit's
`multi`. This assembly was built to k=119 and the reads are 100 bp, so `multi` is a
multiplicity in a de Bruijn graph built from the reads **plus the contigs carried
forward from the smaller k**. It is not a read depth, and we never calibrated it
against one. Four of the six contigs are also mostly backbone - the island is 8%,
0.5%, 22% and 3% of them - so their `multi` is their backbone's, not the island's.
One thing in that column does survive: a contig at twice its neighbours is what a
collapsed two-copy repeat looks like, and that is used again below.

**We had treated the four added contigs as four separate objects.** That is *why*
only depth was left: taken one at a time, each one needs its own provenance, and
there was nothing else to hand. They are not separate. Laid out on the island
reference AM229678.1 (54,791 bp) they tile it, at 99.76-100% identity:

| order | contig | in the bin | island | identity | overlap with previous |
|---|---|---|---|---|---|
| 1 | `k119_70945` | no | 1,842-16,184 | 99.958% | - |
| 2 | `k119_32464` | no | 16,066-17,313 | 99.760% | **119 bp** |
| 3 | `k119_70596` | no | 17,195-21,837 | 100.000% | **119 bp** |
| 4 | `k119_32464` | no | 21,719-22,966 | 99.760% | **119 bp** |
| 5 | `k119_45171` | **yes** | 22,848-51,731 | 99.976% | **119 bp** |
| 6 | `k119_45171` | **yes** | 51,804-54,791 | 99.833% | -72 bp, inside one contig |

**Every junction between two different contigs is an overlap of exactly 119 bp** -
the k-mer length this assembly was built at, the `k119_` in every contig name. That
is not a coordinate artefact: taking the contig sequences themselves, the end of each
contig and the start of the next share **119 bases of identical sequence**, in the
right orientation, at all four junctions. Both numbers are in
`data/kleborate/island_tiling.tsv` and both are asserted. Contigs sharing exactly k
bases are consecutive nodes of one path in the assembly graph. megahit did not
produce four fragments of uncertain origin. It produced one island, 52,878 of the
reference's 54,791 bp, and broke it twice. One of the four, `k119_45171`, is the
contig of item 15: 350 kb of *Klebsiella* chromosome running into the island with
read depth never dropping across the join. **The other three reach the chromosome through it, by
overlap.** That is the step we had missed, and it needs no depth at all.

**Where it broke is the duplication of PITFALLS item 1.** `k119_32464` is in the
tiling twice, at island 16,066 and island 21,719. PITFALLS item 1 measured that
duplication from the reference alone - 1,287 bp starting at **16,066**, the two copies
**5,653 bp** apart - and `code/verify_repeat.py` checks it against ENA. Same start,
same spacing, and it is the one contig at twice its neighbours' multiplicity. The
k-mer scan reads the reference; this reads an assembly of one sample's reads. They
share no input and they land on the same repeat.

**The element has two flanks, and both were in the automatic bin already.** The island
is part of a larger insertion. Measured on the only one of the six references that
carries the whole element - *C. koseri* BAA-895, used **as a ruler, not as a claim
about origin** - the element runs 816,165-953,855, 137,690 bp, tiled by 13 contigs.
Its two ends are reached by `k119_45171` and `k119_25144`, and **neither was added by
hand**. On the *Klebsiella* chromosome, `k119_25144` runs 3,250,687-3,506,454 and
`k119_45171` picks up at 3,506,446: the two intervals meet within 8 bp. One insertion,
one site, both flanks in the automatic bin, the added contigs filling the interior.
That is position - not depth, not similarity.

`k119_25144`'s own end shows the crossing: its last 4 kb runs through 1,358 bp that
aligns to *both* references (an attachment-site repeat - *Klebsiella* 3,433,470-3,434,828
and the ruler at 952,497-953,855), then 1,885 bp of *Klebsiella* chromosome ending at
3,506,454, then 1,037 bp of the element to the end of the contig.

**What this argument assumes**, and section 8's does not: that the two contigs'
positions can be read off HS11286, i.e. that the local chromosome order is the same in
HS11286 (ST11) and in this strain (ST23). A rearrangement between them would break the
"both flanks" reading. The tiling in section 8 uses no reference order at all, only
the contigs' own overlaps, so the two arguments fail in different ways.

It also explains why the similarity argument failed on `k119_70945` specifically. That
contig's backbone aligns to the *Klebsiella* chromosome too, at 3,434,608-3,479,783,
inside a stretch of HS11286 that `k119_25144` skips: HS11286 carries a related element
at a different site. A sequence present in both genera is exactly the sequence whose
species cannot be read off similarity, which is what the section above says, now with
the reason visible.

Seven of the 13 tiling contigs were never added to the bin. The element is larger than
the island and the bin is not the whole genome; nothing in the typing depends on those
seven.

**We also measured the depths properly, and one of them was wrong.** The reads were
mapped back to the whole assembly - all 21,759 contigs in one competitive index, MAPQ
20 (`pipeline/asm_linkage.sh`). On that scale `multi` turns out to track read depth
closely for five of the six contigs (multi/depth 0.736-0.786) and **not at all for
`k119_25170` (0.404)**:

| contig | bp | `multi` | read depth | x the bin | multi/depth |
|---|---|---|---|---|---|
| the bin, 62 contigs | 5,082,743 | 18.08 | **25.13** | 1.000 | - |
| `k119_45171` | 386,079 | 18.01 | 24.25 | 0.965 | 0.743 |
| `k119_25144` | 214,126 | 16.40 | 22.29 | 0.887 | 0.736 |
| `k119_70945` | 63,927 | 18.99 | 25.08 | 0.998 | 0.757 |
| `k119_70596` | 4,643 | 19.83 | 26.58 | 1.058 | 0.746 |
| `k119_25170` | 4,341 | **8.87** | **21.97** | **0.874** | **0.404** |
| `k119_32464` | 1,248 | 36.00 | 45.83 | 1.824 | 0.786 |

The bin's own 46 contigs over 10 kb span 0.76-1.33 times the bin's depth, so all four
added contigs are inside the bin's own spread, and the collapsed repeat is at 1.82 -
near the 2 it has to be. **This corrects one of the two things the section above said
it had found:** "`k119_25170` should probably never have been in the island set. Depth
8.87, about half the genome's" - the depth half of that does not hold. In read depth it
is at 0.874 of the bin, like everything else. The reason to doubt that contig is the
other half, its 137 bp of island, which stands.

**Read-pair linkage.** The same mapping gives pairs whose two mates land on different
contigs, both at MAPQ 20. Partners are counted only among the 160 contigs the
six-genome step assigned at breadth >= 0.5; the other 21,599 are other gut bacteria and
pairs reaching them are shared repeats. The rule was fixed before counting.

| contig | bp | to the bin | to added island | to other species | share bin+island | to bin per kb |
|---|---|---|---|---|---|---|
| `k119_45171` (bin) | 386,079 | 377 | 37 | 27 | 0.939 | 0.98 |
| `k119_25144` (bin) | 214,126 | 216 | 9 | 15 | 0.938 | 1.01 |
| `k119_70945` | 63,927 | **69** | 39 | 7 | **0.939** | 1.08 |
| `k119_70596` | 4,643 | 0 | 83 | 0 | 1.000 | 0.00 |
| `k119_25170` | 4,341 | 11 | 5 | 15 | 0.516 | 2.53 |
| `k119_32464` | 1,248 | **34** | 117 | 0 | 1.000 | 27.24 |

The bin's own contigs over 10 kb span 0.857-1.000 on that share. **`k119_70945` sits at
0.939, inside it**, with 69 pairs landing on bin contigs - and it is the only added
contig long enough for the test to discriminate. That is the question the tiling left
open: a second island copy in another genome would not put its mates on *Klebsiella*
chromosome.

**For the other three the test has no power, and we are saying so rather than reading
the numbers.** Below 10 kb the bin's own contigs span 0.129-1.000 on the same measure,
which covers anything an added contig could score; `k119_70596` (1.000) and
`k119_25170` (0.516) are both inside that. One concrete fact does survive: `k119_32464`
sends 34 pairs to `k119_45171`, and `k119_45171`'s largest single assigned partner is
`k119_32464` - a hand-added contig ahead of every contig in the bin.

**What is still not settled.** Both breaks in the tiling are at the repeat, and a
repeat is the one place where graph adjacency can join paths from different sources.
So it remains possible in principle that `k119_70945` sits on a second copy of the
island in this sample. That is the one question the absence of the alternative is
still carrying, and it is now one narrow question rather than the whole claim: in
this sample the island's read depth is **33 times** the *E. coli* genome's, so an
*E. coli* source would need 33 copies of a 50 kb island per *E. coli* cell. That
ratio is the one in [RESULTS.md](RESULTS.md)'s carrier table and is recomputed by
`verify_numbers.py`.

**What does not change.** The correction itself stands: a typing run on a bin we
assembled by hand is not independent evidence of what the added contigs belong to, and
we had presented ICEKp10 as though it were. What has changed is that the question the
typing begged is now answered by something else. Nothing in the primary result touches
this - it is one sample's strain typing, and the CN-CRC *Klebsiella* count is
unchanged at 2.

**Shipped with this**: `data/kleborate/island_tiling.tsv`,
`data/kleborate/island_locus.tsv`, `data/kleborate/contig_read_depth.tsv` and
`data/kleborate/contig_links.tsv`, generated by `code/make_locus_table.py` and
`code/make_link_table.py` from the measurement in `pipeline/asm_linkage.sh`, and
recomputed with their assertions by `python3 code/verify_extra.py` (sections 8, 9
and 10).

### One thing found while checking this

Kleborate's colibactin scheme types 15 loci and has **no *clbJ* or *clbK* column**.
Those are the two genes holding the duplication of PITFALLS item 1. We do not know the
scheme's reason and do not claim one, but a widely used tool already leaving those two
out is one more reason that item 1 is "known but not applied", not new. Added to
PITFALLS.md item 1.

### What this says about the pattern

Items 13, 15, 16 and 17 are the same mistake four times, and the two amendments above
are a fifth and a sixth. The second amendment is a variant worth naming on its own:
there the claim was right and the evidence we offered for it was real but weak, and we
had stopped at the weak evidence because we had framed the question in a way
(four separate contigs) that made the strong evidence invisible. Each time, **a claim was correct in substance and wrong in what it attributed
to the evidence we had.** Shipping the evidence is what exposed it: items 15 and 16
were both found by someone asking for the underlying file, and the amendment was found
by generating the file that claim depended on rather than re-reading the note that
asserted it.

**2026-11-28 brought a seventh and an eighth, from going looking rather than being
asked.** An audit of the remaining items for this specific failure mode - a number that
is correct and does not support the claim built on it - found two more:

- **Item 8** read a depth ratio as a population fraction ("0.85, i.e. ~85% of *E. coli*
  cells"), which is **item 3 of this same file**, five sections earlier, on a different
  numerator. Item 3 had written down the remedy - use the ratio one-sidedly only - and
  item 8 did not apply it to the sibling quantity. Measuring the scale then showed the
  attribution in that item rests on parsimony rather than on exclusion.
- **Item 10** withdrew a claim about other people's read counts by reasoning about an
  aligner it had not run, and named BWA-MEM, which this project has never used. The
  measurement confirms the conclusion and shows the word "preserved" was about 5%
  too strong.

So the count is eight, and **the fix for the first six did not prevent the seventh and
eighth.** `verify_numbers.py` grew from 84 checks to 209 over those corrections and
could not have caught either, because in both cases the data and the prose agreed
perfectly - what disagreed was the prose and what the measurement was capable of
showing. That is not a thing a number can check.

What came out of 2026-11-28 instead is [CLAIMS.md](CLAIMS.md), which puts one column
next to every load-bearing claim - **what alternative does this measurement exclude** -
and grades the claim on it. Three rows came out **C**: a real measurement that excludes
nothing. One of them is the primary result. `code/verify_claims_ledger.py` checks that
the column is filled and that the citations resolve; it cannot check whether an answer
in it is honest, and nothing can. The point of the file is to make the empty cells
visible, since all eight errors were empty cells that no layout required anyone to
look at.

---

## 17. We audited every number in the prose against the shipped data

The last six items were each found by someone asking for one file. This time we asked
the question in bulk: **extract every number from the five documents and recompute it
from `data/`.** 603 lines carry a number; after setting aside years, accessions, DOIs
and figures quoted from papers, **84 are claims we can recompute here**. Before this
audit **nine of them disagreed with the documents**, from four causes. The rest
reproduced exactly. One of the nine is a measurement error rather than a stale
figure.

### The *K. pneumoniae* abundance column was wrong by 14.93x

`PER_PCT = 313` means "313 marker reads per 1% abundance per million reads". It was
calibrated **once, on the *E. coli* marker set, which is 96,300 bp**. Reads landing on
a marker set scale with abundance **times marker length**, so the constant is specific
to that length. `code/make_table.py` applied the same 313 to *K. pneumoniae*, whose
marker set is **6,450 bp** - understating it by 96,300 / 6,450 = **14.93x**.

The project's own working script `meta/cn_entero.py` had the scaling from the start
(`factor[k] = b / ECOLI_BASES * PER_PCT`), and even printed the caveat that it is an
unvalidated assumption for *Klebsiella*. Only the published column lacked it, and
nothing compared the two.

**That assumption has since been measured and holds** (2026-11-28,
`code/verify_marker_scaling.py`): the six-genome index measures the same species
without reference to marker length, and the multiplier relating marker density to
genome density is the same for both marker sets - *K. pneumoniae* / *E. coli* = 0.929,
95% CI 0.662-1.304, against 14.93 if the scaling were wrong. The interval contains 1
and excludes 14.93 by a factor of 11. It spans a factor of 2.0, so the column is
trustworthy to about that coarseness and no further. `data/marker_counts.tsv` ships
the raw marker read counts this needs; they had never been shipped, which is why the
assumption could not be checked by a reader before now.

| | published until 2026-11-26 | corrected |
|---|---|---|
| median of detected samples | 0.0010% | **0.0143%** |
| maximum | 0.4576% | **6.8318%** |
| samples >= 0.1% (of 221 measured) | 4 | **24** |
| detected (of 221 measured) | 96 (43%) | 96 (43%), unchanged |

**This is what reconciles item 8 with the data.** Item 8 quotes "Kp relative abundance
6.83% against *E. coli* 0.178%, a factor of 38" as the reasoning we had to withdraw.
Until today a reader checking 6.83% against `data/per_sample.tsv` found 0.4576% and
could not reconcile the two: the narrative used the scaled figure and the shipped
column did not. The corrected column now gives 6.8318% for that sample, and the
corrected CN-CRC count at >= 0.1% is 19, matching the working log exactly.

**Effect on the result: none.** The abundance gate uses *E. coli*, whose scale factor
is 96,300 / 96,300 = 1, so no gate decision moves. The carrier analysis uses genome
depth, not marker abundance. Group comparisons use Mann-Whitney, which a positive
scale factor cannot change (DE-CRC *K. pneumoniae* P = 0.502 either way). Regenerating
the tables left `data/depth_hist.tsv` **byte-identical** and changed only the 96
non-zero cells of one column; `code/analyse.py` reports the same 167 judgeable,
Mantel-Haenszel OR 1.91, 95% CI 0.76-4.75, P = 0.1663, I2 = 17%.

### The AT-CRC sensitivity table was stale on every number

RESULTS.md said: "All five rules gave the identical table, 4/20 vs 2/20, P = 0.661."

| | published | recomputed |
|---|---|---|
| judgeable | 20 and 20 | **21 and 21** |
| pre-registered rule | 4/20 vs 2/20 | **5/21 vs 3/21** |
| P | 0.661 | **0.697** |
| do all five rules agree | yes | **no** - *clbB* >= 0.9 gives 4/21 vs 3/21 |

The 20/20 denominators are from a 40-sample pilot stage. The CN-CRC table beside it
was correct, which is why the error survived: the section looked checked. **No
conclusion changes** - no rule reaches P < 0.05 in either cohort - but the claim that
the rules are identical was false, and it was the claim doing the work.

### Three smaller ones

- **The read-rate threshold, 41.7 per million.** Taken from the mapping log's `pks=`
  count, which is not the `island_reads` column we ship (274 against 428 for one
  sample). From the shipped column the threshold is **56.7**. Counts and P are
  unaffected because the threshold is *defined* as the point giving the same number of
  positives as the pre-registered rule - which the table also did not say, and now
  does.
- **"Median reads mapped to the reference: 260."** The two middle values of the 106
  samples are 256 and 260; the median is **258**.
- **"Tataka et al. 2026" had no reference entry.** It is the one citation in these
  documents given by name alone, with no journal and no DOI, while its AUC figures
  carry a sensitivity analysis. Added: *BMC Microbiol* 2026;26:930,
  doi:10.1186/s12866-026-05460-3.

### One class of number cannot be checked at all, and we now say so

PITFALLS.md item 11 reports correlations from Nooij et al.'s Supplementary Figure 1
(*clbJ* vs *clbK* = 0.92; *clb* pairs ~0.79 to 0.99; mobile elements against the
*E. coli* markers 0.47 to 0.76; the 57 *clb*-vs-marker pairs at median +0.02 with 24
negative). Those were read off a figure extracted from the PDF. **Neither the
extracted image nor the 57 transcribed values was kept**, so nothing in this project
can recheck them - including the median and the count of 24, which are our arithmetic
over values we transcribed once. They are now labelled as such.

### Root cause, and the fix

Every published number was computed by a script at some point. **Nothing re-derived
the numbers in the prose from the shipped data**, so a figure went stale whenever a
table was regenerated, and a figure computed from an input we do not ship could not be
caught at all. Items 13, 15, 16 and the five above are all that single gap.

`code/verify_numbers.py` now recomputes **147 of the numbers in these documents
directly from `data/` and checks each one** (84 when it was written; the amendment to
item 16 added the rest), with the claimed value written beside the
computation and a non-zero exit if any disagrees. It runs with no network and no
working tree. A stale figure now fails a check instead of waiting to be noticed.

That does not cover numbers quoted from papers; those are checked against the reading
notes, which are not public. Where a quoted number cannot be rechecked from anything
we hold, it now says so in place.

---

## 18. The 18 genes do not total 49,455 distinct bases

**What we said**, in METHODS, PITFALLS, both READMEs and four scripts: that the 18
*clb* genes total **49,455 bp**, and that this is the denominator of breadth and of the
density ratio R.

**Why it is wrong.** 49,455 is the sum of the 18 interval lengths in
`data/clb_genes.bed`, and four adjacent pairs overlap:

| pair | overlap |
|---|---|
| *clbQ* / *clbP* | 8 bp |
| *clbN* / *clbM* | 4 bp |
| *clbL* / *clbK* | 8 bp |
| *clbG* / *clbF* | 4 bp |

So the sum double-counts **24 bp**, and the number of distinct gene bases is
**49,431**. The project has been carrying both figures for the same object the whole
time: 49,431 is what the six-genome competitive index used as the island target
(`data/carrier_six_genomes.tsv`, `target_bp`), and `code/make_carrier_table.py` says
49,431 in its own docstring on line 9.

**Where it bites.** `hist_by_gene()` in `code/make_table.py` walks the per-base depth
and assigns each position to the **first** interval that contains it, with a `break` -
so each base is counted once, which is right. But it then fills each gene out to its
own length with zeros:

```python
for nm, s, e in regions:
    ln = e - s
    miss = ln - seen[nm]
    if miss > 0:
        h[nm][0] += miss
```

For the later gene of an overlapping pair, `seen[nm]` is short by the overlap, because
those bases went to the earlier gene. The 24 bases are therefore charged to the later
gene **as uncovered**, and `call_from()` divides by `tl`, which sums to 49,455. Every
published breadth is low, and the attainable maximum is 0.99951 rather than 1.

**How much, measured rather than bounded.** Recomputing from `data/depth_hist.tsv`
against the 49,431 denominator, the largest change in any of the 397 samples is
**+0.000952**, or 0.095 of a percentage point:

| | published | corrected |
|---|---|---|
| highest breadth in 397 samples | 0.960692 | 0.961158-0.961644 |
| "everything lost" derived ceiling | 0.9480 | 0.9479 |

The corrected value is a range because the shipped histogram cannot say whether those
24 bases were covered - they are folded indistinguishably into the earlier gene's
histogram - so it is bounded to within 24/49,431 = 0.000486. Recovering the exact
figure needs the raw per-base depth files, which are not in this repository.

**Effect on every result: none, and we checked rather than asserting it.** No sample
crosses the breadth >= 0.5 call threshold (the nearest is nowhere near it). No gene
crosses the 0.5 per-gene threshold: the largest overlap is 8 bp in a 6,501 bp gene, and
the shortest gene, *clbR* at 213 bp, overlaps nothing. The derived "everything lost"
ceiling moves from **0.9480 to 0.9479** - a difference of 0.000026, which does show up
in the fourth place, so item 1 above now carries the right denominator in its own note
rather than being declared unaffected. (This sentence said "0.9480 either way to four
places" when it was first written, which the check in `code/verify_numbers.py`
section 14 rejected on the same day: 0.947927 rounds to 0.9479.) The highest breadth becomes about
0.9616, still far below the 0.99 that item 1 says no sample exceeded. The primary
analysis is unchanged: 167 judgeable, MH OR 1.91, 95% CI 0.76-4.75, P = 0.1663.

**We are not regenerating the data for this, and that is a choice worth stating.** The
exact corrected breadth is not recoverable from what is shipped, the effect is 0.095 of
a percentage point, it moves no call and no conclusion, and regenerating 397 samples to
chase it would mean re-running a pipeline whose inputs we would have to re-download
from ENA. What we have done instead is write down the true figure, keep the number that
was actually used so the published values can be reproduced exactly, and add a check so
the two figures can never silently disagree again.

**The general form of this, which is the transferable part.** Summing the lengths of
BED intervals gives the length of the features, not the number of distinct bases they
cover, and the two differ exactly when features overlap - which adjacent genes in a
compact operon routinely do. A pipeline that takes its numerator from a per-base walk
(where each base is counted once) and its denominator from the interval sum has two
different definitions of the same target, and the ratio is wrong by the overlap. This
is now [PITFALLS.md](PITFALLS.md), in the note under the reference composition.

**How it surfaced.** While adding a check for the new island-ratio table, the island's
`target_bp` in `data/carrier_six_genomes.tsv` read 49,431 where the constant in
`code/verify_extra.py` read 49,455. Two numbers for one object, 24 apart.

---

## What did not change

No pre-specified threshold was altered at any point after results were seen:
`MIN_DEPTH = 1`, `GENE_COVER = 0.5`, `CALL_BREADTH = 0.5`, `CALL_GENES = 14`,
`SUSPECT_BREADTH = 0.1`, `ECOLI_MIN = 0.1`, target depth 20M reads, MAPQ 20.
They are hard-coded at the top of `code/analyse.py` and can be checked against any
version of this history.

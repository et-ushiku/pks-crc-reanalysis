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
and 0.9495; the repeat was later measured at 1,287 bp - see item 9.) Reads that happen to span a
position where the two copies differ (they are 99.840% identical, so 0.16% of
positions) can still be placed, which lifts breadth above that figure by a
sample-dependent amount. **The ceiling is not a fixed number.** The defensible
statements are: 0.9495 is what you get if the duplication is lost completely, and no
sample in 397 exceeded 0.99.

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
2.036 / 0.7638) while their marker columns read 0.000000 and 0.006483. A genome signal
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

**Why it was wrong.** BWA-MEM reports one primary alignment per read. A read from the
repeat is counted **once**, placed arbitrarily at one of the two copies. The total over
*clbJ* + *clbK* is therefore preserved, and nothing is inflated. What is unreliable is
only the split between the two genes.

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

### One thing found while checking this

Kleborate's colibactin scheme types 15 loci and has **no *clbJ* or *clbK* column**.
Those are the two genes holding the duplication of PITFALLS item 1. We do not know the
scheme's reason and do not claim one, but a widely used tool already leaving those two
out is one more reason that item 1 is "known but not applied", not new. Added to
PITFALLS.md item 1.

### What this says about the pattern

Items 13, 15 and 16 are the same mistake three times. Each time, **a claim was correct
in substance and wrong in what it attributed to the evidence we had.** Shipping the
evidence is what exposed it. Items 15 and 16 were both found by someone asking for the
underlying file, not by us re-reading the text.

---

## What did not change

No pre-specified threshold was altered at any point after results were seen:
`MIN_DEPTH = 1`, `GENE_COVER = 0.5`, `CALL_BREADTH = 0.5`, `CALL_GENES = 14`,
`SUSPECT_BREADTH = 0.1`, `ECOLI_MIN = 0.1`, target depth 20M reads, MAPQ 20.
They are hard-coded at the top of `code/analyse.py` and can be checked against any
version of this history.

# What each claim rests on, and what it rules out

`code/verify_numbers.py` checks that every number in these documents can be recomputed
from `data/`. That catches one kind of error - the kind in
[CORRECTIONS.md](CORRECTIONS.md) item 17, where the prose and the data had drifted
apart. It cannot catch the kind in items 13, 15, 16 and the 2026-11-28 amendment to
item 8, which is:

> **the number is correct, and it does not support the claim made from it.**

That error has now been made **eight** times in this project. Six were found by somebody
asking; the last two were found on 2026-11-28 by auditing the other items for this
specific failure mode. None of the eight was found by a script, and none of them could
have been: in every case the prose and the data agreed, and what disagreed was the
prose and what the measurement was capable of showing. This file is the attempt to make it visible
instead. For each claim that something in this repository rests on, it records what
measurement supports it, **what alternative explanation that measurement excludes**,
and how it could be broken.

The middle column is the one that matters. A measurement that excludes nothing is not
evidence for a claim; it is only consistent with it. Every one of the four errors above
was a real measurement in the "evidence" column with an empty "excludes" column, and
nobody noticed because nothing in the layout of the documents made that column exist.

## Grades

| | meaning |
|---|---|
| **A** | A measurement that excludes the stated alternative, without depending on a reference's gene order or on a scale we have not calibrated. |
| **B** | A measurement that excludes the alternative **under an assumption that is stated**. |
| **C** | A real measurement that excludes no alternative. The claim rests on parsimony, or on nobody having offered a competitor. Consistent-with, not evidence-for. |
| **D** | Not recheckable from what this repository ships. |

A **C** is not necessarily wrong. It is a claim whose support is weaker than its
phrasing usually suggests, and the phrasing is what we keep getting wrong.

## The ledger

| # | Claim | Evidence | What it excludes | Grade |
|---|---|---|---|---|
| 1 | *pks* carriage is not associated with CRC in these 397 samples (167 judgeable, MH OR 1.91, 95% CI 0.76-4.75, P = 0.1663, I2 = 17%) | `data/per_sample.tsv`, `code/analyse.py`, thresholds fixed before the cohorts were run, rarefied analysis primary | **Only OR > 4.75.** It does not exclude a real association; the interval contains 1 and also contains 4.7. A null result at this sample size excludes very little, and saying "no association" would be the same overreach as item 4 | **C** |
| 2 | The *pks* reference contains a 1,287 bp duplication 5,653 bp apart | Exhaustive 30-mer scan of AM229678.1 (`code/verify_repeat.py`); 1,198 of 1,206 repeated sequences at exactly that spacing and no other spacing in 5,000-6,500 bp | That it is an artefact of one measurement: an independent assembly produced a contig at double multiplicity starting at the same base with the same spacing (item 16 amendment), sharing no inputs with the scan. Auvray et al. 2021 report the same feature | **A** |
| 3 | The duplication costs ~19% of *clbK* and ~16% of *clbJ* breadth, permanently | Per-gene breadth at >= 10x island depth: 0.8065, 0.8450, against 0.9982 for the other 16 genes | That it is depth-limitation: the loss is 1,920-2,140 bp across 13 samples spanning 10x to 391x, i.e. flat in depth, while everything outside the repeat is depth-dependent | **A** |
| 4 | Breadth-based calling needs rarefaction, and unrarefied group comparisons are unsafe | At >= 10x, 11 of 13 samples have nothing missing outside the repeat; below 2x the median missing runs to thousands of bases. CN-CRC reverses direction when rarefied (item 4) | That rarefaction is a sensitivity check rather than the primary analysis. The direction reversal is in our own data, not a hypothetical | **A** |
| 5 | A read-count (RPKM) method is barely affected by the duplication | Measured on one carrier with bowtie2 (`pipeline/mapq_split.sh`, `data/mapq_split.tsv`): 16,442 primary alignments against 16,442 distinct read names, 0 secondary; the MAPQ 20 filter removes 18.1% of *clbK* and 18.5% of *clbJ* against a 1.2% median elsewhere; the repeat's two copies take 384 and 405 reads, differing by 5.5% | Our own earlier claim that other pipelines' RPKM is inflated (item 10) - no read is counted twice, so the deficit belongs to the filtered pipeline. **Does not exclude a ~5% residual**: unfiltered *clbJ* + *clbK* carry 94.6% of the per-bp rate of the three repeat-free genes over 4 kb, so "the total is preserved" holds to about 5%, not exactly. BWA-MEM, which the withdrawn sentence named, was not run | **A** |
| 6 | *clbJ* vs *clbK* correlation in Nooij et al. Supplementary Figure 1 is not degraded (0.92, mid-range) | One reading of a figure extracted from the PDF | Our own prediction that it would be degraded | **D** - neither the extracted image nor the 57 transcribed values were kept. Irrecoverable without redoing the extraction |
| 7 | `cancer_ERR1018186` carries the island on a *K. pneumoniae* chromosome | Four island-bearing contigs tile the island with **119 bp exact overlaps** at every junction, measured from the contig sequences, 119 being this assembly's *k*; one of the four carries 386 kb of *Klebsiella* chromosome | An *E. coli* source: island depth is **33x** *E. coli*'s genome depth, which would need 33 island copies per *E. coli* cell. Uses no reference gene order | **A** |
| 7b | - the same, from the element's position | The 137,690 bp insertion the island sits in is reached at both ends by contigs **already in the automatic bin**, meeting within 8 bp on the HS11286 chromosome | A second insertion site, under the assumption that HS11286 (ST11) and this strain (ST23) share local chromosome order | **B** |
| 7c | - the same, from read pairs | `k119_70945` (63.9 kb) sends 69 pairs to bin contigs, share 0.939, inside the 0.857-1.000 that the bin's own >= 10 kb contigs occupy | A second island copy in another genome. **No power for the other three added contigs** - below 10 kb the bin's own contigs span 0.129-1.000 | **B** |
| 8 | `cancer_ERR1018282`'s island is carried by *E. coli*, not the 38x more abundant *K. pneumoniae* | Island depth 2.99, *E. coli* 3.52 (ratio 0.850), Kp 109.4 (ratio 0.0273); **plus an assembly** (`data/assembly_282_island.tsv`) in which no contig carries both island genes and chromosome | **Neither species.** 0.85 fits most *E. coli* cells at one copy; 0.0273 fits 1 *Klebsiella* cell in 37, an ordinary frequency for a mobile element. The assembly excludes only the *bulk* *Klebsiella* population (island fragments at multiplicity 3-5 against the Kp chromosome's 83-91) - which the ratio had already said. Two independent measurements, neither excluding either candidate; the verdict is parsimony. **By the stricter standard this is undetermined** | **C** |
| 9 | The island/species depth ratio has no large multiplicative offset; its scale is near 1 | 4 samples where five of six species are excluded at >= 100x give island/*E. coli* of 0.278, 0.854, 0.961, 1.026; selection does not look at the ratio (`data/island_ratio_scale.tsv`) | An offset like R's 1.9x: the ratio never passes 1.03 where *E. coli* must be the carrier. **It does not support reading a value as a cell fraction** - the band spans 3.7x on n = 4 | **B** |
| 10 | The *K. pneumoniae* abundance column needed its 14.93x rescale | The marker/genome multiplier ratio between the two species is 0.929, 95% CI 0.662-1.304 (`code/verify_marker_scaling.py`) | That `PER_PCT = 313` was right for both marker sets: 14.93 is outside the interval by a factor of 11. **Confirmed only to a factor of 2.0**, so this column must not be used to argue a two-fold difference | **B** |
| 11 | Three of the four pitfalls have published antecedents | Six targeted searches and two full-text reads (item 9); Auvray et al. 2021 and Blanco-Miguez et al. 2023 quoted | Our published novelty claim. **Does not exclude** that the fourth is also published - the search was targeted, not systematic | **B** |
| 12 | The strain is ST23 / K1, with ICEKp10 and a complete *clb* 2 | Kleborate 3.2.4 on the species-assigned bin, and on that bin plus the island-bearing contigs | **Nothing about the island's provenance.** The typing is downstream of a hand-made bin, so it cannot be evidence for what the bin contains. METHODS says this. What supports the bin's contents is rows 7, 7b, 7c | **C** |
| 13 | Breadth has no fixed attainable ceiling | A sample came in above the asserted 0.9589, then others to 0.9607 (item 1) | A derived constant. The defensible statements are the "everything lost" figure and that no sample in 397 exceeded 0.99 | **A** |
| 14 | Published breadths are low by at most 0.000952, and no call changes | Recomputed from `data/depth_hist.tsv` against the 49,431 distinct gene bases (item 18) | That it could move a call: no sample crosses the 0.5 threshold, no gene crosses the 0.5 per-gene threshold. **Bounded, not exact** - the 24 overlapping bases' coverage is not recoverable from the shipped histogram, to within 0.000486 | **B** |

## The rows to read first

Row 1, because it is the result, and its "excludes" column is nearly empty. The honest
summary of this project is that it did not find an association and could not have
found a small one.

Row 8, because it is the one claim in RESULTS.md that this ledger downgrades today.
Row 12, because a typing is the kind of output that looks like evidence and is not.
Row 6, because it cannot be rechecked at all and never will be.

Rows 7, 7b and 7c are one claim supported three ways **with different failure modes** -
overlap uses no reference order, position uses HS11286's order, linkage uses no
reference. That is why it is the strongest thing here, and it got that way only after
item 16 was rewritten twice.

## How to use this file when adding a claim

Fill the middle column first. If it is empty, the claim is a **C** and must be written
as one - "consistent with", not "shows". Three of the four errors this file exists
because of would have been caught by that one habit.

`code/verify_claims_ledger.py` checks the mechanical parts: that every `data/` path and
every CORRECTIONS item number cited above exists, and that every grade is one of
A/B/C/D. It cannot check whether a middle column is honest. Nothing can.

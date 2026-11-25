# Four ways *pks* detection from stool metagenomes goes wrong

Each item below is followed by the data in this repository that shows it, so you can
check the claim rather than take it.

## What was already known (checked 2026-11-25)

An earlier version of this file, and of the README, said these four items were "not
documented anywhere". **That was wrong, and we checked it only after publishing.**
Two of the four have their core already in the literature:

| Item | Status |
|---|---|
| 1. The duplication | **The repeat itself is published.** Auvray et al. 2021 report recombination between "two copies of a 1480 bp homologous sequence located in *clbJ* and *clbK*". What we add is its effect on read mapping, and a measurement of its extent that does not match theirs. |
| 2. Read counts are not abundance | **Both the principle and the premise are published.** Breadth-based presence calling over read counts is an established idea in metagenomics methods, and Blanco-Míguez et al. 2023 state that the AM229678 cluster contains "genes not belonging to the island". What we add is the quantification for this reference and the demonstration that *pks* studies still call by read count. |
| 3. The density-ratio offset | We found no report of this specific offset. The general class - mappability and average-genome-size biases in metagenomic gene abundance - is well established. |
| 4. Depth-limited breadth | **The principle is published.** That shallow sequencing cannot support gene presence/absence calls, and that sequencing effort confounded with treatment group is a known source of false positives, are both established. What we add is the threshold for this island (~10x) and a worked example in a real cohort. |

So the honest claim is narrower than the one we first made: **these are four places
where known methodological facts are not being applied to *pks* detection**, plus one
measurement (item 3) we could not find reported. That is still worth writing down, but
it is not four new findings.

The check was targeted, not systematic: six literature searches plus full-text reads
of the two papers named above. It is not a systematic review, and a further item may
yet turn out to be published. See [CORRECTIONS.md](CORRECTIONS.md) item 9.

---

## 1. The reference contains a ~1,290 bp duplication, so breadth cannot reach 1.0

In *E. coli* IHE3034 (GenBank AM229678.1) the same sequence appears twice:

| Copy | Coordinates | Sits inside |
|---|---|---|
| A | 16,066-17,352 | *clbK* (14,619-21,084, 6,465 bp) |
| B | 21,719-23,005 | *clbJ* (21,094-27,595, 6,501 bp) |

That is **1,287 bp at 99.845% identity**. A read from either place aligns equally well
to both, so bowtie2 assigns MAPQ 0-1 and `samtools depth -Q 20` discards it. The
positions are therefore uncovered **no matter how deep you sequence**.

### This repeat is published; its effect on mapping is what we add

Auvray et al. 2021 (*Microbial Genomics* 7:000579) found twelve *E. coli* strains
carrying a 5,651 bp deletion that fuses *clbJ* and *clbK* into one gene, and wrote:

> This deletion presumably resulted from recombination between two copies of a 1480 bp
> homologous sequence located in *clbJ* and *clbK*.

Their deletion size is almost exactly our spacing between the two copies: we measure
**5,653 bp** between equivalent positions, they report **5,651 bp**. It is the same
repeat and the same recombination.

Their only methodological caution is about assembly, not mapping:

> Caution should also be observed during the assembly of sequencing reads as errors
> including deletions may be caused by the presence of tandem repeats.

They give **no percent identity, no coordinates, and nothing about read mapping,
multi-mapping or coverage**. The consequence below - that a MAPQ filter silently caps
breadth near 0.95 - is what this repository adds.

### We cannot reproduce their 1,480 bp

Measured against AM229678.1 three ways, the repeat is ~1,290 bp, not 1,480 bp:

| Method | Extent | Identity |
|---|---|---|
| Exact 30-mer matches | 16,066-17,352 (**1,287 bp**) | 99.845% |
| Windowed identity >= 98% | 16,063-17,354 (1,292 bp) | 99.613% |
| Windowed identity >= 70% | 16,039-17,375 (1,337 bp) | 97.532% |

Even relaxed to 70% windowed identity it does not reach 1,480 bp. An exhaustive 30-mer
scan of the whole island finds 1,206 sequences occurring more than once, of which
**1,198 are at a spacing of exactly 5,653 bp**, and **no other spacing between 5,000
and 6,500 bp exists** - so the repeat is not continued further by an indel-shifted
segment. (The remaining duplicated k-mers are the 8 bp period VNTR in the
*clbB*-*clbR* intergenic region, which Auvray et al. describe separately.)

We do not claim their figure is wrong. They worked from their own hybrid assemblies
and may have measured in a different strain or under a different definition of
"homologous". We report what is measurable in the reference everyone maps to. Anyone
relying on 1,480 bp for AM229678.1 should check it.

### What this costs

If both copies were lost entirely, breadth would be `(49455 - 2574) / 49455 = 0.9480`.
In practice the loss is smaller, because a read overlapping the edge of the repeat has
unique sequence outside it and can still be placed.

In the 13 samples with island depth >= 10x:

```
missing bases inside the duplication:   1,920 - 2,140  (mean 2,000)
samples with zero missing bases outside the duplication:  11 of 13
```

The duplication cost is **nearly constant and independent of depth**; everything else
is depth-limited. Across all 397 samples the highest breadth observed is **0.9607**,
and **no sample exceeded 0.99**.

### Consequences

- **A sample at 0.95 breadth has a complete island.** Do not read it as 95% present.
- A breadth threshold set near 1.0 will reject every true carrier.
- If you see a sample above 0.99, suspect a reference or filtering difference, not a
  better sample.

Reproduce: `data/missing_split.tsv`, columns `missing_in_duplication` and
`missing_outside_duplication`; generated by `code/make_dup_table.py`.

---

## 2. Raw read counts against the island are not a measure of island abundance

The AM229678.1 reference is **55,140 bp**, but the 18 *clb* genes cover only
**49,455 bp**. The remaining **5,685 bp** is intergenic sequence and flanking
backbone, which is *not* island-specific.

Reads from ordinary, island-free *E. coli* land on that 5,685 bp. In 82 samples that
were island-negative (breadth < 0.1) yet still had >= 50 reads mapped to the
reference:

```
median breadth over the clb genes:        0.0038
median ratio of pks reads to E. coli reads: 0.0116
```

So the read count is non-zero while the island is absent.

### This one is not new, and that is the problem

That presence should be called from **breadth of coverage and the distribution of
reads across a reference, not from read counts**, is an established idea in
metagenomics methods work - the standard argument being that mismapped reads pile onto
conserved and mobile regions and inflate apparent presence. And the specific premise
for this reference is published too: Blanco-Míguez et al. 2023 (*mSystems*
8:e00079-23), working on the same AM229678 cluster, state that

> there are genes not belonging to the island distributed in the inner part of the
> cluster

and their Fig. 1 legend labels them. They handle it by requiring 20 of 23 genes to be
detected individually, rather than by counting reads over the cluster.

**The pitfall is therefore not that this is unknown. It is that *pks* metagenomics
papers still report island abundance as reads or RPKM over the cluster**, which
inherits exactly the background measured above. What this repository adds is the size
of that background for this reference: 5,685 bp of non-gene sequence, and a median
`pks`/*E. coli* read ratio of 0.0116 in samples that have no island.

### Consequences

- Never use "reads mapped to the island" as a quantity. Use breadth over the gene
  intervals, or depth restricted to those intervals.
- Any correlation you compute between "pks reads" and a species abundance is
  partly a correlation between that species and itself.

---

## 3. A marker-gene density ratio carries a large systematic offset

Of the four items this is the one we could find no report of. The general class of
bias is well established - mappability differences and average genome size are named
causes of spurious variation in metagenomic gene abundance, and marker-gene based
normalisation exists precisely to correct them - but we found nothing giving the
magnitude for an island-to-marker density ratio. Treat the decomposition below as our
own calibration of our own pipeline, not as a constant.

To ask "could species X account for this island?", we compared the read density on the
island to the read density on that species' marker genes:

```
R = (island reads / 49,455) / (marker reads / marker bp)
```

For an organism carrying one island per genome, R should be about 1. Ours clustered
near **1.9**. The offset decomposes as:

| Factor | Size | Cause |
|---|---|---|
| Length convention | 1.115 | reads counted over 55,140 bp, divided by 49,455 bp |
| Marker edge loss | 1.14 | fragments lose coverage at their ends; the island does not |
| Residual mappability | ~1.50 | marker sequence is not as mappable as island sequence |

The product is ~1.9. It was stable across read lengths of 84, 94 and 100 bp.

### Consequences

- **R is not a population fraction.** Do not say "40% of the *E. coli* carries it".
- The only inference that survives the offset is the one-sided one: **R >> 1 means
  that species cannot account for the island**, because a 700-fold gap is not a 5%
  calibration error. In one sample R was 723 for *E. coli* and 0.79 for
  *K. pneumoniae*; that conclusion is robust to the offset.
- Conversely, R near 1 for a species is *consistent with* carriage, not proof of it.
  In one of our samples *K. pneumoniae* was present at 16.7x with 0.94 genome breadth
  and was still **not** the carrier - *E. coli* was.

---

## 4. Breadth is depth-limited, so depth differences imitate carriage differences

Breadth saturates only once the island is covered several times over. Below that,
"absent" and "not sequenced deeply enough" are the same observation.

```
island depth >= 10x   n =  13   11 of 13 have zero missing bases outside the duplication
island depth 2-10x    n =  18   mean 2,597 bp missing outside the duplication
island depth <  2x    n = 366
```

Only 13 of 397 samples are depth-saturated. 34 of 44 positive calls sat below 10x.

### Why this is dangerous

In our DE-CRC cohort the two groups differ in sequencing depth by a wide margin:

```
CRC     median 69.5M reads   (22.2 - 120.4M)
control median 34.9M reads   (16.3 -  53.7M)
Mann-Whitney P = 1.5e-14
```

Calling carriage on the raw data would find more carriers in the cancer group for a
purely technical reason. In an earlier phase of this work, a logistic regression on a
different cohort showed exactly that: adding read depth as a covariate moved the odds
ratio from 3.11 to 1.31 and P from 0.080 to 0.71.

### This principle is established; the threshold for this island is not

That shallow sequencing cannot support gene presence/absence calls is known - published
estimates put the depth needed to quantify rare genes in fecal metagenomes at around
5 million merged reads - and "sequencing effort confounded with treatment group" is
named in review articles as a principal source of false positives, with rarefaction
the standard remedy. None of that is ours.

Two things here are: the **~10x island depth** figure at which breadth stops being
depth-limited for *this* 49,455 bp target, and a worked case where the confounding is
real rather than hypothetical (the DE-CRC depth split above, P = 1.5e-14). We also note
that the published depth guidance is for *quantifying* genes, not for *presence calls
on a single 49 kb island*, which is a different and more demanding problem: a gene
present at one copy per carrier cell needs the carrier itself to be abundant.

### What to do about it

Rarefy every sample to a common read count **before** calling, and report how many
samples you had to drop for being below it. We used 20M reads, which left 167 of 397
samples judgeable after also applying the abundance gate.

You do not need the FASTQ to do this. See the derivation in [METHODS.md](METHODS.md):
the expected breadth after retaining a fraction `f` is `mean_i [1 - (1-f)^d_i]`, exact
by linearity of expectation, computable from depth histograms alone.

---

## Sources

Read in full and quoted above:

- Auvray F, Perrat A, Arimizu Y, *et al.* Insights into the acquisition of the *pks*
  island and production of colibactin in the *Escherichia coli* population.
  *Microbial Genomics* 2021;7(5):000579. doi:10.1099/mgen.0.000579
- Blanco-Míguez A, *et al.* Targeted depletion of *pks*+ bacteria from a fecal
  microbiota using specific antibodies. *mSystems* 2023;8(3):e00079-23.
  doi:10.1128/msystems.00079-23

Not read in full. These were located by literature search and are given as entry
points for the general principles in items 2, 3 and 4, not as verified support for any
specific number here. **Check them yourself before citing them.**

- Sanguineti *et al.*, Metapresence: species detection from the genome-wide
  distribution of mapping reads. *mSystems* 2024. (breadth vs read counts)
- imGLAD, PMC6216955. (mismapping onto conserved and mobile regions)
- Hakim *et al.*, "genome cover". *mSphere* 2022, PMC9600373.
- Manor O, Borenstein E. MUSiCC. *Genome Biology* 2015;16:53.
  doi:10.1186/s13059-015-0610-8 (mappability as a named sample-level bias)
- Nayfach S, Pollard KS. MicrobeCensus. *Genome Biology* 2015;16:51.
  doi:10.1186/s13059-015-0611-7 (average genome size)

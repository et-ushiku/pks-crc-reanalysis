# Four ways *pks* detection from stool metagenomes goes wrong

Each item below is followed by the data in this repository that shows it, so you can
check the claim rather than take it.

## What was already known (checked 2026-11-25)

An earlier version of this file, and of the README, said these four items were "not
documented anywhere". **That was wrong, and we checked it only after publishing.**
Three of the four have their core already in the literature:

| Item | Status |
|---|---|
| 1. The duplication | **The repeat itself is published.** Auvray et al. 2021 report recombination between "two copies of a 1480 bp homologous sequence located in *clbJ* and *clbK*". What we add is its effect on read mapping, and a measurement of its extent that does not match theirs. **Kleborate's own colibactin typing scheme already excludes both genes**: it types 15 loci, *clbA*-*clbI* and *clbL*-*clbQ*, and has no column for *clbJ* or *clbK* - see `data/kleborate/kp_bin_output.tsv`. It is a pitfall **for breadth-based methods, ours included**, and barely touches read-count methods - see [the correction](#item-1-is-a-pitfall-for-breadth-based-methods-which-includes-ours). |
| 2. Read counts are not abundance | **Both the principle and the premise are published.** Breadth-based presence calling over read counts is an established idea in metagenomics methods, and Blanco-Míguez et al. 2023 state that the AM229678 cluster contains "genes not belonging to the island". What we add is the quantification for this reference and the demonstration that *pks* studies still call by read count. |
| 3. The density-ratio offset | We found no report of this specific offset. The general class - mappability and average-genome-size biases in metagenomic gene abundance - is well established. |
| 4. Depth-limited breadth | **The principle is published, and stated for this island.** Nooij et al. 2021 write that "metagenomics-based assays require high sequencing depths to reliably detect single genes" and subsample to a common 10 million read pairs. What we add is the threshold at which breadth stops being depth-limited (~10x island depth) and a worked example where the confounding is real. |

So the honest claim is narrower than the one we first made: **these are four places
where known methodological facts are not being applied to *pks* detection**, plus one
measurement (item 3) we could not find reported. That is still worth writing down, but
it is not four new findings.

The check was targeted, not systematic: six literature searches plus full-text reads
of the four papers named in [Sources](#sources). It is not a systematic review, and a
further item may yet turn out to be published. See
[CORRECTIONS.md](CORRECTIONS.md) item 9.

For items 2 and 4 the claim "known but not applied" is checkable against a specific
published study, whose methods we read in full: see
[A worked example](#a-worked-example-of-items-2-and-4-in-a-published-study) below.

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

### A typing scheme already leaves these two genes out

Kleborate (v3.2.4) assigns a colibactin sequence type from **15 loci: *clbA* through
*clbI*, and *clbL* through *clbQ*.** There is no *clbJ* or *clbK* column in its output
at all - visible directly in `data/kleborate/kp_bin_output.tsv`, which this repository
ships. Those are exactly the two genes that contain the duplication described above.

We do not know the scheme's stated reason for omitting them and are not claiming one;
the output does not say. But a widely used typing tool treating *clbJ* and *clbK* as
not usable is further evidence that the problem is known, and one more reason to read
our item 1 as "not applied to *pks* detection from stool" rather than as new.

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

**The repeat also turns up on its own in an assembly.** `cancer_ERR1018186` was
assembled de novo and the island came out in four contigs. One of them, `k119_32464`,
is 1,248 bp, aligns to the island at **16,066** and again at **21,719** - the same
start and the same 5,653 bp spacing measured above - and carries twice the k-mer
multiplicity of its neighbours, which is what a collapsed two-copy repeat has to do.
The assembler broke the island at exactly the two copies. That measurement shares no
inputs with the k-mer scan above: see [CORRECTIONS.md](CORRECTIONS.md) item 16 and
`data/kleborate/island_tiling.tsv`.

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

Per gene, in those same 13 samples, the loss is confined to exactly the two genes
Auvray et al. named and nothing else:

| Gene | Length | Mean breadth at >= 10x |
|---|---|---|
| ***clbK*** | 6,465 bp | **0.8065** |
| ***clbJ*** | 6,501 bp | **0.8450** |
| the other 16 *clb* genes | - | **0.9982** (mean) |

**19.3% of *clbK* and 15.5% of *clbJ* are unreachable, against 0.2% for everything
else.** This is the cleanest statement of the problem: at sufficient depth the island
is complete apart from the repeat, and the repeat sits inside those two genes.
Reproduce it from `data/depth_hist.tsv` and `data/clb_genes.bed`.

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

Reads from ordinary, island-free *E. coli* land on that 5,685 bp. Among the **106
samples** that were island-negative (breadth < 0.1) yet still had at least 50 reads
mapped to the reference:

```
median breadth over the clb genes:        0.0018
median reads mapped to the reference:       258
maximum reads mapped to the reference:   10,548
```

So the read count is not merely non-zero: it reaches four figures in samples where the
island is absent. Both columns are in `data/per_sample.tsv` (`island_reads` and
`island_breadth`), so this is checkable - it was not, until this audit, because the
read counts were not shipped. See [CORRECTIONS.md](CORRECTIONS.md) item 13.

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
of that background for this reference: 5,685 bp of non-gene sequence, and a median of
258 reads - up to 10,548 - landing on it in samples that have no island.

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
| Length convention | 1.115 | reads counted over 55,140 bp, divided by 49,455 bp (= 55140/49455, exact) |
| Marker edge loss | 1.14 | fragments lose coverage at their ends; the island does not |
| Residual mappability | ~1.50 | marker sequence is not as mappable as island sequence |

The product is ~1.9. It was stable across read lengths of 84, 94 and 100 bp.

### Consequences

- **R is not a population fraction.** Do not say "40% of the *E. coli* carries it".
- The only inference that survives the offset is the one-sided one: **R >> 1 means
  that species cannot account for the island**, because a 700-fold gap is not a 5%
  calibration error. In one sample R was 764 for *E. coli* and 0.79 for
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
island depth 2-10x    n =  18   mean 2,585 bp missing outside the duplication
island depth <  2x    n = 366
```

Only 13 of 397 samples are depth-saturated. 29 of the 42 positive calls sat below 10x.

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
ratio from **3.11 to 1.78** and P from **0.080 to 0.40**, with depth itself
significant at **P = 0.010**. Run `python3 code/verify_extra.py` to reproduce that,
the depth split above, and the ratio decomposition below; all three were unverifiable
from this repository until `island_reads` was shipped (CORRECTIONS.md items 13-14).
An earlier version of this paragraph gave the adjusted odds ratio as 1.31 and
P as 0.71, from a third model whose specification was not recorded; it is dropped
rather than reconstructed.

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

## Item 1 is a pitfall for breadth-based methods, which includes ours

**Correction to an earlier version of this file**, which said that omitting a MAPQ
filter "inflates RPKM for *clbJ* and *clbK* specifically". **That is wrong.** BWA-MEM
reports one primary alignment per read, so a read from the repeat is counted once,
placed arbitrarily at one of the two copies. The *total* over *clbJ* + *clbK* is
preserved; what is unreliable is the split between them. Nothing is inflated.

The consequence of this is worth stating plainly, because it cuts against us:

- **A read-count method that averages over the genes is barely affected.** The repeat's
  reads are retained and the sum is right, so a mean RPKM across 19 *clb* genes does
  not care that 1,287 bp of it is ambiguous.
- **A breadth method is affected a lot.** Discarding those reads at MAPQ 20 removes
  ~19% of *clbK* and ~16% of *clbJ* permanently, which is why breadth stops at ~0.95.

So item 1 is **not** a defect we found in other people's read-counting pipelines. It is
a property of the reference that anyone calling presence from coverage - us included -
has to know about, and we hit it because of a design choice we made. The per-gene table
above is the evidence either way.

## The upstream source of the one-read rule, and what it gets right

The positive-call rule in the worked example below traces to:

**Nooij S, Ducarmon QR, Laros JFJ, et al. Fecal Microbiota Transplantation Influences
Procarcinogenic *Escherichia coli* in Recipient Recurrent *Clostridioides difficile*
Patients. *Gastroenterology* 2021;161(4):1218-1228.e5.
doi:10.1053/j.gastro.2021.06.009** (CC BY; read in full from the Leiden University
repository, hdl.handle.net/1887/3249076, the publisher's own pages being blocked)

Their rule, verbatim, is weaker than its restatement elsewhere - **one read, not one
gene**:

> Individuals with 0 reads mapped to *clb* genes in their metagenome are considered
> *pks*-negative (*pks*-) and individuals with **1 or more reads mapped to 1 or more
> *clb* genes are *pks*+**.

With that threshold they report 55% of patients and 29% of donor samples positive.

**But this paper anticipates two of our four items, and we should say so before
criticising anything.** It was the upstream method for the example below, and it is
more careful than that example conveys.

It rarefies, in effect, before analysis:

> Raw metagenomic reads were randomly subsampled to 10 million read pairs

(stated as being "to facilitate de novo assembly", and their sequencing target was
10 million reads per sample, so depth is roughly equalised either way - a cap rather
than an exact match, since samples below 10M keep everything.)

It detects homology-driven coverage distortion and **removes the affected genes**,
which is the reasoning behind our items 1 and 3, applied to marker genes:

> Based on the depth of coverage, we found that the MLST genes sometimes showed uneven
> coverage, especially in metagenomes with more non-*E coli* Enterobacteriaceae. This
> indicates that these genes may have regions that are homologous to genes from other
> species, so the MLST genes were excluded.

And it states item 4's principle outright, and combines genes to address it:

> Metagenomics-based assays require high sequencing depths to reliably detect single
> genes. We combined data of all 19 *clb* genes, which we showed are highly
> intercorrelated, to improve the sensitivity of our assay and allow us to detect very
> low levels of *pks*.

They also flag their own detection floor - "Levels of *pks* are generally very low,
near the limit of detection" - which is item 4 by another name.

**What is missing is item 1 only**: no mapping-quality filter and no mention of the
*clbJ*/*clbK* repeat. Per the correction above, that costs their headline quantity
almost nothing.

We predicted it would at least show up in the gene-to-gene correlation matrix they use
to argue the island travels intact, since the repeat's reads are apportioned
arbitrarily between those two genes. **We checked, and it does not.** Their
Supplementary Figure 1 gives ***clbJ* vs *clbK* = 0.92**, squarely inside the ordinary
range for *clb* pairs in that matrix (about 0.79 to 0.99, excluding *clbS*, which is an
outlier throughout at 0.54-0.94). The reason is plain in hindsight: the ambiguous
1,287 bp is only about a fifth of each 6.5 kb gene, the other four fifths are unique,
and arbitrary apportionment averages out over many reads. So the duplication leaves no
trace in their analysis at all. See [CORRECTIONS.md](CORRECTIONS.md) item 11.

### What their Supplementary Figure 1 does show

Reading the same matrix for a different purpose turned up something that bears on the
attribution question above. Correlations of RPKM across their metagenomes:

**These figures cannot be rechecked from anything in this project.** They were read
off the figure after extracting it from the PDF, and neither the extracted image nor
the 57 transcribed cell values was kept - so the median and the count of negatives
below are our arithmetic over a transcription that no longer exists. Anyone relying on
them should read the figure themselves. See [CORRECTIONS.md](CORRECTIONS.md) item 17.

| Pair | Spearman r |
|---|---|
| *clb* gene vs *clb* gene | ~0.79 to 0.99 (excluding *clbS*) |
| The three *E. coli* single-copy markers vs each other | 0.87, 0.87, 0.99 |
| Flanking mobile elements (IS1400 transposases, *intP4*) vs the *E. coli* markers | 0.47 to 0.76 |
| **The 19 *clb* genes vs the three *E. coli* markers** | **median +0.02 over 57 pairs; 24 of the 57 are negative** |

**In their own data, *pks* abundance does not track *E. coli* abundance** - while the
mobile elements sitting in the same reference region do track it, as generic *E. coli*
content should. Their text reports the same thing for one donor: *pks* "correlated
loosely with the relative abundance of *E coli* (Spearman's r = 0.43, P = .11)".

This is **not** proof of a non-*E. coli* carrier. If only one of several co-resident
*E. coli* strains carries the island, *pks* tracks that strain rather than total
*E. coli*, and the correlation collapses just the same - which is the conventional
reading and may well be right. The point is narrower: the figure offered in support of
the island travelling intact is **neutral at best on who carries it**, and the question
is settled by measuring the candidate genomes, which is cheap.

(Numbers in the table above were read off the published figure image, not from
underlying data, so treat the individual cells as approximate.)

### One typo, recorded because it is checkable

Their Figure 1 and Supplementary Figure 1 legends both give the island reference as
**AM227896**. That accession is *Paracentrotus lividus* EST clone MPMGp1174L157Q, a
771 bp sea urchin mRNA. The Methods give the correct **AM229678** (*E. coli* IHE3034
colibactin cluster, 55,140 bp) in all four places it appears, so this is a legend
typo and not a methodological error. We note it only because anyone following the
figure legend to the reference will not find an island there.

### Where our data disagrees with them

Having found *pks* reads in samples with no detectable *E. coli*, they considered the
alternative carriers and ruled them out:

> Although we detected *Klebsiella*, *Enterobacter*, and *Citrobacter* spp using the
> assembly based approach in 6 of 11 donor samples with *pks* but that were negative
> for *E coli* based on the 3 single-copy marker genes, pathogenic isolates of these
> species are unlikely to be present in healthy donors due to thorough screening. We
> therefore consider it more likely that *E coli* was present below the detection limit
> of our assay and is in fact the host of the *pks* island.

**This is the attribution default our carrier analysis argues against**, and it is the
one place where our data speaks directly to theirs. In 6 of 11 relevant samples they
detected a plausible alternative carrier and preferred *E. coli* anyway, on the
grounds that *pathogenic* isolates are screened out of healthy donors. *K. pneumoniae*
is a common gut commensal, and screening for pathogens does not exclude it.

Our [RESULTS.md](RESULTS.md) finds, in 15 positives measured against six complete
genomes in one competitive index: 3 better explained by *K. pneumoniae* than by
*E. coli*, **two of those three being controls rather than cases**, one of them
confirmed by assembly with the island and the *Klebsiella* chromosome on a single
contig - and 4 of 14 that cannot be assigned to any species at all. We also have a
sample where *K. pneumoniae* sits at 109x genome depth and is still not the carrier,
so the inference runs both ways and neither default is safe.

We are not claiming they got their samples wrong; we have not seen their data and
their cohort is FMT donors and *C. difficile* patients, not ours. We are claiming the
reasoning - "the alternative is unlikely in this population, therefore *E. coli*" -
is the step that needs the measurement, and that the measurement is now cheap.

## A worked example of items 2 and 4 in a published study

The claim above is that these are known facts not applied to *pks* detection. Here is
one study where that can be checked, because its methods are explicit. We read the
full text (Europe PMC, PMC10797792).

**Terlouw D, Boot A, Ducarmon QR, et al. Enrichment of colibactin-associated
mutational signatures in unexplained colorectal polyposis patients. *BMC Cancer*
2024;24:102. doi:10.1186/s12885-024-11849-y**

This is a study of colorectal polyposis and mutational signatures, in which fecal
metagenomics on 25 patients is a secondary analysis. It is not a CRC case-control
study and it is not comparable in scope to this reanalysis. We cite it only because
its *pks* detection is documented well enough to compare, and because we could not
find a *pks* metagenomics paper that does apply breadth or a mapping-quality filter.

Their detection, verbatim from Materials and methods:

> The quality-trimmed reads were screened for the presence of the *pks* island by
> mapping to the colibactin gene cluster (accession ID AM229678) using BWA-MEM
> (version 0.7.17). Mapped reads were deduplicated using Picard MarkDuplicates
> (version 2.23.3) to remove technical artifacts and improve quantification.

The positive call, same section (the rule is attributed to an earlier paper by Nooij
et al., which we have not read):

> As previously outlined, fecal samples positive for **at least one *pks* gene** were
> considered *pks*-positive.

The quantity:

> For mean RPKM of the entire *pks* island, RPKMs of all individual genes were summed
> and divided by the total number of 19 *clb* genes.

Against the four items:

| | |
|---|---|
| Item 1, the duplication | Not mentioned, and **mostly does not hurt them** - see the correction below. No MAPQ filter appears anywhere, so the repeat's reads are kept rather than discarded. |
| Item 2, read counts | **This is the pitfall exactly.** Abundance is RPKM over the cluster. The 5,685 bp of AM229678 that is not *clb* gene sequence is excluded by using per-gene coding sequences, which is better than counting over the whole cluster - but the presence call is one gene's worth of reads, with no breadth and no minimum depth. |
| Item 3, density ratios | Not applicable; they do not compute one. |
| Item 4, depth | RPKM normalises *abundance* by library size, but the **presence call does not**. "At least one *pks* gene" is maximally depth-sensitive: a deeper sample is more likely to be called positive for a reason that has nothing to do with carriage. The paper does not report whether sequencing depth differed between its case and control groups. |

Their *pks* result, verbatim from Results:

> 59% (10 out of 17) of the cases were *pks* positive compared to 25% (2 out of 8) of
> controls (p-value = 0.124).

and

> no significant correlation between number of adenomas/carcinomas with *APC* variants
> fitting SBS88 or ID18 and the *pks* RPKM values was observed (Pearson: R = 0.16,
> p-value = 0.45).

**Both are null.** The abstract describes the same comparison as having "revealed
enriched presence of *pks* genes in patients with somatic variants fitting
colibactin-associated signatures", without the P value. We note this because it is the
pattern our own [CORRECTIONS.md](CORRECTIONS.md) item 4 records us falling into: a
headline written from the direction of a non-significant difference.

Their stated limitations are the small sample and the FFPE tissue, plus that "shotgun
fecal metagenomics is a broad analyses but a more sensitive qPCR approach ... could
give more insight". Sensitivity is acknowledged; depth confounding and breadth are not
discussed.

### One further observation, offered as a question rather than a criticism

> None of the bacteria previously associated with colorectal cancer, like *Fusobacterium
> nucleatum*, *Bacteroides fragilis*, *Campylobacter jejuni* and *Clostridioides
> difficile*, or capable of producing colibactin, like *K. pneumoniae*, *E. aerogenes*
> and *C. koseri*, were detected in any of the stool samples (data not shown).

In our 397 samples, *K. pneumoniae* is detected in **96 of the 221 samples where we
measured it (43%)**, mostly at low abundance - **median 0.0143%**, but **maximum
6.83%**, and 19 of the 101 CN-CRC samples are above 0.1%. At the median a detection is
a handful of marker reads, so a profiler with a stricter floor reporting zero is **not
necessarily inconsistent** with our numbers there, and we cannot tell which it was: no
threshold is given and the data are not shown. **At the top of the range that
reconciliation does not work** - a sample at 6.83% is not below anybody's detection
floor - so for those samples the two results genuinely disagree.

(These figures were **0.0010% and 0.4576%** until 2026-11-26, when the published
column was found to be scaled by the *E. coli* marker length rather than
*K. pneumoniae*'s, understating it 14.93-fold. The detection count and rate are
unchanged, because a scale factor cannot change which values are non-zero. See
[CORRECTIONS.md](CORRECTIONS.md) item 17.)

It matters because of what it forecloses. If no *Klebsiella* is detected, no *pks*
island can be attributed to *Klebsiella*, and attribution defaults to *E. coli* by
construction. Across the 23 positives where we looked ([RESULTS.md](RESULTS.md)), the
island is better explained by *K. pneumoniae* than by *E. coli* in 3, and among the 14
of those measured against six competing genomes, 4 cannot be assigned to any species
at all - neither verdict is reachable in a pipeline where the competing species is
absent by definition.

---

## Sources

Read in full and quoted above:

- Auvray F, Perrat A, Arimizu Y, *et al.* Insights into the acquisition of the *pks*
  island and production of colibactin in the *Escherichia coli* population.
  *Microbial Genomics* 2021;7(5):000579. doi:10.1099/mgen.0.000579
- Blanco-Míguez A, *et al.* Targeted depletion of *pks*+ bacteria from a fecal
  microbiota using specific antibodies. *mSystems* 2023;8(3):e00079-23.
  doi:10.1128/msystems.00079-23
- Terlouw D, Boot A, Ducarmon QR, *et al.* Enrichment of colibactin-associated
  mutational signatures in unexplained colorectal polyposis patients. *BMC Cancer*
  2024;24:102. doi:10.1186/s12885-024-11849-y (PMC10797792)
- Nooij S, Ducarmon QR, Laros JFJ, *et al.* Fecal Microbiota Transplantation
  Influences Procarcinogenic *Escherichia coli* in Recipient Recurrent *Clostridioides
  difficile* Patients. *Gastroenterology* 2021;161(4):1218-1228.e5.
  doi:10.1053/j.gastro.2021.06.009 (CC BY; read from hdl.handle.net/1887/3249076)

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

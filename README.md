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

The agreement is one of conclusions, not of methods. Their calls come from Anpan's
gene model, whose sample filter is a per-species k-means boundary fitted to the data
rather than a fixed threshold, so the two analyses cannot be compared on sensitivity
and we do not try. [METHODS.md](METHODS.md) sets out what their criterion actually is
and the three things about it we could not resolve from their text.

**This is not a diagnostic test and must not be used for clinical decisions.**

---

## What this adds

The negative result is a second, independent confirmation of something already
published, so its novelty is modest. Three other things came out of the work, in
descending order of how well they survived checking against the literature.

### 1. Who carries the island, measured instead of assumed

This is the part we would point a reader at first.

All 14 *pks*-positive CN-CRC samples were mapped to **six complete
Enterobacteriaceae genomes in a single competitive index** at MAPQ 20 - *E. coli*
IHE3034, *K. pneumoniae* HS11286, *C. freundii* CFNIH1, *C. koseri* BAA-895,
*E. cloacae* ATCC13047 and *K. aerogenes* KCTC2190 - so that reads from the genomes
compete and ambiguous ones are dropped. The readout is island depth divided by each
species' genome depth, from the same read set. Full table in
[RESULTS.md](RESULTS.md).

| Finding | |
|---|---|
| **Cannot be assigned to any species** | **4 of 14.** *E. coli* and *K. pneumoniae* are both abundant and both ratios fit. In one control, *E. coli* sits at ratio 0.76 and Kp at 0.97, both with genome breadth above 0.7. |
| Better explained by *K. pneumoniae* than *E. coli* | 3 of 15 across both cohorts examined - and **two of the three are controls**, which is the opposite of what a disease association predicts |
| Confirmed by assembly | 1, with the island and the *Klebsiella* chromosome on a single contig |
| **Abundance is not carriage** | Two samples carry *K. pneumoniae* at 16.7x and **109x** genome depth (breadth > 0.88) where the island nonetheless tracks *E. coli* |

**Species attribution of this island from stool reads fails outright in roughly a
quarter to a third of positive samples, and it fails precisely where two
Enterobacteriaceae are co-abundant.** We have not found a paper that reports what
fraction of its positives it could not assign.

What the literature does instead is default to *E. coli*:

- Nooij et al. 2021 found *pks* reads in donor samples with no detectable *E. coli*,
  **detected *Klebsiella*, *Enterobacter* or *Citrobacter* in 6 of 11 such samples**,
  and concluded *E. coli* was "in fact the host" on the grounds that pathogenic
  isolates are screened out of healthy donors. *K. pneumoniae* is a common commensal;
  pathogen screening does not exclude it.
- Their own Supplementary Figure 1 shows the 19 *clb* genes correlating with three
  *E. coli* single-copy markers at a **median of +0.02 over 57 pairs**, 24 of them
  negative, while the markers correlate with each other at 0.87 to 0.99. *pks*
  abundance does not track *E. coli* abundance in their data. (A single *pks*+ strain
  among several co-resident *E. coli* would also produce this, so it is not proof of
  another carrier - it is an absence of support for *E. coli*.)
- Terlouw et al. 2024 report that *K. pneumoniae*, *E. aerogenes* and *C. koseri* were
  "not detected in any of the stool samples (data not shown)" across 25 deep
  metagenomes. A pipeline that detects no *Klebsiella* attributes every island to
  *E. coli* by construction.

Quotations and sources for all three are in [PITFALLS.md](PITFALLS.md).

### 2. The whole analysis reproduces from 753 KB, and the shortcut is checked

The four cohorts are roughly 2.5 TB of FASTQ and none of it is needed, because every
call is a function of the per-gene depth histograms. The substitution is exact by
linearity of expectation, and it is also checked against FASTQ that really was
subsampled and remapped: largest gap **0.0211** over five runs, smaller than the
**0.0182** spread between seeds at the same depth, with **all five calls agreeing**.
See the next section to run it, and [METHODS.md](METHODS.md) for the derivation and
what the check does not settle.

### 3. Four measurements on this reference, mostly with published antecedents

[PITFALLS.md](PITFALLS.md) describes four ways *pks* detection from stool metagenomes
goes wrong. **An earlier version of this README claimed none of the four was
documented anywhere. That was wrong**, and the file now names and quotes the
antecedents; see [CORRECTIONS.md](CORRECTIONS.md) items 9, 10 and 11. Items 3 and 4
are in substance anticipated by Nooij et al. 2021, which already subsamples to a
common depth and already excludes marker genes after finding homology-driven uneven
coverage in them. **Item 1 turned out to be a pitfall for breadth-based calling,
which is our own method, and barely affects read-count methods** - we had that
backwards twice. What is left is a set of measurements:

1. **A ~1,290 bp duplication caps breadth at ~0.95.** In AM229678.1 the same 1,287 bp
   appears at 16,066-17,352 (in *clbK*) and 21,719-23,005 (in *clbJ*), 99.845%
   identical; those reads are discarded at MAPQ 20. At >= 10x island depth the loss is
   confined to those two genes - breadth 0.8065 for *clbK*, 0.8450 for *clbJ*, against
   0.9982 for the other 16. The repeat was published by Auvray et al. 2021 as
   "1480 bp", which we cannot reproduce at any identity threshold.
2. **The 55,140 bp reference contains 5,685 bp that is not *clb* gene sequence**, so
   island-negative samples still accumulate reads on it: median background
   *pks*/*E. coli* read ratio 0.012 across 82 such samples, at breadth ~0.004.
3. **A marker-gene density ratio carries a ~1.9x systematic offset**, decomposed as
   1.115 x 1.14 x ~1.50. Usable for exclusion, not attribution. This is the only one
   of the four for which we found no published antecedent.
4. **Breadth is depth-limited below ~10x island coverage**; 34 of 44 positive samples
   sat below it. In one cohort the between-group depth difference is P = 1.5e-14, so
   calling on raw data would find more carriers in the cancer group for a purely
   technical reason.

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
| [RESULTS.md](RESULTS.md) | **The carrier analysis** (all 14 positives against six genomes), plus every result table per cohort and pooled |
| [METHODS.md](METHODS.md) | Pipeline, reference, pre-specified thresholds, rarefaction derivation |
| [PITFALLS.md](PITFALLS.md) | The four measurements above, what was already published about each, and a worked example in a published study |
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
- **The carrier analysis rests on 14 samples and one assembly.** Of the three
  *K. pneumoniae* attributions, only one is confirmed by linkage on a single contig;
  the other two rest on coverage ratios, which are consistent with carriage but do not
  prove it, and their FASTQ was deleted after measurement so the check cannot be run
  now. The "4 of 14 unassignable" figure is a property of these cohorts at these
  depths, not a constant. We also read one sample as *Klebsiella* on abundance alone
  before measuring it, and were wrong (CORRECTIONS.md item 8).
- **We detect sequence, not function.** Auvray et al. 2021 found that **16 of 109**
  *pks*-positive *E. coli* strains (15%) had lost genotoxicity through mutations in
  *clbB* or *clbD* or a *clbJK* fusion, while mostly still producing small amounts of
  island-derived metabolites. A metagenomic positive call means the genes are present,
  not that colibactin is made. Nothing in this analysis can distinguish a functional
  island from a broken one, and a 15% inactive fraction would dilute any real
  association toward the null.
- **No wet-lab work.** No isolates, no PCR, no colibactin measurement.
- One operator, one analysis, no independent replication of the code.

## Related work by the same author

[et-ushiku/crc-microbiome-reproduction](https://github.com/et-ushiku/crc-microbiome-reproduction)
reanalyses the same kind of public stool metagenomes for a different question:
whether the **community as a whole** separates CRC from control. It does, at a
leave-one-study-out mean AUC of 0.82 over 9 cohorts and 824 people.

Put beside the present result, the pair says something worth stating plainly: **the
community-level signal is real and reproducible, and the carriage of this one island
is not where it lives.** A single mobile element, even one with a demonstrated
mutational mechanism, is not the thing that distinguishes these groups. Neither result
is a screening test, and the limitations in both repositories should be read before
either number is quoted.

## License

Code: MIT (see [LICENSE](LICENSE)).
Documentation and derived data tables: CC BY 4.0.

Underlying sequence data belongs to the original studies and remains under their terms.

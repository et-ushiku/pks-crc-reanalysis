"""Check the three figures that used to be unverifiable from this repository.

Until `island_reads` was added to `data/per_sample.tsv` (CORRECTIONS.md item 13),
three numbers quoted in PITFALLS.md could not be recomputed by a reader, because
their inputs were not shipped:

  1. the DE-CRC between-group sequencing depth difference, P = 1.5e-14
  2. the ~1.9x systematic offset in the marker-gene density ratio R
  3. the logistic regression in which adjusting for depth removed a group effect

All three are now computed here from `data/per_sample.tsv` alone, with no
dependencies beyond the standard library. Run it and compare with what the
documents claim; the expected values are printed alongside.

    python3 code/verify_extra.py
"""
import io
import math
import os
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")

# Constants of the pipeline, documented in METHODS.md.
ISLAND_REF_BP = 55140      # AM229678.1 in full
ISLAND_GENE_BP = 49455     # the 18 clb genes within it
ECOLI_MARKER_BP = 96300    # t__SGB10068 marker set
ECOLI_MARKER_FRAGMENTS = 122
READ_LEN = 100             # nominal; cohorts run 84-100 bp
PER_PCT = 313.0            # marker reads per 1% abundance per million reads


KMER = 119            # megahit の最後の k。contig 名の k119_ がそれ
ISLAND_BP = 54791     # AM229678.1 の長さ
REPEAT_CONTIG = "k119_32464"
REPEAT_START = 16066  # PITFALLS.md 1 番が参照だけから測った重複の開始
REPEAT_GAP = 5653     # 同じ、2 つの複製の間隔


def load(name):
    fh = io.open(os.path.join(DATA, name), encoding="utf-8")
    hd = fh.readline().rstrip("\n").split("\t")
    out = []
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) == len(hd):
            out.append(dict(zip(hd, p)))
    return out


def median(v):
    v = sorted(v)
    n = len(v)
    return v[n // 2] if n % 2 else (v[n // 2 - 1] + v[n // 2]) / 2.0


def norm_sf2(z):
    """Two-sided tail of the standard normal."""
    return math.erfc(abs(z) / math.sqrt(2.0))


# ----------------------------------------------------------------- 1. depth
def depth_difference(rows):
    print("=" * 70)
    print("1. Sequencing depth differs between groups in DE-CRC")
    print("=" * 70)
    de = [r for r in rows if r["cohort"] == "DE-CRC"]
    a = [int(r["total_reads"]) for r in de if r["group"] == "CRC"]
    b = [int(r["total_reads"]) for r in de if r["group"] == "control"]
    print("   CRC      n=%3d  median %5.1fM  range %.1f-%.1fM"
          % (len(a), median(a) / 1e6, min(a) / 1e6, max(a) / 1e6))
    print("   control  n=%3d  median %5.1fM  range %.1f-%.1fM"
          % (len(b), median(b) / 1e6, min(b) / 1e6, max(b) / 1e6))

    # Mann-Whitney U, normal approximation with tie and continuity correction
    allv = sorted(a + b)
    n = len(allv)
    rank = {}
    i = 0
    while i < n:
        j = i
        while j + 1 < n and allv[j + 1] == allv[i]:
            j += 1
        for k in range(i, j + 1):
            rank.setdefault(allv[k], (i + j + 2) / 2.0)
        i = j + 1
    n1, n2 = len(a), len(b)
    u1 = sum(rank[v] for v in a) - n1 * (n1 + 1) / 2.0
    u = min(u1, n1 * n2 - u1)
    tie = sum(t ** 3 - t for t in Counter(allv).values())
    sd = math.sqrt(n1 * n2 / 12.0 * ((n + 1) - tie / (n * (n - 1.0))))
    z = (abs(u - n1 * n2 / 2.0) - 0.5) / sd
    print("   Mann-Whitney U = %.0f   z = %.3f   P = %.3g" % (u, z, norm_sf2(z)))
    print("   PITFALLS.md claims: medians 69.5M and 34.9M, P = 1.5e-14")


# ------------------------------------------------------------------ 2. ratio
def density_ratio(rows):
    print()
    print("=" * 70)
    print("2. The marker-gene density ratio R and its ~1.9x offset")
    print("=" * 70)
    print("   R = (island reads / %d) / (E. coli marker reads / %d)"
          % (ISLAND_GENE_BP, ECOLI_MARKER_BP))
    print("   E. coli marker reads are recovered from the shipped percentage:")
    print("   marker reads = ecoli_pct x %.0f x (total_reads / 1e6)" % PER_PCT)
    print()
    de = [r for r in rows if r["cohort"] == "DE-CRC" and r["call"] == "present"]
    vals = []
    for r in sorted(de, key=lambda x: x["sample"]):
        if r["ecoli_pct"] == "NA":
            continue
        mk = float(r["ecoli_pct"]) * PER_PCT * (int(r["total_reads"]) / 1e6)
        if mk <= 0:
            continue
        rr = (int(r["island_reads"]) / ISLAND_GENE_BP) / (mk / ECOLI_MARKER_BP)
        vals.append((rr, r["sample"], r["group"]))
    vals.sort()
    for rr, s, g in vals:
        tag = "  <- E. coli cannot account for this" if rr > 10 else ""
        print("   %-24s %-8s R = %8.2f%s" % (s, g, rr, tag))
    ec = [v for v, _, _ in vals if v < 10]
    print()
    print("   The %d samples consistent with E. coli: R %.2f-%.2f, median %.2f"
          % (len(ec), min(ec), max(ec), median(ec)))
    print("   PITFALLS.md claims a cluster near 1.9 and one outlier in the hundreds.")
    print("   (The originally published range was 1.37-2.01, computed from")
    print("    MAPQ-filtered island reads, which are not shipped; island_reads is")
    print("    the pre-filter count from the bowtie2 log. See CORRECTIONS.md 14.)")
    print()
    print("   Decomposition of the offset, from constants alone:")
    f1 = ISLAND_REF_BP / float(ISLAND_GENE_BP)
    mean_frag = ECOLI_MARKER_BP / float(ECOLI_MARKER_FRAGMENTS)
    f2 = mean_frag / (mean_frag - READ_LEN + 1.0)
    prod = median(ec)
    f3 = prod / (f1 * f2)
    print("     length convention  %d / %d                  = %.4f"
          % (ISLAND_REF_BP, ISLAND_GENE_BP, f1))
    print("     marker edge loss   %.1f / (%.1f - %d + 1)   = %.4f"
          % (mean_frag, mean_frag, READ_LEN, f2))
    print("     residual           %.2f / (%.4f x %.4f)   = %.4f   (not measured"
          % (prod, f1, f2, f3))
    print("                                                       separately; it is")
    print("                                                       what is left over)")
    print("     product                                       = %.2f" % (f1 * f2 * f3))
    print("   PITFALLS.md claims 1.115 x 1.14 x ~1.50 = ~1.9")


# ------------------------------------------------------------- 3. regression
def _solve(a, b):
    n = len(a)
    m = [row[:] + [b[i]] for i, row in enumerate(a)]
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(m[r][c]))
        if abs(m[p][c]) < 1e-12:
            raise ValueError("singular")
        m[c], m[p] = m[p], m[c]
        pv = m[c][c]
        for r in range(n):
            if r == c:
                continue
            f = m[r][c] / pv
            for k in range(c, n + 1):
                m[r][k] -= f * m[c][k]
    return [m[i][n] / m[i][i] for i in range(n)]


def _logit(x, y, iters=60):
    k = len(x[0])
    beta = [0.0] * k
    for _ in range(iters):
        g = [0.0] * k
        h = [[0.0] * k for _ in range(k)]
        for xi, yi in zip(x, y):
            eta = max(-30.0, min(30.0, sum(b * v for b, v in zip(beta, xi))))
            p = 1.0 / (1.0 + math.exp(-eta))
            w = max(p * (1 - p), 1e-9)
            for a in range(k):
                g[a] += (yi - p) * xi[a]
                for c in range(k):
                    h[a][c] += w * xi[a] * xi[c]
        try:
            d = _solve(h, g)
        except ValueError:
            break
        beta = [b + dd for b, dd in zip(beta, d)]
        if max(abs(v) for v in d) < 1e-10:
            break
    inf = [[0.0] * k for _ in range(k)]
    for xi in x:
        eta = max(-30.0, min(30.0, sum(b * v for b, v in zip(beta, xi))))
        p = 1.0 / (1.0 + math.exp(-eta))
        w = max(p * (1 - p), 1e-9)
        for a in range(k):
            for c in range(k):
                inf[a][c] += w * xi[a] * xi[c]
    se = []
    for a in range(k):
        e = [1.0 if i == a else 0.0 for i in range(k)]
        se.append(math.sqrt(abs(_solve(inf, e)[a])))
    return beta, se


def regression(rows):
    print()
    print("=" * 70)
    print("3. Adjusting for depth removes the group effect in AT-CRC")
    print("=" * 70)
    print("   Logistic regression of 'called present' on group, before and after")
    print("   adding log10(total_reads). Own IRLS implementation, no statsmodels.")
    print()
    print("   %-8s %6s %7s | %-26s %-26s" %
          ("cohort", "n", "pos", "group only", "group + log10(reads)"))
    print("   " + "-" * 78)
    for coh in ("FR-CRC", "AT-CRC", "CN-CRC", "DE-CRC"):
        d = [r for r in rows if r["cohort"] == coh]
        y = [1 if r["call"] == "present" else 0 for r in d]
        grp = [1.0 if r["group"] == "CRC" else 0.0 for r in d]
        dep = [math.log10(int(r["total_reads"])) for r in d]
        mu = sum(dep) / len(dep)
        dep = [v - mu for v in dep]
        out = []
        dp = ""
        for x in ([[1.0, g] for g in grp],
                  [[1.0, g, v] for g, v in zip(grp, dep)]):
            b, se = _logit(x, y)
            out.append("OR %5.2f  P = %.4f" % (math.exp(b[1]), norm_sf2(b[1] / se[1])))
            if len(b) > 2:
                dp = "depth P = %.4f" % norm_sf2(b[2] / se[2])
        print("   %-8s %6d %7d | %-26s %-26s  %s" % (coh, len(d), sum(y), out[0], out[1], dp))
    print()
    print("   PITFALLS.md claims, for AT-CRC: OR 3.11 -> 1.78, P 0.080 -> 0.40,")
    print("   with depth itself significant at P = 0.010.")


# ---------------------------------------------------------------- 4. carrier
def carriers():
    print()
    print("=" * 70)
    print("4. The carrier table: 14 positives against six complete genomes")
    print("=" * 70)
    print("   ratio = island mean depth / that species' genome mean depth,")
    print("   both from the same read set at MAPQ 20.")
    print()
    rows = load("carrier_six_genomes.tsv")
    by = {}
    for r in rows:
        by.setdefault(r["sample"], {})[r["target"]] = r
    print("   %-22s %-8s %8s | %9s %6s %7s | %9s %6s %9s"
          % ("sample", "group", "island", "E.coli", "bread", "ratio",
             "K.pneu", "bread", "ratio"))
    print("   " + "-" * 94)

    def stat(d):
        bp = float(d["target_bp"])
        return int(d["sum_depth"]) / bp, int(d["covered_bp"]) / bp

    out = []
    for s in by:
        isl, _ = stat(by[s]["pks_island_clb_genes"])
        ed, eb = stat(by[s]["Ecoli_IHE3034"])
        kd, kb = stat(by[s]["Kpneumoniae_HS11286"])
        out.append((isl / ed if ed > 0.01 else float("inf"), s, isl, ed, eb, kd, kb))
    for er, s, isl, ed, eb, kd, kb in sorted(out):
        kr = ("%9.2f" % (isl / kd)) if kd > 0.01 else "        -"
        print("   %-22s %-8s %8.2f | %9.2f %6.3f %7.2f | %9.2f %6.3f %s"
              % (s, "CRC" if s.startswith("cancer_") else "control",
                 isl, ed, eb, er, kd, kb, kr))
    print()
    print("   Compare with the table in RESULTS.md. Every species measured is in")
    print("   data/carrier_six_genomes.tsv, including the four not shown here.")


# --------------------------------------------------------------- 5. assembly
def assembly():
    print()
    print("=" * 70)
    print("5. The one carrier confirmed by assembly, and whether it is a chimera")
    print("=" * 70)
    for r in load("assembly_alignment.tsv"):
        print("   %s %s-%s (%s) -> %-22s %s-%s   %s/%s matched   [%s]"
              % (r["contig"], r["contig_start"], r["contig_end"], r["strand"],
                 r["target"], r["target_start"], r["target_end"],
                 r["matches"], r["aln_len"], r["method"]))
    print()
    print("   One contig of 386,079 bp: K. pneumoniae chromosome up to ~354,276,")
    print("   then the pks island to the end. They are contiguous, not two contigs.")
    print()
    d = load("assembly_junction_depth.tsv")
    pos = [int(r["position"]) for r in d]
    dep = [int(r["depth"]) for r in d]
    zero = sum(1 for v in dep if v == 0)
    print("   Read depth across %d-%d (%d bases spanning the junction):"
          % (min(pos), max(pos), len(dep)))
    print("     minimum %d   median %d   maximum %d   bases at zero depth: %d"
          % (min(dep), sorted(dep)[len(dep) // 2], max(dep), zero))
    if zero == 0:
        print("   No gap in coverage, so the join is supported by reads rather than")
        print("   inferred by the assembler. That is what rules out a chimera.")


# -------------------------------------------------------------- 6. kleborate
def _kleborate(name):
    """Read Kleborate's one-row output into a {column: value} dict."""
    fh = io.open(os.path.join(DATA, "kleborate", name), encoding="utf-8")
    hd = fh.readline().rstrip("\n").split("\t")
    val = fh.readline().rstrip("\n").split("\t")
    return dict(zip(hd, val))


def typing():
    print()
    print("=" * 70)
    print("6. The strain typing, and which run each part of it comes from")
    print("=" * 70)
    keys = [
        ("enterobacterales__species__species", "species"),
        ("general__contig_stats__contig_count", "contigs"),
        ("general__contig_stats__total_size", "total bp"),
        ("general__contig_stats__largest_contig", "largest contig"),
        ("klebsiella_pneumo_complex__mlst__ST", "MLST"),
        ("klebsiella_pneumo_complex__kaptive__K_type", "K type"),
        ("klebsiella__ybst__Yersiniabactin", "yersiniabactin"),
        ("klebsiella__cbst__CbST", "CbST"),
        ("klebsiella__cbst__Colibactin", "colibactin"),
    ]
    a, b = _kleborate("kp_bin_output.tsv"), _kleborate("kp_plus_output.tsv")
    print("   %-15s | %-22s | %-24s"
          % ("", "kp_bin (automatic)", "kp_plus (+island contigs)"))
    print("   %-15s-+-%-22s-+-%-24s" % ("-" * 15, "-" * 22, "-" * 24))
    for k, lab in keys:
        x, y = a.get(k, "?"), b.get(k, "?")
        print("   %-15s | %-22s | %-24s%s"
              % (lab, x, y, "" if x == y else "   <- differs"))
    print()
    print("   RESULTS.md quotes ST23/K1/ICEKp10/clb for this sample. The species,")
    print("   the ST and the K type come from the automatic bin; ICEKp10 and a")
    print("   complete clb 2 appear only once the island-bearing contigs are")
    print("   added by hand. See CORRECTIONS.md item 16.")
    print()

    mem = load("kleborate/bin_membership.tsv")
    added = [r for r in mem if r["in_kp_plus"] == "yes" and r["in_kp_bin"] == "no"]
    isl = [r for r in mem if r["carries_island"] == "yes"]
    print("   Bin membership: %d contigs in kp_bin, %d in kp_plus, %d carry island."
          % (sum(1 for r in mem if r["in_kp_bin"] == "yes"),
             sum(1 for r in mem if r["in_kp_plus"] == "yes"), len(isl)))
    ok = all(r["carries_island"] == "yes" for r in added)
    print("   The %d added contigs are %sall island-bearing, so kp_plus is exactly"
          % (len(added), "" if ok else "NOT "))
    print("   kp_bin plus every contig carrying any part of the island.")
    asg = {r["contig"]: r["best_species"]
           for r in load("kleborate/contig_assignment.tsv")}
    for r in added:
        print("     %-14s %7s bp   minimap2 assigned it to %s"
              % (r["contig"], r["contig_bp"], asg.get(r["contig"], "?")))
    print()
    print("   Why none of them could be assigned to Klebsiella:")
    for r in load("clb_in_six_genomes.tsv"):
        n, s = float(r["island_kmers"]), float(r["shared_kmers"])
        print("     %-22s %6.2f%% of the island's 40-mers, %s/%s clb genes"
              % (r["genome"], 100.0 * s / n, r["clb_genes_at_half"],
                 r["clb_genes_total"]))
    print("   K. pneumoniae HS11286 carries no clb island, so a clb-bearing contig")
    print("   can only ever be assigned to E. coli or C. koseri. Adding them was")
    print("   forced by the reference set rather than chosen for this sample - but")
    print("   the typing is then not independent evidence that those contigs are")
    print("   Klebsiella's. The independent evidence is section 5: a single contig")
    print("   carries the chromosome and the island with no gap in read depth.")


def backbone():
    print()
    print("=" * 70)
    print("7. Do the island contigs belong to this Klebsiella genome?")
    print("=" * 70)
    rows = load("kleborate/island_contig_backbone.tsv")
    print("   Alignments split at the island boundary. The last column is matched")
    print("   bases over aligned length from minimap2's PAF, which **includes gaps**")
    print("   and is therefore not an identity: a contig at 99% identity can read as")
    print("   0.81 here. Use it to compare species against each other, nothing more.")
    print("   %-14s %-9s %-22s %9s %9s  %s"
          % ("contig", "region", "species", "matched", "aligned", "mat/aln"))
    for r in rows:
        al = float(r["aligned_bp"])
        print("   %-14s %-9s %-22s %9s %9s  %6.1f%%"
              % (r["contig"], r["region"], r["species"], r["matched_bp"],
                 r["aligned_bp"], 100.0 * int(r["matched_bp"]) / al if al else 0.0))
    print()
    print("   Note k119_70945: its backbone matches C. koseri better than")
    print("   K. pneumoniae, so similarity does NOT place it in Klebsiella.")
    print("   Our working notes claimed the opposite. See CORRECTIONS.md item 16.")
    print()
    d = {r["what"]: r["value"] for r in load("kleborate/contig_depth.tsv")}
    mean = float(d["bin_depth_length_weighted_mean"])
    print("   Coverage depth instead. The bin: %s contigs, %s bp,"
          % (d["bin_contigs"], d["bin_bp"]))
    print("   length-weighted mean depth %.2f, median %s over %s contigs >=10 kb."
          % (mean, d["bin_depth_median_contigs_over_10kb"],
             d["bin_contigs_over_10kb"]))
    seen = []
    for k in d:
        if k.startswith("island_contig_") and k.endswith("_depth"):
            c = k[len("island_contig_"):-len("_depth")]
            seen.append((int(d["island_contig_%s_bp" % c]), c, float(d[k])))
    for bp, c, dep in sorted(seen, reverse=True):
        rel = dep / mean
        note = ""
        if abs(rel - 2.0) < 0.1:
            note = "   <- twice the genome: a collapsed two-copy repeat"
        elif abs(rel - 1.0) < 0.15:
            note = "   <- at the genome's own depth"
        elif rel < 0.7:
            note = "   <- does NOT fit the genome's depth"
        print("     %-14s %8d bp  depth %6.2f  = %.2fx the bin%s"
              % (c, bp, dep, rel, note))
    print()
    print("   C. koseri is not in this sample to be their source: its genome")
    print("   breadth never exceeds 0.028, and the clb island it carries is")
    print("   itself 1.1% of its genome.")
    print()
    print("   **Those depths are megahit `multi` values, not read depths.**")
    print("   This assembly was built to k=119 and the reads are 100 bp, so")
    print("   `multi` is a multiplicity in the de Bruijn graph of reads plus")
    print("   the contigs carried forward from smaller k. It is not")
    print("   proportional to read depth and we never calibrated it. What it")
    print("   does show is one contig at twice its neighbours, which is what a")
    print("   collapsed two-copy repeat looks like. Sections 8 and 9 replace")
    print("   this argument with evidence that is not depth at all.")


def island_tiling():
    print()
    print("=" * 70)
    print("8. The hand-added contigs are not four separate objects")
    print("=" * 70)
    print("   Section 7 argued about them one at a time, which is why only")
    print("   coverage depth was left. They are not independent: laid out on")
    print("   the island reference AM229678.1 they **overlap each other**.")
    print()
    rows = load("kleborate/island_tiling.tsv")
    print("   %-5s %-12s %-6s %10s %10s %9s %8s %8s  %s"
          % ("order", "contig", "in bin", "island", "to", "identity",
             "coords", "sequence", ""))
    for r in rows:
        ov, sq = r["overlap_prev"], r["overlap_seq_bp"]
        note = ""
        if sq != "NA" and int(sq) == KMER:
            note = "   <- exactly k = %d" % KMER
        elif ov != "NA" and int(ov) < 0:
            note = "   <- a gap inside one contig, not a junction"
        print("   %-5s %-12s %-6s %10s %10s %8s%% %8s %8s%s"
              % (r["order"], r["contig"], r["in_kp_bin"], r["island_start"],
                 r["island_end"], r["identity"], ov, sq, note))
    print()
    print("   The last column is not a coordinate difference. It is the length")
    print("   of identical sequence shared between the end of one contig and")
    print("   the start of the next, taken from the contigs themselves.")
    print()

    prev = None
    checked = 0
    for r in rows:
        if prev is not None and r["contig"] != prev:
            assert int(r["overlap_prev"]) == KMER, \
                "junction at %s: coordinate overlap %s, expected %d" % (
                    r["contig"], r["overlap_prev"], KMER)
            assert int(r["overlap_seq_bp"]) == KMER, \
                "junction at %s: shared sequence %s bp, expected %d" % (
                    r["contig"], r["overlap_seq_bp"], KMER)
            checked += 1
        prev = r["contig"]
    print("   All %d junctions between different contigs are overlaps of" % checked)
    print("   exactly %d bp - the k-mer length this assembly was built at." % KMER)
    print("   Contigs sharing exactly k bases are consecutive nodes of one path:")
    print("   megahit did not assemble four fragments. It assembled one island")
    print("   and broke it twice.")
    print()

    rep = [r for r in rows if r["contig"] == REPEAT_CONTIG]
    assert len(rep) == 2, "the repeat contig should appear twice"
    a, b = int(rep[0]["island_start"]), int(rep[1]["island_start"])
    print("   Where it broke: `%s`, which appears **twice** - island %d and"
          % (REPEAT_CONTIG, a))
    print("   island %d, a spacing of %d bp." % (b, b - a))
    print("   PITFALLS.md item 1 measured that duplication independently, from")
    print("   the reference alone: 1,287 bp starting at %d, the two copies %d"
          % (REPEAT_START, REPEAT_GAP))
    print("   bp apart, checked against ENA by verify_repeat.py.")
    assert a == REPEAT_START, "repeat starts at %d, PITFALLS says %d" % (
        a, REPEAT_START)
    assert b - a == REPEAT_GAP, "spacing %d, PITFALLS says %d" % (
        b - a, REPEAT_GAP)
    print("   Same start, same spacing - and in section 7 that is the contig at")
    print("   twice the depth. The assembly broke the island at the duplication,")
    print("   which is what a collapsed two-copy repeat does.")
    print()

    cov = set()
    for r in rows:
        cov |= set(range(int(r["island_start"]), int(r["island_end"]) + 1))
    print("   The four contigs cover %d of AM229678.1's %d bp (%.1f%%), from %d"
          % (len(cov), ISLAND_BP, 100.0 * len(cov) / ISLAND_BP, min(cov)))
    print("   to %d, at 99.76-100%% identity." % max(cov))
    assert len(cov) / ISLAND_BP > 0.95, "the tiling should cover >95%"
    inbin = sorted(set(r["contig"] for r in rows if r["in_kp_bin"] == "yes"))
    assert inbin, "no tiling contig is in the automatic bin"
    print("   One of them, %s, is in the automatic bin, and section 5 shows it"
          % inbin[0])
    print("   carrying 350 kb of Klebsiella chromosome into the island with no")
    print("   drop in read depth across the join. **The other three reach the")
    print("   chromosome through that contig, by overlap - not by depth.**")
    print()
    print("   What this does not settle: both breaks are at the repeat, and a")
    print("   repeat is the one place where graph adjacency can join paths from")
    print("   different sources. So it remains possible in principle that")
    print("   k119_70945 belongs to a second copy of the island in this sample.")
    print("   That is where the absence of the alternative is still doing work -")
    print("   but it now carries one narrow question instead of the whole claim.")


def island_locus():
    print()
    print("=" * 70)
    print("9. The element sits at one site in the Klebsiella chromosome")
    print("=" * 70)
    rows = load("kleborate/island_locus.tsv")
    el = [r for r in rows if r["region"] == "element"]
    kp = [r for r in rows if r["region"] == "kp_chromosome"]
    lo = min(int(r["ruler_start"]) for r in el)
    hi = max(int(r["ruler_end"]) for r in el)
    n = len(set(r["contig"] for r in el))
    print("   The island is part of a larger insertion. Measured on the only one")
    print("   of the six references that carries the whole element - C. koseri")
    print("   BAA-895, used **as a ruler, not as a claim about where it came")
    print("   from** - the element runs %d-%d (%d bp) and is tiled by %d contigs."
          % (lo, hi, hi - lo, n))
    print()
    print("   %-12s %-6s %-6s %-14s %10s %10s %9s"
          % ("contig", "bin", "added", "region", "start", "end", "aln bp"))
    for r in sorted(rows, key=lambda x: (x["region"], int(x["ruler_start"]))):
        print("   %-12s %-6s %-6s %-14s %10s %10s %9s"
              % (r["contig"], r["in_kp_bin"], r["in_kp_plus"], r["region"],
                 r["ruler_start"], r["ruler_end"], r["aln_bp"]))
    print()
    ends = {}
    for r in el:
        c = r["contig"]
        s, e = int(r["ruler_start"]), int(r["ruler_end"])
        o = ends.get(c, (s, e))
        ends[c] = (min(o[0], s), max(o[1], e))
    left = [c for c, (s, e) in ends.items() if s == lo][0]
    right = [c for c, (s, e) in ends.items() if e == hi][0]
    print("   The two ends of the element are reached by %s and %s." % (left, right))
    for c in (left, right):
        assert any(r["contig"] == c and r["in_kp_bin"] == "yes" for r in rows), \
            "%s is not in the automatic bin" % c
    print("   **Both are in the automatic bin.** Neither was added by hand.")
    print()
    print("   Their alignments to the Klebsiella chromosome:")
    for r in sorted(kp, key=lambda x: int(x["ruler_start"])):
        print("     %-12s %9s-%-9s %8s bp"
              % (r["contig"], r["ruler_start"], r["ruler_end"], r["aln_bp"]))
    ka = [int(r["ruler_end"]) for r in kp if r["contig"] == right]
    kb = [int(r["ruler_start"]) for r in kp if r["contig"] == left]
    joint = min(abs(x - y) for x in ka for y in kb)
    print()
    print("   %s leaves the chromosome at %d; %s picks it up at %d. The two"
          % (right, max(ka), left, min(kb)))
    print("   intervals meet within %d bp. **One insertion, one site, both" % joint)
    print("   flanks in the automatic bin**, and the hand-added contigs fill the")
    print("   interior. That is position, not depth and not similarity.")
    assert joint < 100, "the two flanks do not meet: %d bp apart" % joint
    print()
    print("   k119_70945 is worth a line of its own: its backbone also aligns to")
    print("   the Klebsiella chromosome, at 3,434,608-3,479,783, which falls")
    print("   inside a stretch of HS11286 that k119_25144 skips. HS11286 carries")
    print("   a related element at a different site. That is exactly why this")
    print("   contig's similarity is ambiguous between the two genera, and why")
    print("   item 16's similarity argument failed on it.")
    added = set(r["contig"] for r in el if r["in_kp_plus"] == "yes")
    print()
    print("   %d of the %d tiling contigs were never added to the bin: the"
          % (n - len(added), n))
    print("   element is larger than the island, and the bin is not the whole")
    print("   genome. Nothing in the typing depends on those.")


def read_depth_and_links():
    print()
    print("=" * 70)
    print("10. Read depth, measured, and read-pair linkage")
    print("=" * 70)
    print("   Section 7's depths were megahit `multi`. These are read depths from")
    print("   one mapping of the reads back to the whole assembly (21,759 contigs,")
    print("   one competitive index, MAPQ 20) - pipeline/asm_linkage.sh.")
    print()
    rows = load("kleborate/contig_read_depth.tsv")
    byc = dict((r["contig"], r) for r in rows)
    print("   %-12s %9s %8s %9s %9s %9s %10s"
          % ("contig", "bp", "multi", "depth", "breadth", "/bin", "multi/depth"))
    for c in ("kp_bin_all", "k119_45171", "k119_25144", "k119_70945",
              "k119_70596", "k119_25170", "k119_32464"):
        r = byc[c]
        add = "" if r["in_kp_bin"] == "yes" else "   <- added by hand"
        print("   %-12s %9s %8s %9s %9s %9s %10s%s"
              % (c, r["contig_bp"], r["megahit_multi"], r["read_depth"],
                 r["breadth"], r["depth_over_bin"], r["multi_over_depth"], add))
    print()
    big = [r for r in rows if r["in_kp_bin"] == "yes"
           and r["contig"] != "kp_bin_all" and int(r["contig_bp"]) >= 10000]
    lo = min(float(r["depth_over_bin"]) for r in big)
    hi = max(float(r["depth_over_bin"]) for r in big)
    print("   The bin's own %d contigs over 10 kb span %.2f-%.2f times the bin's"
          % (len(big), lo, hi))
    print("   depth. Three of the four added contigs are at %s, %s and %s;"
          % tuple(byc[c]["depth_over_bin"]
                  for c in ("k119_70945", "k119_70596", "k119_25170")))
    print("   all inside that range. The fourth, k119_32464, is at %s - the"
          % byc["k119_32464"]["depth_over_bin"])
    print("   collapsed two-copy repeat, which has to be near 2.")
    for c in ("k119_70945", "k119_70596", "k119_25170"):
        v = float(byc[c]["depth_over_bin"])
        assert lo <= v <= hi, "%s at %.3f is outside the bin's own range" % (c, v)
    assert 1.7 < float(byc["k119_32464"]["depth_over_bin"]) < 2.3, \
        "the collapsed repeat should sit near twice the bin"
    print()
    print("   **And this is where `multi` misled us.** multi/depth is %s-%s for"
          % (min(byc[c]["multi_over_depth"] for c in
                 ("k119_45171", "k119_25144", "k119_70945", "k119_70596",
                  "k119_32464")),
             max(byc[c]["multi_over_depth"] for c in
                 ("k119_45171", "k119_25144", "k119_70945", "k119_70596",
                  "k119_32464"))))
    print("   five of the six contigs and %s for k119_25170. In read depth"
          % byc["k119_25170"]["multi_over_depth"])
    print("   k119_25170 is at %s of the bin, not the half that section 7's"
          % byc["k119_25170"]["depth_over_bin"])
    print("   `multi` column reported. See the amendment to CORRECTIONS item 16:")
    print("   the reason to doubt that contig is its 137 bp of island, not depth.")
    print()
    print("   Read-pair linkage. Pairs whose two mates land on different contigs,")
    print("   both at MAPQ 20. Partners are counted only when they are one of the")
    print("   160 contigs the six-genome step assigned at breadth >= 0.5; the other")
    print("   21,599 are other gut bacteria and pairs reaching them are shared")
    print("   repeats.")
    lk = load("kleborate/contig_links.tsv")
    bylk = dict((r["contig"], r) for r in lk)
    print()
    print("   %-12s %9s %7s %7s %7s %7s %9s %9s"
          % ("contig", "bp", "bin", "island", "other", "total", "share",
             "bin/kb"))
    for c in ("k119_45171", "k119_25144", "k119_70945", "k119_70596",
              "k119_25170", "k119_32464"):
        r = bylk[c]
        add = "" if r["in_kp_bin"] == "yes" else "   <-"
        print("   %-12s %9s %7s %7s %7s %7s %9.3f %9.2f%s"
              % (c, r["contig_bp"], r["pairs_to_kp_bin"],
                 r["pairs_to_added_island"], r["pairs_to_other_species"],
                 r["pairs_assigned_total"],
                 float(r["share_kp_bin_or_island"]),
                 float(r["pairs_to_kp_bin_per_kb"]), add))
    print()

    def band(lo_bp, hi_bp):
        v = [float(r["share_kp_bin_or_island"]) for r in lk
             if r["in_kp_bin"] == "yes" and lo_bp <= int(r["contig_bp"]) < hi_bp
             and int(r["pairs_assigned_total"]) > 0]
        return (len(v), min(v), max(v)) if v else (0, 0, 0)

    for a, b, nm in ((10000, 10 ** 9, "over 10 kb"), (2000, 10000, "2-10 kb"),
                     (0, 2000, "under 2 kb")):
        n, x, y = band(a, b)
        print("   bin contigs %-11s %2d  share %.3f-%.3f" % (nm, n, x, y))
    print()
    n, x, y = band(10000, 10 ** 9)
    v = float(bylk["k119_70945"]["share_kp_bin_or_island"])
    print("   **k119_70945, at %.3f, is inside the range the bin's own contigs" % v)
    print("   of that length occupy (%.3f-%.3f).** It is 63.9 kb, the only added"
          % (x, y))
    print("   contig long enough for this test to discriminate, and it does not")
    print("   behave differently from a bin contig. Its mates reach %s pairs onto"
          % bylk["k119_70945"]["pairs_to_kp_bin"])
    print("   bin contigs. That is the question section 8 left open.")
    assert x <= v <= y, "k119_70945 at %.3f is outside %.3f-%.3f" % (v, x, y)
    print()
    print("   **For the other three the test has no power and we say so.** Below")
    print("   10 kb the bin's own contigs span %.3f-1.000, which covers every"
          % band(2000, 10000)[1])
    print("   value an added contig could take. k119_70596 (1.000) and")
    print("   k119_25170 (%.3f) are both inside it. One concrete fact survives:"
          % float(bylk["k119_25170"]["share_kp_bin_or_island"]))
    print("   k119_32464 sends %s pairs to k119_45171, and k119_45171's single"
          % "34")
    print("   largest assigned partner is k119_32464 - a hand-added contig ahead")
    print("   of every contig in the bin.")


def ratio_scale():
    print()
    print("=" * 70)
    print("11. The scale of the island/species depth ratio, measured")
    print("=" * 70)
    print("   CORRECTIONS item 8 read this ratio as a population fraction")
    print("   (\"0.85, i.e. ~85% of E. coli cells\"), which is item 3's mistake")
    print("   on a different numerator. Item 3 measured R's offset; this")
    print("   measures this ratio's. The selection rule was fixed first and")
    print("   puts NO condition on the ratio - conditioning on it would pick")
    print("   values near 1 and make the calibration circular.")
    print()
    rows = load("island_ratio_scale.tsv")
    cal = [r for r in rows if r["calibrator_for"] == "Ecoli_IHE3034"]
    kcal = [r for r in rows if r["calibrator_for"] == "Kpneumoniae_HS11286"]
    print("   %-22s %9s %9s %9s %9s"
          % ("calibrator for E. coli", "island", "E. coli", "isl/Ec", "rep.bias"))
    for r in sorted(cal, key=lambda x: float(x["island_over_ecoli"])):
        off = (float(r["island_depth_covered_only"]) / float(r["island_depth"]))
        print("   %-22s %9.3f %9.3f %9.3f %9.4f"
              % (r["sample"], float(r["island_depth"]),
                 float(r["Ecoli_IHE3034"]), float(r["island_over_ecoli"]), off))
    v = sorted(float(r["island_over_ecoli"]) for r in cal)
    offs = [float(r["island_depth_covered_only"]) / float(r["island_depth"])
            for r in cal]
    print()
    print("   %d samples calibrate E. coli; the band is %.3f-%.3f, a spread of"
          % (len(cal), v[0], v[-1]))
    print("   %.1fx. The ratio never passes %.3f where E. coli must be the"
          % (v[-1] / v[0], v[-1]))
    print("   carrier, so the scale is near 1 - R's offset was 1.9.")
    print("   The repeat interior biases the island's depth down by %.3f-%.3fx,"
          % (min(offs), max(offs)))
    print("   so the one-copy-in-every-cell ceiling sits near %.2f."
          % (v[-1] * sum(offs) / len(offs)))
    assert len(cal) == 4, "expected 4 E. coli calibrators, got %d" % len(cal)
    assert v[-1] < 1.05, "a calibrator above 1.05 would mean a real offset"
    assert v[-1] / v[0] > 3.0, \
        "a narrow band would mean the selection had been circular"
    for o in offs:
        assert 1.03 < o < 1.06, "repeat bias out of the measured range: %.4f" % o
    print()
    print("   **There is no calibrator for Klebsiella: %d samples.** The scale"
          % len(kcal))
    print("   of this ratio is measured for E. coli only, so reading a")
    print("   Klebsiella ratio as copies per cell - the 1.30 in the item 16")
    print("   amendment included - assumes the two genomes behave alike here.")
    print("   Item 16 does not rest on it; its depth evidence was replaced.")
    assert not kcal, "a Klebsiella calibrator would change what item 8 can say"
    print()
    e = [r for r in rows if r["sample"] == "cancer_ERR1018282"][0]
    rk = float(e["island_over_kp"])
    print("   cancer_ERR1018282: island %.3f, E. coli %.3f (ratio %.3f),"
          % (float(e["island_depth"]), float(e["Ecoli_IHE3034"]),
             float(e["island_over_ecoli"])))
    print("   Kp %.3f (ratio %.4f = 1 cell in %.1f)."
          % (float(e["Kpneumoniae_HS11286"]), rk, 1.0 / rk))
    print("   Both fit: ~%.0f%% of a %.1fx E. coli population at one copy, or"
          % (100 * float(e["island_over_ecoli"]), float(e["Ecoli_IHE3034"])))
    print("   %.1f%% of a %.0fx Klebsiella population at one copy. **The ratio"
          % (100 * rk, float(e["Kpneumoniae_HS11286"])))
    print("   excludes neither.** 1 in %.0f is an ordinary frequency for a" % (1 / rk))
    print("   mobile element, not an exclusion, and item 8 said \"1 in 33\" by")
    print("   inverting the rounded 0.03.")
    assert 1 / rk > 33.5, \
        "if 1/ratio really were 33 the published figure would have been right"
    assert float(e["island_over_ecoli"]) <= v[-1], \
        "282's E. coli ratio should sit inside the calibrated band"


def main():
    rows = load("per_sample.tsv")
    depth_difference(rows)
    density_ratio(rows)
    regression(rows)
    carriers()
    assembly()
    typing()
    backbone()
    island_tiling()
    island_locus()
    read_depth_and_links()
    ratio_scale()


if __name__ == "__main__":
    main()

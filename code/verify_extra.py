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
    print("   itself 1.1% of its genome. So the argument is coverage depth plus")
    print("   the absence of the alternative - weaker than section 5's single")
    print("   contig, and not presented as equivalent to it.")


def main():
    rows = load("per_sample.tsv")
    depth_difference(rows)
    density_ratio(rows)
    regression(rows)
    carriers()
    assembly()
    typing()
    backbone()


if __name__ == "__main__":
    main()

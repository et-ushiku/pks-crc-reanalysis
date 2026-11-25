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


def main():
    rows = load("per_sample.tsv")
    depth_difference(rows)
    density_ratio(rows)
    regression(rows)


if __name__ == "__main__":
    main()

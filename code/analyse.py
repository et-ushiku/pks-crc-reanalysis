"""公開データだけで、主結果を最初から再現する。

このスクリプトは data/depth_hist.tsv と data/per_sample.tsv しか読まない。
fastq も bam も要らない。標準ライブラリだけで動く（scipy 不要）。

事前に決めた判定（結果を見てから一度も変えていない）
  MIN_DEPTH     = 1     深さ 1 以上で「覆われた」とする
  GENE_COVER    = 0.5   遺伝子の半分以上が覆われたら「その遺伝子はある」
  CALL_BREADTH  = 0.5   島の覆いが 0.5 以上
  CALL_GENES    = 14    かつ 18 遺伝子のうち 14 以上  → 「あり」
  SUSPECT_BREADTH = 0.1 覆い 0.1 以上 0.5 未満        → 「疑い」
  ECOLI_MIN     = 0.1   大腸菌が 0.1% 未満の検体は「判定不能」として外す
  TARGET        = 2000万読み取りに揃えてから判定する（主解析）

間引きの式（fastq が無くてもできる理由）
  位置 i の深さを d_i、残す割合を f とすると
      その位置が覆われたまま残る確率 = 1 - (1-f)^d_i
      覆いの期待値 = 平均_i [ 1 - (1-f)^d_i ]
  期待値の線形性から、位置どうしの相関によらず厳密。
  深さの度数分布さえあれば同じ値になるので、度数分布を配っている。

使い方: python3 analyse.py [目標の読み取り数]   既定 20000000
"""
import math
import os
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")

MIN_DEPTH, GENE_COVER, CALL_BREADTH, CALL_GENES, SUSPECT_BREADTH = 1, 0.5, 0.5, 14, 0.1
ECOLI_MIN = 0.1
ORDER = ["FR-CRC", "AT-CRC", "CN-CRC", "DE-CRC"]


def load_bed():
    reg = []
    for line in open(os.path.join(DATA, "clb_genes.bed"), encoding="utf-8"):
        q = line.split()
        if len(q) >= 4:
            reg.append((q[3], int(q[1]), int(q[2])))
    return reg


def load_hist():
    h = defaultdict(lambda: defaultdict(dict))
    fh = open(os.path.join(DATA, "depth_hist.tsv"), encoding="utf-8")
    fh.readline()
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) < 4:
            continue
        h[p[0]][p[1]][int(p[2])] = int(p[3])
    return h


def load_samples():
    rows = []
    fh = open(os.path.join(DATA, "per_sample.tsv"), encoding="utf-8")
    hd = fh.readline().rstrip("\n").split("\t")
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) != len(hd):
            continue
        rows.append(dict(zip(hd, p)))
    return rows


def call_from(h, reg, f):
    tl = tc = ng = 0
    for nm, s, e in reg:
        g = h.get(nm, {})
        ln = e - s
        if f >= 1.0:
            c = sum(n for d, n in g.items() if d >= MIN_DEPTH)
        else:
            c = sum(n * (1.0 - (1.0 - f) ** d) for d, n in g.items())
        tl += ln
        tc += c
        if c / ln >= GENE_COVER:
            ng += 1
    br = tc / tl if tl else 0.0
    cl = "present" if (br >= CALL_BREADTH and ng >= CALL_GENES) \
        else ("suspect" if br >= SUSPECT_BREADTH else "absent")
    return br, ng, cl


def fisher(a, b, c, d):
    def lc(n, k):
        return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)
    n = a + b + c + d
    r1, r2, c1 = a + b, c + d, a + c
    base = lc(n, c1)
    obs = lc(r1, a) + lc(r2, c) - base
    tot = 0.0
    for i in range(max(0, c1 - r2), min(r1, c1) + 1):
        p = lc(r1, i) + lc(r2, c1 - i) - base
        if p <= obs + 1e-9:
            tot += math.exp(p)
    return min(1.0, tot)


def norm_sf(z):
    return 0.5 * math.erfc(z / math.sqrt(2))


def meta(tbl):
    """Mantel-Haenszel（層＝コホート）と Robins-Breslow-Greenland の分散"""
    num = den = 0.0
    for _, a, b, c, d in tbl:
        n = a + b + c + d
        if n:
            num += a * d / n
            den += b * c / n
    if den == 0:
        return None
    mh = num / den
    s = sp = sq = 0.0
    for _, a, b, c, d in tbl:
        n = a + b + c + d
        if not n:
            continue
        P, Q = (a + d) / n, (b + c) / n
        R, S = a * d / n, b * c / n
        s += P * R
        sp += P * S + Q * R
        sq += Q * S
    var = s / (2 * num * num) + sp / (2 * num * den) + sq / (2 * den * den)
    se = math.sqrt(var)
    z = abs(math.log(mh)) / se
    # 不均一性 I2
    ws, ls = [], []
    for _, a, b, c, d in tbl:
        aa, bb, cc, dd = a + .5, b + .5, c + .5, d + .5
        ls.append(math.log(aa * dd / (bb * cc)))
        ws.append(1.0 / (1 / aa + 1 / bb + 1 / cc + 1 / dd))
    W = sum(ws)
    mu = sum(w * l for w, l in zip(ws, ls)) / W
    Q2 = sum(w * (l - mu) ** 2 for w, l in zip(ws, ls))
    k = len(ls)
    I2 = max(0.0, (Q2 - (k - 1)) / Q2 * 100) if Q2 > 0 else 0.0
    return (mh, math.exp(math.log(mh) - 1.96 * se),
            math.exp(math.log(mh) + 1.96 * se), 2 * norm_sf(z), Q2, I2)


def table(rows, key):
    out = []
    for coh in ORDER:
        r = [x for x in rows if x["cohort"] == coh]
        a = sum(1 for x in r if x["group"] == "CRC" and x[key] == "present")
        na = sum(1 for x in r if x["group"] == "CRC")
        c = sum(1 for x in r if x["group"] == "control" and x[key] == "present")
        nc = sum(1 for x in r if x["group"] == "control")
        if na and nc:
            out.append((coh, a, na - a, c, nc - c))
    return out


def show(tbl, title):
    print("=" * 70)
    print(title)
    print("=" * 70)
    print("%-10s %12s %12s %9s %8s" % ("cohort", "CRC", "control", "OR", "P"))
    for coh, a, b, c, d in tbl:
        orr = ((a + .5) * (d + .5)) / ((b + .5) * (c + .5))
        print("%-10s %5d/%-6d %5d/%-6d %9.2f %8.4f"
              % (coh, a, a + b, c, c + d, orr, fisher(a, b, c, d)))
    m = meta(tbl)
    if m:
        print("\nMantel-Haenszel  OR %.2f  95%% CI %.2f-%.2f  P = %.4f"
              % (m[0], m[1], m[2], m[3]))
        print("heterogeneity    Q = %.2f  I2 = %.0f%%" % (m[4], m[5]))
        ps = [fisher(a, b, c, d) for _, a, b, c, d in tbl]
        print("smallest cohort P = %.4f  ->  Bonferroni x%d = %.4f"
              % (min(ps), len(ps), min(1.0, min(ps) * len(ps))))
    print()


def main():
    M = int(sys.argv[1]) if len(sys.argv) > 1 else 20000000
    reg = load_bed()
    h = load_hist()
    rows = load_samples()
    print("pks island carriage and colorectal cancer -- reanalysis of 4 public cohorts")
    print("samples: %d   target depth: %dM reads\n" % (len(rows), M // 1000000))

    # 1) 配った度数分布から、表の覆いを作り直せるか確かめる
    worst = 0.0
    for r in rows:
        br, ng, cl = call_from(h[r["sample"]], reg, 1.0)
        worst = max(worst, abs(br - float(r["island_breadth"])))
        r["call_raw"] = cl
        assert cl == r["call"], (r["sample"], cl, r["call"])
    print("check: breadth recomputed from depth_hist.tsv matches per_sample.tsv")
    print("       largest difference over %d samples = %.2e\n" % (len(rows), worst))

    # 2) 共通の読み取り数に間引く
    for r in rows:
        N = r["total_reads"]
        if N == "NA" or int(N) < M:
            r["call_rar"] = "excluded"
            continue
        r["call_rar"] = call_from(h[r["sample"]], reg, M / int(N))[2]

    # 3) 事前に決めた門（大腸菌 >= 0.1%）
    def gated(r):
        return r["ecoli_pct"] != "NA" and float(r["ecoli_pct"]) >= ECOLI_MIN

    prim = [r for r in rows if r["call_rar"] != "excluded" and gated(r)]
    print("judgeable after rarefying and gating: %d of %d\n" % (len(prim), len(rows)))

    show(table(prim, "call_rar"),
         "PRIMARY  rarefied to %dM reads, E. coli >= %.1f%% (both pre-specified)"
         % (M // 1000000, ECOLI_MIN))
    show(table([r for r in rows if r["call_rar"] != "excluded"], "call_rar"),
         "reference  rarefied, NO abundance gate (not the pre-specified analysis)")
    show(table([r for r in rows if gated(r)], "call_raw"),
         "reference  gated but NOT rarefied (not the pre-specified analysis)")


if __name__ == "__main__":
    main()

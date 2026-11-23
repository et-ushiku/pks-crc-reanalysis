"""【後付けの感度分析】Tataka 2026 の別の判定ルールで、結論が変わるかを見る。

これは事前登録の判定ではない
  本番の判定は `覆い >= 0.5 かつ 遺伝子 >= 14/18`（集計前に決め、変更していない）。
  Tataka 2026（BMC Microbiol 26:930、日本人 223 名）は、同じ問題について
      clbB 単独          AUC 0.826
      19 遺伝子のリード総和 AUC 0.838
      「>= 8 遺伝子」ルール AUC 0.734〜0.773
  と報告しており、**遺伝子数ルールが劣る**と示している。

  だから「遺伝子数ルールをやめて clbB 単独にすればよい」のではない。
  結果を見た後に閾値を動かせば事前登録の意味が消える。
  やるべきことは**事前登録の判定を主とし、別法を後付けと明記して併記する**こと。
  ここで見たいのは一点だけ ―
      **ルールを替えたら結論（差があるか）が変わるのか。**
  変わらなければ、結論はルールの選び方に依存していないと言える。
  変わるなら、それも正直に書く。どちらでも閾値は動かさない。

使い方: python3 sens_tataka.py [コホート]   既定 cn
"""
import glob
import os
import re
import sys
from collections import defaultdict

B = "/mnt/c/crc_strain_pilot"
BED = B + "/pks/ref/clb_genes.bed"

COH = {
    "cn": ("CN-CRC（Yu 2017・香港）", B + "/pks/out_cn", B + "/meta/cn_map.log"),
    "at": ("AT-CRC（Feng 2015・墺）", B + "/pks/out_at", B + "/meta/at_map.log"),
    "de": ("DE-CRC（独・途中）", B + "/pks/out_de", B + "/meta/de_map.log"),
}

MIN_DEPTH = 1
GENE_COVER = 0.5
CALL_BREADTH = 0.5
CALL_GENES = 14
PER_PCT = 313.0
ECOLI_MIN = 0.1


def load_bed(path):
    out = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            p = line.split()
            if len(p) >= 4:
                out.append((p[3], int(p[1]), int(p[2])))
    return out


def gene_cov(path, regions):
    cov = defaultdict(int)
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            p = line.split("\t")
            if len(p) < 3:
                continue
            if int(p[2]) >= MIN_DEPTH:
                pos = int(p[1])
                for name, s, e in regions:
                    if s < pos <= e:
                        cov[name] += 1
                        break
    per = {}
    tl = tc = 0
    for name, s, e in regions:
        ln = e - s
        c = cov.get(name, 0)
        per[name] = c / ln
        tl += ln
        tc += c
    return per, (tc / tl if tl else 0.0)


def read_tot(log):
    if not os.path.exists(log):
        return None
    with open(log, encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            m = re.match(r"\s*(\d+) reads; of these:", line)
            if m:
                return int(m.group(1))
    return None


def fisher(a, b, c, d):
    """2x2 の両側 Fisher。a=がん陽性 b=がん陰性 c=対照陽性 d=対照陰性"""
    from math import lgamma

    def lc(n, k):
        return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)

    n = a + b + c + d
    r1, r2 = a + b, c + d
    c1 = a + c
    base = lc(n, c1)
    obs = lc(r1, a) + lc(r2, c) - base
    tot = 0.0
    for i in range(max(0, c1 - r2), min(r1, c1) + 1):
        p = lc(r1, i) + lc(r2, c1 - i) - base
        if p <= obs + 1e-9:
            tot += pow(2.718281828459045, p)
    return min(1.0, tot)


def main():
    key = sys.argv[1] if len(sys.argv) > 1 else "cn"
    label, out, maplog = COH[key]
    regions = load_bed(BED)

    eco_hit = {}
    if os.path.exists(maplog):
        with open(maplog, encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                if line.startswith("--- 終了"):
                    m = re.search(r"大腸菌=(\d+)", line)
                    if m:
                        eco_hit[line.split()[2]] = int(m.group(1))
    pks_hit = {}
    if os.path.exists(maplog):
        with open(maplog, encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                if line.startswith("--- 終了"):
                    m = re.search(r"pks=(\d+)", line)
                    if m:
                        pks_hit[line.split()[2]] = int(m.group(1))

    rows = []
    for dep in sorted(glob.glob(out + "/*.depth.tsv")):
        s = os.path.basename(dep)[:-len(".depth.tsv")]
        per, br = gene_cov(dep, regions)
        ng = sum(1 for v in per.values() if v >= GENE_COVER)
        tot = read_tot(out + "/" + s + ".bowtie2.log")
        eh = eco_hit.get(s)
        eco = ((eh / (tot / 1e6)) / PER_PCT) if (eh is not None and tot) else None
        pk = pks_hit.get(s, 0)
        # リード総和は読み取り量で割って揃える（100 万本あたり）
        rate = (pk / (tot / 1e6)) if tot else None
        grp = "cancer" if s.startswith("cancer") else "ctrl"
        rows.append(dict(s=s, grp=grp, br=br, ng=ng, clbB=per.get("clbB", 0.0),
                         rate=rate, eco=eco, tot=tot))

    # 大腸菌 >= 0.1% の検体だけを使う（事前に決めた判定可能の条件）
    use = [r for r in rows if r["eco"] is not None and r["eco"] >= ECOLI_MIN]

    print("【後付けの感度分析】%s" % label)
    print("  事前登録の判定は変更していない。ここはルールを替えて結論が動くかの確認。")
    print("  測れた検体 %d / 大腸菌 >= %.1f%% で判定可能 %d"
          % (len(rows), ECOLI_MIN, len(use)))
    nc = sum(1 for r in use if r["grp"] == "cancer")
    print("    うち がん %d / 対照 %d" % (nc, len(use) - nc))
    print()

    # リード総和ルールの閾値は、事前登録の判定と同じ陽性数になる点に置く
    # （閾値を外から持ち込めないので、「同じ陽性率で較べる」形にする）
    def call_pre(r):
        return r["br"] >= CALL_BREADTH and r["ng"] >= CALL_GENES

    npos = sum(1 for r in use if call_pre(r))
    rates = sorted((r["rate"] for r in use if r["rate"] is not None), reverse=True)
    thr = rates[npos - 1] if 0 < npos <= len(rates) else None

    rules = [
        ("事前登録: 覆い>=0.5 かつ 遺伝子>=14/18", call_pre, "（これが主）"),
        ("後付け: clbB 単独の覆い >= 0.5", lambda r: r["clbB"] >= 0.5,
         "Tataka 2026 で最良のひとつ（AUC 0.826）"),
        ("後付け: clbB 単独の覆い >= 0.9", lambda r: r["clbB"] >= 0.9, "より厳しく"),
        ("後付け: 遺伝子 >= 8/18（覆いの条件なし）", lambda r: r["ng"] >= 8,
         "Tataka 2026 が劣ると示したルール"),
    ]
    if thr is not None:
        rules.append(("後付け: pks リード率 >= %.1f 本/百万（陽性数を揃えた点）" % thr,
                      lambda r: r["rate"] is not None and r["rate"] >= thr,
                      "Tataka 2026 で最良（AUC 0.838）の代用"))

    print("%-44s %12s %12s %9s %7s"
          % ("ルール", "がん陽性", "対照陽性", "オッズ比", "P"))
    for name, fn, note in rules:
        a = sum(1 for r in use if r["grp"] == "cancer" and fn(r))
        b = sum(1 for r in use if r["grp"] == "cancer" and not fn(r))
        c = sum(1 for r in use if r["grp"] == "ctrl" and fn(r))
        d = sum(1 for r in use if r["grp"] == "ctrl" and not fn(r))
        orr = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
        p = fisher(a, b, c, d)
        print("%-44s %4d/%-7d %4d/%-7d %9.2f %7.3f"
              % (name, a, a + b, c, c + d, orr, p))
        print("%-44s   %s" % ("", note))

    print()
    print("読み方")
    print("  どのルールでも P が 0.05 を下回らないなら、**結論はルールの選び方に")
    print("  依存していない**と言える。これが確認したかった一点である。")
    print("  いずれかで有意になるなら、それは『そのルールなら有意になった』という")
    print("  後付けの所見にすぎず、事前登録の結論を置き換えるものではない。")
    print("  （複数のルールを試して有意なものを拾えば、それが多重比較そのものになる。）")


if __name__ == "__main__":
    main()

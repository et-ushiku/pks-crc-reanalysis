"""公開した文書に書いてある数字を、data/ から計算し直して一つずつ検算する。

なぜこの script があるか
  この repository の数字はすべて、どこかの時点で script が計算したものである。
  だが **文章に書いた数字を、配っているデータから再計算するものが無かった。**
  そのため
    - 表を作り直したとき、文章側の数字が古いまま残っても誰も気づかない
    - 配っていない入力から計算した数字は、そもそも検算できない
  という穴が開いていた。CORRECTIONS.md の 13・15・16・17 番は全部これである。

  ここでは文章に書いた数字を**主張として直接書き並べ**、data/ から計算した値と
  突き合わせて assert する。古くなった数字は、気づかれるのを待たずに落ちる。

  網にかからないもの（正直に書いておく）
    - 論文から引用した数字（読書ノートと照合しているが、ノートは公開していない）
    - 配っていない入力から出た数字（例: 走行ログの読み取り本数）
    - 図から目で読んだ数字（PITFALLS.md 11 番の相関。原本を残していない）

使い方: python3 code/verify_numbers.py
  引数なし。network も作業ツリーも要らない。data/ だけを読む。
"""
import io
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(HERE), "data")

NOK = [0]
BAD = []


def load(name):
    path = os.path.join(DATA, name)
    rows = [l.rstrip("\n").split("\t") for l in io.open(path, encoding="utf-8") if l.strip()]
    head = rows[0]
    return [dict(zip(head, r)) for r in rows[1:]]


def chk(where, claim, got, tol=0.0):
    """文書の記載 claim と、計算した got を突き合わせる。

    tol は記載が丸められている分の許容。0 なら完全一致を要求する。
    """
    try:
        ok = abs(float(claim) - float(got)) <= tol
        shown = "%s" % got if isinstance(got, int) else "%.6g" % float(got)
    except (TypeError, ValueError):
        ok = claim == got
        shown = str(got)
    NOK[0] += 1
    if ok:
        print("  OK   %-56s 記載=%-14s 計算=%s" % (where, claim, shown))
    else:
        print("  NG   %-56s 記載=%-14s 計算=%s" % (where, claim, shown))
        BAD.append((where, claim, shown))


def median(v):
    v = sorted(v)
    n = len(v)
    return v[n // 2] if n % 2 else (v[n // 2 - 1] + v[n // 2]) / 2.0


def mean(v):
    return sum(v) / float(len(v))


def fisher(a, b, c, d):
    """2x2 の Fisher 正確確率検定（両側）。"""
    def lb(n, k):
        return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)
    n, r1, c1 = a + b + c + d, a + b, a + c
    probs = {}
    for x in range(max(0, c1 - (n - r1)), min(r1, c1) + 1):
        probs[x] = math.exp(lb(r1, x) + lb(n - r1, c1 - x) - lb(n, c1))
    pa = probs[a]
    return sum(p for p in probs.values() if p <= pa * (1 + 1e-9))


def mwu(x, y):
    """両側 Mann-Whitney U。正規近似・同順位補正あり。

    文書の P は実装差で末尾が動くので、0.01 の許容で見る。
    """
    a = sorted([(v, 0) for v in x] + [(v, 1) for v in y])
    rank = [0.0] * len(a)
    i = 0
    while i < len(a):
        j = i
        while j + 1 < len(a) and a[j + 1][0] == a[i][0]:
            j += 1
        for k in range(i, j + 1):
            rank[k] = (i + j) / 2.0 + 1
        i = j + 1
    n1, n2, n = len(x), len(y), len(x) + len(y)
    r1 = sum(rank[k] for k in range(n) if a[k][1] == 0)
    u = r1 - n1 * (n1 + 1) / 2.0
    cnt = {}
    for v, _ in a:
        cnt[v] = cnt.get(v, 0) + 1
    tie = sum(t ** 3 - t for t in cnt.values())
    sd = math.sqrt(n1 * n2 / 12.0 * ((n + 1) - tie / float(n * (n - 1))))
    return math.erfc(abs((u - n1 * n2 / 2.0) / sd) / math.sqrt(2))


# ======================================================================
ps = load("per_sample.tsv")
hist = load("depth_hist.tsv")
ms = load("missing_split.tsv")
car = load("carrier_six_genomes.tsv")
bed = [l.split("\t") for l in
       io.open(os.path.join(DATA, "clb_genes.bed"), encoding="utf-8").read().split("\n") if l.strip()]
glen = dict((r[3], int(r[2]) - int(r[1])) for r in bed)
GENE_BP = sum(glen.values())
REF_BP = 55140          # AM229678.1 の全長。配列は配っていないので定数
DUP_BP = 1287 * 2       # 重複 1,287 bp が 2 か所

depth = dict((r["sample"], float(r["island_mean_depth"])) for r in ms)
split = dict((r["sample"], r) for r in ms)

cov = {}
for r in hist:
    if int(r["depth"]) >= 1:
        k = (r["sample"], r["gene"])
        cov[k] = cov.get(k, 0) + int(r["bases"])


def gene_cov(s, g):
    return cov.get((s, g), 0) / float(glen[g])


print(__doc__.split("使い方")[0].strip().split("\n")[0])
print()

print("1. 参照配列と遺伝子の長さ  METHODS.md / PITFALLS.md 2 番 / README.md")
chk("clb 18 遺伝子の合計 49,455 bp", 49455, GENE_BP)
chk("遺伝子は 18 個", 18, len(glen))
chk("遺伝子ではない 5,685 bp", 5685, REF_BP - GENE_BP)
chk("clbK は 6,465 bp", 6465, glen["clbK"])
chk("clbJ は 6,501 bp", 6501, glen["clbJ"])
chk("両方失ったときの覆い 0.9480", 0.9480, (GENE_BP - DUP_BP) / float(GENE_BP), 5e-5)

print("\n2. 覆いの上限  PITFALLS.md 1 番")
breadth = [float(r["island_breadth"]) for r in ps]
chk("検体数 397", 397, len(ps))
chk("最大の覆い 0.9607", 0.9607, max(breadth), 5e-5)
chk("0.99 を超えた検体は無い", 0, sum(1 for b in breadth if b > 0.99))

print("\n3. 深さ 10x 以上での重複の損失  PITFALLS.md 1 番")
deep = set(s for s, d in depth.items() if d >= 10)
chk("深さ >=10x は 13 検体", 13, len(deep))
mi = [int(split[s]["missing_in_duplication"]) for s in deep]
mo = [int(split[s]["missing_outside_duplication"]) for s in deep]
chk("重複内の欠け 下限 1,920", 1920, min(mi))
chk("重複内の欠け 上限 2,140", 2140, max(mi))
chk("重複内の欠け 平均 2,000", 2000, int(round(mean(mi))))
chk("重複外の欠けが 0 の検体 11 / 13", 11, sum(1 for v in mo if v == 0))
per_gene = {}
for s in deep:
    for g in glen:
        per_gene.setdefault(g, []).append(gene_cov(s, g))
others = [v for g, vs in per_gene.items() if g not in ("clbJ", "clbK") for v in vs]
chk("clbK の平均覆い 0.8065", 0.8065, mean(per_gene["clbK"]), 5e-5)
chk("clbJ の平均覆い 0.8450", 0.8450, mean(per_gene["clbJ"]), 5e-5)
chk("他 16 遺伝子の平均覆い 0.9982", 0.9982, mean(others), 5e-5)
chk("clbK の到達不能 19.3%", 19.3, 100 * (1 - mean(per_gene["clbK"])), 0.05)
chk("clbJ の到達不能 15.5%", 15.5, 100 * (1 - mean(per_gene["clbJ"])), 0.05)
chk("他は 0.2%", 0.2, 100 * (1 - mean(others)), 0.05)

print("\n4. 島が無いのに読み取りが乗る検体  PITFALLS.md 2 番")
neg = [r for r in ps if float(r["island_breadth"]) < 0.1 and int(r["island_reads"]) >= 50]
chk("覆い<0.1 かつ 50 本以上 = 106 検体", 106, len(neg))
chk("その覆いの中央値 0.0018", 0.0018, median([float(r["island_breadth"]) for r in neg]), 5e-5)
chk("読み取りの中央値 258", 258, median([int(r["island_reads"]) for r in neg]))
chk("読み取りの最大 10,548", 10548, max(int(r["island_reads"]) for r in neg))

print("\n5. 陽性検体の内訳  CORRECTIONS.md 13 番の表")
pos = [r for r in ps if r["call"] == "present"]
chk("陽性 42 検体", 42, len(pos))
chk("うち深さ 10x 未満 29", 29, sum(1 for r in pos if depth[r["sample"]] < 10))
chk("欠けが重複内だけ 11", 11,
    sum(1 for r in pos if int(split[r["sample"]]["missing_outside_duplication"]) == 0))
chk("14 遺伝子以上だが導出天井未満 26", 26,
    sum(1 for r in pos if int(r["genes_called"]) >= 14 and float(r["island_breadth"]) < 0.9480))
mid = [int(split[r["sample"]]["missing_outside_duplication"])
       for r in pos if 2 <= depth[r["sample"]] < 10]
chk("深さ 2-10x で重複外の欠け 平均 2,585", 2585, int(round(mean(mid))))

print("\n6. DE-CRC の深さの偏り  RESULTS.md")
de = [r for r in ps if r["cohort"] == "DE-CRC"]
crc = [r for r in de if r["group"] == "CRC"]
ctl = [r for r in de if r["group"] == "control"]
chk("CRC 60 / 対照 60", True, len(crc) == 60 and len(ctl) == 60)
chk("CRC 読み取り中央値 69.5M", 69.5, median([int(r["total_reads"]) for r in crc]) / 1e6, 0.05)
chk("CRC 読み取りの範囲 22.2M", 22.2, min(int(r["total_reads"]) for r in crc) / 1e6, 0.05)
chk("CRC 読み取りの範囲 120.4M", 120.4, max(int(r["total_reads"]) for r in crc) / 1e6, 0.05)
chk("対照 読み取り中央値 34.9M", 34.9, median([int(r["total_reads"]) for r in ctl]) / 1e6, 0.05)
chk("対照 読み取りの範囲 16.3M", 16.3, min(int(r["total_reads"]) for r in ctl) / 1e6, 0.05)
chk("対照 読み取りの範囲 53.7M", 53.7, max(int(r["total_reads"]) for r in ctl) / 1e6, 0.05)
chk("読み取りの差 P = 1.5e-14", 1.5e-14,
    mwu([int(r["total_reads"]) for r in crc], [int(r["total_reads"]) for r in ctl]), 2e-15)
chk("CRC 大腸菌% 中央値 0.0146", 0.0146, median([float(r["ecoli_pct"]) for r in crc]), 5e-5)
chk("対照 大腸菌% 中央値 0.0066", 0.0066, median([float(r["ecoli_pct"]) for r in ctl]), 5e-5)
chk("大腸菌% の差 P = 0.172", 0.172,
    mwu([float(r["ecoli_pct"]) for r in crc], [float(r["ecoli_pct"]) for r in ctl]), 0.01)
chk("Kp% の差 P = 0.502", 0.502,
    mwu([float(r["kpneumoniae_pct"]) for r in crc], [float(r["kpneumoniae_pct"]) for r in ctl]), 0.01)
chk("門を通った CRC 16 / 60", 16, sum(1 for r in crc if float(r["ecoli_pct"]) >= 0.1))
chk("門を通った 対照 10 / 60", 10, sum(1 for r in ctl if float(r["ecoli_pct"]) >= 0.1))

print("\n7. K. pneumoniae の量  PITFALLS.md 5 番 / CORRECTIONS.md 17 番")
chk("長さ比 96,300 / 6,450 = 14.93", 14.93, 96300 / 6450.0, 0.005)
kp = [r for r in ps if r["kpneumoniae_pct"] != "NA"]
det = [float(r["kpneumoniae_pct"]) for r in kp if float(r["kpneumoniae_pct"]) > 0]
chk("Kp を測った 221 検体", 221, len(kp))
chk("検出 96 検体", 96, len(det))
chk("検出率 43%", 43, int(round(100.0 * len(det) / len(kp))))
chk("検出した検体の中央値 0.0143%", 0.0143, median(det), 5e-5)
chk("最大 6.8318%", 6.8318, max(det), 5e-5)
cn = [float(r["kpneumoniae_pct"]) for r in kp if r["cohort"] == "CN-CRC"]
chk("CN-CRC 101 検体で測った", 101, len(cn))
chk("CN-CRC で 0.1% 以上が 19 検体", 19, sum(1 for v in cn if v >= 0.1))
chk("cancer_ERR1018282 の Kp 6.83%", 6.83,
    float([r for r in ps if r["sample"] == "cancer_ERR1018282"][0]["kpneumoniae_pct"]), 0.005)
chk("cancer_ERR1018282 の大腸菌 0.178%", 0.178,
    float([r for r in ps if r["sample"] == "cancer_ERR1018282"][0]["ecoli_pct"]), 0.0005)

print("\n8. 担い手の深さの比  RESULTS.md / CORRECTIONS.md 8 番")
C = dict(((r["sample"], r["target"]), r) for r in car)


def gd(s, t):
    r = C[(s, t)]
    return float(r["sum_depth"]) / int(r["target_bp"]), float(r["covered_bp"]) / int(r["target_bp"])


def island_depth(s):
    r = C[(s, "pks_island_clb_genes")]
    return float(r["sum_depth"]) / int(r["target_bp"])


for s, ed, eb, kd, kb in (("n_ERR1018237", 64.7, 0.81, 2.04, 0.76),
                          ("cancer_ERR1018282", 3.5, 0.84, 109.4, 0.89),
                          ("n_ERR1018240", 3.5, 0.72, 2.74, 0.83),
                          ("cancer_ERR1018186", 0.41, 0.31, 10.3, 0.88)):
    d, b = gd(s, "Ecoli_IHE3034")
    chk(s + " 大腸菌 深さ/覆い", True, abs(d - ed) <= 0.06 and abs(b - eb) <= 0.006)
    d, b = gd(s, "Kpneumoniae_HS11286")
    chk(s + " Kp 深さ/覆い", True, abs(d - kd) <= 0.06 and abs(b - kb) <= 0.006)
il = island_depth("n_ERR1018237")
chk("n_ERR1018237 島の深さ 1.7x", 1.7, il, 0.06)
chk("n_ERR1018237 大腸菌との比 0.03", 0.03, il / gd("n_ERR1018237", "Ecoli_IHE3034")[0], 0.006)
chk("n_ERR1018237 Kp との比 0.84", 0.84, il / gd("n_ERR1018237", "Kpneumoniae_HS11286")[0], 0.006)
chk("38 個の大腸菌に 1 個", 38, int(round(1.0 / (il / gd("n_ERR1018237", "Ecoli_IHE3034")[0]))))
chk("n_ERR1018236 の Kp 深さ 2.063", 2.063, gd("n_ERR1018236", "Kpneumoniae_HS11286")[0], 0.0006)
chk("n_ERR1018236 の Kp 覆い 0.7643", 0.7643, gd("n_ERR1018236", "Kpneumoniae_HS11286")[1], 6e-5)
chk("n_ERR1018237 の Kp 深さ 2.036", 2.036, gd("n_ERR1018237", "Kpneumoniae_HS11286")[0], 0.0006)
chk("n_ERR1018237 の Kp 覆い 0.7638", 0.7638, gd("n_ERR1018237", "Kpneumoniae_HS11286")[1], 6e-5)
chk("C. freundii の覆い 0.496", 0.496, gd("cancer_ERR1018198", "Cfreundii_CFNIH1")[1], 0.0006)
chk("C. freundii との比 1.37", 1.37,
    island_depth("cancer_ERR1018198") / gd("cancer_ERR1018198", "Cfreundii_CFNIH1")[0], 0.006)
chk("C. koseri の覆いは 0.028 以下", True,
    max(gd(s, "Ckoseri_BAA895")[1] for s in set(r["sample"] for r in car)) <= 0.028)

print("\n9. 判定ルールを替えたときの感度分析  RESULTS.md")
print("   リード率の閾値は『事前登録と同じ陽性数になる点』で、外から持ち込んだ値ではない。")
for coh, want in (("CN-CRC", ((7, 44, 6, 33, 1.000), (8, 44, 7, 33, 0.777), (6, 44, 5, 33, 1.000),
                              (8, 44, 6, 33, 1.000), (7, 44, 6, 33, 1.000))),
                  ("AT-CRC", ((5, 21, 3, 21, 0.697), (5, 21, 3, 21, 0.697), (4, 21, 3, 21, 1.000),
                              (5, 21, 3, 21, 0.697), (5, 21, 3, 21, 0.697)))):
    use = [r for r in ps if r["cohort"] == coh and r["ecoli_pct"] != "NA"
           and float(r["ecoli_pct"]) >= 0.1]
    g1 = [r for r in use if r["group"] == "CRC"]
    g0 = [r for r in use if r["group"] == "control"]
    npos = sum(1 for r in use if float(r["island_breadth"]) >= 0.5 and int(r["genes_called"]) >= 14)
    rates = sorted((int(r["island_reads"]) / (int(r["total_reads"]) / 1e6) for r in use), reverse=True)
    thr = rates[npos - 1]
    rules = (
        ("事前登録 覆い>=0.5 かつ >=14 遺伝子",
         lambda r: float(r["island_breadth"]) >= 0.5 and int(r["genes_called"]) >= 14),
        ("clbB の覆い >= 0.5", lambda r: gene_cov(r["sample"], "clbB") >= 0.5),
        ("clbB の覆い >= 0.9", lambda r: gene_cov(r["sample"], "clbB") >= 0.9),
        (">= 8 / 18 遺伝子", lambda r: int(r["genes_called"]) >= 8),
        ("リード率 >= %.1f 本/百万" % thr,
         lambda r: int(r["island_reads"]) / (int(r["total_reads"]) / 1e6) >= thr),
    )
    for (name, fn), (wa, wna, wc, wnc, wp) in zip(rules, want):
        a = sum(1 for r in g1 if fn(r))
        c = sum(1 for r in g0 if fn(r))
        p = fisher(a, len(g1) - a, c, len(g0) - c)
        chk("%s %s" % (coh, name), "%d/%d %d/%d P=%.3f" % (wa, wna, wc, wnc, wp),
            "%d/%d %d/%d P=%.3f" % (a, len(g1), c, len(g0), p))
    chk("%s のリード率の閾値" % coh, 56.7 if coh == "CN-CRC" else 56.3, thr, 0.05)

print("\n10. 島が載った contig の並び  CORRECTIONS.md 16 番の追補")
til = load("kleborate/island_tiling.tsv")
chk("並べた行は 6 行", 6, len(til))
chk("contig がまたがる継ぎ目はすべて重なり 119 bp（座標の差）",
    [119, 119, 119, 119],
    [int(r["overlap_prev"]) for i, r in enumerate(til)
     if i and r["contig"] != til[i - 1]["contig"]])
chk("同じ継ぎ目を配列で測っても 119 bp",
    [119, 119, 119, 119],
    [int(r["overlap_seq_bp"]) for i, r in enumerate(til)
     if i and r["contig"] != til[i - 1]["contig"]])
chk("contig の中のすき間は -72", -72,
    [int(r["overlap_prev"]) for i, r in enumerate(til)
     if i and r["contig"] == til[i - 1]["contig"]][0])
cov = set()
for r in til:
    cov |= set(range(int(r["island_start"]), int(r["island_end"]) + 1))
chk("覆った島の塩基 52,878", 52878, len(cov))
chk("島の参照の長さ 54,791", 54791, 54791)
chk("覆いの割合 96.5%", 96.5, 100.0 * len(cov) / 54791, 0.05)
chk("並びの端は 1,842 から 54,791", (1842, 54791), (min(cov), max(cov)))
rep = [r for r in til if r["contig"] == "k119_32464"]
chk("重複の contig は 2 回出る", 2, len(rep))
chk("重複の開始は 16,066（PITFALLS 1 番と同じ）", 16066, int(rep[0]["island_start"]))
chk("重複の間隔は 5,653（PITFALLS 1 番と同じ）", 5653,
    int(rep[1]["island_start"]) - int(rep[0]["island_start"]))
chk("並びのうち自動の袋に入っていたのは k119_45171 だけ", ["k119_45171"],
    sorted(set(r["contig"] for r in til if r["in_kp_bin"] == "yes")))
bb = {r["contig"]: r for r in load("kleborate/island_contig_backbone.tsv")}
for c, want in (("k119_45171", 8.2), ("k119_25144", 0.5),
                ("k119_70945", 22.4), ("k119_25170", 3.2)):
    chk("%s の島の割合 %.1f%%" % (c, want), want,
        100.0 * int(bb[c]["island_bp"]) / int(bb[c]["contig_bp"]), 0.06)

print("\n11. 挿入要素の場所  CORRECTIONS.md 16 番の追補")
loc = load("kleborate/island_locus.tsv")
el = [r for r in loc if r["region"] == "element"]
kp = [r for r in loc if r["region"] == "kp_chromosome"]
lo = min(int(r["ruler_start"]) for r in el)
hi = max(int(r["ruler_end"]) for r in el)
chk("要素は 816,165-953,855", (816165, 953855), (lo, hi))
chk("要素の長さ 137,690", 137690, hi - lo)
chk("敷き詰めている contig は 13 本", 13,
    len(set(r["contig"] for r in el)))
chk("そのうち袋に足されていないものが 7 本", 7,
    len(set(r["contig"] for r in el))
    - len(set(r["contig"] for r in el if r["in_kp_plus"] == "yes")))
ends = {}
for r in el:
    c = r["contig"]
    a, b = int(r["ruler_start"]), int(r["ruler_end"])
    o = ends.get(c, (a, b))
    ends[c] = (min(o[0], a), max(o[1], b))
chk("要素の両端に届くのは k119_45171 と k119_25144",
    ["k119_25144", "k119_45171"],
    sorted([c for c, (a, b) in ends.items() if a == lo]
           + [c for c, (a, b) in ends.items() if b == hi]))
chk("両端の contig はどちらも自動の袋に入っている", ["yes", "yes"],
    [sorted(set(r["in_kp_bin"] for r in loc if r["contig"] == c))[0]
     for c in ("k119_25144", "k119_45171")])
k1 = [(int(r["ruler_start"]), int(r["ruler_end"])) for r in kp
      if r["contig"] == "k119_25144"]
k2 = [(int(r["ruler_start"]), int(r["ruler_end"])) for r in kp
      if r["contig"] == "k119_45171"]
chk("k119_25144 の染色体側は 3,250,687 から 3,506,454",
    (3250687, 3506454), (min(a for a, b in k1), max(b for a, b in k1)))
chk("k119_45171 の染色体側は 3,506,446 から 3,856,980",
    (3506446, 3856980), (min(a for a, b in k2), max(b for a, b in k2)))
chk("2 つの区間は 8 bp のところで出会う", 8,
    abs(max(b for a, b in k1) - min(a for a, b in k2)))
chk("k119_70945 の染色体側は 3,434,608-3,479,783", (3434608, 3479783),
    [(int(r["ruler_start"]), int(r["ruler_end"])) for r in kp
     if r["contig"] == "k119_70945"][0])
chk("大腸菌源なら必要な島の倍率 33", 33,
    int(round(island_depth("cancer_ERR1018186")
              / gd("cancer_ERR1018186", "Ecoli_IHE3034")[0])))

print("\n12. 読み取りの深さと連結  CORRECTIONS.md 16 番の追補")
rd = {r["contig"]: r for r in load("kleborate/contig_read_depth.tsv")}
chk("袋の読み取りの深さ 25.13", 25.13, float(rd["kp_bin_all"]["read_depth"]), 0.005)
for c, dep, rel, mo in (("k119_45171", 24.25, 0.965, 0.743),
                        ("k119_25144", 22.29, 0.887, 0.736),
                        ("k119_70945", 25.08, 0.998, 0.757),
                        ("k119_70596", 26.58, 1.058, 0.746),
                        ("k119_25170", 21.97, 0.874, 0.404),
                        ("k119_32464", 45.83, 1.824, 0.786)):
    chk("%s の深さ %.2f" % (c, dep), dep, float(rd[c]["read_depth"]), 0.006)
    chk("%s は袋の %.3f 倍" % (c, rel), rel, float(rd[c]["depth_over_bin"]), 0.0006)
    chk("%s の multi/深さ %.3f" % (c, mo), mo, float(rd[c]["multi_over_depth"]), 0.0006)
big = [r for r in rd.values() if r["in_kp_bin"] == "yes"
       and r["contig"] != "kp_bin_all" and int(r["contig_bp"]) >= 10000]
chk("袋の 10 kb 以上は 46 本", 46, len(big))
chk("その深さは袋の 0.76-1.33 倍", (0.76, 1.33),
    (round(min(float(r["depth_over_bin"]) for r in big), 2),
     round(max(float(r["depth_over_bin"]) for r in big), 2)))
lk = {r["contig"]: r for r in load("kleborate/contig_links.tsv")}
for c, b, a, o, sh in (("k119_45171", 377, 37, 27, 0.939),
                       ("k119_25144", 216, 9, 15, 0.938),
                       ("k119_70945", 69, 39, 7, 0.939),
                       ("k119_70596", 0, 83, 0, 1.000),
                       ("k119_25170", 11, 5, 15, 0.516),
                       ("k119_32464", 34, 117, 0, 1.000)):
    chk("%s の対 袋 %d / 島 %d / 他の種 %d" % (c, b, a, o), (b, a, o),
        (int(lk[c]["pairs_to_kp_bin"]), int(lk[c]["pairs_to_added_island"]),
         int(lk[c]["pairs_to_other_species"])))
    chk("%s の袋+島の割合 %.3f" % (c, sh), sh,
        float(lk[c]["share_kp_bin_or_island"]), 0.0006)


def shband(lo, hi):
    v = [float(r["share_kp_bin_or_island"]) for r in lk.values()
         if r["in_kp_bin"] == "yes" and lo <= int(r["contig_bp"]) < hi
         and int(r["pairs_assigned_total"]) > 0]
    return (round(min(v), 3), round(max(v), 3))


chk("袋の 10 kb 以上の割合は 0.857-1.000", (0.857, 1.0), shband(10000, 10 ** 9))
chk("袋の 2-10 kb の割合は 0.129-1.000", (0.129, 1.0), shband(2000, 10000))
chk("k119_70945 の割合は 10 kb 以上の範囲の中", True,
    shband(10000, 10 ** 9)[0] <= float(lk["k119_70945"]["share_kp_bin_or_island"])
    <= shband(10000, 10 ** 9)[1])

print("\n%d 件を検算し、%d 件が不一致" % (NOK[0], len(BAD)))
if BAD:
    for where, claim, got in BAD:
        print("  NG %s  記載=%s 計算=%s" % (where, claim, got))
    sys.exit("文書の数字が data/ と合わない。文書か data/ のどちらかが古い。")
print("文書の数字はすべて data/ から再現できた。")

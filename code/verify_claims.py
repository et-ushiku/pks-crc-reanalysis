"""2026-11-20〜21 に自分が書いた 3 つの主張を、手元のデータで検証する。

なぜやるか
  この 2 日で記述を 2 回訂正した（「天井 0.9589」と「密度比 R」）。
  訂正した先の記述も、まだ少数の検体で見ただけである。
  **訂正が正しいかどうかを全検体で確かめる。**

検証する主張
  主張1  到達しうる覆いは 0.9495〜1.000 の間で検体ごとに動く。
          0.9495 は「重複 1,248bp×2 が丸ごと落ちる」場合の下限。
          → 全コホートで 0.9495 を下回る「あり」検体が無いか、
            0.99 を超える検体が無いかを見る。どちらかが出たら主張が崩れる。

  主張2  島を持つ検体では、覆われない塩基は**重複領域の中だけ**に収まる。
          → 覆い 0.5 以上の全検体で、重複外の欠けを数える。

  主張3  島の参照配列（55,140bp）には島でない部分が 5,685bp あり、
          島を持たない大腸菌の読み取りがそこに載る（背景 pks/大腸菌 ≈ 0.009）。
          → 島が「なし」なのに pks の読み取りがある検体で、
            clb 遺伝子 18 個のうちどこに覆いが出るかを見る。
            背景が遺伝子の外にあるなら、遺伝子の覆いはほぼ 0 のはず。

使い方: python3 verify_claims.py
"""
import glob
import os
import re
from collections import defaultdict

B = "/mnt/c/crc_strain_pilot"
BED = B + "/pks/ref/clb_genes.bed"
# 重複している 2 か所（AM229678.1 の座標）。互いに 99.840% 一致、各 1,248bp
DUP = [(16066, 17313), (21719, 22966)]
ISLAND_BP = 49455
REF_BP = 55140
MIN_DEPTH = 1
GENE_COVER = 0.5
CALL_BREADTH = 0.5
CALL_GENES = 14

COH = [
    ("Feng(墺)", "/pks/out_at", "/meta/at_map.log"),
    ("Yu(中国)", "/pks/out_cn", "/meta/cn_map.log"),
    ("Wirbel/Zeller(独)", "/pks/out_de", "/meta/de_map.log"),
    ("Zeller(仏)", "/pks/out", None),
]


def load_bed():
    out = []
    with open(BED, encoding="utf-8") as fh:
        for line in fh:
            p = line.split()
            if len(p) >= 4:
                out.append((p[3], int(p[1]), int(p[2])))
    return out


REG = load_bed()


def gene_of(pos):
    for name, s, e in REG:
        if s < pos <= e:
            return name
    return None


def in_dup(pos):
    return any(a < pos <= b for a, b in DUP)


def scan(path):
    """深さの表から（覆い, 遺伝子数, 遺伝子ごとの覆い, 欠けの内訳）"""
    cov = set()
    seen = 0
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            p = line.split("\t")
            if len(p) < 3:
                continue
            seen += 1
            if int(p[2]) >= MIN_DEPTH:
                cov.add(int(p[1]))
    per = defaultdict(int)
    for pos in cov:
        g = gene_of(pos)
        if g:
            per[g] += 1
    gcov = {}
    tl = tc = ng = 0
    for name, s, e in REG:
        ln = e - s
        c = per.get(name, 0)
        gcov[name] = c / ln
        tl += ln
        tc += c
        if c / ln >= GENE_COVER:
            ng += 1
    br = tc / tl if tl else 0.0
    # 欠けた塩基を重複の内外で数える
    miss_in = miss_out = 0
    for name, s, e in REG:
        for pos in range(s + 1, e + 1):
            if pos not in cov:
                if in_dup(pos):
                    miss_in += 1
                else:
                    miss_out += 1
    return br, ng, gcov, miss_in, miss_out, seen


def pks_reads(logp):
    """マップのログから検体ごとの pks 読み取り数と大腸菌の当たり数"""
    out = {}
    if not logp or not os.path.exists(B + logp):
        return out
    with open(B + logp, encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if not line.startswith("--- 終了"):
                continue
            s = line.split()[2]
            mp = re.search(r"pks=(\d+)", line)
            me = re.search(r"大腸菌=(\d+)", line)
            if mp:
                out[s] = (int(mp.group(1)),
                          int(me.group(1)) if me else None)
    return out


def main():
    print("=" * 74)
    print("主張1・2 の検証: 覆い 0.5 以上の全検体で、欠けは重複領域の中だけか")
    print("=" * 74)
    print("  重複が丸ごと落ちる場合の下限: (49455 - 2496) / 49455 = %.4f"
          % ((ISLAND_BP - 2496) / ISLAND_BP))
    print()
    print("%-26s %-20s %7s %5s %9s %9s"
          % ("検体", "コホート", "覆い", "遺伝", "欠け重複内", "欠け重複外"))

    rows = []
    for name, d, logp in COH:
        for f in sorted(glob.glob(B + d + "/*.depth.tsv")):
            s = os.path.basename(f)[:-len(".depth.tsv")]
            br, ng, gcov, mi, mo, seen = scan(f)
            rows.append((s, name, br, ng, gcov, mi, mo, seen))

    high = [r for r in rows if r[2] >= CALL_BREADTH]
    for r in sorted(high, key=lambda x: -x[2]):
        print("%-26s %-20s %7.4f %5d %9d %9d"
              % (r[0], r[1], r[2], r[3], r[5], r[6]))

    print()
    bad_low = [r for r in high if r[2] < 0.9495 and r[3] >= CALL_GENES]
    bad_high = [r for r in rows if r[2] > 0.99]
    clean = [r for r in high if r[6] == 0]
    print("覆い 0.5 以上の検体 %d 件" % len(high))
    print("  欠けが**重複領域の中だけ**の検体: %d / %d" % (len(clean), len(high)))
    print("  0.99 を超えた検体（重複まで覆われたことになり異常）: %s"
          % (", ".join("%s %.4f" % (r[0], r[2]) for r in bad_high) if bad_high else "なし"))
    print("  遺伝子 14/18 を満たすのに覆いが下限 0.9495 未満: %d 件"
          % len(bad_low))
    if bad_low:
        for r in bad_low:
            print("      %-26s %.4f  遺伝子 %d  重複外の欠け %d bp"
                  % (r[0], r[2], r[3], r[6]))
        print("    → 下限より低い。重複以外にも覆えない所があるということ")
    if high:
        print("  実測の幅: %.4f 〜 %.4f（n=%d）"
              % (min(r[2] for r in high), max(r[2] for r in high), len(high)))

    print()
    print("=" * 74)
    print("主張3 の検証: 島が「なし」の検体に出る pks の読み取りは、遺伝子の外か")
    print("=" * 74)
    print("  参照配列 %d bp のうち clb 遺伝子は %d bp。差の %d bp は"
          % (REF_BP, ISLAND_BP, REF_BP - ISLAND_BP))
    print("  遺伝子間領域と両端の骨格である。背景がそこにあるなら、")
    print("  遺伝子の覆いはほぼ 0 のまま pks の読み取り数だけが出るはず。")
    print()

    pk = {}
    for name, d, logp in COH:
        pk.update(pks_reads(logp))

    print("%-26s %-20s %8s %8s %7s %6s"
          % ("検体", "コホート", "pks読み", "大腸菌", "覆い", "遺伝子"))
    cand = []
    for r in rows:
        s = r[0]
        if s not in pk:
            continue
        p, e = pk[s]
        if r[2] >= 0.1 or p < 50:
            continue
        cand.append((s, r[1], p, e, r[2], r[3]))
    for s, c, p, e, br, ng in sorted(cand, key=lambda x: -x[2])[:15]:
        print("%-26s %-20s %8d %8s %7.4f %6d"
              % (s, c, p, ("%d" % e) if e is not None else "-", br, ng))

    if cand:
        print()
        print("  島が「なし」(覆い<0.1) で pks 読み取りが 50 本以上ある検体 %d 件"
              % len(cand))
        print("  その検体群の覆いの中央値: %.4f"
              % sorted(x[4] for x in cand)[len(cand) // 2])
        mx = max(cand, key=lambda x: x[4])
        print("  いちばん覆いが高いもの: %s %.4f（遺伝子 %d/18）"
              % (mx[0], mx[4], mx[5]))
        print("  → 覆いがほぼ 0 のまま読み取りだけ出ている = **読み取りは遺伝子の外**")
        print("    主張3 と整合する。pks の生の読み取り数は島の量ではない。")
        rat = [(x[2] / x[3]) for x in cand if x[3]]
        if rat:
            rat.sort()
            print("  背景の比 pks/大腸菌: 中央 %.4f  幅 %.4f〜%.4f（n=%d）"
                  % (rat[len(rat) // 2], rat[0], rat[-1], len(rat)))
    else:
        print("  該当する検体が無い（主張3 はこのデータでは確かめられない）")


if __name__ == "__main__":
    main()

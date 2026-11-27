"""目印の長さに比例する、という換算の仮定を、独立な測定で検証する。

検証する仮定
  目印に当たる本数は「量 x 目印の総塩基」に比例する。したがって
  「1% あたり 100 万本中 313 本」という較正（大腸菌の目印 96,300 塩基で
  一度だけ求めた値）は、長さの違う目印には `marker_bp / 96,300` で伸縮させる。
  K. pneumoniae の目印は 6,450 塩基しかないので 20.96 本になる。
  CORRECTIONS.md 17 番でこの補正を入れた。そこでは「検証できていない仮定」と
  書いた。これはその検証である。

どう検証するか
  同じ検体を、**長さに依存しない別の測定**でも測っている。六ゲノム競合索引
  （完全ゲノム 6 個を 1 つの索引に入れ MAPQ 20）である。そこでは種の量は
  ゲノム全体に当たった塩基の総和であり、目印の長さは関係しない。

  種ごとに
      目印側 = 目印に当たった本数 / 目印の実効塩基
      ゲノム側 = ゲノムに当たった塩基 / ゲノム塩基 / 読み取り長
  を出す。どちらも「その種のゲノム 1 塩基あたりに何本乗ったか」である。
  その比（倍率）は、目印が普通のゲノム配列と同じ率で読み取りを集めるなら 1 に
  なる。1 からずれても構わない（ずれは 313 の較正に吸収されている）。

  **長さ比例の仮定が要求しているのは、倍率が種によらず同じであること**である。
  6,450 塩基の目印と 96,300 塩基の目印で倍率が同じなら、長さで伸縮させる扱いは
  正しい。倍率が 14.93 倍（= 96,300 / 6,450）ずれるなら、伸縮させないほうが
  正しかったことになる。そこを測る。

  端の損: bowtie2 は end-to-end なので、断片の端から外に出る読み取りは当たらない。
  断片 1 つあたり (読み取り長 - 1) 塩基が実質的に死ぬ。大腸菌は 122 断片で
  96,300 塩基、K. pneumoniae は 11 断片で 6,450 塩基なので、損の割合は種で違う。
  入れた場合と入れない場合の両方を出す。

読み取り長
  比の検定には読み取り長は効かない（両種で同じ数で割るので消える）。
  端の損の大きさだけに効くので、84 / 94 / 100 塩基の 3 通りで出し、
  結論がそこに依存しないことを示す。CN-CRC は ENA の base_count/read_count から
  約 94 塩基。

使い方: python3 code/verify_marker_scaling.py
  引数なし。network も作業ツリーも要らない。data/ だけを読む。
"""
import io
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(HERE), "data")

SGB2GENOME = {
    "t__SGB10068": "Ecoli_IHE3034",
    "t__SGB10115": "Kpneumoniae_HS11286",
    "t__SGB10079": "Ckoseri_BAA895",
    "t__SGB10116": "Kaerogenes_KCTC2190",
}
EC, KP = "t__SGB10068", "t__SGB10115"
MIN_BREADTH = 0.7          # 両種ともこれ以上覆われている検体だけを主に使う
# t 分布の両側 95% 点。df 1-12 だけ持つ（検体数がこの範囲なので）
T95 = {1: 12.706, 2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571, 6: 2.447,
       7: 2.365, 8: 2.306, 9: 2.262, 10: 2.228, 11: 2.201, 12: 2.179}


def load(name):
    rows = [l.rstrip("\n").split("\t") for l in
            io.open(os.path.join(DATA, name), encoding="utf-8") if l.strip()]
    return [dict(zip(rows[0], r)) for r in rows[1:]]


def geomean(v):
    return math.exp(sum(math.log(x) for x in v) / len(v))


def logci(v):
    """対数で平均と両側 95% 区間を出し、元の尺度に戻す。"""
    n = len(v)
    lg = [math.log(x) for x in v]
    m = sum(lg) / n
    if n < 2:
        return math.exp(m), None, None
    sd = math.sqrt(sum((x - m) ** 2 for x in lg) / (n - 1))
    half = T95.get(n - 1, 1.96) * sd / math.sqrt(n)
    return math.exp(m), math.exp(m - half), math.exp(m + half)


mark = load("marker_counts.tsv")
car = load("carrier_six_genomes.tsv")
C = dict(((r["sample"], r["target"]), r) for r in car)
M = dict(((r["sample"], r["sgb"]), r) for r in mark)
samples = []
for r in mark:
    if r["sample"] not in samples:
        samples.append(r["sample"])


def genome(s, sgb):
    """その種のゲノムに当たった塩基 / ゲノム塩基 と、覆い。"""
    g = C[(s, SGB2GENOME[sgb])]
    return (float(g["sum_depth"]) / int(g["target_bp"]),
            float(g["covered_bp"]) / int(g["target_bp"]))


def multiplier(s, sgb, col, read_len, edge):
    """目印側の密度 / ゲノム側の密度。"""
    r = M[(s, sgb)]
    bp = float(r["marker_bp"])
    if edge and r["marker_fragments"] != "NA":
        bp -= int(r["marker_fragments"]) * (read_len - 1)
    rho = int(r[col]) / bp
    depth, _ = genome(s, sgb)
    return rho / (depth / float(read_len))


print(__doc__.split("\n")[0])
print()

bp_ec = float(M[(samples[0], EC)]["marker_bp"])
bp_kp = float(M[(samples[0], KP)]["marker_bp"])
UNSCALED = bp_ec / bp_kp
print("目印の長さ: E. coli %d bp / %s 断片, K. pneumoniae %d bp / %s 断片"
      % (bp_ec, M[(samples[0], EC)]["marker_fragments"],
         bp_kp, M[(samples[0], KP)]["marker_fragments"]))
print("長さ比 %.0f / %.0f = %.2f" % (bp_ec, bp_kp, UNSCALED))
print("  長さで伸縮させるのが正しければ   倍率の kp/ec は 1 付近")
print("  伸縮させないのが正しければ       倍率の kp/ec は %.2f 付近" % UNSCALED)
print()

core = [s for s in samples
        if genome(s, EC)[1] >= MIN_BREADTH and genome(s, KP)[1] >= MIN_BREADTH]
print("検体 %d 件のうち、両種ともゲノムの覆いが %.1f 以上の %d 件を使う。"
      % (len(samples), MIN_BREADTH, len(core)))
print("覆いが低い種は、参照株と検体の株が違う分だけ深さを取りこぼすため。")
print()

print("%-20s %9s %9s %9s %9s %8s" %
      ("sample", "ec 覆い", "kp 覆い", "ec 倍率", "kp 倍率", "kp/ec"))
rows = {}
for s in samples:
    be, bk = genome(s, EC)[1], genome(s, KP)[1]
    me = multiplier(s, EC, "reads_mapq20", 94, True)
    mk = multiplier(s, KP, "reads_mapq20", 94, True)
    flag = "" if s in core else "   <- 覆いが足りないので主な計算から外す"
    print("%-20s %9.3f %9.3f %9.3f %9.3f %8.3f%s" % (s, be, bk, me, mk, mk / me, flag))
    rows[s] = mk / me
print()

print("読み取り長と端の損を変えても結論が動かないかを見る（kp/ec の幾何平均）")
print("%-26s %10s %10s %10s" % ("", "84 bp", "94 bp", "100 bp"))
for col, cname in (("reads_no_mapq", "MAPQ なし"), ("reads_mapq20", "MAPQ 20")):
    for edge, ename in ((False, "端の損なし"), (True, "端の損あり")):
        out = []
        for rl in (84, 94, 100):
            v = [multiplier(s, KP, col, rl, edge) / multiplier(s, EC, col, rl, edge)
                 for s in core]
            out.append(geomean(v))
        print("%-26s %10.3f %10.3f %10.3f" % (cname + " / " + ename, out[0], out[1], out[2]))
print("  端の損を入れないと kp/ec が下がる。11 断片しかない K. pneumoniae の目印は")
print("  端の損が相対的に小さいので、入れないと K. pneumoniae を過小に見積もる。")
print("  bowtie2 は end-to-end なので、入れるほうが正しい。")
print()

main = [multiplier(s, KP, "reads_mapq20", 94, True) / multiplier(s, EC, "reads_mapq20", 94, True)
        for s in core]
gm, lo, hi = logci(main)
print("=" * 70)
print("主な計算: MAPQ 20・端の損あり・読み取り長 94 bp・覆い %.1f 以上の %d 検体"
      % (MIN_BREADTH, len(core)))
print("  kp/ec の幾何平均 = %.3f   両側 95%% 区間 %.3f - %.3f" % (gm, lo, hi))
print("  伸縮させない場合に必要な値 = %.2f" % UNSCALED)
print("=" * 70)

assert lo <= 1.0 <= hi, "kp/ec の 95% 区間が 1 を含まない: %.3f - %.3f" % (lo, hi)
assert hi < UNSCALED, "伸縮させない説を排除できていない"
assert lo > 1.0 / UNSCALED, "逆向きの伸縮も排除できていない"
print()
print("結論")
print("  1. 区間は 1 を含む。**長さに比例させる扱いは、この測定と矛盾しない。**")
print("  2. 区間の上端 %.3f は %.2f の %.0f 分の 1 でしかない。"
      % (hi, UNSCALED, UNSCALED / hi))
print("     **伸縮させない扱いは、桁で排除される。**")
print("  3. ただし区間は上端と下端で %.1f 倍の開きがある。**長さ比例が確かめられた"
      % (hi / lo))
print("     のは「2 倍以内」という粗さまで**で、それより細かい議論には使えない。")
print("     この列で 2 倍の差を論じてはいけない。")
ec_m = [multiplier(s, EC, "reads_mapq20", 94, True) for s in core]
kp_m = [multiplier(s, KP, "reads_mapq20", 94, True) for s in core]
print("  4. 倍率そのものも、両種とも 1 の近くに出ている"
      "（E. coli %.2f-%.2f / K. pneumoniae %.2f-%.2f）。"
      % (min(ec_m), max(ec_m), min(kp_m), max(kp_m)))
print("     端の損を入れた目印は、**自分のゲノムとほぼ同じ率で読み取りを集めている**。")
print("     長さ比例が成り立つ物理的な理由はこれで、仮定の中身が直接見えている。")
print()
print("ばらつきの出どころ（検証できていない部分）")
print("  - 参照株と検体の株が違う。覆いは %.2f - %.2f で、欠けた領域の深さは取りこぼす。"
      % (min(min(genome(s, EC)[1], genome(s, KP)[1]) for s in core),
         max(min(genome(s, EC)[1], genome(s, KP)[1]) for s in core)))
print("  - 六ゲノム索引は競合なので、6 種の外にいる腸内細菌科の読み取りが")
print("    どれかに寄ることがある。目印側にも同じことが起きる。")
print("  - プラスミドと可動性要素。参照に無いものは落ち、共有されるものは寄る。")
print("  - 検体が %d 件しかなく、どれも CN-CRC である。他のコホートでは測っていない。"
      % len(core))
print("  - C. koseri と K. aerogenes の目印（279,700 / 138,200 / 282,850 塩基）も")
print("    同じ索引で測ってあるが、この 11 検体ではどちらも不在（覆い 0.06 以下）")
print("    なので、長さの範囲を 44 倍に広げる検証には使えなかった。")
print("    検証できたのは 6,450 対 96,300 の 2 点、つまり問題の 14.93 倍そのものである。")

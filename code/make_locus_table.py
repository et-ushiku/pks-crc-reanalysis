"""島が載った contig どうしの並びを出す。CORRECTIONS.md 16 番の作り直し。

なぜ要るか
  16 番は「手で足した 4 本の contig がこの Klebsiella のものか」を
  **被覆深さと代案の不在**で論じていた。自分でそれを「15 番の接合部より弱い」と
  書いた。弱い理由は 2 つある。
    1. 使った深さは megahit の `multi` である。この組み立ては k=119 で、
       読み取りは 100 塩基しかない。つまり multi は読み取りの深さではなく、
       反復のあいだ持ち越した contig を含む de Bruijn グラフ上の重複度である。
       読み取りの深さに比例する保証がない。
    2. 4 本を**独立な 4 個**として扱っていた。だから 1 本ずつ出どころを
       論じる必要が生じ、深さしか残らなかった。

  しかし 4 本は独立ではない。**島の座標の上で重なって並んでいる。**
  それを測るのがこの script である。深さも配列の似かたも使わない。

出すもの
  data/kleborate/island_tiling.tsv
      pks 島の参照 AM229678.1（54,791 bp）の上に contig を並べたもの。
      隣との重なりを 2 通りで出す。1 つは島の座標の差、もう 1 つは
      **contig の配列そのものを突き合わせた長さ**（前の contig の末尾と
      次の contig の先頭が何塩基一致するか）。座標の差は当たりの境目に
      左右されるので、配列から測ったほうを根拠にする。
      重なりが k = 119 ちょうどなら、その 2 本は組み立てグラフの上で
      **連続した 1 本の道**である。
  data/kleborate/island_locus.tsv
      その島を含む挿入要素ぜんたいを、要素を丸ごと持つ唯一の参照
      （C. koseri BAA-895）を**物差しとして**並べたもの。
      あわせて各 contig の K. pneumoniae 染色体側の当たりも出す。
      要素の範囲は、島を載せた contig から**重なりでつながる範囲**として
      決める（すき間 JOIN_GAP 塩基まで許す）。窓を手で決めない。

  物差しに C. koseri を使うのは、6 種のうちこの要素を丸ごと持つのが
  それだけだからで、**出どころの主張ではない**。16 番が記録したとおり、
  この要素の配列は 2 つの属に相同のものがあり、似かたで種は決まらない。

使い方: python3 make_locus_table.py <出力先ディレクトリ>
"""
import os
import sys

B = "/mnt/c/crc_strain_pilot"
A = B + "/asm/cancer_ERR1018186"
HITS = A + "/megahit/pks_hits.tsv"       # blastn: 問い合わせ = AM229678.1
PAF = A + "/kleborate/vs_g6.paf"         # minimap2: contig を 6 種に当てたもの
CTG = A + "/megahit/final.contigs.fa"
BIN = A + "/kleborate/kp_contigs.txt"
PLUS = A + "/kleborate/kp_plus.txt"

KMER = 119            # 組み立ての k。contig 名の k119_ がそれ
MIN_ID = 99.0         # 島の並べ方に使う当たりの同一性の下限
MIN_LEN = 1000        # 同じ、長さの下限
MIN_ALN = 500         # 要素の並べ方に使う当たりの長さの下限
JOIN_GAP = 200        # これ以下のすき間なら「つながっている」とみなす
RULER = "Ckoseri_BAA895"
KPREF = "Kpneumoniae_HS11286"


def read_lines(p):
    return [l.strip() for l in open(p, encoding="utf-8") if l.strip()]


def contig_len():
    out = {}
    for line in open(CTG, encoding="utf-8"):
        if line.startswith(">"):
            q = line[1:].split()
            out[q[0]] = int(q[-1].split("=")[1])
    return out


def contig_seq(want):
    """要る contig だけ配列を読む。全体は 122 MB あるので選って読む。"""
    out, cur, buf = {}, None, []
    for line in open(CTG, encoding="utf-8"):
        if line.startswith(">"):
            if cur:
                out[cur] = "".join(buf)
            n = line[1:].split()[0]
            cur, buf = (n if n in want else None), []
        elif cur:
            buf.append(line.strip())
    if cur:
        out[cur] = "".join(buf)
    return out


def rc(x):
    return x[::-1].translate(str.maketrans("ACGTacgtN", "TGCAtgcaN"))


def oriented(seq, a, b):
    """blastn の sstart/send を島の向きに揃えて切り出す。"""
    x = seq[min(a, b) - 1:max(a, b)]
    return rc(x) if a > b else x


def shared(x, y):
    """x の末尾と y の先頭が一致する最長の長さ。"""
    for k in range(min(len(x), len(y)), 0, -1):
        if x[-k:] == y[:k]:
            return k
    return 0


def paf():
    """contig ごと・参照ごとの当たり。(contig, 参照, 参照開始, 参照終わり,
    contig 開始, contig 終わり, 向き, 当たった長さ)"""
    out = []
    for line in open(PAF, encoding="utf-8"):
        q = line.rstrip("\n").split("\t")
        if len(q) < 11:
            continue
        ref = q[5].split("|")[0]
        aln = int(q[3]) - int(q[2])
        if aln < MIN_ALN:
            continue
        out.append((q[0], ref, int(q[7]), int(q[8]), int(q[2]), int(q[3]),
                    q[4], aln))
    return out


def main():
    outdir = sys.argv[1]
    ln = contig_len()
    binc = set(read_lines(BIN))
    plus = set(read_lines(PLUS))

    # ---------------------------------------------------------- 1. 島の並び
    tile = []
    for line in open(HITS, encoding="utf-8"):
        q = line.rstrip("\n").split("\t")
        if len(q) < 8:
            continue
        pid, aln = float(q[2]), int(q[3])
        if pid < MIN_ID or aln < MIN_LEN:
            continue
        tile.append((int(q[4]), int(q[5]), q[1], pid, aln, int(q[6]), int(q[7])))
    tile.sort()
    d = os.path.join(outdir, "kleborate")
    if not os.path.isdir(d):
        os.makedirs(d)
    seqs = contig_seq(set(c for _, _, c, _, _, _, _ in tile))
    with open(d + "/island_tiling.tsv", "w", encoding="utf-8", newline="") as fh:
        fh.write("order\tcontig\tcontig_bp\tin_kp_bin\tisland_start\tisland_end"
                 "\taln_bp\tidentity\tcontig_start\tcontig_end\toverlap_prev"
                 "\toverlap_seq_bp\n")
        prev = None
        prevseg = None
        prevc = None
        for i, (s, e, c, pid, aln, cs, ce) in enumerate(tile):
            ov = "NA" if prev is None else str(prev - s + 1)
            seg = oriented(seqs[c], cs, ce)
            # 同じ contig どうしの行は継ぎ目ではないので配列は突き合わせない
            osq = "NA" if (prevseg is None or c == prevc) else str(
                shared(prevseg, seg))
            fh.write("%d\t%s\t%d\t%s\t%d\t%d\t%d\t%.3f\t%d\t%d\t%s\t%s\n"
                     % (i + 1, c, ln[c], "yes" if c in binc else "no",
                        s, e, aln, pid, cs, ce, ov, osq))
            prev = max(prev or 0, e)
            prevseg, prevc = seg, c
    print("島の並び %d 行（同一性 %.1f%% 以上・%d bp 以上）"
          % (len(tile), MIN_ID, MIN_LEN))
    isl_contigs = set(c for _, _, c, _, _, _, _ in tile)

    # ------------------------------------------------- 2. 挿入要素ぜんたい
    rows = paf()
    ck = [r for r in rows if r[1] == RULER]
    # 島を載せた contig の物差し上の区間から始めて、重なりでつなげて広げる
    seed = [(r[2], r[3]) for r in ck if r[0] in isl_contigs]
    lo, hi = min(s for s, _ in seed), max(e for _, e in seed)
    grew = True
    while grew:
        grew = False
        for s, e in [(r[2], r[3]) for r in ck]:
            if e >= lo - JOIN_GAP and s <= hi + JOIN_GAP:
                if s < lo - 0:
                    lo, grew = min(lo, s), True
                if e > hi:
                    hi, grew = max(hi, e), True
    inside = sorted((r for r in ck if r[3] >= lo and r[2] <= hi),
                    key=lambda r: r[2])
    want = set(r[0] for r in inside)
    with open(d + "/island_locus.tsv", "w", encoding="utf-8", newline="") as fh:
        fh.write("contig\tcontig_bp\tin_kp_bin\tin_kp_plus\tcarries_island"
                 "\tregion\truler\truler_start\truler_end\taln_bp"
                 "\tcontig_start\tcontig_end\tstrand\n")
        out = []
        for r in inside:
            out.append((r[2], r[0], "element", r))
        for r in rows:
            if r[0] in want and r[1] == KPREF:
                out.append((r[2], r[0], "kp_chromosome", r))
        for _, _, region, r in sorted(out, key=lambda x: (x[2], x[0])):
            c = r[0]
            fh.write("%s\t%d\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%s\n"
                     % (c, ln[c], "yes" if c in binc else "no",
                        "yes" if c in plus else "no",
                        "yes" if c in isl_contigs else "no",
                        region, r[1], r[2], r[3], r[7], r[4], r[5], r[6]))
    print("挿入要素 %s %d-%d（%d bp）を contig %d 本が敷き詰めている"
          % (RULER, lo, hi, hi - lo, len(want)))
    print("  うち自動の袋に入っていたもの %d 本、手で足したもの %d 本"
          % (len(want & binc), len(want & (plus - binc))))


if __name__ == "__main__":
    main()

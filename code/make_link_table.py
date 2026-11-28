"""contig ごとの読み取りの深さと、contig をまたぐ読み取り対を表にする。

元になる測定は pipeline/asm_linkage.sh（アセンブル 21,759 本すべてを 1 つの
索引にして貼り直す）。ここはその出力を公開用の 2 つの表にするだけ。

出すもの
  data/kleborate/contig_read_depth.tsv
      contig ごとの **読み取りの深さ**。CORRECTIONS.md 16 番が使っていた
      megahit の `multi` も並べて入れる。`multi` は k=119 のグラフ上の
      重複度で読み取り長は 100 なので、読み取りの深さではない。
      並べて配るのは、どれだけずれるかを読者が見られるようにするため。
  data/kleborate/contig_links.tsv
      contig の組ごとの対の数。相手は **6 種のどれかに被覆 0.5 以上で
      割り振られた contig だけ**に絞る。アセンブル 21,759 本のうち
      割り振りがあるのは 160 本で、残りは便の中の別の細菌である。
      そこへ渡る対は寄せ集めの反復配列なので、分母に入れない。
      絞る条件を先に決めてから数えている。

使い方: python3 make_link_table.py <出力先ディレクトリ>
"""
import gzip
import os
import sys

B = "/mnt/c/crc_strain_pilot"
A = B + "/asm/cancer_ERR1018186"
LINKDIR = B + "/meta/link186"
CTG = A + "/megahit/final.contigs.fa"
ASSIGN = A + "/kleborate/assign.tsv"
BIN = A + "/kleborate/kp_contigs.txt"
PLUS = A + "/kleborate/kp_plus.txt"
MIN_COV = 0.5        # 「割り振りがある」とみなす被覆。asm_kleborate.sh と同じ
CTRL_BP = 50000      # 比べるための袋の contig の長さの下限
CTRL_N = 8           # 同じ、本数


def read_lines(p):
    return [l.strip() for l in open(p, encoding="utf-8") if l.strip()]


def main():
    outdir = sys.argv[1]
    d = os.path.join(outdir, "kleborate")
    if not os.path.isdir(d):
        os.makedirs(d)

    ln, multi = {}, {}
    for line in open(CTG, encoding="utf-8"):
        if line.startswith(">"):
            q = line[1:].split()
            ln[q[0]] = int(q[-1].split("=")[1])
            # 見出しは ">名前 flag=1 multi=18.0100 len=386079"。
            # 列の位置ではなく名前で取る（最初は flag を multi として読んでいた）
            multi[q[0]] = float([t for t in q if t.startswith("multi=")][0]
                                .split("=")[1])
    sp, cov = {}, {}
    for line in open(ASSIGN, encoding="utf-8"):
        q = line.rstrip("\n").split("\t")
        sp[q[0]], cov[q[0]] = q[1], float(q[3])
    binc, plus = set(read_lines(BIN)), set(read_lines(PLUS))

    # ------------------------------------------------------ 読み取りの深さ
    dsum, dcov, dlen = {}, {}, {}
    with gzip.open(LINKDIR + "/keep.depth.tsv.gz", "rt") as fh:
        for line in fh:
            q = line.split("\t")
            c, v = q[0], int(q[2])
            dlen[c] = dlen.get(c, 0) + 1
            if v > 0:
                dcov[c] = dcov.get(c, 0) + 1
                dsum[c] = dsum.get(c, 0) + v
    bin_bp = sum(dlen.get(c, 0) for c in binc)
    bin_depth = sum(dsum.get(c, 0) for c in binc) / bin_bp
    want = sorted(plus, key=lambda c: -ln[c])
    with open(d + "/contig_read_depth.tsv", "w", encoding="utf-8",
              newline="") as fh:
        fh.write("contig\tcontig_bp\tin_kp_bin\tin_kp_plus\tmegahit_multi"
                 "\tread_depth\tbreadth\tdepth_over_bin\tmulti_over_depth\n")
        fh.write("kp_bin_all\t%d\tyes\tyes\tNA\t%.3f\tNA\t1.000\tNA\n"
                 % (bin_bp, bin_depth))
        for c in want:
            dep = dsum.get(c, 0) / dlen[c]
            fh.write("%s\t%d\t%s\tyes\t%.2f\t%.3f\t%.4f\t%.3f\t%.3f\n"
                     % (c, ln[c], "yes" if c in binc else "no", multi[c], dep,
                        dcov.get(c, 0) / dlen[c], dep / bin_depth,
                        multi[c] / dep))
    print("読み取りの深さ: 袋 %d bp の深さ %.2f、contig %d 本"
          % (bin_bp, bin_depth, len(want)))

    # 長さを揃えて比べる。短い contig は反復がたたんであって深さが跳ねるので、
    # 「袋の contig のうち同じくらいの長さのもの」と比べないと意味がない。
    def span(lo, hi, only_bin=True):
        v = sorted(dsum.get(c, 0) / dlen[c] / bin_depth
                   for c in (binc if only_bin else plus)
                   if lo <= ln[c] < hi and c in dlen)
        return (len(v), v[0], v[-1]) if v else (0, 0, 0)

    for lo, hi, nm in ((10000, 10 ** 9, "10 kb 以上"),
                       (2000, 10000, "2-10 kb"),
                       (0, 2000, "2 kb 未満")):
        n, a, b = span(lo, hi)
        print("  袋の contig %-10s %2d 本  深さは袋の %.2f-%.2f 倍"
              % (nm, n, a, b))

    # -------------------------------------------------------------- 連結
    def assigned(c):
        return c in binc or c in plus or (c in sp and cov.get(c, 0) >= MIN_COV)

    def group(c):
        if c in binc:
            return "kp_bin"
        if c in plus:
            return "added_island"
        return sp[c]

    part = {}
    for line in open(LINKDIR + "/links.tsv", encoding="utf-8"):
        q = line.rstrip("\n").split("\t")
        if len(q) < 3:
            continue
        n = int(q[2])
        part.setdefault(q[0], {})[q[1]] = n
        part.setdefault(q[1], {})[q[0]] = n
    ctrl = [c for c in sorted(binc, key=lambda x: -ln[x])
            if ln[c] >= CTRL_BP][:CTRL_N]
    rows = []
    for c in want + [x for x in ctrl if x not in want]:
        p = part.get(c, {})
        g = {}
        for q, n in p.items():
            if assigned(q) and q != c:
                g[group(q)] = g.get(group(q), 0) + n
        tot = sum(g.values())
        rows.append((c, g, tot))
    with open(d + "/contig_links.tsv", "w", encoding="utf-8", newline="") as fh:
        fh.write("contig\tcontig_bp\tin_kp_bin\tin_kp_plus\tpairs_to_kp_bin"
                 "\tpairs_to_added_island\tpairs_to_other_species"
                 "\tpairs_assigned_total\tshare_kp_bin_or_island"
                 "\tpairs_to_kp_bin_per_kb\n")
        for c, g, tot in rows:
            b = g.get("kp_bin", 0)
            a = g.get("added_island", 0)
            o = tot - b - a
            fh.write("%s\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%.4f\t%.3f\n"
                     % (c, ln[c], "yes" if c in binc else "no",
                        "yes" if c in plus else "no", b, a, o, tot,
                        (b + a) / tot if tot else 0.0, 1000.0 * b / ln[c]))
    print("連結: contig %d 本について、割り振りのある相手との対を数えた"
          % len(rows))


if __name__ == "__main__":
    main()

"""アセンブルによる担い手確認の証拠を data/ に出す。

RESULTS.md は cancer_ERR1018186 について「島と Klebsiella の染色体が
ひと続きのコンティグに載る」と書いているが、**その証拠は配布していなかった。**
この研究でいちばん強い所見が検算できない状態だった（CORRECTIONS.md の 15 番）。

出すもの
  data/assembly_alignment.tsv       … コンティグのどこが何に当たるか
  data/assembly_junction_depth.tsv  … 接合部をまたぐ 1 塩基ごとの読み取りの深さ

後者が要るのは、**アセンブラのつなぎ間違い（キメラ）でないこと**を示すため。
接合部に読み取りが切れ目なく載っていれば、つなぎ間違いではない。

使い方: python3 make_junction_table.py <出力先ディレクトリ>
"""
import io
import os
import sys

A = "/mnt/c/crc_strain_pilot/asm/cancer_ERR1018186"
CONTIG = "k119_45171"
ISLAND_FA = "/mnt/c/crc_strain_pilot/pks/ref/pks_island.fasta"
K = 50


def read_fasta(path):
    name, buf = None, []
    for line in io.open(path, encoding="utf-8", errors="replace"):
        if line.startswith(">"):
            if name:
                yield name, "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        yield name, "".join(buf)


def main():
    outdir = sys.argv[1]
    out = []

    # 1. minimap2 の結果から、このコンティグの行だけ取る
    paf = A + "/kleborate/vs_g6.paf"
    for line in io.open(paf, encoding="utf-8"):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 11 and p[0] == CONTIG:
            out.append((CONTIG, int(p[1]), int(p[2]), int(p[3]), p[4],
                        p[5].split("|")[0], int(p[7]), int(p[8]),
                        int(p[9]), int(p[10]), "minimap2"))

    # 2. 島の参照との完全一致 50-mer で、島がどこに載っているかを出す
    island = list(read_fasta(ISLAND_FA))[0][1].upper()
    target = None
    for n, s in read_fasta(A + "/kleborate/island_contigs.fa"):
        if n == CONTIG:
            target = s.upper()
            break
    idx = {}
    for i in range(len(island) - K + 1):
        idx.setdefault(island[i:i + K], i + 1)
    comp = str.maketrans("ACGT", "TGCA")
    hits = []
    for i in range(len(target) - K + 1):
        km = target[i:i + K]
        pos, strand = idx.get(km), "+"
        if pos is None:
            pos, strand = idx.get(km.translate(comp)[::-1]), "-"
        if pos:
            hits.append((i + 1, pos, strand))
    if hits:
        cs = [h[0] for h in hits]
        rs = [h[1] for h in hits]
        out.append((CONTIG, len(target), min(cs), max(cs) + K - 1,
                    hits[0][2], "AM229678.1_pks_island",
                    min(rs), max(rs) + K - 1, len(hits), len(hits),
                    "exact_%dmer" % K))

    path = outdir + "/assembly_alignment.tsv"
    with open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write("contig\tcontig_bp\tcontig_start\tcontig_end\tstrand\ttarget"
                 "\ttarget_start\ttarget_end\tmatches\taln_len\tmethod\n")
        for r in out:
            fh.write("\t".join(str(v) for v in r) + "\n")
    print("%s に %d 行" % (path, len(out)))

    # 3. 接合部の深さ
    src = A + "/junction/depth_junction.tsv"
    n = 0
    path = outdir + "/assembly_junction_depth.tsv"
    with open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write("contig\tposition\tdepth\n")
        for line in io.open(src, encoding="utf-8"):
            p = line.rstrip("\n").split("\t")
            if len(p) >= 3:
                fh.write("%s\t%s\t%s\n" % (p[0], p[1], p[2]))
                n += 1
    print("%s に %d 行" % (path, n))


if __name__ == "__main__":
    main()

"""MAPQ で絞る前と後の当たり本数を、公開用の 1 つの表にする。

元になる測定は pipeline/mapq_split.sh。CORRECTIONS.md 10 番が
「整列器は 1 本の読み取りに primary を 1 つしか出すから水増しは起きない」と
**理屈で** 直したところを、実測に置き換えるためのもの。

出すもの
  data/mapq_split.tsv
    kind=alignment   貼り付けの総数と読み取り名の種類数。等しければ
                     1 本が 2 回数えられることは無い（= 水増しは無い）
    kind=gene        clb 18 遺伝子それぞれの、絞る前と MAPQ 20 の本数
    kind=repeat_copy 重複の 2 つの複製それぞれの本数

検体は cancer_ERR1018186（読み取りが手元に残っている唯一の担体）。
整列器は bowtie2。**BWA-MEM は走らせていない**（この環境に無い）。

使い方: python3 make_mapq_table.py <出力先ディレクトリ>
"""
import io
import os
import re
import sys

SRC = "/mnt/c/crc_strain_pilot/meta/mapq186"


def main():
    outdir = sys.argv[1]
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    txt = io.open(os.path.join(SRC, "counts.txt"), encoding="utf-8").read()

    def one(pat):
        m = re.search(pat + r"\s+(\d+)", txt)
        if not m:
            raise SystemExit("counts.txt に %s が無い" % pat)
        return int(m.group(1))

    prim = one("primary")
    names = one("読み取り名の種類数")
    sec = one("secondary")
    sup = one("supplementary")

    rows = [("alignment", "primary_alignments", "NA", prim, "NA"),
            ("alignment", "distinct_read_names", "NA", names, "NA"),
            ("alignment", "secondary", "NA", sec, "NA"),
            ("alignment", "supplementary", "NA", sup, "NA")]

    for name, path, kind in (("gene", "per_gene.tsv", "gene"),
                             ("copy", "repeat_split.tsv", "repeat_copy")):
        with io.open(os.path.join(SRC, path), encoding="utf-8") as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                q = line.rstrip("\n").split("\t")
                r = dict(zip(hdr, q))
                if kind == "gene":
                    rows.append((kind, r["gene"], r["gene_bp"],
                                 int(r["reads_unfiltered"]),
                                 int(r["reads_mapq20"])))
                else:
                    rows.append((kind, "%s_%s-%s" % (r["copy"], r["start"],
                                                     r["end"]),
                                 str(int(r["end"]) - int(r["start"])),
                                 int(r["reads_unfiltered"]),
                                 int(r["reads_mapq20"])))

    path = os.path.join(outdir, "mapq_split.tsv")
    with io.open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write("kind\tname\tbp\treads_unfiltered\treads_mapq20"
                 "\tmapq20_loss\n")
        for kind, name, bp, a, b in rows:
            loss = "NA" if b == "NA" or not a else "%.4f" % (1 - b / float(a))
            fh.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (kind, name, bp, a, b, loss))

    print("書いた: %s（%d 行）" % (path, len(rows)))
    print("primary %d / 読み取り名 %d / secondary %d / supplementary %d"
          % (prim, names, sec, sup))
    if prim == names and sec == 0 and sup == 0:
        print("  → 1 本の読み取りが 2 回数えられることは無い。水増しは起きない")


if __name__ == "__main__":
    main()

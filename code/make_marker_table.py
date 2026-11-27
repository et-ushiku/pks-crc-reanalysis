"""目印に当たった読み取りの表を書き出す。data/marker_counts.tsv を作る。

なぜ配るか
  `data/per_sample.tsv` の `ecoli_pct` / `kpneumoniae_pct` は、目印に当たった
  本数を「1% あたり 100 万本中 313 本」で割って出している。この 313 は
  **大腸菌の目印 96,300 塩基で一度だけ較正した値**で、長さの違う目印には
  長さ比で伸縮させる必要がある（CORRECTIONS.md 17 番）。

  その長さ比例という仮定を検証するには、**目印に当たった生の本数**が必要で、
  これまで配っていなかった。per_sample.tsv の % からは逆算できない
  （総読み取り数で割った後の値なので、較正係数を仮定しないと戻せない）。

  混合索引なので 5 つの SGB が同時に測られている。目印の長さは 6,450 から
  282,850 塩基まで 44 倍ちがう。長さ比例が成り立っているかを見るための材料は
  この表にある。

書き出すもの
  data/marker_counts.tsv
    sample            検体
    sgb               MetaPhlAn の SGB
    species           その SGB の種
    reads_no_mapq     目印に当たった本数（MAPQ の門なし。313 を較正した条件）
    reads_mapq20      同じく MAPQ >= 20（六ゲノム索引と同じ条件）
    marker_bp         目印の総塩基（counts.tsv の 4 列目。覚えている数ではない）
    marker_fragments  目印の断片数。bowtie2-inspect で数えた分だけ。他は NA

  担い手の表（data/carrier_six_genomes.tsv）に載っている検体のうち、
  混合索引の counts.tsv があるものだけを書く。両方そろっていないと
  突き合わせられないため。

使い方: python3 make_marker_table.py <出力先ディレクトリ>
"""
import os
import sys

B = "/mnt/c/crc_strain_pilot"
NOUT = B + "/pks/entero_cn"
CARRIER = B + "/release/pks-crc-reanalysis/data/carrier_six_genomes.tsv"

SPECIES = {
    "t__SGB10068": "Escherichia coli",
    "t__SGB10115": "Klebsiella pneumoniae",
    "t__SGB10079": "Citrobacter koseri",
    "t__SGB10116": "Klebsiella aerogenes",
    "t__SGB165717": "Klebsiella aerogenes (2nd SGB)",
}
# bowtie2-inspect で数えた断片数。数えていないものは入れない。
FRAGMENTS = {"t__SGB10068": 122, "t__SGB10115": 11}
ORDER = ["t__SGB10068", "t__SGB10115", "t__SGB10079", "t__SGB10116", "t__SGB165717"]


def carrier_samples(path):
    out = []
    with open(path, encoding="utf-8") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        i = head.index("sample")
        for line in fh:
            q = line.rstrip("\n").split("\t")
            if len(q) > i and q[i] not in out:
                out.append(q[i])
    return out


def read_counts(path):
    out = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            q = line.rstrip("\n").split("\t")
            if len(q) >= 4:
                out[q[0]] = (int(q[1]), int(q[2]), int(q[3]))
    return out


def main():
    outdir = sys.argv[1]
    samples = carrier_samples(CARRIER)
    rows = []
    skipped = []
    for s in samples:
        fc = "%s/%s.counts.tsv" % (NOUT, s)
        if not os.path.exists(fc):
            skipped.append(s)
            continue
        c = read_counts(fc)
        for sgb in ORDER:
            if sgb not in c:
                continue
            n0, n20, bp = c[sgb]
            rows.append((s, sgb, SPECIES[sgb], n0, n20, bp,
                         FRAGMENTS.get(sgb, "NA")))

    path = outdir + "/marker_counts.tsv"
    with open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write("sample\tsgb\tspecies\treads_no_mapq\treads_mapq20"
                 "\tmarker_bp\tmarker_fragments\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")

    print("検体 %d 件 / 行 %d 件を書き出した" % (len(set(r[0] for r in rows)), len(rows)))
    if skipped:
        print("混合索引の counts.tsv が無く外した検体 %d 件: %s"
              % (len(skipped), ", ".join(skipped)))
        print("  （CN-CRC の先行分は大腸菌専用の索引だけで走らせたため。")
        print("    CORRECTIONS.md 7 番と同じ事情）")


if __name__ == "__main__":
    main()

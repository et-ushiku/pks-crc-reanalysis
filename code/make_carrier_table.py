"""data/carrier_six_genomes.tsv を作る。

RESULTS.md の担い手の表は、14 検体を 6 種の完全ゲノムに競わせて測った結果に
乗っているが、**その測定結果は配布データに入っていなかった**。
公開資料の目玉が検算できない状態だったので、ここで出す（CORRECTIONS.md の 15 番）。

元になるファイル（作業ツリー側）
  <検体>.g6.depth.tsv  … 6 種それぞれの 参照長 / 覆われた塩基 / 深さの総和
  <検体>.depth.tsv     … pks 島（clb 遺伝子 49,431 塩基）の 1 塩基ごとの深さ
どちらも **同じ読み取り集合**を、同じ MAPQ 20 で測ったもの。だから
「島の深さ ÷ その種のゲノムの深さ」という比が意味を持つ。

出力の列
  sample, group, target, target_bp, covered_bp, sum_depth
平均深さ = sum_depth / target_bp、覆い = covered_bp / target_bp、
比 = 島の平均深さ ÷ その種の平均深さ。すべて読み手が計算し直せる。

使い方: python3 make_carrier_table.py <出力先ディレクトリ>
"""
import glob
import os
import sys

OUT = "/mnt/c/crc_strain_pilot/pks/asm_out"
ISLAND = "pks_island_clb_genes"


def main():
    outdir = sys.argv[1]
    rows = []
    for f in sorted(glob.glob(OUT + "/*.g6.depth.tsv")):
        s = os.path.basename(f)[:-len(".g6.depth.tsv")]
        grp = "CRC" if s.startswith("cancer_") else "control"
        for line in open(f, encoding="utf-8"):
            q = line.rstrip("\n").split("\t")
            if len(q) >= 4:
                rows.append((s, grp, q[0], int(q[1]), int(q[2]), int(q[3])))
        # 島のほう（同じ読み取り集合）
        dep = OUT + "/" + s + ".depth.tsv"
        n = cov = tot = 0
        if os.path.exists(dep):
            for line in open(dep, encoding="utf-8"):
                p = line.split("\t")
                if len(p) >= 3:
                    v = int(p[2])
                    n += 1
                    tot += v
                    if v >= 1:
                        cov += 1
        rows.append((s, grp, ISLAND, n, cov, tot))

    path = outdir + "/carrier_six_genomes.tsv"
    with open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write("sample\tgroup\ttarget\ttarget_bp\tcovered_bp\tsum_depth\n")
        for r in rows:
            fh.write("%s\t%s\t%s\t%d\t%d\t%d\n" % r)
    print("%d 行を書き出した（検体 %d）"
          % (len(rows), len(set(r[0] for r in rows))))


if __name__ == "__main__":
    main()

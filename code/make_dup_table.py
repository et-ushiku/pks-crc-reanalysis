"""重複領域の検証用の表を書き出す（位置が必要なので度数分布では足りない）。

AM229678.1 には同じ 1,248 塩基が 2 か所ある。
    16,066-17,313  （clbK の中、clbK は 14,619-21,084 で 6,465 塩基）
    21,719-22,966  （clbJ の中、clbJ は 21,094-27,595 で 6,501 塩基）
互いに 99.840% 一致する。読み取りがどちらの場所のものか決まらないため
MAPQ 20 で落ち、**どれだけ深く読んでも覆いはここで止まる**。

この表は「覆われなかった塩基」を重複の内と外で分けて数えたもの。
  - 重複の内の欠けは深さにほとんどよらない（上限の性質）
  - 重複の外の欠けは深さで決まる（読み取り不足）
を確かめられる。

使い方: python3 make_dup_table.py <出力先ディレクトリ>
"""
import glob
import os
import re
import sys

B = "/mnt/c/crc_strain_pilot/pks"
BED = B + "/ref/clb_genes.bed"
DUP = [(16066, 17313), (21719, 22966)]
MIN_DEPTH = 1
COHORTS = {
    "FR-CRC": [B + "/out", B + "/out_desktop"],
    "AT-CRC": [B + "/out_at"],
    "CN-CRC": [B + "/out_cn"],
    "DE-CRC": [B + "/out_de"],
}


def main():
    outdir = sys.argv[1]
    reg = []
    for line in open(BED, encoding="utf-8"):
        q = line.split()
        if len(q) >= 4:
            reg.append((q[3], int(q[1]), int(q[2])))

    fh = open(outdir + "/missing_split.tsv", "w", encoding="utf-8", newline="")
    fh.write("cohort\tsample\tisland_mean_depth\tmissing_in_duplication"
             "\tmissing_outside_duplication\n")
    n = 0
    for coh, dirs in COHORTS.items():
        for d in dirs:
            for f in sorted(glob.glob(d + "/*.depth.tsv")):
                s = os.path.basename(f)[:-len(".depth.tsv")]
                cov, tot, nb = set(), 0, 0
                for line in open(f, encoding="utf-8"):
                    p = line.split("\t")
                    if len(p) < 3:
                        continue
                    v = int(p[2])
                    pos = int(p[1])
                    for nm, a, b in reg:
                        if a < pos <= b:
                            tot += v
                            nb += 1
                            if v >= MIN_DEPTH:
                                cov.add(pos)
                            break
                mi = mo = 0
                for nm, a, b in reg:
                    for pos in range(a + 1, b + 1):
                        if pos not in cov:
                            if any(x < pos <= y for x, y in DUP):
                                mi += 1
                            else:
                                mo += 1
                fh.write("%s\t%s\t%.3f\t%d\t%d\n"
                         % (coh, s, (tot / nb) if nb else 0.0, mi, mo))
                n += 1
    fh.close()
    print("検体 %d 件" % n)


if __name__ == "__main__":
    main()

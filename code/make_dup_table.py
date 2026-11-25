"""重複領域の検証用の表を書き出す（位置が必要なので度数分布では足りない）。

AM229678.1 には同じ配列が 2 か所ある。広がりは 2026-11-25 に測り直した。
    16,066-17,352  （clbK の中、clbK は 14,619-21,084 で 6,465 塩基）
    21,719-23,005  （clbJ の中、clbJ は 21,094-27,595 で 6,501 塩基）
長さ 1,287 塩基、互いに 99.845% 一致する。読み取りがどちらの場所のものか
決まらないため MAPQ 20 で落ち、**どれだけ深く読んでも覆いはここで止まる**。

測り方: 2 コピーは間隔 5,653 塩基で同じ向きに並ぶので、位置 i と i+5653 を
直接比べれば広がりが出る。さらに島全体の 30-mer を総当たりしたところ、
2 回以上現れる 30-mer 1,206 種のうち **1,198 本が間隔 5,653 ちょうど**で、
5,000〜6,500 に別の位置差は無い（= 挿入欠失でずれた続きは無い）。
残りは clbB-clbR 間の 8 塩基周期の VNTR である。

**以前は 1,248 塩基（16,066-17,313）と書いていた。約 39 塩基短かった。**
Auvray ら 2021 はこの反復配列を「1480 bp の相同配列」と報告しているが、
AM229678.1 では再現できない（窓一致率 70% まで緩めても 1,337 塩基）。
CORRECTIONS.md の 9 番を見ること。

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
DUP = [(16066, 17352), (21719, 23005)]
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

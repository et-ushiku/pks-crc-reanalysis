"""島の深さ ÷ 種のゲノム深さ、という比の **目盛り** を測る。

なぜ要るか
  CORRECTIONS.md 3 番で、目印遺伝子から作った比 R を「集団のうち何割が
  持つか」と読んでいたのが間違いだった、と直した。R の目盛りは測って
  いなかった（1.9 倍ずれていた）。R は以後「R≫1 ならその種では説明が
  つかない」という **片側** の推論にしか使わないことにした。

  ところが 8 番は、六ゲノムの競合索引から出した **別の** 比について
  「0.85、つまり E. coli の細胞のおよそ 85% が 1 コピー持っている」と
  書いた。これは 3 番で捨てたのと同じ読み方である。目盛りを測っていない
  比を割合として読んでいる。3 番の教訓が、兄弟分の量に適用されていない。

  ここでその目盛りを測る。R と違い、この比は 1 塩基あたりの深さの比なので
  長さの約束事や断片の端の損失は入らない。入るのは別のずれである。

測り方（先に決めてある。結果を見てから変えない）
  EXCLUDE = 100 倍。島の深さがある種のゲノム深さの 100 倍以上なら、
  その種は担体として排除する。
  **較正検体** = 候補のうち 1 種以外がすべて排除される検体。その検体では
  残った 1 種が担体であることに疑いが無いので、そこでの比の値が
  「その種のほぼ全部の細胞が 1 コピー持つ」場合の目盛りになる。
  比そのものには条件を付けない。付けると 1 の近くを選んでしまい循環する。

ずれの成分を 1 つだけ数字で出す
  島の深さは 49,431 bp 全体の平均で、そこには MAPQ 20 では決して覆われない
  重複の内側（約 2,000 bp）が含まれる。だから島の深さは数 % 低く出る。
  覆われた塩基だけで割った値も併せて出して、その成分の大きさを見せる。

出すもの
  data/island_ratio_scale.tsv

使い方: python3 make_ratio_scale.py <出力先ディレクトリ>
"""
import csv
import io
import os
import sys

ISLAND = "pks_island_clb_genes"
ECOLI = "Ecoli_IHE3034"
KP = "Kpneumoniae_HS11286"
SPECIES = (ECOLI, KP, "Ckoseri_BAA895", "Cfreundii_CFNIH1",
           "Ecloacae_ATCC13047", "Kaerogenes_KCTC2190")
EXCLUDE = 100.0   # この倍率以上なら、その種は担体として排除する


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    src = os.path.join(here, "..", "data", "carrier_six_genomes.tsv")
    outdir = sys.argv[1]
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    d, grp = {}, {}
    with io.open(src, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            tb, cb = float(r["target_bp"]), float(r["covered_bp"])
            d.setdefault(r["sample"], {})[r["target"]] = (
                float(r["sum_depth"]) / tb,
                float(r["sum_depth"]) / cb if cb else 0.0)
            grp[r["sample"]] = r["group"]

    rows = []
    for s in sorted(d):
        t = d[s]
        isl, isl_cov = t[ISLAND]
        rec = {"sample": s, "group": grp[s],
               "island_depth": isl, "island_depth_covered_only": isl_cov}
        excluded = []
        for sp in SPECIES:
            v = t[sp][0]
            rec[sp] = v
            if v > 0 and isl / v >= EXCLUDE:
                excluded.append(sp)
        rec["island_over_ecoli"] = isl / t[ECOLI][0] if t[ECOLI][0] else -1
        rec["island_over_kp"] = isl / t[KP][0] if t[KP][0] else -1
        rec["n_excluded"] = len(excluded)
        left = [sp for sp in SPECIES if sp not in excluded]
        # 較正検体は「残りが 1 種だけ」の検体。比の値は条件に使わない
        rec["calibrator_for"] = left[0] if len(left) == 1 else "NA"
        rows.append(rec)

    cols = (["sample", "group", "island_depth", "island_depth_covered_only"]
            + list(SPECIES)
            + ["island_over_ecoli", "island_over_kp", "n_excluded",
               "calibrator_for"])
    path = os.path.join(outdir, "island_ratio_scale.tsv")
    with io.open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(
                r[c] if isinstance(r[c], str) else "%.4f" % r[c]
                for c in cols) + "\n")

    cal = [r for r in rows if r["calibrator_for"] == ECOLI]
    v = sorted(r["island_over_ecoli"] for r in cal)
    print("検体 %d 件を書いた: %s" % (len(rows), path))
    print("E. coli の較正検体 %d 件（他の 5 種がすべて %.0f 倍以上で排除される）"
          % (len(cal), EXCLUDE))
    if v:
        print("  そこでの 島/E.coli は %.3f-%.3f（幅 %.1f 倍）、上端 %.3f"
              % (v[0], v[-1], v[-1] / v[0], v[-1]))
    calkp = [r for r in rows if r["calibrator_for"] == KP]
    print("K. pneumoniae の較正検体: %d 件" % len(calkp))
    if not calkp:
        print("  **ゼロ**。この比の目盛りは Klebsiella については測れていない。")
    ratio = [r["island_depth_covered_only"] / r["island_depth"] for r in rows]
    print("覆われた塩基だけで割ると島の深さは %.3f-%.3f 倍になる"
          % (min(ratio), max(ratio)))


if __name__ == "__main__":
    main()

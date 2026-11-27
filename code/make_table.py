"""手元の深さの表から、公開用の 2 つのデータファイルを書き出す。

なぜこの形にするか
  元データは 4 コホート 397 検体の fastq で、合計およそ 2.5TB ある。
  これを配ることはできないし、ENA から取り直すには数日かかる。
  一方、この研究の判定はすべて
      「clb 遺伝子 18 個の、位置ごとの深さ」
  だけから決まる。さらに間引き後の覆いの期待値は
      平均_i [ 1 - (1-f)^d_i ]
  なので、**位置の順番は要らず、遺伝子ごとの深さの度数分布さえあれば厳密に同じ値が出る**。

  だから深さの度数分布を配る。これで
    - 覆い、遺伝子数、判定
    - 任意の読み取り数に間引いた後の覆いと判定
    - 重複領域を除いた欠けの数え直し（別ファイル）
  が、fastq 無しで完全に再現できる。

書き出すもの
  data/depth_hist.tsv   検体 / 遺伝子 / 深さ / その深さだった塩基数
  data/per_sample.tsv   検体ごとの要約（覆い・遺伝子数・判定・読み取り数・大腸菌%）

使い方: python3 make_table.py <出力先ディレクトリ>
"""
import glob
import os
import re
import sys
from collections import defaultdict

B = "/mnt/c/crc_strain_pilot/pks"
BED = B + "/ref/clb_genes.bed"
MIN_DEPTH, GENE_COVER, CALL_BREADTH, CALL_GENES, SUSPECT_BREADTH = 1, 0.5, 0.5, 14, 0.1
PER_PCT = 313.0
ECOLI_MARK_BP = 96300.0   # 313 を較正したときの目印の総塩基。marker_pct を参照

COHORTS = {
    "FR-CRC": [B + "/out", B + "/out_desktop"],
    "AT-CRC": [B + "/out_at"],
    "CN-CRC": [B + "/out_cn"],
    "DE-CRC": [B + "/out_de"],
}
STUDY = {
    "FR-CRC": "Zeller 2014 (France)",
    "AT-CRC": "Feng 2015 (Austria)",
    "CN-CRC": "Yu 2017 (China)",
    "DE-CRC": "Wirbel 2019 + Zeller 2014 (Germany)",
}
ZELLER_REF = B + "/pks_vs_ecoli.txt"
EXTRA_N = "/mnt/c/crc_strain_pilot/meta/zeller_desktop_reads.tsv"


def load_bed(path):
    r = []
    for line in open(path, encoding="utf-8"):
        q = line.split()
        if len(q) >= 4:
            r.append((q[3], int(q[1]), int(q[2])))
    return r


def hist_by_gene(path, regions):
    """遺伝子ごとに {深さ: 塩基数}。深さ 0 の位置も数える（表に出ない位置は後で足す）"""
    h = defaultdict(lambda: defaultdict(int))
    seen = defaultdict(int)
    for line in open(path, encoding="utf-8"):
        q = line.split("\t")
        if len(q) < 3:
            continue
        pos, v = int(q[1]), int(q[2])
        for nm, s, e in regions:
            if s < pos <= e:
                h[nm][v] += 1
                seen[nm] += 1
                break
    for nm, s, e in regions:
        ln = e - s
        miss = ln - seen[nm]
        if miss > 0:
            h[nm][0] += miss
    return h


def call_from(h, regions, f):
    tl = tc = ng = 0
    for nm, s, e in regions:
        g = h.get(nm, {})
        ln = e - s
        if f >= 1.0:
            c = sum(n for d, n in g.items() if d >= MIN_DEPTH)
        else:
            c = sum(n * (1.0 - (1.0 - f) ** d) for d, n in g.items())
        tl += ln
        tc += c
        if c / ln >= GENE_COVER:
            ng += 1
    br = tc / tl if tl else 0.0
    cl = "present" if (br >= CALL_BREADTH and ng >= CALL_GENES) \
        else ("suspect" if br >= SUSPECT_BREADTH else "absent")
    return br, ng, cl


def total_reads(log):
    if not os.path.exists(log):
        return None
    for line in open(log, encoding="utf-8", errors="ignore"):
        m = re.match(r"\s*(\d+) reads; of these:", line)
        if m:
            return int(m.group(1))
    return None


def aligned(log):
    if not os.path.exists(log):
        return 0
    h = 0
    for line in open(log, encoding="utf-8", errors="ignore"):
        for pat in (r"(\d+) \([\d.]+%\) aligned exactly 1 time",
                    r"(\d+) \([\d.]+%\) aligned >1 times"):
            m = re.search(pat, line)
            if m:
                h += int(m.group(1))
    return h


def marker_pct(d, s, sgb):
    """目印の表から種の % を出す。

    表そのものが無い検体では None（= NA）を返す。2026-11-22 に直した。
    それ以前は hit=0 のまま計算して **0.000000 を返していた**。
    CN-CRC の先行 27 検体は大腸菌専用の索引だけで走らせていて、混合索引の
    counts.tsv が無い。そのため「測っていない」が「Klebsiella はいない」として
    表に載っていた。測っていないものは NA と書く。

    **換算係数は目印の総塩基数に比例させる。2026-11-26 に直した。**
    PER_PCT = 313 は「大腸菌 1% あたり 100 万本中 313 本」で、*大腸菌の目印
    96,300 塩基*で一度だけ較正した値である。目印に当たる本数は
        量(%) × 目印の総塩基
    に比例するので、長さの違う目印にそのまま 313 を当ててはいけない。
    K. pneumoniae の目印は **6,450 塩基しかなく**、同じ 313 を使うと
    96,300 / 6,450 = **14.93 倍小さい値**が出る。長さ比で較正を伸ばす。

    同じ長さ比の補正を meta/cn_entero.py が最初から入れていて、そちらが出した
    6.83% が CORRECTIONS.md の 8 番に載っている。公開した per_sample.tsv の列
    だけが補正されていなかった。

    注意: これは「読み取りの集まり方が目印の長さに比例する」という仮定である。
    cn_entero.py は「Klebsiella では検証できていない仮定」と明記していた。
    **2026-11-28 に検証した**（code/verify_marker_scaling.py）。六ゲノム競合索引で
    同じ検体を測り、目印の密度とゲノムの密度の倍率が種で同じかを見ると
    kp/ec = 0.929（95% 区間 0.662-1.304）で、伸縮させない説に必要な 14.93 を
    11 倍の差で排除する。区間の幅が 2 倍あるので、この列は倍率 2 倍以内の
    粗さまでしか信用できない。
    """
    fc = d + "/" + s + ".counts.tsv"
    if not os.path.exists(fc):
        return None
    lg = d + "/" + s + ".bowtie2.log"
    t = total_reads(lg)
    if not t:
        return None
    hit = 0
    mark_bp = None
    for line in open(fc, encoding="utf-8"):
        q = line.rstrip().split("	")
        if len(q) >= 2 and q[0] == sgb:
            hit = int(q[1])
            # 4 列目が目印の総塩基。覚えている数ではなく表から読む。
            if len(q) >= 4:
                mark_bp = float(q[3])
    if mark_bp is None:
        raise SystemExit("目印の塩基数が表に無い: %s %s" % (fc, sgb))
    return (hit / t * 1e6) / (PER_PCT * mark_bp / ECOLI_MARK_BP)


def zeller_ecoli():
    out = {}
    for line in open(ZELLER_REF, encoding="utf-8"):
        p = line.split()
        if len(p) >= 5 and re.match(r"^(cancer|n)_", p[0]):
            try:
                out[p[0]] = float(p[-1])
            except ValueError:
                pass
    return out


def ecoli_of(coh, s, zref):
    if coh == "FR-CRC":
        return zref.get(s)
    if coh == "DE-CRC":
        return marker_pct(B + "/entero_de", s, "t__SGB10068")
    if coh == "AT-CRC":
        lg = B + "/ecoli_at/" + s + ".bowtie2.log"
        t = total_reads(lg)
        return (aligned(lg) / t * 1e6) / PER_PCT if t else None
    lg = B + "/ecoli_cn/" + s + ".bowtie2.log"
    t = total_reads(lg)
    if t:
        return (aligned(lg) / t * 1e6) / PER_PCT
    return marker_pct(B + "/entero_cn", s, "t__SGB10068")


def kp_of(coh, s):
    if coh == "DE-CRC":
        return marker_pct(B + "/entero_de", s, "t__SGB10115")
    if coh == "CN-CRC":
        return marker_pct(B + "/entero_cn", s, "t__SGB10115")
    return None


def main():
    outdir = sys.argv[1]
    regions = load_bed(BED)
    zref = zeller_ecoli()
    extra = {}
    if os.path.exists(EXTRA_N):
        for line in open(EXTRA_N, encoding="utf-8"):
            p = line.split("\t")
            if len(p) >= 2:
                extra[p[0]] = int(p[1])

    fh_h = open(outdir + "/depth_hist.tsv", "w", encoding="utf-8", newline="")
    fh_h.write("sample\tgene\tdepth\tbases\n")
    fh_s = open(outdir + "/per_sample.tsv", "w", encoding="utf-8", newline="")
    fh_s.write("cohort\tstudy\tsample\tgroup\ttotal_reads\tisland_reads"
               "\tisland_breadth\tgenes_called\tcall\tecoli_pct"
               "\tkpneumoniae_pct\n")
    n = 0
    for coh, dirs in COHORTS.items():
        for d in dirs:
            for f in sorted(glob.glob(d + "/*.depth.tsv")):
                s = os.path.basename(f)[:-len(".depth.tsv")]
                h = hist_by_gene(f, regions)
                br, ng, cl = call_from(h, regions, 1.0)
                N = total_reads(d + "/" + s + ".bowtie2.log") or extra.get(s)
                # 参照配列に貼り付いた本数（MAPQ で絞る前）。PITFALLS 項目 2 用
                nr = aligned(d + "/" + s + ".bowtie2.log")
                eco = ecoli_of(coh, s, zref)
                kp = kp_of(coh, s)
                for nm, a, b in regions:
                    for dep in sorted(h.get(nm, {})):
                        fh_h.write("%s\t%s\t%d\t%d\n" % (s, nm, dep, h[nm][dep]))
                fh_s.write("%s\t%s\t%s\t%s\t%s\t%d\t%.6f\t%d\t%s\t%s\t%s\n"
                           % (coh, STUDY[coh], s,
                              "CRC" if s.startswith("cancer_") else "control",
                              N if N else "NA", nr, br, ng, cl,
                              ("%.6f" % eco) if eco is not None else "NA",
                              ("%.6f" % kp) if kp is not None else "NA"))
                n += 1
    fh_h.close()
    fh_s.close()
    print("検体 %d 件を書き出した" % n)


if __name__ == "__main__":
    main()

"""DE-CRC の検体ごとに「島を担いでいるのはどの種か」を密度比で見る。

考え方
  島は 49,455 bp。目印配列は種ごとに長さが分かっている。
  同じ検体の中で
      島の読み取り密度 = 島の読み取り数 / 49,455
      種の読み取り密度 = 目印の読み取り数 / 目印の長さ
  の比 R を取る。
      R ≈ 1  → その種の集団のほぼ全員が島を 1 本持っている（担い手として整合）
      R < 1  → その種の一部の株だけが持っている
      R ≫ 1  → **その種では島の量を説明できない**（担い手は別）

限界（これは証明ではない）
  密度比は「量の整合」でしかない。ERR1018186 で学んだとおり、
  担い手の確定には**島と骨格が同じ分子に載っていること**を示す必要がある
  （アセンブル＋接合部の読み取り）。ここで出るのは候補の絞り込みまで。
  目印 t__SGB10115 は 6,450 bp しかなく、読み取り数が少ないと比は荒れる。

使い方: python3 de_carrier.py
"""
import glob
import os
import re
from collections import defaultdict

B = "/mnt/c/crc_strain_pilot"
BED = B + "/pks/ref/clb_genes.bed"
OUT = B + "/pks/out_de"
NOUT = B + "/pks/entero_de"
ISLAND_BP = 49455
MIN_DEPTH = 1
GENE_COVER = 0.5
CALL_BREADTH = 0.5
CALL_GENES = 14
SUSPECT_BREADTH = 0.1

NAME = {
    "t__SGB10068": "大腸菌",
    "t__SGB10079": "大腸菌の近縁",
    "t__SGB10115": "K. pneumoniae",
    "t__SGB10116": "Klebsiella の近縁",
    "t__SGB165717": "Enterobacter 近縁",
}


def load_bed(path):
    out = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            p = line.split()
            if len(p) >= 4:
                out.append((p[3], int(p[1]), int(p[2])))
    return out


def gene_of(pos, regions):
    for name, s, e in regions:
        if s < pos <= e:
            return name
    return None


def breadth_of(path, regions):
    cov = defaultdict(int)
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            p = line.split("\t")
            if len(p) < 3:
                continue
            if int(p[2]) >= MIN_DEPTH:
                g = gene_of(int(p[1]), regions)
                if g:
                    cov[g] += 1
    tl = tc = ng = 0
    for name, s, e in regions:
        ln = e - s
        c = cov.get(name, 0)
        tl += ln
        tc += c
        if c / ln >= GENE_COVER:
            ng += 1
    return (tc / tl if tl else 0.0), ng


def call(b, n):
    if b >= CALL_BREADTH and n >= CALL_GENES:
        return "あり"
    if b >= SUSPECT_BREADTH:
        return "疑い"
    return "なし"


def markers(s):
    """目印ごとに (MAPQ無しの数, 目印の長さ)"""
    f = NOUT + "/" + s + ".counts.tsv"
    out = {}
    if not os.path.exists(f):
        return out
    with open(f, encoding="utf-8") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 4:
                out[p[0]] = (int(p[1]), int(p[3]))
    return out


def pks_reads(s):
    m = None
    lg = B + "/meta/de_map.log"
    with open(lg, encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith("--- 終了") and line.split()[2] == s:
                mm = re.search(r"pks=(\d+)", line)
                if mm:
                    m = int(mm.group(1))
    return m


def main():
    regions = load_bed(BED)
    print("DE-CRC: 島の担い手の候補を密度比で絞る")
    print("  R = (島の読み取り/49,455) / (目印の読み取り/目印の長さ)")
    print("  R≈1 でその種が担い手として整合、R≫1 でその種では説明できない")
    print()

    rows = []
    for dep in sorted(glob.glob(OUT + "/*.depth.tsv")):
        s = os.path.basename(dep)[:-len(".depth.tsv")]
        br, ng = breadth_of(dep, regions)
        c = call(br, ng)
        if c == "なし":
            continue
        pk = pks_reads(s)
        if not pk:
            continue
        mk = markers(s)
        dens_i = pk / ISLAND_BP
        best = None
        rs = {}
        for k, (n, L) in mk.items():
            if n <= 0:
                continue
            r = dens_i / (n / L)
            rs[k] = (r, n)
            if best is None or abs(r - 1) < abs(best[1] - 1):
                best = (k, r, n)
        rows.append((s, c, br, ng, pk, rs, best))

    print("%-24s %-4s %6s %5s %8s  %s"
          % ("検体", "判定", "覆い", "遺伝", "島の読み", "いちばん整合する種（R / 目印の読み）"))
    for s, c, br, ng, pk, rs, best in sorted(rows, key=lambda x: -x[2]):
        b = ("%s R=%.2f (%d本)" % (NAME.get(best[0], best[0]), best[1], best[2])
             ) if best else "目印が無い"
        print("%-24s %-4s %6.3f %5d %8d  %s" % (s, c, br, ng, pk, b))

    print()
    print("大腸菌では説明できない検体（大腸菌の R が 10 を超える）")
    hit = []
    for s, c, br, ng, pk, rs, best in rows:
        if "t__SGB10068" in rs and rs["t__SGB10068"][0] > 10:
            hit.append((s, c, br, pk, rs))
    if not hit:
        print("  なし")
    for s, c, br, pk, rs in sorted(hit, key=lambda x: -x[2]):
        print("  %-24s 判定 %s  覆い %.3f  島 %d 本" % (s, c, br, pk))
        for k in ("t__SGB10068", "t__SGB10115", "t__SGB10116"):
            if k in rs:
                print("      %-18s R=%8.2f  目印 %5d 本"
                      % (NAME.get(k, k), rs[k][0], rs[k][1]))
            else:
                print("      %-18s 0 本（この種は居ない）" % NAME.get(k, k))
        print("      → 大腸菌では量が合わない。候補は Klebsiella 側")

    print()
    print("次にやること: 候補になった検体をアセンブルして、")
    print("  島と骨格が同じ contig に載っているか（接合部の読み取り）を確かめる。")
    print("  密度比だけでは担い手は決まらない ― ERR1018186 で確認した手順を踏む。")


if __name__ == "__main__":
    main()

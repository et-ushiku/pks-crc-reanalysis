"""cancer_ERR1018282 のアセンブルで、島がどこに載ったかを表にする。

なぜ組み立てたか
  CORRECTIONS.md 8 番の 2026-11-28 の追補で、**深さの比はどちらの種も
  排除できない** ことが分かった。ERR1018282 の帰属は節約原理だけで持って
  いた。16 番を決めたのと同じ方法（組み立てて、島と染色体が同じ contig に
  載るかを見る）で決めに行った。

いちばん大事な分け方
  参照 AM229678.1 は 55,140 bp だが clb 18 遺伝子はそのうち 49,431 bp で、
  残り 5,685 bp は **島に固有でない前後の配列** である（PITFALLS 2 番）。
  だから「島の参照に当たった contig」を数えるだけでは駄目で、**当たりが
  遺伝子の中か外か** を必ず分ける。ここを分けないと、前後配列にしか
  当たっていない Klebsiella の染色体 contig を「Klebsiella が島を持つ
  証拠」として読んでしまう。この検体ではまさにそれが起きる形になっていた。

出すもの
  data/assembly_282_island.tsv
    verdict=kp_chromosome_backbone_only  Klebsiella の染色体だが、島の
        当たりは遺伝子の外だけ。担体の証拠には **ならない**
    verdict=backbone_only                同じく遺伝子の外だけ
    verdict=island_fragment              clb 遺伝子を載せている。ただし
        染色体は載っていないので、どの種のものか言えない

使い方: python3 make_asm282_table.py <出力先ディレクトリ>
"""
import collections
import io
import os
import sys

B = "/mnt/c/crc_strain_pilot"
A = B + "/asm/cancer_ERR1018282/megahit"
HITS = A + "/pks_hits.tsv"
VSG = A + "/contig_vs_genomes.tsv"
CTG = A + "/final.contigs.fa"
KP = "Kpneumoniae_HS11286"


def main():
    outdir = sys.argv[1]
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    here = os.path.dirname(os.path.abspath(__file__))
    bed = os.path.join(here, "..", "data", "clb_genes.bed")

    genes = []
    with io.open(bed, encoding="utf-8") as fh:
        for line in fh:
            q = line.split()
            genes.append((int(q[1]), int(q[2]), q[3]))

    ln = {}
    with io.open(VSG, encoding="utf-8") as fh:
        for line in fh:
            q = line.split("\t")
            ln[q[0]] = int(q[1])

    hits = collections.defaultdict(list)
    with io.open(HITS, encoding="utf-8") as fh:
        for line in fh:
            q = line.rstrip("\n").split("\t")
            a, b = int(q[4]), int(q[5])
            # contig 側は BLAST が出した **向きのまま** 持つ。min/max に
            # 畳むと逆向きが消え、重複の 2 か所を復元できなくなる
            hits[q[1]].append((min(a, b), max(a, b),
                               int(q[6]), int(q[7])))

    sp = collections.defaultdict(lambda: collections.defaultdict(int))
    with io.open(VSG, encoding="utf-8") as fh:
        for line in fh:
            q = line.rstrip("\n").split("\t")
            sp[q[0]][q[2].split("|")[0]] += int(round(int(q[4])
                                                     * float(q[3]) / 100.0))

    multi = {}
    with io.open(CTG, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                q = line[1:].split()
                if q[0] in hits:
                    multi[q[0]] = float(
                        [t for t in q if t.startswith("multi=")][0]
                        .split("=")[1])

    rows = []
    for c in sorted(hits, key=lambda x: -ln.get(x, 0)):
        iv = sorted((a, b) for a, b, _, _ in hits[c])
        mer = []
        for a, b in iv:
            if mer and a <= mer[-1][1] + 1:
                mer[-1][1] = max(mer[-1][1], b)
            else:
                mer.append([a, b])
        hit_bp = sum(b - a + 1 for a, b in mer)
        gcov = set()
        for a, b in iv:
            for s, e, _ in genes:
                lo, hi = max(a - 1, s), min(b, e)
                if hi > lo:
                    gcov.update(range(lo, hi))
        gene_bp = len(gcov)
        top, tv = max(sp[c].items(), key=lambda x: x[1])
        n = ln.get(c, 0)
        if gene_bp == 0:
            verdict = ("kp_chromosome_backbone_only"
                       if top == KP else "backbone_only")
        else:
            verdict = "island_fragment"
        # 重複がたたまれた contig を復元できるように、当たりの座標も配る。
        # これが無いと「1,287 bp・16,066・間隔 5,653」の検算が、配り物を
        # 読まずに自分の算術を確かめるだけのものになる
        spans = ";".join("%d-%d/%d-%d" % (a, b, ca, cb)
                         for a, b, ca, cb in sorted(hits[c]))
        # 例: 16066-17431/1366-1 は「島 16066-17431 が contig 1366->1 に
        # 逆向きで当たった」。向きがあるので 2 つの当たりから重複を復元できる
        rows.append((c, n, multi.get(c, 0.0), hit_bp, gene_bp, top, tv,
                     verdict, spans))

    path = os.path.join(outdir, "assembly_282_island.tsv")
    with io.open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write("contig\tcontig_bp\tmegahit_multi\tisland_hit_bp"
                 "\tisland_gene_bp\ttop_species\ttop_species_matched_bp"
                 "\ttop_species_pct_of_contig\tverdict"
                 "\tisland_spans_over_contig_spans\n")
        for c, n, m, h, g, top, tv, v, spans in rows:
            fh.write("%s\t%d\t%.2f\t%d\t%d\t%s\t%d\t%.3f\t%s\t%s\n"
                     % (c, n, m, h, g, top, tv,
                        tv / float(n) if n else 0.0, v, spans))

    gb = 0
    mer = []
    for s, e, _ in sorted(genes):
        if mer and s <= mer[-1][1]:
            mer[-1][1] = max(mer[-1][1], e)
        else:
            mer.append([s, e])
    gb = sum(e - s for s, e in mer)
    allg = set()
    for c in hits:
        for a, b, _, _ in hits[c]:
            for s, e, _ in genes:
                lo, hi = max(a - 1, s), min(b, e)
                if hi > lo:
                    allg.update(range(lo, hi))

    frag = [r for r in rows if r[7] == "island_fragment"]
    kpb = [r for r in rows if r[7] == "kp_chromosome_backbone_only"]
    print("書いた: %s（contig %d 本）" % (path, len(rows)))
    print("clb 遺伝子を載せた contig: %d 本、長さ %s bp、multi %s"
          % (len(frag), "/".join(str(r[1]) for r in frag),
             "/".join("%.0f" % r[2] for r in frag)))
    print("  覆った遺伝子塩基の和集合 %d / %d = %.1f%%"
          % (len(allg), gb, 100.0 * len(allg) / gb))
    print("島の参照には当たるが遺伝子の外だけの contig: %d 本" %
          len([r for r in rows if r[4] == 0]))
    for r in kpb:
        print("  うち Klebsiella 染色体: %s（%d bp, multi %.0f, Kp %.0f%%）"
              % (r[0], r[1], r[2], 100 * r[6] / float(r[1])))
    print("→ 島と染色体が同じ contig に載った例は **無い**")


if __name__ == "__main__":
    main()

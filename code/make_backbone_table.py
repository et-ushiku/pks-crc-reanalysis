"""島を載せた contig について、島以外の部分がどの種に当たるか、と contig の深さ。

なぜ要るか
  contig の種への割り振りは 6 種への一致塩基の多寡で決めている。ところが
  *C. koseri* BAA-895 自身が clb 島を持つ（既報。Putze 2009 の代表株）ので、
  島の配列が *C. koseri* の点数を押し上げる。**可動性要素を含めて種を決めては
  いけない**という話になる。だから島の内と外を分けて数える。

  分けて数えた結果は、こちらの作業ノートに書いてあった結論と**合わなかった**。
  ノートは「島以外の部分は *K. pneumoniae* しか当たらない（*C. koseri* は 0）」と
  書いていたが、k119_70945 の骨格は *C. koseri* が全長を 92.8% で覆い、
  *K. pneumoniae* は 45.6 kb を 81.8% でしか覆わない。
  **似ている度合いでは、この contig は Klebsiella のものだと言えない。**
  CORRECTIONS.md の 16 番に書いた。

  言えるのは別の筋からである。だからこの表は contig の**深さ**も出す。
  megahit が contig ごとに出す `multi=` は k-mer の被覆深さで、同じゲノムから
  来た contig は同じ値のあたりに揃う。*C. koseri* はこの検体にほぼ居ない
  （ゲノムの覆い 0.014）。したがって、Klebsiella のゲノムと同じ深さに居て、
  *C. koseri* が居ないなら、担ぎ手は Klebsiella 以外に残らない。

やり方
  1. 島の参照を 40-mer に刻み、contig 上で完全一致する位置を拾う。
     **一致位置の最初と最後で挟んではいけない。** そうすると、島とほとんど
     関係のない contig（k119_25144 は一致が数百 bp ずつ 4 箇所に散るだけ）で
     「島 39,448 bp」と出てしまう。隙間が 200 bp を超えたら切り、
     **連続した塊**として扱う。
  2. PAF の各行の当たった区間が島の塊と重なる割合で、一致塩基を按分する。

出すもの
  data/kleborate/island_contig_backbone.tsv
    contig, contig_bp, kmer_depth, island_bp, island_blocks, largest_block,
    region, species, matched_bp, aligned_bp
  data/kleborate/contig_depth.tsv
    bin 全体の深さの分布（島の contig をこれと比べるため）

使い方:
  python3 make_backbone_table.py <島 FASTA> <contig FASTA> <PAF> \
      <megahit final.contigs.fa> <kp_contigs.txt> <出力先 data/>
"""
import collections
import io
import os
import re
import sys

K = 40
GAP = 200
COMP = str.maketrans("ACGT", "TGCA")


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


def megahit_depths(path):
    """megahit のヘッダから contig ごとの len と multi を読む。"""
    out = {}
    for line in io.open(path, encoding="utf-8", errors="replace"):
        if not line.startswith(">"):
            continue
        m = re.match(r">(\S+)\s+flag=\d+\s+multi=([\d.]+)\s+len=(\d+)", line)
        if m:
            out[m.group(1)] = (int(m.group(3)), float(m.group(2)))
    return out


def island_blocks(seq, iset):
    hits = [i for i in range(len(seq) - K + 1) if seq[i:i + K] in iset]
    if not hits:
        return []
    runs, a, prev = [], hits[0], hits[0]
    for i in hits[1:]:
        if i - prev > GAP:
            runs.append((a, prev + K))
            a = i
        prev = i
    runs.append((a, prev + K))
    return runs


def main():
    island_fa, contig_fa, paf, allctg, binlist, data = sys.argv[1:7]
    out = os.path.join(data, "kleborate")
    if not os.path.isdir(out):
        os.makedirs(out)

    island = list(read_fasta(island_fa))[0][1].upper()
    iset = set()
    for i in range(len(island) - K + 1):
        km = island[i:i + K]
        if set(km) <= set("ACGT"):
            iset.add(km)
            iset.add(km.translate(COMP)[::-1])

    depth = megahit_depths(allctg)

    blocks, length = {}, {}
    for name, seq in read_fasta(contig_fa):
        seq = seq.upper()
        length[name] = len(seq)
        runs = island_blocks(seq, iset)
        if not runs:
            print("  %-14s %7d bp  島の 40-mer は一つも無い" % (name, len(seq)))
            continue
        blocks[name] = runs
        big = max(runs, key=lambda r: r[1] - r[0])
        print("  %-14s %7d bp  深さ %6.2f  島 %6d bp（塊 %d、最大 %d-%d）"
              % (name, len(seq), depth.get(name, (0, 0.0))[1],
                 sum(b - a for a, b in runs), len(runs), big[0] + 1, big[1]))

    agg = collections.defaultdict(lambda: [0, 0])
    for line in io.open(paf, encoding="utf-8"):
        p = line.rstrip("\n").split("\t")
        if len(p) < 11 or p[0] not in blocks:
            continue
        c, a, b = p[0], int(p[2]), int(p[3])
        sp = p[5].split("|")[0]
        mat, aln = int(p[9]), b - a
        ov = sum(max(0, min(b, hi) - max(a, lo)) for lo, hi in blocks[c])
        frac = (ov / float(aln)) if aln else 0.0
        for region, share in (("island", frac), ("backbone", 1.0 - frac)):
            if share <= 0:
                continue
            k = (c, region, sp)
            agg[k][0] += int(round(mat * share))
            agg[k][1] += int(round(aln * share))

    path = os.path.join(out, "island_contig_backbone.tsv")
    with open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write("contig\tcontig_bp\tkmer_depth\tisland_bp\tisland_blocks"
                 "\tlargest_block\tregion\tspecies\tmatched_bp\taligned_bp\n")
        for (c, region, sp) in sorted(agg, key=lambda k: (-length[k[0]], k[1], k[2])):
            r = blocks[c]
            big = max(r, key=lambda x: x[1] - x[0])
            m, al = agg[(c, region, sp)]
            fh.write("%s\t%d\t%.2f\t%d\t%d\t%d-%d\t%s\t%s\t%d\t%d\n"
                     % (c, length[c], depth.get(c, (0, 0.0))[1],
                        sum(b - x for x, b in r), len(r), big[0] + 1, big[1],
                        region, sp, m, al))
    print("%s に %d 行" % (path, len(agg)))

    # bin 全体の深さ。島の contig をこれと比べるために出す。
    names = set(l.split()[0] for l in io.open(binlist, encoding="utf-8")
                if l.strip())
    tl = tw = 0
    big = []
    for n in names:
        if n not in depth:
            continue
        ln, mu = depth[n]
        tl += ln
        tw += ln * mu
        if ln >= 10000:
            big.append(mu)
    big.sort()
    path = os.path.join(out, "contig_depth.tsv")
    with open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write("what\tvalue\n")
        fh.write("bin_contigs\t%d\n" % len(names))
        fh.write("bin_bp\t%d\n" % tl)
        fh.write("bin_depth_length_weighted_mean\t%.2f\n" % (tw / tl))
        fh.write("bin_depth_median_contigs_over_10kb\t%.2f\n" % big[len(big) // 2])
        fh.write("bin_contigs_over_10kb\t%d\n" % len(big))
        for c in sorted(blocks, key=lambda x: -length[x]):
            fh.write("island_contig_%s_depth\t%.2f\n" % (c, depth[c][1]))
            fh.write("island_contig_%s_bp\t%d\n" % (c, length[c]))
    print("%s を書いた（bin の深さ 加重平均 %.2f、10kb 以上の中央値 %.2f）"
          % (path, tw / tl, big[len(big) // 2]))


if __name__ == "__main__":
    main()

"""6 種の参照ゲノムのうち、どれが clb 島を自分で持っているかを測る。

なぜ要るか
  アセンブルした contig を 6 種に当てて種へ割り振ると、clb の載った contig は
  *C. koseri* BAA-895 に行き、*K. pneumoniae* HS11286 には一本も行かなかった。
  これを「種同定が外れた」と読むか「HS11286 に clb 島が無いのだから
  他へ行くしかない」と読むかで、意味が正反対になる。
  だから**参照ゲノムの側**を測る。contig は関係ない。

やり方
  島の参照（AM229678.1、55,140 bp）を 40-mer に刻み、各ゲノムの 40-mer 集合に
  あるかどうかを見る。両鎖をまとめるため正規化（自分と相補鎖の小さいほう）する。
  完全一致なので、閾値も当てはめ方も無い。

出すもの
  data/clb_in_six_genomes.tsv
    genome, genome_bp, island_kmers, shared_kmers, clb_genes_at_half
  最後の列は「その遺伝子の 40-mer の半分以上が共有されている」遺伝子の数。

使い方: python3 make_clb_in_genomes.py <島 FASTA> <clb_genes.bed> <6種 FASTA> <出力先 data/>
"""
import io
import os
import sys

K = 40
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


def canon(s):
    r = s.translate(COMP)[::-1]
    return s if s < r else r


def kmers(seq):
    ok = set("ACGT")
    for i in range(len(seq) - K + 1):
        km = seq[i:i + K]
        if set(km) <= ok:
            yield i, canon(km)


def main():
    island_fa, bed, g6, data = sys.argv[1:5]

    island = list(read_fasta(island_fa))[0][1].upper()
    n_pos = len(island) - K + 1
    pos_of = {}
    for i, km in kmers(island):
        pos_of.setdefault(km, []).append(i)
    print("島 %d bp、%d-mer の位置 %d 箇所" % (len(island), K, n_pos))

    genes = []
    for line in io.open(bed, encoding="utf-8"):
        p = line.split()
        if len(p) >= 4:
            genes.append((int(p[1]), int(p[2]), p[3]))   # 0-based 開始, 終わり
    genes.sort()

    rows = []
    for name, seq in read_fasta(g6):
        seq = seq.upper()
        sp = name.split("|")[0]
        gset = set(km for _, km in kmers(seq))
        found = bytearray(n_pos)
        shared = 0
        for km, ps in pos_of.items():
            if km in gset:
                shared += len(ps)
                for i in ps:
                    found[i] = 1
        half = 0
        for a, b, g in genes:
            n = (b - a) - K + 1
            if n <= 0:
                continue
            c = sum(found[a:a + n])
            if c >= n / 2.0:
                half += 1
        rows.append((sp, len(seq), n_pos, shared, half, len(genes)))
        print("  %-22s %7.2f%% 共有、clb 遺伝子 %d/%d"
              % (sp, 100.0 * shared / n_pos, half, len(genes)))

    path = os.path.join(data, "clb_in_six_genomes.tsv")
    with open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write("genome\tgenome_bp\tisland_kmers\tshared_kmers"
                 "\tclb_genes_at_half\tclb_genes_total\n")
        for r in sorted(rows, key=lambda x: -x[3]):
            fh.write("%s\t%d\t%d\t%d\t%d\t%d\n" % r)
    print("%s に %d 行（k=%d）" % (path, len(rows), K))


if __name__ == "__main__":
    main()

"""Kleborate の出力と、その判定が成り立つ根拠を data/kleborate/ に出す。

RESULTS.md は cancer_ERR1018186 の担い手を「*K. pneumoniae*、ST23/K1/ICEKp10/*clb*」と
書いていたが、**Kleborate の出力を配っていなかった**（CORRECTIONS.md の 15 番）。
出してみると、その 1 行が**2 回の別々の実行を混ぜたもの**だと分かった
（CORRECTIONS.md の 16 番）。だから生の出力だけでなく、

  - どの contig がどの袋に入ったか
  - なぜそう入ったか（contig ごと・種ごとの一致塩基）
  - 6 種のうちどれが clb 島を自分で持っているか

も一緒に出す。最後のものが要るのは、clb の載った contig が自動では
*C. koseri* に割り振られた理由が、**K. pneumoniae HS11286 が clb 島を
まったく持っていないから**であって、種を間違えたからではないことを
読み手が自分で確かめられるようにするため。

出すもの
  data/kleborate/kp_bin_output.tsv          … 自動で選んだ袋の Kleborate 出力
  data/kleborate/kp_plus_output.tsv         … そこに島の contig を足した袋の出力
  data/kleborate/kp_bin_hAMRonization.tsv   … 同・薬剤耐性の内訳
  data/kleborate/kp_plus_hAMRonization.tsv
  data/kleborate/contig_assignment.tsv      … contig ごとの勝った種と被覆
  data/kleborate/contig_species_matches.tsv … contig ごと・種ごとの一致塩基
  data/kleborate/bin_membership.tsv         … どの contig がどの袋に入ったか

使い方: python3 make_kleborate_table.py <kleborate ディレクトリ> <出力先 data/>
"""
import collections
import io
import os
import shutil
import sys


def read_names(path):
    """FASTA の名前だけ、または 1 行 1 名のファイルを読む。

    FASTA では `>` の行だけを拾う。megahit の出力は配列を折り返さないので、
    行を区別せずに拾うと配列の行まで名前として数えてしまう。
    """
    lines = [l.strip() for l in io.open(path, encoding="utf-8", errors="replace")]
    lines = [l for l in lines if l]
    if any(l.startswith(">") for l in lines):
        return [l[1:].split()[0] for l in lines if l.startswith(">")]
    return [l.split()[0] for l in lines]


def main():
    src, data = sys.argv[1], sys.argv[2]
    out = os.path.join(data, "kleborate")
    if not os.path.isdir(out):
        os.makedirs(out)

    # 1. Kleborate の生の出力をそのまま置く（列名も値も一切加工しない）
    for tgt in ("kp_bin", "kp_plus"):
        for kind, dst in (("klebsiella_pneumo_complex_output.txt",
                           "%s_output.tsv" % tgt),
                          ("klebsiella_pneumo_complex_hAMRonization_output.txt",
                           "%s_hAMRonization.tsv" % tgt)):
            p = os.path.join(src, "out_" + tgt, kind)
            if os.path.exists(p):
                shutil.copyfile(p, os.path.join(out, dst))
                print("  置いた: kleborate/%s" % dst)

    # 2. contig ごとの割り振り（勝った種と、その種に対する被覆）
    rows = []
    for line in io.open(os.path.join(src, "assign.tsv"), encoding="utf-8"):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 4:
            rows.append(p[:4])
    p = os.path.join(out, "contig_assignment.tsv")
    with open(p, "w", encoding="utf-8", newline="") as fh:
        fh.write("contig\tbest_species\tcontig_bp\tcoverage_of_contig\n")
        for r in sorted(rows, key=lambda x: -int(x[2])):
            fh.write("\t".join(r) + "\n")
    print("  kleborate/contig_assignment.tsv に %d 行" % len(rows))

    # 3. その割り振りの根拠。contig ごと・種ごとに一致塩基と当たった長さを合計する。
    #    1 本の contig は複数の断片に分かれて当たるので、合計してからでないと
    #    比べられない（asm_kleborate.sh の注記と同じ理由）。
    mat = collections.defaultdict(int)
    aln = collections.defaultdict(int)
    length = {}
    for line in io.open(os.path.join(src, "vs_g6.paf"), encoding="utf-8"):
        p = line.rstrip("\n").split("\t")
        if len(p) < 11:
            continue
        sp = p[5].split("|")[0]
        length[p[0]] = int(p[1])
        mat[(p[0], sp)] += int(p[9])
        aln[(p[0], sp)] += int(p[3]) - int(p[2])
    path = os.path.join(out, "contig_species_matches.tsv")
    with open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write("contig\tcontig_bp\tspecies\tmatched_bp\taligned_bp\n")
        for (c, sp) in sorted(mat, key=lambda k: (-length[k[0]], k[1])):
            fh.write("%s\t%d\t%s\t%d\t%d\n"
                     % (c, length[c], sp, mat[(c, sp)], aln[(c, sp)]))
    print("  kleborate/contig_species_matches.tsv に %d 行" % len(mat))

    # 4. どの contig がどの袋に入ったか。kp_plus のつくり方そのものが
    #    これで検算できる（= kp_bin ∪ 島が載った contig）。
    bin_c = set(read_names(os.path.join(src, "kp_contigs.txt")))
    plus_c = set(read_names(os.path.join(src, "kp_plus.txt")))
    isl_c = set(read_names(os.path.join(src, "island_contigs.fa")))
    path = os.path.join(out, "bin_membership.tsv")
    with open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write("contig\tcontig_bp\tin_kp_bin\tin_kp_plus\tcarries_island\n")
        for c in sorted(bin_c | plus_c | isl_c,
                        key=lambda x: -length.get(x, 0)):
            fh.write("%s\t%d\t%s\t%s\t%s\n"
                     % (c, length.get(c, 0),
                        "yes" if c in bin_c else "no",
                        "yes" if c in plus_c else "no",
                        "yes" if c in isl_c else "no"))
    print("  kleborate/bin_membership.tsv に %d 行（袋 %d / 足した後 %d / 島 %d）"
          % (len(bin_c | plus_c | isl_c), len(bin_c), len(plus_c), len(isl_c)))
    assert plus_c == bin_c | isl_c, "kp_plus が「袋＋島の contig」になっていない"
    print("  確認: kp_plus = kp_bin ∪ 島の contig（ちょうど一致）")


if __name__ == "__main__":
    main()

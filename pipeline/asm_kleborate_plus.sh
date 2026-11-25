#!/bin/bash
# asm_kleborate.sh が作った袋に、島を載せた contig を足して Kleborate をもう一度当てる。
#
# なぜこれを別の script にして配るか
#   RESULTS.md に書いた「ST23/K1/ICEKp10」は、**2 回の実行を混ぜたもの**だった
#   （CORRECTIONS.md の 16 番）。最初の実行で何が出て、何を足したら何が増えたのかが
#   分からないと、読み手は型別を評価できない。だから足す規則を script にして出す。
#   もとは手元で 1 回打った命令で、記録が残っていなかった。
#
# 足す規則（ちょうどこれだけ）
#   kp_plus = kp_bin（自動で K. pneumoniae に割り振られた contig）
#             ∪ 島を載せた contig（pks 島が当たった contig の全部）
#   data/kleborate/bin_membership.tsv がこの関係を満たすことを
#   code/verify_extra.py が assert で確かめている。
#
# なぜ足す必要があったか
#   contig の種への割り振りは 6 種の完全ゲノムに当てて勝ち負けで決めている。
#   ところが **K. pneumoniae HS11286 は clb 島を持っていない**
#   （島の 40-mer のうち共有は 0.11%、clb 遺伝子 0/18）。一方
#   **C. koseri BAA-895 は持っている**（94.88%、18/18）。
#   だから clb でできた contig は、この参照の組では K. pneumoniae に
#   絶対に割り振られない。足したのは参照の組のせいで、検体を見てからの選択ではない。
#   data/clb_in_six_genomes.tsv がこの数字で、code/make_clb_in_genomes.py が作る。
#
# それでも残る限界
#   手で組んだ袋から出した型別は、**その contig がこの Klebsiella のものだという
#   独立した証拠にはならない**。独立した証拠は asm_junction.sh のほう
#   （染色体と島がひと続きで、接合部の読み取りの深さが 0 にならない）であって、
#   そちらに使う contig は自動の袋に最初から入っている。
#
# 使い方: bash asm_kleborate_plus.sh [検体名]
set -u
S=${1:-cancer_ERR1018186}
B=/mnt/c/crc_strain_pilot
W=$B/asm/$S/kleborate
KLEB=/root/kleb_venv/bin/kleborate   # Kleborate 3.2.4
CTG=$B/asm/$S/megahit/final.contigs.fa

[ -s "$W/kp_contigs.txt" ]    || { echo "先に asm_kleborate.sh を走らせること"; exit 1; }
[ -s "$W/island_contigs.fa" ] || { echo "島の contig が無い"; exit 1; }

# 袋の名前 ∪ 島の contig の名前
grep '^>' "$W/island_contigs.fa" | sed 's/^>//; s/[[:space:]].*//' \
  | cat - "$W/kp_contigs.txt" | sort -u > "$W/kp_plus.txt"
echo "袋 $(wc -l < "$W/kp_contigs.txt") 本 -> 足した後 $(wc -l < "$W/kp_plus.txt") 本"

# その名前の contig だけを取り出す（megahit の出力は折り返さないので 1 本 2 行）
awk 'NR==FNR{k[$1]=1; next}
     /^>/{n=substr($1,2); f=(n in k)}
     f' "$W/kp_plus.txt" "$CTG" > "$W/kp_plus.fa"

"$KLEB" -a "$W/kp_plus.fa" -o "$W/out_kp_plus" -p kpsc

# 2 回の実行の違いだけを並べる
for T in kp_bin kp_plus; do
  F="$W/out_$T/klebsiella_pneumo_complex_output.txt"
  [ -s "$F" ] || continue
  echo "--- $T"
  paste <(head -1 "$F" | tr '\t' '\n') <(sed -n 2p "$F" | tr '\t' '\n') \
    | grep -iE 'species__species$|mlst__ST$|K_type|Yersiniabactin|CbST|Colibactin|total_size|contig_count' \
    | sed 's/^/    /'
done

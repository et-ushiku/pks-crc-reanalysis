#!/bin/bash
# 読み取りを**アセンブル全体に貼り直し**、contig ごとの読み取りの深さと、
# contig をまたぐ読み取り対（連結）を数える。
#
# なぜ要るか
#   CORRECTIONS.md 16 番は、手で足した 4 本の contig の出どころを
#   megahit の `multi` で論じていた。この組み立ては k=119 で読み取りは 100 塩基
#   なので、`multi` は「読み取り + 小さい k から持ち越した contig」で作った
#   de Bruijn グラフ上の重複度であり、**読み取りの深さではない**。
#   較正もしていなかった。だから読み取りの深さを直接測る。
#
#   同時に、対の片方が別の contig に落ちた組を数える。手で足した contig が
#   自動の袋の contig と**物理的につながっている**なら、そこに対が渡る。
#
# なぜ索引をアセンブル全体にするか（袋だけにしてはいけない）
#   袋だけに当てると、便の中の別の細菌の読み取りも行き場が無くて袋に寄る。
#   21,759 本すべてを 1 つの索引に入れれば競合になり、読み取りは
#   アセンブルの中で一番合う場所に落ちる。MAPQ 20 で絞れば寄せ集めも除ける。
#
# なぜ -f 2（proper pair）で絞らないか
#   contig をまたぐ対は bowtie2 の定義では concordant ではない。
#   -f 2 を付けると**数えたい信号そのものが消える**。
#   代わりに「相手が別の参照に落ちていて、相手も当たっている」で選ぶ。
#
# 出すもの（$B/meta/link186/ の下）
#   keep.bam            対象 contig に落ちた読み取りだけ（深さ用）
#   keep.depth.tsv.gz   その位置ごとの深さ（MAPQ 20）
#   links.tsv           contig の組ごとの対の数
#
# 使い方: bash asm_linkage.sh [検体名] [スレッド数]
#   索引作りに数分、貼り付けに 20 分ほど（26,427,396 対・16 スレッド）。
set -u
S=${1:-cancer_ERR1018186}
NP=${2:-14}
B=/mnt/c/crc_strain_pilot
A=$B/asm/$S
CTG=$A/megahit/final.contigs.fa
W=/tmp/link_$S
OUT=$B/meta/link186
MAPQ=20
mkdir -p "$W" "$OUT"
LOG=$OUT/run.log
echo "=== 開始 $S $(date)" > "$LOG"

[ -s "$CTG" ] || { echo "  contig が無い $CTG" >> "$LOG"; exit 1; }

if [ ! -s "$W/asm.1.bt2" ] && [ ! -s "$W/asm.1.bt2l" ]; then
  bowtie2-build --threads "$NP" "$CTG" "$W/asm" >> "$LOG" 2>&1 \
    || { echo "  索引の作成に失敗" >> "$LOG"; exit 1; }
fi

# 対象 contig: Kleborate に渡した袋 + 手で足した島の contig
sort -u "$A/kleborate/kp_plus.txt" > "$OUT/keep.txt"

bowtie2 -x "$W/asm" -1 "$A"/*.1.fq.gz -2 "$A"/*.2.fq.gz \
    -p "$NP" --no-unal -X 1000 2>> "$LOG" \
 | samtools view -h -q "$MAPQ" - \
 | awk -F'\t' -v keep="$OUT/keep.txt" -v links="$OUT/links.tsv" -v sam="$W/keep.sam" '
     BEGIN { while ((getline l < keep) > 0) K[l] = 1 }
     /^@/ { print > sam; next }
     {
       if ($3 in K) print > sam
       # 7 列目 = 相手の参照。"=" でなければ別の contig。
       # and($2,8)==0 は「相手も当たっている」。片側が対象 contig の組だけ数える。
       if ($7 != "=" && $7 != "*" && and($2, 8) == 0 && (($3 in K) || ($7 in K))) {
         a = $3; b = $7; if (a > b) { t = a; a = b; b = t }
         c[a "\t" b]++
       }
     }
     END { for (k in c) printf "%s\t%d\n", k, c[k] > links }'

samtools sort -@ 6 -o "$OUT/keep.bam" "$W/keep.sam" >> "$LOG" 2>&1
samtools index "$OUT/keep.bam" >> "$LOG" 2>&1
samtools depth -a -Q "$MAPQ" "$OUT/keep.bam" | gzip -c > "$OUT/keep.depth.tsv.gz"
echo "=== おわり $S $(date)" >> "$LOG"

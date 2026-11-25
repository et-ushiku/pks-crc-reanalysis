#!/bin/bash
# 接合部の検証 ― 島と K. pneumoniae 骨格が「本当に同じ分子」に載っているか。
#
# なぜ要るか
#   アセンブルで 386,079 bp の contig（k119_45171）が組み上がり、
#     contig       1 〜 354,278  … K. pneumoniae HS11286 に 99.0〜99.5% で連続一致
#     contig 354,225 〜 386,079  … pks 島（AM229678.1 に 99.98%）
#   と並んだ。しかしこれは**アセンブラが繋いだ結果**であって、
#   繋ぎ間違い（キメラ）の可能性が原理上まだ残っている。
#
# どう潰すか
#   元の読み取りを contig に貼り直し、**接合部をまたぐ証拠**を数える。
#     (a) またぐペア   片方が骨格側、もう片方が島側に、正しい向きと距離で載るペア
#     (b) またぐ単読   1本の読み取りが接合部の両側に 20bp 以上かかっている
#   キメラなら (a) も (b) も 0 になる。実在なら両方が周辺の深さなりに出る。
#   深さの落ち込みも見る（接合部だけ深さが落ちるならそこが継ぎ目）。
#
# 使い方: bash asm_junction.sh [スレッド数]
set -u
S=cancer_ERR1018186
C=k119_45171
J=354250          # 接合部（骨格の終わり 354,278 / 島の始まり 354,225 の間）
B=/mnt/c/crc_strain_pilot
A=$B/asm/$S
W=$A/junction
LOG=$B/meta/asm_junction.log
NP=${1:-4}
mkdir -p "$W"
echo "=== 接合部の検証 開始 $S $C 接合部=$J $(date)" >> "$LOG"

# 1) 島が載った contig だけ取り出す
if [ ! -s "$W/target.fa" ]; then
  awk -v c="$C" '$0 ~ "^>"c" " {f=1; print ">"c; next} /^>/{f=0} f' \
      "$A/megahit/final.contigs.fa" > "$W/target.fa"
fi
LEN=$(grep -v '^>' "$W/target.fa" | tr -d '\n' | wc -c)
echo "  取り出した contig の長さ $LEN bp" >> "$LOG"
if [ "$LEN" -lt 300000 ]; then echo "  contig が取れていない" >> "$LOG"; exit 1; fi

# 2) 索引
if [ ! -f "$W/t.1.bt2" ]; then
  bowtie2-build -q --threads "$NP" "$W/target.fa" "$W/t" >> "$LOG" 2>&1 \
    || { echo "  索引の作成が失敗" >> "$LOG"; exit 1; }
fi

R1=$(ls "$A"/*.1.fq.gz 2>/dev/null | head -1)
R2=$(ls "$A"/*.2.fq.gz 2>/dev/null | head -1)
if [ -z "$R1" ] || [ -z "$R2" ]; then echo "  読み取りが無い" >> "$LOG"; exit 1; fi
echo "  読み取り $(basename "$R1") / $(basename "$R2")" >> "$LOG"

# 3) ペアのまま貼る（向きと距離を使うので -1/-2。-U では判定できない）
if [ ! -s "$W/j.bam" ]; then
  bowtie2 -x "$W/t" -1 "$R1" -2 "$R2" --no-unal -p "$NP" 2>> "$LOG" \
    | samtools sort -@ 2 -m 1G -o "$W/j.bam" - >> "$LOG" 2>&1 \
    || { echo "  貼り付けが失敗" >> "$LOG"; exit 1; }
  samtools index "$W/j.bam"
fi

# 4) 接合部まわりの深さ（MAPQ>=20）
samtools depth -a -Q 20 -r "$C:352000-357000" "$W/j.bam" > "$W/depth_junction.tsv"
awk -v j=$J -F'\t' '
  { d[$2]=$3 }
  END{
    for (p in d) {
      if (p>=352000 && p<j-500)      { a+=d[p]; na++ }
      else if (p>=j-200 && p<=j+200) { b+=d[p]; nb++ }
      else if (p>j+500 && p<=357000) { c+=d[p]; nc++ }
    }
    printf "  平均深さ  骨格側(352000-%d) %.1f / 接合部±200 %.1f / 島側(%d-357000) %.1f\n",
           j-500, (na?a/na:0), (nb?b/nb:0), j+500, (nc?c/nc:0)
  }' "$W/depth_junction.tsv" >> "$LOG"
awk -v j=$J -F'\t' '$2>=j-200 && $2<=j+200 {if($3<mn||mn==""){mn=$3; mp=$2}} \
  END{printf "  接合部±200 の最小深さ %s（位置 %s）\n", mn, mp}' \
  "$W/depth_junction.tsv" >> "$LOG"

# 5) またぐペア（-f 2 = 正しく向き合って載ったペアのみ）
samtools view -q 20 -f 2 "$W/j.bam" "$C:350000-359000" 2>/dev/null \
  | awk -v j=$J '{ p=$4; m=$8; l=length($10);
      if ((p+l-1 < j-25 && m > j+25) || (m+l-1 < j-25 && p > j+25)) n++ }
      END{ printf "  (a) またぐペア（正しく向き合ったもの） %d 本\n", n+0 }' >> "$LOG"

# 6) またぐ単読（1本が接合部の両側に20bp以上かかる）
samtools view -q 20 "$W/j.bam" "$C:$((J-150))-$((J+150))" 2>/dev/null \
  | awk -v j=$J '{ p=$4; l=length($10);
      if (p <= j-20 && p+l-1 >= j+20) n++ }
      END{ printf "  (b) またぐ単読（両側20bp以上） %d 本\n", n+0 }' >> "$LOG"

# 7) 対照 ― 骨格の内側（島から遠い場所）で同じ数え方をする。
#    ここは繋ぎ目ではないので、(a)(b) の「当たり前の値」が分かる。
for REF in 200000 300000; do
  samtools view -q 20 -f 2 "$W/j.bam" "$C:$((REF-5000))-$((REF+5000))" 2>/dev/null \
    | awk -v j=$REF '{ p=$4; m=$8; l=length($10);
        if ((p+l-1 < j-25 && m > j+25) || (m+l-1 < j-25 && p > j+25)) n++ }
        END{ printf "  対照 位置%d: またぐペア %d 本\n", j, n+0 }' >> "$LOG"
done

echo "=== 接合部の検証 おわり $(date)" >> "$LOG"

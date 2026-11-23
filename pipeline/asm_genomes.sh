#!/bin/bash
# 腸内細菌科6種の**完全ゲノム**に当てて、種ごとのゲノムの深さを測る。
#
# なぜ目印配列ではなくゲノムか
#   MetaPhlAn の目印は K. pneumoniae で 6,450 塩基しかなく、深さの推定が粗い。
#   また目印は5系統しか用意していないので、それ以外の担い手を見落とす。
#   完全ゲノム（1種あたり約500万塩基）なら深さを精密に測れるうえ、
#   覆われ方（breadth）も見えるので「本当にいるのか」を判定できる。
#
#   6種を1本の索引にまとめるのが肝心。近縁種にまたがる読み取りは
#   最も似た1か所に落ち、MAPQ で絞れば種をまたぐ染み出しを除ける。
#
# 使い方: bash asm_genomes.sh 検体名 [スレッド数]
set -u
S=$1 ; NPROC=${2:-10}
B=/mnt/c/crc_strain_pilot
D=$B/asm/$S
IDX=$B/asm/ref/combined/entero6_idx
LEN=$B/asm/ref/entero6_len.tsv
OUT=$B/pks/asm_out
LOG=$B/meta/asm_measure.log
MAPQ=20
mkdir -p "$OUT"

FQ=$(ls "$D"/*.fq.gz 2>/dev/null)
[ -z "$FQ" ] && FQ=$(ls "$D"/*.fastq.gz 2>/dev/null)
echo "=== 6種ゲノムに当てる $S $(date)" >> "$LOG"
echo "  使うファイル: $(echo $FQ | tr ' ' '\n' | xargs -n1 basename | tr '\n' ' ')" >> "$LOG"

bowtie2 -x "$IDX" -U $(echo $FQ | tr ' ' ',') --very-sensitive --no-unal \
    -p "$NPROC" -S "$OUT/$S.g6.sam" 2> "$OUT/$S.g6.bowtie2.log"
samtools sort -@ 4 -o "$OUT/$S.g6.bam" "$OUT/$S.g6.sam" 2>> "$LOG"
samtools index "$OUT/$S.g6.bam"
rm -f "$OUT/$S.g6.sam"

# 種ごとの 当たり本数（MAPQ20）と 覆われた塩基数
samtools idxstats "$OUT/$S.g6.bam" | awk '$1!="*"{print $1"\t"$2"\t"$3}' > "$OUT/$S.g6.idxstats.tsv"
samtools depth -a -Q "$MAPQ" "$OUT/$S.g6.bam" \
  | awk -F'\t' '{split($1,a,"|"); k=a[1]; n[k]++; if($3>0){c[k]++; s[k]+=$3}}
      END{for(k in n) printf "%s\t%d\t%d\t%d\n", k, n[k], c[k]+0, s[k]+0}' \
  | sort > "$OUT/$S.g6.depth.tsv"
echo "--- 終了 $S $(date)" >> "$LOG"
cat "$OUT/$S.g6.depth.tsv" >> "$LOG"

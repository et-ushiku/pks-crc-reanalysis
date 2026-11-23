#!/bin/bash
# 落とし直した検体を測り直す（pks 島 ＋ 腸内細菌科5系統）。
#
#   cn_map.sh と同じ条件で当てる。違うのは
#     - fastq を消さない（このあとアセンブルに使うため）
#     - 出力先が pks/asm_out
#
# 使い方: bash asm_measure.sh 検体名 [スレッド数]
set -u
S=$1 ; NPROC=${2:-10}
B=/mnt/c/crc_strain_pilot
D=$B/asm/$S
IDX=$B/pks/ref/pks_idx
NIDX=$B/pks/ref/entero/entero_idx
BED=$B/pks/ref/clb_genes.bed
NLEN=$B/pks/ref/entero/entero_len.tsv
OUT=$B/pks/asm_out
LOG=$B/meta/asm_measure.log
MAPQ=20
mkdir -p "$OUT"
echo "=== 開始 $S $(date)" >> "$LOG"

FQ=$(ls "$D"/*.fq.gz 2>/dev/null)
[ -z "$FQ" ] && FQ=$(ls "$D"/*.fastq.gz 2>/dev/null)
echo "  使うファイル: $(echo $FQ | tr ' ' '\n' | xargs -n1 basename | tr '\n' ' ')" >> "$LOG"

# 1パス目: pks 島
bowtie2 -x "$IDX" -U $(echo $FQ | tr ' ' ',') --very-sensitive --no-unal \
    -p "$NPROC" -S "$OUT/$S.sam" 2> "$OUT/$S.bowtie2.log"
samtools sort -@ 4 -o "$OUT/$S.bam" "$OUT/$S.sam" 2>> "$LOG"
samtools index "$OUT/$S.bam"
rm -f "$OUT/$S.sam"
samtools depth -a -Q "$MAPQ" -b "$BED" "$OUT/$S.bam" > "$OUT/$S.depth.tsv"
npks=$(samtools view -c -q "$MAPQ" "$OUT/$S.bam")

# 2パス目: 腸内細菌科5系統（種別に数える）
bowtie2 -x "$NIDX" -U $(echo $FQ | tr ' ' ',') --very-sensitive --no-unal \
    -p "$NPROC" -S "$OUT/$S.entero.sam" 2> "$OUT/$S.entero.bowtie2.log"
samtools view "$OUT/$S.entero.sam" \
  | awk -F'\t' -v LEN="$NLEN" -v Q="$MAPQ" '
      BEGIN{while ((getline l < LEN) > 0) {split(l,b,"\t"); L[b[1]]=b[2]}}
      {split($3,a,"__"); k=a[1]"__"a[2]; n[k]++; if ($5+0 >= Q) m[k]++}
      END{for(k in L) printf "%s\t%d\t%d\t%d\n", k, n[k]+0, m[k]+0, L[k]}' \
  | sort > "$OUT/$S.counts.tsv"
rm -f "$OUT/$S.entero.sam"

echo "--- 終了 $S $(date)  pks(MAPQ20)=$npks" >> "$LOG"
cat "$OUT/$S.counts.tsv" >> "$LOG"

#!/bin/bash
# pks 島（コリバクチンを作る遺伝子群）が各検体にあるかを調べる。
#   参照: AM229678（E. coli IHE3034 の pks 島 55,140 bp。Nougayrede 2006）
#   手順: 生の読み取りを pks 島に当て、clb 遺伝子ごとに「どれだけ覆われたか」を数える
# 偽陽性対策: 当たりの質（MAPQ）で絞り、18個の clb 遺伝子に満遍なく当たるかを見る。
#            よそから来た読み取りは、保存された一部分にだけ集まる。
B=/mnt/c/crc_strain_pilot
DATA=$B/data_full
REF=$B/pks/ref/pks_island.fasta
BED=$B/pks/ref/clb_genes.bed
IDX=$B/pks/ref/pks_idx
OUT=$B/pks/out
LOG=$B/meta/pks_map.log
NPROC=${1:-12}
MAPQ=20

mkdir -p "$OUT"
if [ ! -s "${IDX}.1.bt2" ]; then
  echo "=== 索引を作る $(date)" >> "$LOG"
  bowtie2-build --threads "$NPROC" "$REF" "$IDX" >> "$LOG" 2>&1
fi

echo "=== 開始 $(date) nproc=$NPROC MAPQ>=$MAPQ" >> "$LOG"
for d in "$DATA"/*/; do
  S=$(basename "$d")
  [ -f "$d/.complete" ] || { echo "飛ばす（未完成） $S" >> "$LOG"; continue; }
  if [ -s "$OUT/$S.depth.tsv" ]; then echo "済み $S" >> "$LOG"; continue; fi
  echo "--- 開始 $S $(date)" >> "$LOG"

  # 総読み取り数も数えておく（正規化に使う）
  bowtie2 -x "$IDX" -U "$d"/*.fastq.gz --very-sensitive --no-unal \
      -p "$NPROC" -S "$OUT/$S.sam" 2> "$OUT/$S.bowtie2.log"
  ex=$?
  if [ $ex -ne 0 ]; then echo "--- 失敗 $S EXIT=$ex" >> "$LOG"; rm -f "$OUT/$S.sam"; continue; fi

  samtools sort -@ 4 -o "$OUT/$S.bam" "$OUT/$S.sam" 2>> "$LOG"
  samtools index "$OUT/$S.bam"
  rm -f "$OUT/$S.sam"
  # 質の低い当たりを捨てて、位置ごとの深さを出す
  samtools depth -a -Q "$MAPQ" -b "$BED" "$OUT/$S.bam" > "$OUT/$S.depth.tsv"
  n=$(samtools view -c -q "$MAPQ" "$OUT/$S.bam")
  echo "--- 終了 $S $(date) 当たった読み取り(MAPQ>=$MAPQ)=$n" >> "$LOG"
  rm -f "$OUT/$S.bam" "$OUT/$S.bam.bai"
done
echo "=== 全部終了 $(date)" >> "$LOG"

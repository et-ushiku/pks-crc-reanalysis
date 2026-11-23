#!/bin/bash
# AT-CRC 109検体を pks 島に当てる。at_download.py と並行して動かす。
#   .complete が置かれた検体を順に処理し、終わったら fastq を消して .mapped を置く
#   （落とす→当てる→消す の流し込み。ディスクを使いすぎない）
#
# 設定は既存40検体（pks_map.sh）と完全に同一にしている。変えると比較できなくなる。
B=/mnt/c/crc_strain_pilot
DATA=$B/data_at
REF=$B/pks/ref/pks_island.fasta
BED=$B/pks/ref/clb_genes.bed
IDX=$B/pks/ref/pks_idx
OUT=$B/pks/out_at
LOG=$B/meta/at_map.log
NPROC=${1:-12}
MAPQ=20

mkdir -p "$OUT"
echo "=== 開始 $(date) nproc=$NPROC MAPQ>=$MAPQ" >> "$LOG"

while true; do
  did=0
  for d in "$DATA"/*/; do
    S=$(basename "$d")
    [ -f "$d/.complete" ] || continue
    [ -f "$d/.mapped" ] && continue
    echo "--- 開始 $S $(date)" >> "$LOG"

    bowtie2 -x "$IDX" -U "$d"/*.fastq.gz --very-sensitive --no-unal \
        -p "$NPROC" -S "$OUT/$S.sam" 2> "$OUT/$S.bowtie2.log"
    ex=$?
    if [ $ex -ne 0 ]; then
      echo "--- 失敗 $S EXIT=$ex" >> "$LOG"; rm -f "$OUT/$S.sam"; continue
    fi

    samtools sort -@ 4 -o "$OUT/$S.bam" "$OUT/$S.sam" 2>> "$LOG"
    samtools index "$OUT/$S.bam"
    rm -f "$OUT/$S.sam"
    samtools depth -a -Q "$MAPQ" -b "$BED" "$OUT/$S.bam" > "$OUT/$S.depth.tsv"
    n=$(samtools view -c -q "$MAPQ" "$OUT/$S.bam")
    rm -f "$OUT/$S.bam" "$OUT/$S.bam.bai"

    # 当て終わったので生の読み取りを消す（ディスクを空ける）。印は残す。
    sz=$(du -sm "$d" 2>/dev/null | cut -f1)
    rm -f "$d"/*.fastq.gz
    touch "$d/.mapped"
    echo "--- 終了 $S $(date) 当たった読み取り(MAPQ>=$MAPQ)=$n fastq ${sz}MB を削除" >> "$LOG"
    did=1
  done

  if [ $did -eq 0 ]; then
    if [ -f "$DATA/.download_finished" ]; then
      echo "=== 全部終了 $(date)" >> "$LOG"; exit 0
    fi
    sleep 120
  fi
done

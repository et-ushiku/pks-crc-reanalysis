#!/bin/bash
# 検体をアセンブルして、pks 島がどの種のゲノムに繋がっているかを見る。
#
# 考え方
#   深さの比（cn_carrier.py）は「大腸菌では説明がつかない」までは言えるが、
#   「では誰か」には答えられない。島とその前後が**同じ連続配列（contig）に載る**ことを
#   示して初めて担い手が決まる。そのためにアセンブルする。
#
# 手順
#   1. megahit で組み立てる
#   2. 組み上がった contig に pks 島（AM229678.1、56kb）を当て、島を載せた contig を拾う
#   3. その contig を腸内細菌科6種の完全ゲノムに当て、**島以外の部分**がどの種に一致するかを見る
#      島そのものは6種のどれにも載りうる（可動性要素なので）。効くのは前後の染色体部分。
#
# 使い方: bash asm_assemble.sh 検体名 [スレッド数]
set -u
S=$1 ; NPROC=${2:-10}
B=/mnt/c/crc_strain_pilot
D=$B/asm/$S
OUT=$D/megahit
REF=$B/asm/ref/combined/entero6.fasta
PKS=$B/pks/ref/pks_island.fasta
LOG=$B/meta/asm_assemble.log

R1=$(ls "$D"/*.1.fq.gz 2>/dev/null | head -1)
R2=$(ls "$D"/*.2.fq.gz 2>/dev/null | head -1)
if [ -z "$R1" ] || [ -z "$R2" ]; then
  echo "対になる fastq が見つからない（$D）" >> "$LOG"; exit 1
fi
echo "=== アセンブル開始 $S $(date)" >> "$LOG"
echo "  $R1" >> "$LOG"; echo "  $R2" >> "$LOG"

# megahit は作業用に名前付きパイプ（FIFO）を作るが、Windows 側の /mnt/c では
#   OSError: [Errno 95] Operation not supported
# になる。DrvFs が FIFO に対応していないため。
# そこで**組み立ては WSL 側の領域（ext4）で行い**、結果だけ /mnt/c に戻す。
WORK=/root/megahit_work/$S
rm -rf "$OUT" "$WORK"
mkdir -p /root/megahit_work
megahit -1 "$R1" -2 "$R2" -o "$WORK" -t "$NPROC" -m 0.7 \
        --min-contig-len 1000 >> "$LOG" 2>&1 || { echo "megahit が失敗" >> "$LOG"; exit 1; }
mkdir -p "$OUT"
cp "$WORK/final.contigs.fa" "$OUT/" || { echo "contig が無い" >> "$LOG"; exit 1; }
cp "$WORK/log" "$OUT/megahit.log" 2>/dev/null
rm -rf "$WORK"      # 中間ファイルは大きいので消す
C=$OUT/final.contigs.fa
echo "--- 組み上がり $(grep -c '^>' "$C") 本  $(date)" >> "$LOG"

# 島を載せた contig を拾う
makeblastdb -in "$C" -dbtype nucl -out "$OUT/contigs" >/dev/null
blastn -query "$PKS" -db "$OUT/contigs" -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore' \
       -evalue 1e-20 -num_threads "$NPROC" > "$OUT/pks_hits.tsv"
cut -f2 "$OUT/pks_hits.tsv" | sort -u > "$OUT/pks_contigs.txt"
echo "--- 島が載った contig $(wc -l < "$OUT/pks_contigs.txt") 本" >> "$LOG"

# その contig を取り出す
awk 'NR==FNR{k[$1]=1; next} /^>/{nm=substr($1,2); p=(nm in k)} p' \
    "$OUT/pks_contigs.txt" "$C" > "$OUT/pks_contigs.fa"

# 6種の完全ゲノムに当てる
makeblastdb -in "$REF" -dbtype nucl -out "$OUT/entero6" >/dev/null
blastn -query "$OUT/pks_contigs.fa" -db "$OUT/entero6" \
       -outfmt '6 qseqid qlen sseqid pident length qstart qend sstart send evalue bitscore' \
       -evalue 1e-20 -num_threads "$NPROC" > "$OUT/contig_vs_genomes.tsv"
echo "=== 終了 $S $(date)" >> "$LOG"

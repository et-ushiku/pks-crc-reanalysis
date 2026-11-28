#!/usr/bin/env bash
# CORRECTIONS.md 10 番を、理屈ではなく実測で確かめる。
#
# 10 番はこう直した: 「MAPQ で絞らないと clbJ と clbK の RPKM が水増しされる」は
# 間違いで、整列器は 1 本の読み取りに primary を 1 つしか出すから、重複領域の
# 読み取りは **1 回しか数えられない**。だから 2 遺伝子の合計は保たれ、
# 壊れるのは 2 遺伝子の間の **振り分け** だけである、と。
#
# ところがこの直し方自体が「整列器はこう振る舞うはずだ」という理屈で、
# 15・16 番と同じ型（測っていないことを根拠として書く）だった。ここで測る。
#
# 測るのは 3 つ。
#   (1) primary の本数 == 読み取り名の種類数 か。
#       等しければ、1 本の読み取りが 2 回数えられることは無い。これが
#       「水増しされない」の中身そのもの。secondary/supplementary も数える。
#   (2) clbJ+clbK の当たり本数を、絞らない場合と MAPQ 20 で絞る場合で比べる。
#       絞らない側が真の本数で、絞る側はそこから重複ぶんだけ **減る**。
#       「水増し」ではなく「欠損」であることを数字で見せる。
#   (3) 振り分けが恣意的か。重複の 2 つの複製それぞれに、MAPQ 20 未満の
#       読み取りがどう積まれるかを見る。
#
# 使っている整列器は bowtie2 である（この研究の全工程がそう）。既定では
# 1 本につき最良の 1 か所しか報告しないので、10 番の論の土台はここで測れる。
# **BWA-MEM そのものは走らせていない**（この環境に入っていない）。10 番の
# 文は BWA-MEM の名前で書いてあるので、そこは直す。
#
# 検体は cancer_ERR1018186。読み取りが手元に残っている唯一の担体である
# （島の読み取り 16,513 本、覆い 0.9578）。ENA から取り直さない。
set -u

B=/mnt/c/crc_strain_pilot
REF=$B/pks/ref/pks_island.fasta
A=$B/asm/cancer_ERR1018186
BED=$B/release/pks-crc-reanalysis/data/clb_genes.bed
OUT=$B/meta/mapq186
W=/tmp/mapq186
NP=${NPROC:-14}
MAPQ=20
# 重複の 2 つの複製。PITFALLS.md 1 番 / CORRECTIONS.md 9 番の実測値
REP1_S=16066; REP1_E=17352
REP2_S=21719; REP2_E=23005

mkdir -p "$OUT" "$W"
LOG=$OUT/run.log
echo "=== 開始 $(date) ===" > "$LOG"

if [ ! -f "$W/ref.1.bt2" ]; then
  bowtie2-build -q --threads "$NP" "$REF" "$W/ref" >> "$LOG" 2>&1 \
    || { echo "bowtie2-build が失敗" >> "$LOG"; exit 1; }
fi

# 本番の工程（asm_measure.sh）と同じ呼び方。両側を別々の単独読み取りとして渡す
FQ=$(ls "$A"/*.rmhost.1.fq.gz "$A"/*.rmhost.2.fq.gz | tr '\n' ',' | sed 's/,$//')
echo "reads: $FQ" >> "$LOG"

if [ ! -f "$OUT/all.bam" ]; then
  bowtie2 -x "$W/ref" -U "$FQ" --very-sensitive --no-unal -p "$NP" \
      2>> "$LOG" \
    | samtools sort -@ 4 -o "$OUT/all.bam" - 2>> "$LOG" \
    || { echo "bowtie2 が失敗" >> "$LOG"; exit 1; }
  samtools index "$OUT/all.bam"
fi

# ------------------------------------------------- (1) 二重計上は起きるか
{
echo "## (1) 1 本の読み取りが 2 回数えられるか"
prim=$(samtools view -c -F 0x900 "$OUT/all.bam")
sec=$(samtools view -c -f 0x100 "$OUT/all.bam")
sup=$(samtools view -c -f 0x800 "$OUT/all.bam")
names=$(samtools view -F 0x900 "$OUT/all.bam" | cut -f1 | sort -u | wc -l)
echo "primary            $prim"
echo "読み取り名の種類数 $names"
echo "secondary          $sec"
echo "supplementary      $sup"
} | tee "$OUT/counts.txt"

# ------------------------------------------- (2) 遺伝子ごと、絞る前と後
echo "## (2) 遺伝子ごとの当たり本数" >> "$OUT/counts.txt"
printf "gene\tgene_bp\treads_unfiltered\treads_mapq20\n" > "$OUT/per_gene.tsv"
while IFS=$'\t' read -r c s e g; do
  [ -z "${g:-}" ] && continue
  n0=$(samtools view -c "$OUT/all.bam" "$c:$((s+1))-$e")
  n1=$(samtools view -c -q "$MAPQ" "$OUT/all.bam" "$c:$((s+1))-$e")
  printf "%s\t%d\t%s\t%s\n" "$g" "$((e-s))" "$n0" "$n1" >> "$OUT/per_gene.tsv"
done < "$BED"

# ------------------------------------------- (3) 振り分けは恣意的か
echo "## (3) 重複の 2 つの複製への積まれ方" >> "$OUT/counts.txt"
C=$(cut -f1 "$BED" | head -1)
printf "copy\tstart\tend\treads_unfiltered\treads_mapq20\n" > "$OUT/repeat_split.tsv"
for p in "copy1 $REP1_S $REP1_E" "copy2 $REP2_S $REP2_E"; do
  set -- $p
  n0=$(samtools view -c "$OUT/all.bam" "$C:$2-$3")
  n1=$(samtools view -c -q "$MAPQ" "$OUT/all.bam" "$C:$2-$3")
  printf "%s\t%s\t%s\t%s\t%s\n" "$1" "$2" "$3" "$n0" "$n1" >> "$OUT/repeat_split.tsv"
done

# MAPQ の分布（重複の中と外）
samtools view "$OUT/all.bam" "$C:$REP1_S-$REP1_E" | awk '{print $5}' | sort -n \
  | uniq -c > "$OUT/mapq_hist_repeat.txt"
samtools view "$OUT/all.bam" "$C:25000-35000" | awk '{print $5}' | sort -n \
  | uniq -c > "$OUT/mapq_hist_unique.txt"

echo "=== 終了 $(date) ===" >> "$LOG"
echo "出力: $OUT"

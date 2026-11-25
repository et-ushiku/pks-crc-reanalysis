#!/bin/bash
# アセンブルした contig を Kleborate に当て、clb の型（CbST）と ICEKp の構造を判定する。
#
# なぜ要るか
#   「*K. pneumoniae* が clb 島を持つ」こと自体は 2009 年から既知（Putze 2009）で、
#   Lam 2018 は 2,498 ゲノムの 14% が clb を持つと報告し、**65 種類の CbST** に
#   型別するスキームまで作っている。だから「持っていた」だけでは何も新しくない。
#
#   意味があるのは**どの型か**である。Wami 2021 によれば clb 系統は概ね種特異的で、
#   *Klebsiella* 型（clb1 / clb2A / clb2B）と *E. coli* 型が区別できる。
#     - Klebsiella 型なら → 既知の ICEKp10 保有株。ありふれた所見
#     - E. coli 型なら    → Wami 2021 が報告した**種間伝達の実例**に該当しうる
#   また ICEKp10 が完全か（KpZM モジュール・T4SS・int/xis）も型別で分かる。
#
# なぜ contig を選り分けるか
#   Kleborate は 1 株のゲノムを前提にしている。こちらは便のメタゲノムなので、
#   21,759 本の contig をそのまま渡すと種同定も MLST も混ざって無意味になる。
#   そこで **K. pneumoniae に当たる contig だけを集めた袋**を作って渡す。
#   minimap2 で当てて、一致率と当たった長さの割合で選ぶ。
#
# 使い方: bash asm_kleborate.sh [検体名] [スレッド数]
set -u
S=${1:-cancer_ERR1018186}
NP=${2:-4}
B=/mnt/c/crc_strain_pilot
A=$B/asm/$S
CTG=$A/megahit/final.contigs.fa
KPREF=$B/asm/ref/Kpneumoniae_HS11286.fa
W=$A/kleborate
LOG=$B/meta/asm_kleborate.log
KLEB=/root/kleb_venv/bin/kleborate
mkdir -p "$W"
echo "=== Kleborate 開始 $S $(date)" >> "$LOG"

[ -s "$CTG" ]   || { echo "  contig が無い $CTG"   >> "$LOG"; exit 1; }
[ -s "$KPREF" ] || { echo "  Kp 参照が無い $KPREF" >> "$LOG"; exit 1; }

# 1) contig を **6 種のゲノムまとめて**に当てる
#
#    なぜ Kp 単独ではだめか（最初に間違えたところ）
#      Kp だけに当てると、腸内細菌科の近縁種の contig も「当たった」ことになる。
#      さらに PAF の一致率（列10/列11）は**ギャップを含む**ので、
#      実際に 99% 一致している contig でも 0.81 などと出る。これに 0.95 の門を
#      かけたので、ほぼ全部落ちた。
#    どう直すか
#      6 種に同時に当て、contig ごとに**一致塩基がいちばん多い種**へ割り振る。
#      絶対的な一致率ではなく**相対的な勝ち負け**で決めるので、
#      ギャップの扱いに左右されない。
G6=$B/asm/ref/combined/entero6.fasta
[ -s "$G6" ] || { echo "  6種ゲノムが無い $G6" >> "$LOG"; exit 1; }
if [ ! -s "$W/vs_g6.paf" ]; then
  minimap2 -x asm10 -t "$NP" "$G6" "$CTG" > "$W/vs_g6.paf" 2>> "$LOG" \
    || { echo "  minimap2 が失敗" >> "$LOG"; exit 1; }
fi
echo "  当たった行 $(wc -l < "$W/vs_g6.paf") / 異なる contig $(cut -f1 "$W/vs_g6.paf" | sort -u | wc -l)" >> "$LOG"

# 2) Kp らしい contig を選ぶ
#    PAF: 1=名前 2=長さ 3=開始 4=終わり 10=一致塩基 11=当たった長さ
#
#    注意（最初に間違えたところ）
#      1 本の contig は**複数の断片に分かれて**当たる。だから 1 行ごとに
#      (終わり-開始)/長さ を見ると、どの行も 0.5 を超えず全部落ちる。
#      contig ごとに**合計してから**割らなければならない。
#    割り振り: contig ごと・種ごとに一致塩基を合計し、最多の種を勝ちとする。
#    さらに「勝った種に対する被覆が 0.5 以上」を要求して、
#    ごく一部だけ当たった contig を除く。
awk -F'\t' '{
    split($6, a, "|"); sp = a[1]
    len[$1] = $2
    k = $1 SUBSEP sp
    mat[k] += $10
    aln[k] += $4 - $3
    seen[$1] = 1
  }
  END{
    for (c in seen) {
      bs = ""; bm = -1; ba = 0
      for (k in mat) {
        split(k, p, SUBSEP)
        if (p[1] != c) continue
        if (mat[k] > bm) { bm = mat[k]; bs = p[2]; ba = aln[k] }
      }
      cov = (len[c] > 0) ? ba / len[c] : 0
      printf "%s\t%s\t%d\t%.3f\n", c, bs, len[c], cov
    }
  }' "$W/vs_g6.paf" | sort > "$W/assign.tsv"

echo "  種ごとの contig 本数と塩基数（被覆 0.5 以上のもの）" >> "$LOG"
awk -F'\t' '$4>=0.5 {n[$2]++; b[$2]+=$3} END{for(s in n) printf "    %-24s %4d 本  %9d bp\n", s, n[s], b[s]}' \
    "$W/assign.tsv" | sort -k3,3nr >> "$LOG"

awk -F'\t' '$2=="Kpneumoniae_HS11286" && $4>=0.5 {print $1}' "$W/assign.tsv" \
    | sort -u > "$W/kp_contigs.txt"
NK=$(wc -l < "$W/kp_contigs.txt")
echo "  Kp らしい contig $NK 本" >> "$LOG"
if [ "$NK" -lt 10 ]; then echo "  少なすぎる。条件を見直すこと" >> "$LOG"; fi

# 3) その contig だけを取り出して袋にする
awk 'NR==FNR{k[$1]=1; next}
     /^>/{n=substr($1,2); f=(n in k)}
     f' "$W/kp_contigs.txt" "$CTG" > "$W/kp_bin.fa"
BP=$(grep -v '^>' "$W/kp_bin.fa" | tr -d '\n' | wc -c)
echo "  袋の大きさ $BP bp（K. pneumoniae のゲノムは約 5,300,000 bp）" >> "$LOG"

# 4) 島が載った contig も単独で用意する（袋に入らなかった場合の保険）
cp "$A/megahit/pks_contigs.fa" "$W/island_contigs.fa" 2>/dev/null

# 5) Kleborate を当てる
#    袋のほう: ST・ybST・CbST・ICEKp をまとめて見る
for TGT in kp_bin island_contigs; do
  [ -s "$W/$TGT.fa" ] || continue
  echo "  --- Kleborate: $TGT" >> "$LOG"
  "$KLEB" -a "$W/$TGT.fa" -o "$W/out_$TGT" -p kpsc >> "$LOG" 2>&1 \
    || echo "    失敗（$TGT）" >> "$LOG"
done

# 6) 要点を抜き出す
for TGT in kp_bin island_contigs; do
  R="$W/out_$TGT/klebsiella_pneumo_complex_output.txt"
  [ -s "$R" ] || R=$(ls "$W/out_$TGT"/*output.txt 2>/dev/null | head -1)
  [ -s "$R" ] || continue
  echo "  ### $TGT の結果" >> "$LOG"
  head -1 "$R" | tr '\t' '\n' | nl > "$W/$TGT.hdr"
  awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) h[i]=$i; next}
    {for(i=1;i<=NF;i++) if($i!="" && $i!="-" && $i!="0")
       printf "    %s = %s\n", h[i], $i}' "$R" >> "$LOG"
done

echo "=== Kleborate おわり $S $(date)" >> "$LOG"

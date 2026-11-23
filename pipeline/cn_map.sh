#!/bin/bash
# CN-CRC（Yu 2017、中国、128検体）を落ちてきた順に処理する。
#
# at_map2.sh からの変更点
#   3パス目を追加した。pks 島を持ちうる腸内細菌科4種（大腸菌・Klebsiella pneumoniae・
#   Citrobacter koseri・Klebsiella aerogenes）の目印配列に当て、種ごとの量を出す。
#   Putze 2009 のとおり pks 島は大腸菌専用ではないので、大腸菌だけを見ていると
#   「大腸菌は少ないが Klebsiella で pks を持つ」検体を取りこぼす。
#
# 大腸菌専用のパス（2パス目）は**残す**。
#   換算係数 313 は大腸菌だけの索引で求めた値で、既存149検体もそれで測っている。
#   索引を混ぜると MAPQ が変わりうるので、比較可能性のために別に測る。
#   3パス目の大腸菌の値と突き合わせれば、混ぜても大丈夫かが分かる。
#
# 1検体あたり: pks 約2分 + 大腸菌 約4分 + 腸内細菌科 約5分 ≒ 11〜14分
#   ダウンロードは1検体約17分なので、これで追いつく。
#
# 使い方: bash cn_map.sh 12     （引数は bowtie2 のスレッド数）
set -u
B=/mnt/c/crc_strain_pilot
DATA=$B/data_cn
IDX=$B/pks/ref/pks_idx                  # pks 島（AM229678.1、IHE3034）
EIDX=$B/pks/ref/ecoli/ecoli_idx         # 大腸菌の目印（t__SGB10068）
NIDX=$B/pks/ref/entero/entero_idx       # 腸内細菌科4種の目印
BED=$B/pks/ref/clb_genes.bed            # clb 遺伝子18個の位置
NLEN=$B/pks/ref/entero/entero_len.tsv   # SGB ごとの目印の総塩基数（固定表）
OUT=$B/pks/out_cn                       # pks の深さ
EOUT=$B/pks/ecoli_cn                    # 大腸菌の深さ
NOUT=$B/pks/entero_cn                   # 腸内細菌科の種別の数
LOG=$B/meta/cn_map.log
NPROC=${1:-12}
MAPQ=20

mkdir -p "$OUT" "$EOUT" "$NOUT"
echo "=== 開始 $(date) nproc=$NPROC" >> "$LOG"

# 索引が揃うまで待つ（entero_ref.py が作り終わるのを待つ場合がある）
while [ ! -s "$NIDX.1.bt2" ] && [ ! -s "$NIDX.1.bt2l" ]; do
  echo "$(date '+%H:%M:%S') 腸内細菌科の索引を待っている" >> "$LOG"
  sleep 120
done

idle=0
while true; do
  did=0
  for d in "$DATA"/*/; do
    S=$(basename "$d")
    [ -f "$d/.complete" ] || continue
    [ -f "$d/.mapped" ] && continue
    ls "$d"*.fastq.gz >/dev/null 2>&1 || continue
    echo "--- 開始 $S $(date)" >> "$LOG"
    t0=$(date +%s)

    # 1パス目: pks 島
    bowtie2 -x "$IDX" -U "$d"/*.fastq.gz --very-sensitive --no-unal \
        -p "$NPROC" -S "$OUT/$S.sam" 2> "$OUT/$S.bowtie2.log"
    if [ $? -ne 0 ]; then
      echo "--- 失敗(pks) $S" >> "$LOG"; rm -f "$OUT/$S.sam"; continue
    fi
    samtools sort -@ 4 -o "$OUT/$S.bam" "$OUT/$S.sam" 2>> "$LOG"
    samtools index "$OUT/$S.bam"
    rm -f "$OUT/$S.sam"
    samtools depth -a -Q "$MAPQ" -b "$BED" "$OUT/$S.bam" > "$OUT/$S.depth.tsv"
    npks=$(samtools view -c -q "$MAPQ" "$OUT/$S.bam")
    rm -f "$OUT/$S.bam" "$OUT/$S.bam.bai"

    # 2パス目（大腸菌専用）は 2026-11-18 に外した。
    #   混合索引（3パス目）の中の大腸菌が、専用索引の当たり数を**完全に再現する**ことを
    #   25検体で確認した（比 1.000、最大のずれが 740,329本中3本。`cn_entero.py`）。
    #   MetaPhlAn の目印は系統固有なので、近縁種にまたがって落ちる読み取りがほぼない。
    #   → 大腸菌の量は3パス目の MAPQ 無しの列から取る。換算係数 313 はそのまま使える。
    #   1検体あたり約4分の節約。ダウンロード（14分）とマッピング（18分→14分）が釣り合う。
    #   ※ ecoli_cn/ に残っている先行52検体ぶんの結果はそのまま使える（同じ値なので）。

    # 3パス目: 腸内細菌科4種。配列名の先頭が SGB 番号なので種別に数えられる
    bowtie2 -x "$NIDX" -U "$d"/*.fastq.gz --very-sensitive --no-unal \
        -p "$NPROC" -S "$NOUT/$S.sam" 2> "$NOUT/$S.bowtie2.log"
    if [ $? -ne 0 ]; then
      echo "--- 失敗(腸内細菌科) $S" >> "$LOG"; rm -f "$NOUT/$S.sam"; continue
    fi
    # 配列名は t__SGB10115__UniRef90_... の形。"__" で切った**最初の2片**が SGB になる。
    #   前の版は1片目だけを鍵にしていたので、5種が全部 "t" の1行に潰れていた（2026-11-18 修正）。
    # 当たり数は2通り出す。使い分けがあるため。
    #   all : MAPQ で絞らない全部。大腸菌の量(%)の換算係数 313 はこの条件で求めたので、
    #         既存149検体の大腸菌%と比べるならこちら。
    #   q20 : MAPQ 20 以上。近縁種にまたがって当たる読み取りを落とすので、
    #         「Klebsiella が本当にいるか」を見るならこちら（大腸菌の読み取りの染み出しを防ぐ）。
    # 並べ替えも索引も要らないので SAM をそのまま流す。種ごとの塩基数は固定表から引く。
    samtools view "$NOUT/$S.sam" \
      | awk -F'\t' -v LEN="$NLEN" -v Q="$MAPQ" '
          BEGIN{while ((getline l < LEN) > 0) {split(l,b,"\t"); L[b[1]]=b[2]}}
          {split($3,a,"__"); k=a[1]"__"a[2]; n[k]++; if ($5+0 >= Q) m[k]++}
          END{for(k in L) printf "%s\t%d\t%d\t%d\n", k, n[k]+0, m[k]+0, L[k]}' \
      | sort > "$NOUT/$S.counts.tsv"
    rm -f "$NOUT/$S.sam"

    sz=$(du -sm "$d" | cut -f1)
    rm -f "$d"/*.fastq.gz
    touch "$d/.mapped"
    t1=$(date +%s)
    # 大腸菌は MAPQ 無しの列（換算係数313を求めた条件）、Klebsiella は MAPQ20 の列
    neco=$(awk -F'\t' '$1=="t__SGB10068"{print $2}' "$NOUT/$S.counts.tsv")
    kp=$(awk -F'\t' '$1=="t__SGB10115"{print $3}' "$NOUT/$S.counts.tsv")
    oth=$(awk -F'\t' '$1!="t__SGB10068"{s+=$3} END{print s+0}' "$NOUT/$S.counts.tsv")
    echo "--- 終了 $S $(date) $(( (t1-t0)/60 ))分  pks=$npks 大腸菌=$neco Kp=${kp:-0} 大腸菌以外=${oth:-0}  fastq ${sz}MB を削除" >> "$LOG"
    did=1
  done

  if [ "$did" -eq 0 ]; then
    # ダウンロードが全部終わっていて、処理するものがなければ抜ける
    if [ -f "$DATA/.download_finished" ]; then
      echo "=== 全部終了 $(date)" >> "$LOG"; exit 0
    fi
    idle=$((idle+1))
    [ $((idle % 10)) -eq 1 ] && echo "$(date '+%H:%M:%S') 待機中（落ちてくるのを待っている）" >> "$LOG"
    sleep 60
  else
    idle=0
  fi
done

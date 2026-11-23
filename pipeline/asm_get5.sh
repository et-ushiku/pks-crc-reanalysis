#!/bin/bash
# 検体の fastq を落とす（md5 で検証し、ダメなら投稿元ファイルに切り替える版）。
#
# なぜ 4 つ目を作るのか（2026-11-22）
#   asm_get4.sh で分割をやめたが、219 と 224 はまだ落ちた。しかもログを読むと
#     1 回目 3014MB / 3062MB → 2 回目 3062MB / 3062MB → gzip 検査で失敗
#   **満量まで届いてから壊れている。**「途中で切れた」ではない。
#   つまり分割は原因ではなかった。asm_get4.sh の前提は間違っていた。
#
#   原因として残るのは「中断したあと -C - で継ぎ足した境目」である。
#   1 回目が SSL eof で切れたとき、最後に書いたところが中途半端だと、
#   そこから継ぎ足しても**その 1 か所だけ壊れたファイルが満量で完成する**。
#   大きさは合うので、大きさでは絶対に気づけない。
#
# 直したところ
#   1. **md5 が合わなければ継ぎ足さない。消して最初から取る。**
#      これが本当の直しである。真ん中の 1 バイトが違うファイルは、
#      続きから取っても永久に直らない。asm_get4.sh は「大きさが足りていれば
#      次へ」だったので、壊れたまま満量になったファイルを何度でも作り直さずに
#      放置していた。
#   2. **md5 で検証する。** ENA の API が fastq_md5 を返す。
#      ただし gzip -t を過小評価するのは間違いなので、はっきり書いておく。
#      **gzip -t は弱い検査ではない。** 展開した中身全体の CRC32 と長さを
#      照合するので、真ん中が壊れていればほぼ確実に落ちる。実際そうやって
#      219・224・245 を捕まえた。すでに測った 11 検体は gzip -t を通って
#      いるので、**中身が壊れている心配はしなくてよい。**
#      md5 が足すのはそこではなく、「完全に正しい gzip だが**別の中身**」を
#      捕まえられる点である（サーバが違うオブジェクトを返した、など）。
#      外部が宣言した値と照合するので、検査の基準が自分の手元から離れる。
#      CORRECTIONS.md の 6 番に書いた教訓と同じ向きの直しだが、
#      **今回の失敗の原因ではない。**
#   3. **それでもダメなら submitted_ftp に切り替える。**
#      ENA が生成した fastq とは別の実体（投稿者が出した元ファイル）で、
#      サーバ上の別オブジェクトである。中身の読み取りは同じもの。
#
# 接続は 1 本だけ。
#
# 使い方: bash asm_get5.sh 検体名 ERR番号
set -u
S=$1 ; E=$2
B=/mnt/c/crc_strain_pilot
OUT=$B/asm/$S
LOG=$B/meta/asm_get.log
API="https://www.ebi.ac.uk/ena/portal/api/filereport"
TRIES=6
mkdir -p "$OUT"
echo "=== 開始 $S ($E) $(date) [asm_get5 md5で検証]" >> "$LOG"

TSV=""
for t in $(seq 1 5); do
  TSV=$(curl -sS --max-time 120 \
    "$API?accession=$E&result=read_run&fields=fastq_ftp,fastq_bytes,fastq_md5,submitted_ftp,submitted_bytes,submitted_md5&format=tsv")
  [ -n "$TSV" ] && break
  echo "  ENA の応答が空（$t 回目）" >> "$LOG"; sleep 20
done
[ -z "$TSV" ] && { echo "=== 失敗 $S（API の応答が空）" >> "$LOG"; exit 1; }

col() { printf '%s\n' "$TSV" | awk -F'\t' -v want="$1" '
  NR==1 { for (i = 1; i <= NF; i++) if ($i == want) c = i; next }
  NR==2 { if (c) print $c; exit }'; }

# ひとつの系列（生成 fastq か投稿ファイル）を試す。
# 成功したら 0、ダメなら 1 を返す。
try_series() {
  local URLS=$1 BYTS=$2 MD5S=$3 NAME=$4
  local n i u w m F have got
  [ -z "$URLS" ] && { echo "  $NAME は URL が無い" >> "$LOG"; return 1; }
  [ -z "$MD5S" ] && { echo "  $NAME は md5 が無い。検証できないので使わない" >> "$LOG"; return 1; }
  n=$(printf '%s\n' "$URLS" | tr ';' '\n' | grep -c .)
  echo "  $NAME を試す（$n ファイル・1本ずつ・md5 で検証）" >> "$LOG"
  i=1
  while [ "$i" -le "$n" ]; do
    u=$(printf '%s\n' "$URLS" | cut -d';' -f"$i")
    w=$(printf '%s\n' "$BYTS" | cut -d';' -f"$i")
    m=$(printf '%s\n' "$MD5S" | cut -d';' -f"$i")
    F=$OUT/$(basename "$u")
    case "$m" in ''|*[!0-9a-f]*) echo "  md5 が読めない（$m）" >> "$LOG"; return 1 ;; esac
    if [ -s "$F" ]; then
      got=$(md5sum "$F" | cut -d' ' -f1)
      if [ "$got" = "$m" ]; then
        echo "  $(basename "$F") はすでにある（md5 一致）" >> "$LOG"
        i=$((i+1)); continue
      fi
    fi
    echo "  $(date '+%H:%M:%S') $(basename "$F") $((w/1048576))MB  md5 $m" >> "$LOG"
    local ok=0
    for t in $(seq 1 $TRIES); do
      have=0 ; [ -f "$F" ] && have=$(stat -c%s "$F")
      # 大きさが足りていないあいだだけ継ぎ足す。満量なら一度 md5 を見る
      if [ "$have" -lt "$w" ]; then
        curl -sS -L -C - --connect-timeout 30 --max-time 10800 -o "$F" "https://$u" >> "$LOG" 2>&1
        have=0 ; [ -f "$F" ] && have=$(stat -c%s "$F")
      fi
      if [ "$have" -ge "$w" ]; then
        got=$(md5sum "$F" | cut -d' ' -f1)
        if [ "$got" = "$m" ]; then
          echo "    $t 回目 $((have/1048576))MB md5 一致" >> "$LOG"
          ok=1; break
        fi
        # 満量なのに md5 が違う = どこかが壊れている。継ぎ足しでは直らない
        echo "    $t 回目 満量だが md5 が違う（$got）。消して最初から取る" >> "$LOG"
        rm -f "$F"
      else
        echo "    $t 回目 $((have/1048576))MB / $((w/1048576))MB（途中）" >> "$LOG"
      fi
      sleep 15
    done
    [ "$ok" -ne 1 ] && { echo "  **$(basename "$F") は md5 が合わない**" >> "$LOG"; return 1; }
    i=$((i+1))
  done
  return 0
}

if try_series "$(col fastq_ftp)" "$(col fastq_bytes)" "$(col fastq_md5)" "ENA生成(fastq_ftp)"; then
  echo "=== 完了 $S $(date) ENA生成 md5 検証OK" >> "$LOG"
  touch "$OUT/.complete"; exit 0
fi

echo "  ENA生成が md5 で通らない。**投稿元ファイルに切り替える**" >> "$LOG"
rm -f "$OUT"/*.fastq.gz
if try_series "$(col submitted_ftp)" "$(col submitted_bytes)" "$(col submitted_md5)" "投稿ファイル(submitted_ftp)"; then
  echo "=== 完了 $S $(date) 投稿ファイル md5 検証OK" >> "$LOG"
  touch "$OUT/.complete"; exit 0
fi

echo "=== 失敗 $S（両方の系列が md5 で通らない）" >> "$LOG"
exit 1

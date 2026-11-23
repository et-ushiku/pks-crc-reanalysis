#!/bin/bash
# 検体の fastq を落とす（asm_get2.sh の弱点を 1 つ直した版）。
#
# asm_get2.sh で 12 検体のうち 7 検体が落ちた理由（2026-11-21）
#   URL は直っていた（ENA の API に聞いているので 00<下1桁> で正しい）。
#   落ちたのは **大きさの取り方**である。grab() は
#       curl -sSI -L "$U" | awk '…content-length…'
#   の 1 回だけで大きさを決め、取れなければ即 return 1 していた。
#   ENA は HEAD を返さないことがあり（同じ時間帯に OpenSSL の
#   "unexpected eof while reading" が多発している）、
#   **cancer_ERR1018219 は _1 を 3,063MB 落としきったあと、_2 の HEAD が
#   返らなかっただけで検体ごと捨てられた。** 3GB が無駄になっている。
#
#   直し方は 2 つある。
#     1. HEAD を数回やり直す
#     2. **そもそも大きさは ENA の API（fastq_bytes）で既に分かっている**
#   2 を主、1 を予備にする。API の数字があれば HEAD は要らない。
#
# それ以外は asm_get2.sh と同じ（見出し名で列を引く / URL は 1 行 1 件の
# ファイルに展開して数えるのは grep -c . / 接続は 4 本まで）。
#
# 使い方: bash asm_get3.sh 検体名 ERR番号
set -u
S=$1 ; E=$2
B=/mnt/c/crc_strain_pilot
OUT=$B/asm/$S
LOG=$B/meta/asm_get.log
PARTS=4
API="https://www.ebi.ac.uk/ena/portal/api/filereport"
mkdir -p "$OUT"
URLFILE=$OUT/.urls          # 1 行 = URL<TAB>バイト数
echo "=== 開始 $S ($E) $(date) [asm_get3]" >> "$LOG"

TSV=""
for t in 1 2 3 4 5; do
  TSV=$(curl -sS --max-time 120 \
    "$API?accession=$E&result=read_run&fields=fastq_ftp,fastq_bytes,submitted_ftp,submitted_bytes&format=tsv")
  [ -n "$TSV" ] && break
  echo "  ENA の応答が空（$t 回目）。やり直す" >> "$LOG"
  sleep 20
done
if [ -z "$TSV" ]; then
  echo "=== 失敗 $S（ENA の応答が 5 回とも空）" >> "$LOG"; exit 1
fi

col() {   # col <見出し名>
  printf '%s\n' "$TSV" | awk -F'\t' -v want="$1" '
    NR==1 { for (i = 1; i <= NF; i++) if ($i == want) c = i; next }
    NR==2 { if (c) print $c; exit }'
}
F_FTP=$(col fastq_ftp) ; F_BYT=$(col fastq_bytes)
S_FTP=$(col submitted_ftp) ; S_BYT=$(col submitted_bytes)

pick_ok() {   # pick_ok <ftp列> <bytes列>  …全ファイルが 1MB 以上なら採用
  local ftp=$1 byt=$2 n i v
  [ -z "$ftp" ] && return 1
  n=$(printf '%s\n' "$ftp" | tr ';' '\n' | grep -c .)
  [ "$n" -lt 1 ] && return 1
  i=1
  while [ "$i" -le "$n" ]; do
    v=$(printf '%s\n' "$byt" | cut -d';' -f"$i")
    case "$v" in ''|*[!0-9]*) return 1 ;; esac
    [ "$v" -lt 1048576 ] && return 1
    i=$((i+1))
  done
  return 0
}

if pick_ok "$F_FTP" "$F_BYT"; then
  URLS=$F_FTP ; BYTS=$F_BYT ; WHICH="ENA生成(fastq_ftp)"
elif pick_ok "$S_FTP" "$S_BYT"; then
  URLS=$S_FTP ; BYTS=$S_BYT ; WHICH="投稿ファイル(submitted_ftp)"
  echo "  ENA 生成ファイルが怪しいので投稿ファイルに切り替えた" >> "$LOG"
else
  echo "  どちらの系列も健全でない fastq_bytes=[$F_BYT] submitted_bytes=[$S_BYT]" >> "$LOG"
  echo "=== 失敗 $S" >> "$LOG"; exit 1
fi

# URL とバイト数を 1 行 1 件に組む
: > "$URLFILE"
NU=$(printf '%s\n' "$URLS" | tr ';' '\n' | grep -c .)
i=1
while [ "$i" -le "$NU" ]; do
  u=$(printf '%s\n' "$URLS" | cut -d';' -f"$i")
  w=$(printf '%s\n' "$BYTS" | cut -d';' -f"$i")
  printf '%s\t%s\n' "$u" "$w" >> "$URLFILE"
  i=$((i+1))
done
NFILE=$(grep -c . "$URLFILE")
echo "  $WHICH を使う（$NFILE ファイル、大きさは API の fastq_bytes）" >> "$LOG"
[ "$NFILE" -lt 1 ] && { echo "=== 失敗 $S（URL が 0 件）" >> "$LOG"; exit 1; }

fetch_part() {   # fetch_part <URL> <出力> <開始> <終了>
  local U=$1 F=$2 a=$3 b=$4 want have t
  want=$((b - a + 1))
  for t in $(seq 1 12); do
    have=0
    [ -f "$F" ] && have=$(stat -c%s "$F")
    [ "$have" -ge "$want" ] && return 0
    curl -sS -L -r "$((a + have))-$b" --connect-timeout 30 --max-time 7200 "$U" >> "$F"
    sleep 5
  done
  have=0
  [ -f "$F" ] && have=$(stat -c%s "$F")
  [ "$have" -ge "$want" ]
}

head_len() {   # 予備: HEAD で大きさを聞く（3 回まで）
  local U=$1 L t
  for t in 1 2 3; do
    L=$(curl -sSI -L --connect-timeout 30 --max-time 120 "$U" \
        | tr -d '\r' | awk 'tolower($1)=="content-length:"{print $2}' | tail -1)
    case "$L" in ''|*[!0-9]*) sleep 10 ;; *) printf '%s' "$L"; return 0 ;; esac
  done
  return 1
}

grab() {   # grab <URL> <出力> <APIのバイト数>
  local U=$1 F=$2 LEN=$3 CH p a b i fail
  local pids=() ps=()
  case "$LEN" in ''|*[!0-9]*) LEN=$(head_len "$U") || {
      echo "  大きさが取れない（API も HEAD も） $U" >> "$LOG"; return 1; } ;;
  esac
  if [ "$LEN" -lt 1048576 ]; then
    echo "  小さすぎる（$LEN バイト） $U" >> "$LOG"; return 1
  fi
  echo "  $(date '+%H:%M:%S') $(basename "$F") $((LEN/1048576))MB を${PARTS}分割で取る" >> "$LOG"
  CH=$(( (LEN + PARTS - 1) / PARTS ))
  for p in $(seq 0 $((PARTS-1))); do
    a=$((p*CH)); b=$(( (p+1)*CH - 1 ))
    [ "$b" -ge "$LEN" ] && b=$((LEN-1))
    fetch_part "$U" "$F.part$p" "$a" "$b" &
    pids+=($!); ps+=("$p")
  done
  i=0; fail=0
  for pid in "${pids[@]}"; do
    wait "$pid" || { echo "  ${ps[$i]} 番目の分割が失敗" >> "$LOG"; fail=1; }
    i=$((i+1))
  done
  [ "$fail" -eq 1 ] && { rm -f "$F".part*; return 1; }
  cat "$F".part* > "$F" && rm -f "$F".part*
  if gzip -t "$F" 2>/dev/null; then
    echo "  $(date '+%H:%M:%S') $(basename "$F") 完了 $(du -m "$F"|cut -f1)MB gzip 検査OK" >> "$LOG"
    return 0
  fi
  echo "  $(date '+%H:%M:%S') $(basename "$F") gzip が壊れている" >> "$LOG"
  rm -f "$F"; return 1
}

FAIL=0
while IFS=$'\t' read -r u w; do
  [ -z "$u" ] && continue
  bn=$(basename "$u")
  if [ -s "$OUT/$bn" ] && gzip -t "$OUT/$bn" 2>/dev/null; then
    echo "  $bn はすでにある" >> "$LOG"
    continue
  fi
  grab "https://$u" "$OUT/$bn" "$w" || { FAIL=1; break; }
done < "$URLFILE"

[ "$FAIL" -ne 0 ] && { echo "=== 失敗 $S" >> "$LOG"; exit 1; }

NG=0
while IFS=$'\t' read -r u w; do
  [ -z "$u" ] && continue
  bn=$(basename "$u")
  if [ -s "$OUT/$bn" ] && gzip -t "$OUT/$bn" 2>/dev/null; then
    NG=$((NG+1))
  else
    echo "  揃っていない: $bn" >> "$LOG"
  fi
done < "$URLFILE"

if [ "$NG" -ne "$NFILE" ]; then
  echo "=== 失敗 $S（$NG/$NFILE ファイルしか揃わなかった）" >> "$LOG"; exit 1
fi
touch "$OUT/.complete"
echo "=== 完了 $S $(date) $NG/$NFILE ファイル" >> "$LOG"

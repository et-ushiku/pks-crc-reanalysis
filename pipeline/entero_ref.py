"""pks 島を持ちうる腸内細菌科4種の「目印配列」を MetaPhlAn の索引から抜き出す。

なぜ必要か
  pks 島は大腸菌だけのものではない。Putze 2009（Infect Immun、PMID 19720753）は
  1,565株の PCR 調査で、Klebsiella pneumoniae 5/141、Enterobacter aerogenes 3/11、
  Citrobacter koseri 1/1 に島を検出した。島は ICE 様の可動性要素に乗って
  種を越えて水平伝達している。
  → 「大腸菌 ≥0.1%」で判定を条件づけると、Klebsiella 等で pks を持つ検体を落とす。

なぜ extract_markers.py を使わないか
  あれは1クレードごとに索引(.pkl)を読み直すので、E. coli 1種で1時間20分かかった。
  4種で5時間超になる。ここでは索引を1回読み、目印配列の fasta を1回流すだけにする。
  （中身は同じことをしている。E. coli で前回と同じ結果になるかを検算に使う）

出力
  <OUT>/entero.fna     4種の目印配列。名前の先頭に SGB 番号を付け、後で種別に数えられるようにする
  <OUT>/entero_idx.*   bowtie2 の索引
  <OUT>/entero_map.tsv 配列名 → SGB の対応表
"""
import bz2
import os
import pickle
import subprocess
import sys
import time

DB = "/root/metaphlan_db/mpa_vJan26_CHOCOPhlAnSGB_202605.pkl"
FNA = "/root/metaphlan_db/mpa_vJan26_CHOCOPhlAnSGB_202605.fna"
OUT = "/mnt/c/crc_strain_pilot/pks/ref/entero"

# Putze 2009 が pks 島を検出した種 + 大腸菌（検算用）
TARGETS = {
    "t__SGB10068": "Escherichia_coli",
    "t__SGB10115": "Klebsiella_pneumoniae",
    "t__SGB10079": "Citrobacter_koseri",
    "t__SGB10116": "Klebsiella_aerogenes",
    "t__SGB165717": "Klebsiella_aerogenes_2",
}
# 前回 extract_markers.py で得た値（検算の基準）
ECOLI_EXPECT = (122, 96300)


def log(m):
    print(f'{time.strftime("%H:%M:%S")} {m}', flush=True)


def load(path):
    try:
        with bz2.open(path, "rb") as fh:
            return pickle.load(fh)
    except OSError:
        with open(path, "rb") as fh:
            return pickle.load(fh)


def main():
    os.makedirs(OUT, exist_ok=True)
    log("索引(.pkl)を読む")
    db = load(DB)
    markers = db.get("markers", {})
    log(f"目印の総数 {len(markers):,}")

    want = {}          # 配列名 -> SGB
    per = {k: 0 for k in TARGETS}
    for name, info in markers.items():
        cl = info.get("clade") if isinstance(info, dict) else None
        if cl in TARGETS:
            want[name] = cl
            per[cl] += 1
    log("種ごとの目印の数:")
    for sgb, sp in TARGETS.items():
        log(f"  {sgb:<14} {sp:<24} {per[sgb]:>5} 本")
    if not want:
        log("目印が見つからない。索引の構造が想定と違う")
        sys.exit(1)

    log(f"目印配列の fasta を流して {len(want):,} 本を抜き出す（21.5GB、数分）")
    seen, wrote, bases = set(), 0, {k: 0 for k in TARGETS}
    keep = False
    cur_sgb = None
    with open(FNA, encoding="utf-8", errors="ignore") as fin, \
            open(f"{OUT}/entero.fna", "w", encoding="utf-8") as fout, \
            open(f"{OUT}/entero_map.tsv", "w", encoding="utf-8") as fmap:
        fmap.write("seq_name\tsgb\tspecies\n")
        for line in fin:
            if line.startswith(">"):
                nm = line[1:].strip().split()[0]
                cur_sgb = want.get(nm)
                keep = cur_sgb is not None
                if keep:
                    seen.add(nm)
                    wrote += 1
                    # 名前の先頭に SGB を付ける（後で種別に数えるため）
                    fout.write(f">{cur_sgb}__{nm}\n")
                    fmap.write(f"{cur_sgb}__{nm}\t{cur_sgb}\t{TARGETS[cur_sgb]}\n")
            elif keep:
                fout.write(line)
                bases[cur_sgb] += len(line.strip())

    log(f"抜き出した配列 {wrote:,} 本")
    miss = set(want) - seen
    if miss:
        log(f"※ fasta に見つからなかった目印 {len(miss)} 本（索引の版ずれ）")
    log("種ごとの塩基数:")
    for sgb, sp in TARGETS.items():
        log(f"  {sgb:<14} {sp:<24} {per[sgb]:>5} 本  {bases[sgb]:>10,} 塩基")

    # 検算: 大腸菌が前回の extract_markers.py と一致するか
    n_e, b_e = per["t__SGB10068"], bases["t__SGB10068"]
    if (n_e, b_e) == ECOLI_EXPECT:
        log(f"検算OK: 大腸菌 {n_e}本 / {b_e:,}塩基 は前回と一致")
    else:
        log(f"検算: 大腸菌 {n_e}本 / {b_e:,}塩基。前回は "
            f"{ECOLI_EXPECT[0]}本 / {ECOLI_EXPECT[1]:,}塩基 → **不一致**。"
            "換算係数313をそのまま使えるか確認が必要")

    log("bowtie2 の索引を作る")
    r = subprocess.run(["bowtie2-build", "--threads", "4",
                        f"{OUT}/entero.fna", f"{OUT}/entero_idx"],
                       capture_output=True, text=True)
    if r.returncode != 0:
        log(f"索引作成に失敗: {r.stderr[-500:]}")
        sys.exit(1)
    log(f"完了 → {OUT}/entero_idx")


if __name__ == "__main__":
    main()

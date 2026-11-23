"""6種の完全ゲノムに当てた結果から、pks 島の担い手を絞る。

読み方
  目印配列ではなく完全ゲノム（1種あたり約500万塩基）に当てているので、
  深さと覆われ方の両方が見える。

    平均深さ  その種のゲノムが何重に読まれているか
    覆い      ゲノムのどれだけが1本以上に覆われたか
              ほんとうにいる種なら、深さ1でも覆い0.6前後にはなる。
              覆いが極端に低いのに深さだけある＝近縁種からの染み出し

  pks 島の深さと**同じ深さ**の種が担い手の候補。
  島は可動性要素なので、持っている株と持っていない株が同じ種の中に混ざる。
  したがって「島の深さ <= その種の深さ」が条件で、比が1に近いほど
  その種の集団のほとんどが島を持っていることになる。

使い方: python3 asm_report.py 検体名
"""
import os
import sys

B = "/mnt/c/crc_strain_pilot"
OUT = B + "/pks/asm_out"
PKS_BP = 49455

NAME = {
    "Ecoli_IHE3034": "大腸菌 (E. coli IHE3034)",
    "Kpneumoniae_HS11286": "K. pneumoniae HS11286",
    "Ckoseri_BAA895": "C. koseri BAA-895",
    "Kaerogenes_KCTC2190": "K. aerogenes KCTC2190",
    "Ecloacae_ATCC13047": "E. cloacae ATCC13047",
    "Cfreundii_CFNIH1": "C. freundii CFNIH1",
}


def read_total(log):
    if not os.path.exists(log):
        return None
    with open(log, encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if "reads; of these:" in line:
                return int(line.split()[0])
    return None


def main():
    s = sys.argv[1]
    g6 = OUT + "/" + s + ".g6.depth.tsv"
    if not os.path.exists(g6):
        print("まだ結果がない: " + g6)
        return
    ntot = read_total(OUT + "/" + s + ".g6.bowtie2.log")

    # pks 島のほう（同じ読み取り集合で測ったもの）
    pks_depth = pks_reads = None
    dep = OUT + "/" + s + ".depth.tsv"
    if os.path.exists(dep):
        v = []
        with open(dep, encoding="utf-8") as fh:
            for line in fh:
                p = line.split("\t")
                if len(p) >= 3:
                    v.append(int(p[2]))
        if v:
            pks_depth = sum(v) / len(v)

    rows = []
    with open(g6, encoding="utf-8") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 4:
                continue
            key = p[0]
            tot, cov, ssum = int(p[1]), int(p[2]), int(p[3])
            rows.append((key, tot, cov, ssum))

    print("検体 %s   貼り付けた読み取り %s" % (s, "%d" % ntot if ntot else "不明"))
    if pks_depth is not None:
        print("pks 島の平均深さ %.2f（同じ読み取り集合、MAPQ>=20）" % pks_depth)
    print()
    print("%-26s %10s %8s %10s %12s"
          % ("種", "平均深さ", "覆い", "覆われた塩基", "島の深さ/この種"))
    rows.sort(key=lambda r: -(r[3] / r[1] if r[1] else 0))
    for key, tot, cov, ssum in rows:
        d = ssum / tot if tot else 0.0
        br = cov / tot if tot else 0.0
        ratio = (pks_depth / d) if (pks_depth and d > 0.01) else None
        print("%-26s %10.3f %8.3f %10s %12s"
              % (NAME.get(key, key), d, br, "%d" % cov,
                 ("%.2f" % ratio) if ratio is not None else "-"))

    print()
    print("判定の目安")
    print("  比が 1 に近い種 → その種の集団のほとんどが島を持つ。担い手の第一候補")
    print("  比が 1 より大きい → その種だけでは島の深さを説明できない")
    print("  比が 1 より小さい → その種はいるが、島を持つのは一部の株だけ")
    print("  覆いが 0.1 未満の種は、近縁種からの染み出しを疑う（本当にはいない）")


if __name__ == "__main__":
    main()

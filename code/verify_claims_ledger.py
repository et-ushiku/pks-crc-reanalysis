"""CLAIMS.md の、機械で確かめられるところだけを確かめる。

確かめられるのは
  - 台帳のどの行にも「何を排除するか」の欄が埋まっていること
  - 等級が A/B/C/D のどれかであること
  - 引いている data/ code/ pipeline/ のファイルが実在すること
  - 引いている CORRECTIONS.md の項目番号が実在すること

確かめられないのは
  **「何を排除するか」の欄が正直かどうか**。これは人が読むしかない。
  この研究が 4 回間違えたのはそこで、機械化できないと書いてあるとおりである。

使い方: python3 verify_claims_ledger.py
"""
import io
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
GRADES = ("A", "B", "C", "D")
BAD = []


def chk(name, ok, detail=""):
    print("  %-4s %s%s" % ("OK" if ok else "NG", name,
                           "" if ok else "   <- " + detail))
    if not ok:
        BAD.append(name)


def main():
    path = os.path.join(ROOT, "CLAIMS.md")
    s = io.open(path, encoding="utf-8").read()

    print("1. 台帳の行")
    rows = []
    for line in s.split("\n"):
        if not line.startswith("| "):
            continue
        c = [x.strip() for x in line.strip().strip("|").split("|")]
        # 台帳の行は 5 列。等級の表（2 列）や区切り行は除く
        if len(c) != 5 or set(c[0]) <= set("-: ") or c[0] == "#":
            continue
        rows.append(c)
    chk("行が 10 行以上ある", len(rows) >= 10, "%d 行" % len(rows))
    print("     %d 行を読んだ" % len(rows))

    for num, claim, ev, excl, grade in rows:
        tag = "行 %s" % num
        chk("%s に主張が書いてある" % tag, len(claim) > 10, claim)
        chk("%s に根拠が書いてある" % tag, len(ev) > 10, ev)
        # ここが本題。空でないだけでなく、中身が無い定型句も弾く
        bare = re.sub(r"[\s*`.\-]", "", excl).lower()
        chk("%s の『何を排除するか』が埋まっている" % tag,
            len(bare) > 15 and bare not in ("na", "none", "nothing", "tbd"),
            repr(excl))
        m = re.match(r"\*\*([A-D])\*\*", grade)
        chk("%s の等級が A/B/C/D" % tag, bool(m), grade)
        if m:
            chk("%s は等級 %s" % (tag, m.group(1)), m.group(1) in GRADES)

    print()
    print("2. 引いているファイルが実在するか")
    paths = sorted(set(re.findall(
        r"`((?:data|code|pipeline)/[A-Za-z0-9_./-]+)`", s)))
    for p in paths:
        chk(p, os.path.exists(os.path.join(ROOT, p)))
    chk("引いたファイルが 5 本以上ある", len(paths) >= 5, "%d 本" % len(paths))

    print()
    print("3. 引いている CORRECTIONS.md の項目番号が実在するか")
    corr = io.open(os.path.join(ROOT, "CORRECTIONS.md"),
                   encoding="utf-8").read()
    have = set(re.findall(r"^## (\d+)\.", corr, re.M))
    cited = sorted(set(re.findall(r"item (\d+)", s)), key=int)
    for n in cited:
        chk("item %s" % n, n in have, "CORRECTIONS.md に %s 番が無い" % n)
    chk("項目を 4 つ以上引いている", len(cited) >= 4, "%d 個" % len(cited))

    print()
    print("4. 等級ごとの数")
    g = {}
    for r in rows:
        m = re.match(r"\*\*([A-D])\*\*", r[4])
        if m:
            g[m.group(1)] = g.get(m.group(1), 0) + 1
    for k in GRADES:
        print("     %s: %d 行" % (k, g.get(k, 0)))
    # C と D が 1 行も無い台帳は、採点していないのと同じ
    chk("C か D の行がある（全部 A の台帳は採点になっていない）",
        g.get("C", 0) + g.get("D", 0) > 0)

    print()
    if BAD:
        print("%d 件が駄目" % len(BAD))
        for b in BAD:
            print("  NG %s" % b)
        sys.exit(1)
    print("台帳の機械で見られるところはすべて通った。"
          "『何を排除するか』が正直かどうかは人が読むこと。")


if __name__ == "__main__":
    main()

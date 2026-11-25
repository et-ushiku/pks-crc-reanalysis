"""Recompute the duplication in AM229678.1 from the public sequence.

PITFALLS.md item 1 rests on facts about the reference itself - a 1,287 bp repeat at
99.845% identity, the 5,653 bp spacing between its two copies, and the absence of any
other large repeat. The reference is not shipped here (it belongs to GenBank/ENA and
this repository redistributes no sequence), so those numbers could not be checked from
`data/` alone. This script fetches AM229678.1 and recomputes them.

    python3 code/verify_repeat.py            # downloads ~55 kb from ENA
    python3 code/verify_repeat.py file.fa    # or point it at a local copy

**This is the only script here that needs network access.** Everything else runs from
`data/`. The sequence is fetched from the ENA browser API over https; nothing is sent
but the accession.

Method, and its limits: the two copies are co-linear at a fixed offset, so position i
is compared with position i+offset directly. That is exact for substitutions and
blind to indels, so an exhaustive 30-mer scan is run as well - it assumes nothing
about co-linearity and would reveal a repeat continued at a shifted offset.
"""
import os
import sys
from collections import defaultdict

ACC = "AM229678.1"
URL = "https://www.ebi.ac.uk/ena/browser/api/fasta/" + ACC
OFFSET = 5653        # spacing between the two copies, measured below
K = 30


def get_sequence(path=None):
    if path:
        txt = open(path, encoding="utf-8").read()
    else:
        import urllib.request
        print("fetching %s from ENA ..." % ACC)
        with urllib.request.urlopen(URL, timeout=120) as fh:
            txt = fh.read().decode("ascii", "replace")
    seq = "".join(l.strip() for l in txt.splitlines() if l and not l.startswith(">"))
    return seq.upper()


def runs_at(seq, offset, thresh, win=50, minlen=200):
    n = len(seq) - offset
    match = [1 if seq[i] == seq[i + offset] else 0 for i in range(n)]
    ok = [sum(match[i:i + win]) / float(win) >= thresh for i in range(n - win + 1)]
    out, start = [], None
    for i, v in enumerate(ok):
        if v and start is None:
            start = i
        elif not v and start is not None:
            out.append((start + 1, i + win))
            start = None
    if start is not None:
        out.append((start + 1, len(ok) + win))
    res = []
    for a, b in out:
        if b - a + 1 < minlen:
            continue
        m = sum(match[a - 1:b])
        res.append((a, b, b - a + 1, 100.0 * m / (b - a + 1)))
    return res


def main():
    seq = get_sequence(sys.argv[1] if len(sys.argv) > 1 else None)
    print("sequence length %d   (PITFALLS.md: 55,140)" % len(seq))
    print()

    print("Windowed identity at offset %d:" % OFFSET)
    for th in (0.98, 0.95, 0.70):
        for a, b, ln, idp in runs_at(seq, OFFSET, th):
            print("   >= %.0f%%  copy A %d-%d / copy B %d-%d   %d bp   %.3f%% identical"
                  % (th * 100, a, b, a + OFFSET, b + OFFSET, ln, idp))
    print("   PITFALLS.md: 1,292 bp at >=98%, 1,337 bp at >=70%")
    print()

    print("Exhaustive %d-mer scan (assumes nothing about co-linearity):" % K)
    idx = defaultdict(list)
    for i in range(len(seq) - K + 1):
        idx[seq[i:i + K]].append(i + 1)
    dup = {k: v for k, v in idx.items() if len(v) > 1}
    spacing = defaultdict(int)
    at_off = []
    for v in dup.values():
        for a in range(len(v)):
            for b in range(a + 1, len(v)):
                d = v[b] - v[a]
                spacing[d] += 1
                if d == OFFSET:
                    at_off.append(v[a])
    print("   %d-mers occurring more than once: %d" % (K, len(dup)))
    top = sorted(spacing.items(), key=lambda x: -x[1])[:4]
    for d, n in top:
        note = "  <- the clbJ/clbK repeat" if d == OFFSET else \
               ("  <- the 8 bp VNTR between clbB and clbR" if d % 8 == 0 else "")
        print("   spacing %6d : %5d pairs%s" % (d, n, note))
    if at_off:
        lo, hi = min(at_off), max(at_off) + K - 1
        print("   repeat extent from exact %d-mers: %d-%d  (%d bp)"
              % (K, lo, hi, hi - lo + 1))
        m = sum(1 for i in range(lo, hi + 1) if seq[i - 1] == seq[i - 1 + OFFSET])
        print("   identity over that span: %.3f%%" % (100.0 * m / (hi - lo + 1)))
    other = sorted(d for d in spacing if d != OFFSET and 5000 < d < 6500)
    print("   other spacings between 5,000 and 6,500 bp: %s"
          % (", ".join(str(d) for d in other) if other else "none"))
    print("   PITFALLS.md: 1,287 bp at 99.845%, 1,198 of 1,206 pairs at this spacing,")
    print("   no other spacing in that range, and Auvray et al.'s 1,480 bp not reproducible")


if __name__ == "__main__":
    main()

"""Check the analytic rarefaction against actually rarefied reads.

The whole repository rests on one claim: that rarefying can be done on the depth
tables instead of on the reads, because

    E[breadth | f] = mean_i [ 1 - (1-f)^d_i ]

is exact by linearity of expectation. That makes the per-gene depth histogram
sufficient, which is why 753 KB of `data/depth_hist.tsv` replaces 2.5 TB of FASTQ.

A derivation is not a check. For two samples we also did it the slow way: subsampled
the FASTQ with `seqtk sample`, remapped, and recomputed breadth from the result. Five
such runs are in `data/rarefy_check_hist.tsv`, as per-gene depth histograms in the
same form as `data/depth_hist.tsv`. This script puts the two side by side.

What to look for: the analytic value is an *expectation*, so it should land near the
actual value, not on it. The useful comparison is the gap against the spread between
seeds for the same target - three runs at 40M reads on the same sample are included
for exactly that reason.

    python3 code/verify_rarefy.py

No arguments, no dependencies beyond the standard library.
"""
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import analyse as A  # noqa: E402


def load_actual():
    """data/rarefy_check_hist.tsv -> {(sample, target, seed): {gene: {depth: bases}}}"""
    out = defaultdict(lambda: defaultdict(dict))
    fh = open(os.path.join(A.DATA, "rarefy_check_hist.tsv"), encoding="utf-8")
    fh.readline()
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) < 6:
            continue
        out[(p[0], int(p[1]), int(p[2]))][p[3]][int(p[4])] = int(p[5])
    return out


def main():
    reg = A.load_bed()
    hist = A.load_hist()
    total = {r["sample"]: int(r["total_reads"]) for r in A.load_samples()}
    actual = load_actual()

    print("analytic rarefaction vs actually rarefied reads")
    print("thresholds are the pre-specified ones in analyse.py\n")
    print("%-24s %10s %5s %7s | %8s %8s %8s | %6s %6s | %s"
          % ("sample", "target", "seed", "f", "full", "analytic", "actual",
             "genes", "genes", "call agrees"))
    print("-" * 118)

    gaps = []
    by_target = defaultdict(list)
    for (s, target, seed) in sorted(actual):
        h = hist.get(s)
        if h is None:
            print("%-24s  no full-depth histogram for this sample" % s)
            continue
        n = total[s]
        f = target / n
        fb, fg, _ = A.call_from(h, reg, 1.0)
        pb, pg, pc = A.call_from(h, reg, f)
        ab, ag, ac = A.call_from(actual[(s, target, seed)], reg, 1.0)
        gaps.append(pb - ab)
        by_target[(s, target)].append(ab)
        print("%-24s %10d %5d %7.4f | %8.4f %8.4f %8.4f | %6d %6d | %s %s"
              % (s, target, seed, f, fb, pb, ab, pg, ag,
                 "yes" if pc == ac else "NO",
                 "(%s vs %s)" % (pc, ac) if pc != ac else ""))

    print()
    print("largest gap  |analytic - actual|  = %.4f over %d runs"
          % (max(abs(g) for g in gaps), len(gaps)))
    print("sign of gap: analytic below actual in %d of %d runs"
          % (sum(1 for g in gaps if g < 0), len(gaps)))
    for (s, target), vals in sorted(by_target.items()):
        if len(vals) > 1:
            print("seed-to-seed spread of the ACTUAL breadth at %s reads, %s:"
                  % (target, s))
            print("   %s   spread %.4f  mean %.4f"
                  % (" ".join("%.4f" % v for v in sorted(vals)),
                     max(vals) - min(vals), sum(vals) / len(vals)))
    print()
    print("Read the gap against that spread. A gap smaller than the spread between")
    print("seeds means the analytic value is inside the range that actually rarefying")
    print("would have produced, which is all an expectation can be asked to do.")
    print("The calls are what the analysis uses, and those are compared too.")


if __name__ == "__main__":
    main()

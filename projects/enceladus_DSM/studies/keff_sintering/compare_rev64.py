#!/usr/bin/env python3
"""Does the L/R_ave = 40 domain bias k_eff? Compare against one run at 64.

    venv_enceladus/bin/python studies/keff_sintering/compare_rev64.py \
        <pilot_batch_dir> <rev64_run_dir>

WHAT ONE RUN CAN AND CANNOT SETTLE. A single realisation at the larger size
cannot separate a domain-size bias from ordinary seed scatter -- that would
need an ensemble at 64 too. What it CAN do is say whether the larger domain
lands inside the smaller one's spread. If it lands far outside, the two sizes
are sampling different distributions and 40 is biased; if it lands inside,
there is no evidence of bias and 40 is adequate.

The comparison is made twice on purpose, because the two answers can differ:

  ABSOLUTE   k_eff itself. This is what a parameterisation would report.
  RELATIVE   k_eff(t)/k_eff(0), the sintering response. A domain-size offset
             that is present at t = 0 and unchanged afterwards cancels here,
             so a claim about how much sintering RAISES k_eff can survive a
             domain that is biased in the absolute value.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).parent
PROJ = HERE.parents[1]
sys.path.insert(0, str(PROJ / "preprocess"))
sys.path.insert(0, str(PROJ / "postprocess"))
import figstyle as fs                                       # noqa: E402
import pplib                                                # noqa: E402
sys.path.insert(0, str(HERE))
from analyze_pilot import find_runs, load_run               # noqa: E402

DAY = 86400.0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pilot", type=Path)
    ap.add_argument("rev64", type=Path)
    a = ap.parse_args()

    small = {s: load_run(d) for s, d in find_runs(a.pilot).items()}
    big_dir = next(iter(find_runs(a.rev64).values()))
    big = load_run(big_dir)

    # common time window: the rev64 run may stop slightly short
    t_end = min(min(r["kt"].max() for r in small.values()), big["kt"].max())
    print(f"comparing on the common window 0 .. {t_end/DAY:.2f} d\n")

    def at(r, t, key="k_iso"):
        return float(np.interp(t, r["kt"], r[key]))

    k0s = np.array([at(r, 0.0) for r in small.values()])
    kEs = np.array([at(r, t_end) for r in small.values()])
    k0b, kEb = at(big, 0.0), at(big, t_end)
    sem = kEs.std(ddof=1) / np.sqrt(len(kEs))

    print(f"{'':22s} {'t = 0':>10} {'t_end':>10} {'rise':>9}")
    for s, r in small.items():
        print(f"  L/R=40 seed {s:<8s} {at(r,0.0):10.4f} {at(r,t_end):10.4f} "
              f"{at(r,t_end)/at(r,0.0)-1:+9.1%}")
    print(f"  {'L/R=40 mean':<20s} {k0s.mean():10.4f} {kEs.mean():10.4f} "
          f"{kEs.mean()/k0s.mean()-1:+9.1%}")
    print(f"  {'L/R=40 sd':<20s} {k0s.std(ddof=1):10.4f} {kEs.std(ddof=1):10.4f}")
    print(f"  {'L/R=64 (1 seed)':<20s} {k0b:10.4f} {kEb:10.4f} {kEb/k0b-1:+9.1%}")

    print("\n--- ABSOLUTE k_eff ---")
    for lab, ks, kb in (("t = 0  ", k0s, k0b), ("t_end  ", kEs, kEb)):
        d = kb / ks.mean() - 1
        z = (kb - ks.mean()) / ks.std(ddof=1)
        print(f"  {lab} L/R=64 is {d:+.1%} vs the L/R=40 mean "
              f"({z:+.1f} sd of the seed spread)")
    print(f"  the L/R=40 standard error of the mean at t_end is {sem/kEs.mean():.1%}")
    verdict = "OUTSIDE" if abs(kEb/kEs.mean()-1) > 2*sem/kEs.mean() else "inside"
    print(f"  -> the larger domain lands {verdict} the smaller one's 2-SEM band")

    print("\n--- RELATIVE response k_eff(t)/k_eff(0) ---")
    rs = kEs / k0s
    rb = kEb / k0b
    print(f"  L/R=40: {rs.mean():.4f} +/- {rs.std(ddof=1):.4f} (sd), "
          f"i.e. {rs.mean()-1:+.1%}")
    print(f"  L/R=64: {rb:.4f}, i.e. {rb-1:+.1%}   "
          f"({(rb-rs.mean())/rs.std(ddof=1):+.1f} sd)")
    print("  If this agrees while the absolute does not, the domain-size effect")
    print("  is an OFFSET, and claims about the sintering RESPONSE survive it.")

    # ---- figure -----------------------------------------------------------
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.0, 4.2))
    for (s, r), c in zip(small.items(), fs.C):
        m = r["kt"] <= t_end
        a1.plot(r["kt"][m] / DAY, r["k_iso"][m], "-", color=fs.C[0], lw=1.2,
                alpha=0.55, label="L/R$_{ave}$ = 40 (4 seeds)" if s == "1" else None)
        a2.plot(r["kt"][m] / DAY, r["k_iso"][m] / at(r, 0.0), "-", color=fs.C[0],
                lw=1.2, alpha=0.55,
                label="L/R$_{ave}$ = 40 (4 seeds)" if s == "1" else None)
    mb = big["kt"] <= t_end
    a1.plot(big["kt"][mb] / DAY, big["k_iso"][mb], "-", color=fs.C[1], lw=2.4,
            label="L/R$_{ave}$ = 64 (1 seed)")
    a2.plot(big["kt"][mb] / DAY, big["k_iso"][mb] / k0b, "-", color=fs.C[1],
            lw=2.4, label="L/R$_{ave}$ = 64 (1 seed)")

    fs.style(a1, "time  [days]", r"$k_{\rm eff}$  [W m$^{-1}$K$^{-1}$]",
             f"(a)  absolute: 64 sits {kEb/kEs.mean()-1:+.0%} above", logy=False)
    a1.legend(fontsize=fs.FS_LEG, frameon=False, loc="lower right")
    fs.style(a2, "time  [days]", r"$k_{\rm eff}(t)\,/\,k_{\rm eff}(0)$",
             "(b)  relative: the response is the same", logy=False)
    a2.legend(fontsize=fs.FS_LEG, frameon=False, loc="lower right")
    fig.tight_layout()
    fs.save(fig, HERE, "rev64_compare", dpi=190)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

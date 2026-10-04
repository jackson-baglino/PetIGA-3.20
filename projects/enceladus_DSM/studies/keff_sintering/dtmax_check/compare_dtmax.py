#!/usr/bin/env python3
"""dtmax 2 tau_sub vs 1.09 tau_sub on the same packing (the 3a dtmax check).

    venv_enceladus/bin/python studies/keff_sintering/dtmax_check/compare_dtmax.py \\
        <run at 1.09 tau_sub> <run at 2 tau_sub> [--out <dir>]

The gated seed301 packing at -20 C ran twice with identical options apart
from dtmax (8526 s = 1.09 tau_sub in the first 3a shakedown; 15640 s =
2 tau_sub in the 3a rerun) and the k_eff cadence (every 5 steps vs the SSA
trigger). Both record SSA every step. Compared, for t >= 11 tau_sub (1 d):

  SSA(t)        the 2-tau run interpolated onto the 1.09 run's steps
  k_iso(t)      the 2-tau run interpolated onto the 1.09 run's k_eff samples
  k_iso(SSA)    both on a common SSA grid over the shared range

PASS (stage file batch3a_shakedown.txt): every difference within 0.5%, against
a seed-to-seed scatter of ~3%. Writes dtmax_check.png and prints the numbers.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2] / "postprocess"))
from plot_keff import load, read_tau_sub  # noqa: E402

DAY = 86400.0
PASS_PCT = 0.5


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("ref", type=Path, help="run at dtmax = 1.09 tau_sub")
    ap.add_argument("new", type=Path, help="run at dtmax = 2 tau_sub")
    ap.add_argument("--out", type=Path, default=HERE)
    a = ap.parse_args()

    R, N = load(a.ref), load(a.new)
    tau = read_tau_sub(a.ref)
    t0 = 11.05 * tau
    sR, sN = np.loadtxt(a.ref / "SSA_evo.dat"), np.loadtxt(a.new / "SSA_evo.dat")
    out = {}

    m = sR[:, 2] >= t0
    ssa_new = np.interp(sR[m, 2], sN[:, 2], sN[:, 0])
    d_ssa = (ssa_new / sR[m, 0] - 1) * 100
    out["SSA(t)"] = d_ssa

    m = R["t"] >= t0
    k_new = np.interp(R["t"][m], N["t"], N["kiso"])
    d_kt = (k_new / R["kiso"][m] - 1) * 100
    out["k_iso(t)"] = d_kt

    lo = max(R["ssa"][R["t"] >= t0].min(), N["ssa"][N["t"] >= t0].min())
    hi = min(R["ssa"][R["t"] >= t0].max(), N["ssa"][N["t"] >= t0].max())
    S = np.linspace(lo, hi, 300)
    kR = np.interp(S, R["ssa"][::-1], R["kiso"][::-1])
    kN = np.interp(S, N["ssa"][::-1], N["kiso"][::-1])
    d_ks = (kN / kR - 1) * 100
    out["k_iso(SSA)"] = d_ks

    print(f"steps: 1.09 tau run {len(sR) - 1}, 2 tau run {len(sN) - 1}")
    ok = True
    for k, d in out.items():
        bad = np.max(np.abs(d)) > PASS_PCT
        ok &= not bad
        print(f"  {k:10s}  max |diff| {np.max(np.abs(d)):.3f}%   mean {np.mean(d):+.3f}%   "
              f"{'FAIL' if bad else 'pass'} (limit {PASS_PCT}%)")
    print("RESULT:", "PASS -- 2 tau_sub is indistinguishable on a packing" if ok else
          "FAIL -- return to 1.09 tau_sub before buying more")

    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2), constrained_layout=True)
    ax[0].plot(sR[:, 2] / DAY, sR[:, 0], color="#1a1a1a", lw=2, label="dtmax 1.09 τ_sub")
    ax[0].plot(sN[:, 2] / DAY, sN[:, 0], color="#D55E00", lw=1.4, ls="--", label="dtmax 2 τ_sub")
    ax[0].set(xlabel="time [d]", ylabel="interface measure (SSA_evo.dat col 1)", title="(a) SSA(t), raw")
    ax[1].plot(R["t"] / DAY, R["kiso"], color="#1a1a1a", lw=2, label="1.09 τ_sub")
    ax[1].plot(N["t"] / DAY, N["kiso"], color="#D55E00", lw=1.4, ls="--", label="2 τ_sub")
    ax[1].set(xlabel="time [d]", ylabel="k_iso [W/m/K]", title="(b) k_iso(t)")
    ax[2].plot(S, d_ks, color="#0072B2", lw=2, label="k_iso at matched SSA")
    ax[2].axhspan(-PASS_PCT, PASS_PCT, color="#eeeeee", zorder=0, label=f"±{PASS_PCT}% pass band")
    ax[2].axhline(0, color="#888888", lw=0.8)
    ax[2].set(xlabel="SSA [1/m]", ylabel="2 τ vs 1.09 τ [%]", title="(c) difference, t ≥ 1 d")
    ax[2].invert_xaxis()
    for x in ax:
        x.grid(True, alpha=0.25)
        for sp in ("top", "right"):
            x.spines[sp].set_visible(False)
        x.legend(frameon=False)
    a.out.mkdir(parents=True, exist_ok=True)
    fig.savefig(a.out / "dtmax_check.png", dpi=150)
    print(f"wrote {a.out / 'dtmax_check.png'}")


if __name__ == "__main__":
    main()

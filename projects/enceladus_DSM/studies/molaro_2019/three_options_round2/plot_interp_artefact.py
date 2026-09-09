#!/usr/bin/env python3
"""
plot_interp_artefact.py — the sinusoidal second mode, and where it came from.

The neck-growth curves carried an oscillation riding on the power law. It is
not physics. `.vts` snapshots are written on a coarser grid than the solve
(1080 x 541 against 5394 x 2697), so the sample spacing is dy = 3.5 eps, and
the OLD sub-cell crossing interpolated phi linearly across it. The error of
that interpolation depends on where the interface falls inside a cell, so a
neck growing steadily outward sweeps through sub-cell offsets at a steady rate
and the error appears as a sinusoid. Interpolating in logit(phi) instead is
exact for the model's own equilibrium profile at any spacing.

This reads each snapshot ONCE and measures the neck both ways, so it is
reproducible from the run directories alone.

The analytic version of the same statement, with no simulation involved, is
`studies/molaro_2019/verification/verify_neck_interp.py`.

Usage:  python studies/molaro_2019/three_options_round2/plot_interp_artefact.py
"""

import glob
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent.parent
sys.path.insert(0, str(REPO / "postprocess"))

from neck_width import read_vts_phi, refine_min, _cross          # noqa: E402
from pplib import step_times                                     # noqa: E402

BATCH = Path("/Users/jacksonbaglino/SimulationResults/HPC_results/enceladus_DSM/"
             "batch_2026-09-08__17.20.46_molaro_T-20_round2")
GEOM = "molaro_2D_L450x225um_eps0.12um_axisym_T-20pair_r14um_dom2"
ARMS = [("untuned", "h0.99715_2h_a1e-1_dirichlet", "#B23A48"),
        (r"$D_v\times30$", "h0.99928_2h_a1e-1_Dv30", "#2E86AB"),
        (r"$D_v\times100$", "h0.99923_2h_a1e-1_Dv100", "#1B4965")]


def linear_cross(ya, yb, pa, pb, level):
    """The OLD estimator: linear in phi."""
    if pa == pb:
        return ya
    return ya + (pa - level) / (pa - pb) * (yb - ya)


def chord(col, y, level, cross):
    above = col >= level
    if not above.any():
        return 0.0
    idx = np.flatnonzero(above)
    lo_i, hi_i = idx[0], idx[-1]
    y_lo = y[lo_i] if lo_i == 0 else cross(y[lo_i - 1], y[lo_i],
                                           col[lo_i - 1], col[lo_i], level)
    y_hi = (y[hi_i] if hi_i >= len(y) - 1
            else cross(y[hi_i], y[hi_i + 1], col[hi_i], col[hi_i + 1], level))
    return y_hi - y_lo


def series(run):
    """(t, w_linear, w_logit) in seconds and metres, one pass over the files."""
    files = sorted(glob.glob(str(run / "vtkOut" / "solV_*.vts")))
    tmap = step_times(run / "outp.txt")
    ts, wl, wg = [], [], []
    centers = None
    for fn in files:
        step = int(Path(fn).stem.split("_")[1])
        phi, x, y = read_vts_phi(fn)
        wlin = 2.0 * np.array([chord(phi[:, j], y, 0.5, linear_cross)
                               for j in range(len(x))])
        wlog = 2.0 * np.array([chord(phi[:, j], y, 0.5, _cross)
                               for j in range(len(x))])
        if centers is None:
            thr = 0.5 * wlog.max()
            peaks = [j for j in range(5, len(wlog) - 5)
                     if wlog[j] >= thr and wlog[j] == wlog[j - 5:j + 6].max()]
            centers = (peaks[0], peaks[-1])
        lo, hi = centers
        interior = np.arange(lo + 1, hi)
        out = []
        for w in (wlin, wlog):
            sel = interior[w[interior] > 0]
            jn = sel[np.argmin(w[sel])]
            out.append(refine_min(w, x, jn)[0])
        ts.append(tmap.get(step, np.nan)); wl.append(out[0]); wg.append(out[1])
    o = np.argsort(ts)
    return np.array(ts)[o], np.array(wl)[o], np.array(wg)[o]


def fast(w, k=9):
    """High-frequency part: the curve minus its own 9-point moving mean."""
    return (w - np.convolve(w, np.ones(k) / k, mode="same")) * 1e6


def main():
    k = 9
    fig, axes = plt.subplots(1, 2, figsize=(12.2, 4.9), sharey=True)
    fig.subplots_adjust(wspace=0.08, left=0.075, right=0.985, top=0.80, bottom=0.13)
    stats = {}
    for name, suf, col in ARMS:
        run = BATCH / f"{GEOM}__molaro_T-20_{suf}"
        if not run.is_dir():
            sys.exit(f"missing run directory: {run}")
        t, wl, wg = series(run)
        a, b = fast(wl)[k:-k], fast(wg)[k:-k]
        tt = t[k:-k] / 60.0
        axes[0].plot(tt, a, color=col, lw=1.4, alpha=0.9, label=name)
        axes[1].plot(tt, b, color=col, lw=1.4, alpha=0.9, label=name)
        stats[name] = (a.std(), b.std())
        print(f"{name:16s} fast-residual RMS  linear {a.std():.4f} → "
              f"logit {b.std():.4f} µm   ({100*(b.std()/a.std()-1):+.0f} %)")

    lo = min(v[0] for v in stats.values()); hi = max(v[0] for v in stats.values())
    lo2 = min(v[1] for v in stats.values()); hi2 = max(v[1] for v in stats.values())
    for ax, title, sub in (
            (axes[0], r"A. Before — linear in $\phi$", f"RMS {lo:.3f}–{hi:.3f} µm"),
            (axes[1], r"B. After — linear in $\mathrm{logit}(\phi)$",
             f"RMS {lo2:.3f}–{hi2:.3f} µm")):
        ax.axhline(0, color="0.6", lw=0.8)
        ax.set_ylim(-0.20, 0.20)
        ax.grid(alpha=0.25)
        ax.set_xlabel("time  [min]")
        ax.set_title(f"{title}\n{sub}", fontsize=10.5, loc="left")
        ax.legend(fontsize=8.5, frameon=False, ncol=3, loc="lower center")
    axes[0].set_ylabel("neck width − its own 9-point moving mean  [µm]")

    fig.suptitle("The sinusoid riding on the neck-growth curve is a measurement "
                 "artefact, not physics:\nsub-cell interpolation across a .vts grid "
                 "sampled at 3.5·ε", fontsize=12, y=0.97)
    out = HERE / "neck_interp_effect.png"
    fig.savefig(out, dpi=200)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()

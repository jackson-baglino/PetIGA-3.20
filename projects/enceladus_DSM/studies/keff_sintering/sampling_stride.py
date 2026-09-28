#!/usr/bin/env python3
"""How often does k_eff need sampling? Thin every-step data and look.

    venv_enceladus/bin/python studies/keff_sintering/sampling_stride.py <batch> [--seed 1]
        [--strides 5 10] [--out <dir>]

The warm-end batch sampled k_eff at every accepted step, and the k_eff solves
were 67-78% of its wall time. This keeps every Nth sample (plus the last, as a
thinned run would) and draws them over the full-resolution curve, against time
and against SSA, one column per temperature.

It also measures what thinning loses: k_iso rebuilt from the thinned samples by
linear interpolation (in time, and in SSA), against the every-step value. That
error is printed per stride and in each panel's corner. Only t >= 1 d counts --
the first day is IC relaxation and never enters a result.

Writes one figure per stride, <out>/keff_every<N>_seed<S>.png, lines only and
on identical axes so they can be flipped through; N = 1 is the every-step
reference (default out: <batch>/compare/sampling_test/).
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

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "postprocess"))
from plot_keff import DAY, load          # noqa: E402

# Okabe-Ito, as in plot_keff.py; stride 1 is the reference, drawn in grey ink.
C_STRIDE = ["#0072B2", "#D55E00", "#009E73"]
MARKERS = ["o", "s", "^"]


def thin(n, stride):
    idx = np.arange(0, n, stride)
    return idx if idx[-1] == n - 1 else np.append(idx, n - 1)


def max_err(x, y, idx, live):
    """Max |relative error| of y rebuilt from y[idx] by linear interpolation in x."""
    o = np.argsort(x[idx])
    yr = np.interp(x, x[idx][o], y[idx][o])
    return float(np.max(np.abs(yr[live] / y[live] - 1)) * 100)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("batch", type=Path)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--strides", type=int, nargs="+", default=[5, 10])
    p.add_argument("--out", type=Path, default=None)
    a = p.parse_args(argv)

    runs = {}
    for d in sorted(a.batch.glob(f"*seed{a.seed}_*")):
        m = re.search(r"_T(-?\d+)__", d.name)
        data = load(d) if m else None
        if data is not None:
            runs[int(m.group(1))] = data
    if not runs:
        print(f"no seed-{a.seed} runs with k_eff under {a.batch}")
        return 1
    temps = sorted(runs)

    out = a.out or a.batch / "compare" / "sampling_test"
    out.mkdir(parents=True, exist_ok=True)
    print(f"  seed {a.seed}: max |error| in k_iso after 1 d, linear interpolation "
          "between kept samples")
    # One figure per stride, lines only, identical axes across figures so they
    # can be flipped through. Stride 1 is the every-step reference.
    lims = {}
    for T, d in runs.items():
        lims[T] = ((d["t"][0] / DAY, d["t"][-1] / DAY), (d["ssa"].min(), d["ssa"].max()),
                   (d["kiso"].min(), d["kiso"].max()))
    for s in [1] + [x for x in a.strides if x != 1]:
        fig, axes = plt.subplots(2, len(temps), figsize=(5.2 * len(temps), 9.0),
                                 squeeze=False)
        for j, T in enumerate(temps):
            d = runs[T]
            n = len(d["t"])
            live = d["t"] >= DAY
            idx = thin(n, s)
            tday = d["t"] / DAY
            e_t = max_err(tday, d["kiso"], idx, live)
            e_s = max_err(d["ssa"], d["kiso"], idx, live)
            if s > 1:
                print(f"    T {T:>4}  every {s:>2} ({len(idx):>3} of {n} samples): "
                      f"{e_t:.3f}% vs time, {e_s:.3f}% vs SSA")
            for row, (x, xlabel, err) in enumerate(((tday, "Time [d]", e_t),
                                                    (d["ssa"], r"SSA  [m$^{-1}$]", e_s))):
                ax = axes[row, j]
                ax.plot(x[idx], d["kiso"][idx], "-", color=C_STRIDE[0], lw=1.8)
                if row == 0:
                    ax.axvspan(0, 1, color="#f0f0f0", zorder=0, lw=0)
                    ax.set_title(f"T = {T} °C  ({len(idx)} of {n} samples)", fontsize=14)
                    ax.set_xlim(*lims[T][0])
                else:
                    ax.set_xlim(*lims[T][1])
                ax.set_ylim(*lims[T][2])
                if s > 1:
                    ax.text(0.97, 0.04 if row == 0 else 0.96,
                            f"max err after 1 d: {err:.2f}%", transform=ax.transAxes,
                            ha="right", va="bottom" if row == 0 else "top",
                            fontsize=10, color="#444444")
                ax.set_xlabel(xlabel, fontsize=12)
                if j == 0:
                    ax.set_ylabel(r"$k_\mathrm{iso}$  [W m$^{-1}$ K$^{-1}$]", fontsize=12)
                ax.grid(True, alpha=0.25, lw=0.6)
                for sp in ("top", "right"):
                    ax.spines[sp].set_visible(False)
        what = "every step (reference)" if s == 1 else f"every {s} steps"
        fig.suptitle(f"k_eff sampled {what}, seed {a.seed}", fontsize=15)
        fig.tight_layout()
        f = out / f"keff_every{s}_seed{a.seed}.png"
        fig.savefig(f, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"  wrote {f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

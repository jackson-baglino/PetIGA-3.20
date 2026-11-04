#!/usr/bin/env python3
"""
plot_scalars.py  —  Plot scalar time-series quantities from SSA_evo.dat.

Works for 1D, 2D, and 3D runs. The SSA_evo.dat file is written each monitor
step and has 5 columns (4-column legacy files are also accepted):

    sub_interf/eps   tot_ice   t[s]   step   dt[s]

Derived quantities computed here:
  - Ice volume (= tot_ice in 1D; tot_ice / (Lx*Ly) in 2D)
  - Interface density = sub_interf / eps  (already the first column)
  - Change in ice volume (sublimation proxy)

Usage
-----
  # Single run in current directory
  python plot_scalars.py

  # Specify path
  python plot_scalars.py --file /path/to/SSA_evo.dat --save scalars.png

  # Normalise time to hours
  python plot_scalars.py --time-unit h
"""

import argparse
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from pplib import load_ssa


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Time unit helpers
# ---------------------------------------------------------------------------

TIME_SCALES = {
    "s":   (1.0,    "Time  [s]"),
    "min": (60.0,   "Time  [min]"),
    "h":   (3600.0, "Time  [h]"),
    "d":   (86400., "Time  [days]"),
}


def convert_time(t_sec: np.ndarray, unit: str):
    scale, label = TIME_SCALES.get(unit, (1.0, "Time  [s]"))
    return t_sec / scale, label


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_scalars(data_path: str, time_unit: str = "h", save_path: str = None,
                 title: str = None):
    """
    Plot scalar quantities from SSA_evo.dat in a 5-panel figure:
      1. Ice volume vs time
      2. Change in ice volume vs time
      3. Interface density vs time
      4. dt vs time
      5. dt vs step number
    """
    data = load_ssa(data_path)
    if data is None:
        sys.exit(f"ERROR: no usable SSA data in: {data_path}")

    intf   = data[:, 0]            # sub_interf / eps
    ice    = data[:, 1]            # tot_ice
    t_sec  = data[:, 2]
    steps  = data[:, 3].astype(int)
    dt_col = data[:, 4]            # dt from TSGetTimeStep (may be NaN for legacy)

    t, xlabel = convert_time(t_sec, time_unit)

    d_ice = ice - ice[0]

    # Decide dt source: prefer stored column; fall back to diff
    have_dt = not np.all(np.isnan(dt_col))
    if have_dt:
        dt_vals  = dt_col            # per-step dt (all steps)
        dt_steps = steps
        dt_times = t
    else:
        # Legacy: derive from time differences (one fewer point)
        dt_vals  = np.diff(t_sec)
        dt_steps = steps[1:]
        dt_times = t[1:]

    fig, axes = plt.subplots(5, 1, figsize=(9, 14), sharex=False)

    # ── shared x (time) for panels 0-3 ──────────────────────────────────────
    for ax in axes[:4]:
        ax.set_xlabel(xlabel, fontsize=11)

    # Panel 0: Ice volume
    axes[0].plot(t, ice, color="#1f77b4", lw=2)
    axes[0].set_ylabel(r"$\int \phi_i \, dV$  [m or m²]", fontsize=12)
    axes[0].set_title(
        (title or f"Scalar evolution — {os.path.dirname(os.path.abspath(data_path))}"),
        fontsize=13,
    )

    # Panel 1: Change in ice volume
    axes[1].plot(t, d_ice, color="#d62728", lw=2)
    axes[1].axhline(0, color="gray", lw=0.8, ls="--")
    axes[1].set_ylabel(r"$\Delta\int\phi_i\,dV$  (from $t_0$)", fontsize=12)

    # Panel 2: Interface density
    axes[2].plot(t, intf, color="#2ca02c", lw=2)
    axes[2].set_ylabel(r"Interface density  $\Sigma/\varepsilon$", fontsize=12)

    # Panel 3: dt vs time
    if len(dt_vals) > 0:
        axes[3].semilogy(dt_times, dt_vals, color="#9467bd", lw=1.5,
                         label="from file" if have_dt else "np.diff(t)")
        if not have_dt:
            axes[3].set_title("dt (legacy: computed from Δt — upgrade to 5-col SSA_evo.dat)",
                               fontsize=9, color="gray")
    else:
        axes[3].text(0.5, 0.5, "Single step — no dt info",
                     ha="center", va="center", transform=axes[3].transAxes)
    axes[3].set_ylabel(r"$\Delta t$  [s]", fontsize=12)

    # Panel 4: dt vs step number (independent x-axis)
    if len(dt_vals) > 0:
        axes[4].semilogy(dt_steps, dt_vals, color="#8c564b", lw=1.5, marker=".",
                         markersize=3)
    else:
        axes[4].text(0.5, 0.5, "Single step — no dt info",
                     ha="center", va="center", transform=axes[4].transAxes)
    axes[4].set_xlabel("Time step  #", fontsize=11)
    axes[4].set_ylabel(r"$\Delta t$  [s]", fontsize=12)

    for ax in axes:
        ax.grid(True, alpha=0.3)
        ax.tick_params(labelsize=10)

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Figure saved to: {save_path}")
    else:
        plt.show()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Plot scalar time series from SSA_evo.dat (1D/2D/3D)."
    )
    p.add_argument("--file",      default="SSA_evo.dat",
                   help="Path to SSA_evo.dat (default: ./SSA_evo.dat)")
    p.add_argument("--time-unit", default="h",
                   choices=["s", "min", "h", "d"],
                   help="Time unit for x-axis (default: h)")
    p.add_argument("--save",      default=None,
                   help="Save figure to this path (omit to display)")
    p.add_argument("--title",     default=None,
                   help="Figure title")
    return p.parse_args()


def main():
    args = parse_args()
    plot_scalars(args.file, time_unit=args.time_unit,
                 save_path=args.save, title=args.title)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
compare_runs.py  —  Overlay scalar time-series from multiple simulation runs.

Reads SSA_evo.dat from two or more output directories and plots all runs on
the same axes so you can compare different parameters (temperature, humidity,
mesh resolution, etc.).

Usage
-----
  python compare_runs.py  dir1  dir2  [dir3 ...]

  # With custom labels
  python compare_runs.py  dir1  dir2  --labels "T=-20C" "T=-10C"

  # Save to file
  python compare_runs.py  dir1  dir2  --save comparison.png

  # Change time unit
  python compare_runs.py  dir1  dir2  --time-unit d

  # Normalise quantities by their t=0 value
  python compare_runs.py  dir1  dir2  --normalise
"""

import argparse
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

SSA_COLUMNS = {
    0: ("Interface density  $\\Sigma/\\varepsilon$", "Interface density"),
    1: ("$\\int \\phi_i \\, dV$", "Ice volume"),
}

TIME_SCALES = {
    "s":   (1.0,    "Time  [s]"),
    "min": (60.0,   "Time  [min]"),
    "h":   (3600.0, "Time  [h]"),
    "d":   (86400., "Time  [days]"),
}


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_ssa(run_dir: str) -> np.ndarray:
    path = os.path.join(run_dir, "SSA_evo.dat")
    if not os.path.isfile(path):
        print(f"  WARNING: SSA_evo.dat not found in '{run_dir}' — skipping.")
        return None
    try:
        data = np.genfromtxt(path, dtype=float, comments="#",
                              invalid_raise=False)
        if data.ndim == 1:
            data = data[np.newaxis, :]
        data = data[~np.isnan(data).any(axis=1)]
        return data if len(data) > 0 else None
    except Exception as e:
        print(f"  WARNING: Could not read '{path}': {e}")
        return None


# ---------------------------------------------------------------------------
# Main comparison plot
# ---------------------------------------------------------------------------

def compare(run_dirs: list, labels: list = None, time_unit: str = "h",
            normalise: bool = False, save_path: str = None):
    """
    Overlay ice volume, interface density, and Δice from multiple runs.
    """
    if labels is None or len(labels) != len(run_dirs):
        labels = [os.path.basename(d.rstrip("/")) or d for d in run_dirs]

    scale, xlabel = TIME_SCALES.get(time_unit, (1.0, "Time  [s]"))
    colors = [cm.tab10(i / 10.0) for i in range(len(run_dirs))]

    fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)

    any_data = False
    x_max    = 0.0

    for i, (run_dir, label) in enumerate(zip(run_dirs, labels)):
        data = load_ssa(run_dir)
        if data is None:
            continue

        any_data = True
        t_sec = data[:, 2]
        intf  = data[:, 0]
        ice   = data[:, 1]
        t     = t_sec / scale
        x_max = max(x_max, t[-1])

        c  = colors[i]
        lw = 2.0

        if normalise and ice[0] != 0:
            ice_plot  = ice  / ice[0]
            intf_plot = intf / intf[0] if intf[0] != 0 else intf
            ice_label = "$\\phi_i / \\phi_i(t_0)$"
            intf_label = "$\\Sigma / \\Sigma(t_0)$"
        else:
            ice_plot  = ice
            intf_plot = intf
            ice_label = "$\\int \\phi_i \\, dV$"
            intf_label = "$\\Sigma / \\varepsilon$"

        axes[0].plot(t, ice_plot,         color=c, lw=lw, label=label)
        axes[1].plot(t, intf_plot,        color=c, lw=lw, label=label)
        axes[2].plot(t, ice - ice[0],     color=c, lw=lw, label=label)

    if not any_data:
        sys.exit("No valid SSA_evo.dat files found in any of the specified directories.")

    axes[0].set_ylabel(ice_label,            fontsize=12)
    axes[1].set_ylabel(intf_label,           fontsize=12)
    axes[2].set_ylabel(r"$\Delta\int\phi_i\,dV$ (from $t_0$)", fontsize=12)
    axes[2].set_xlabel(xlabel,               fontsize=12)
    axes[2].axhline(0, color="gray", lw=0.8, ls="--")

    axes[0].set_title("Multi-run comparison — scalar evolution", fontsize=13)

    for ax in axes:
        ax.grid(True, alpha=0.3)
        ax.tick_params(labelsize=11)
        ax.legend(fontsize=9, loc="best")

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Figure saved to: {save_path}")
    else:
        plt.show()


# ---------------------------------------------------------------------------
# 1D profile comparison (overlay final snapshot from each run)
# ---------------------------------------------------------------------------

def compare_1D_profiles(run_dirs: list, labels: list = None,
                        snapshot: str = "last", save_path: str = None,
                        iga_file: str = "igasol.dat"):
    """
    Overlay the ice profile φ_i(x) from each run at the same snapshot.

    snapshot : 'last' (final sol_*.dat), 'first', or an integer step number.
    """
    try:
        from igakit.io import PetIGA
    except ImportError:
        sys.exit("ERROR: igakit is required. Install with: pip install igakit")

    import glob

    if labels is None or len(labels) != len(run_dirs):
        labels = [os.path.basename(d.rstrip("/")) or d for d in run_dirs]

    colors = [cm.tab10(i / 10.0) for i in range(len(run_dirs))]

    fig, axes = plt.subplots(3, 1, figsize=(9, 9), sharex=True)

    any_data = False
    for i, (run_dir, label) in enumerate(zip(run_dirs, labels)):
        iga_path = os.path.join(run_dir, iga_file)
        if not os.path.isfile(iga_path):
            print(f"  WARNING: {iga_path} not found — skipping.")
            continue

        nrb = PetIGA().read(iga_path)
        ctrl = nrb.control
        x    = ctrl[:, 0] * 1e3  # convert m → mm

        # Find snapshot file
        all_sols = sorted(
            glob.glob(os.path.join(run_dir, "sol_*.dat")),
            key=lambda p: int(
                os.path.splitext(os.path.basename(p))[0]
                .lstrip("abcdefghijklmnopqrstuvwxyz_")
            )
        )
        if not all_sols:
            print(f"  WARNING: no sol_*.dat in '{run_dir}' — skipping.")
            continue

        if snapshot == "last":
            sol_file = all_sols[-1]
        elif snapshot == "first":
            sol_file = all_sols[0]
        else:
            target = int(snapshot)
            matches = [s for s in all_sols
                       if int(os.path.splitext(os.path.basename(s))[0]
                               .lstrip("abcdefghijklmnopqrstuvwxyz_")) == target]
            sol_file = matches[0] if matches else all_sols[-1]

        try:
            sol = PetIGA().read_vec(sol_file, nrb)
        except Exception as e:
            print(f"  WARNING: Could not read {sol_file}: {e}")
            continue

        any_data = True
        c = colors[i]
        axes[0].plot(x, sol[:, 0], color=c, lw=2.0, label=label)
        axes[1].plot(x, sol[:, 1], color=c, lw=2.0)
        axes[2].plot(x, sol[:, 2], color=c, lw=2.0)

    if not any_data:
        sys.exit("No valid 1D profiles could be loaded.")

    axes[0].set_ylabel(r"$\phi_i$ (ice)",         fontsize=12)
    axes[1].set_ylabel("Temperature  [°C]",        fontsize=12)
    axes[2].set_ylabel(r"Vapor density  [kg/m³]",  fontsize=12)
    axes[2].set_xlabel("x  [mm]",                  fontsize=12)

    snap_label = f"snapshot: {snapshot}"
    axes[0].set_title(f"1D profile comparison  ({snap_label})", fontsize=13)
    axes[0].legend(fontsize=9)
    axes[0].set_ylim(-0.05, 1.1)

    for ax in axes:
        ax.grid(True, alpha=0.3)
        ax.tick_params(labelsize=11)

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
        description="Compare scalar evolution or 1D profiles across multiple runs."
    )
    p.add_argument("dirs", nargs="+",
                   help="Output directories containing SSA_evo.dat and/or sol_*.dat")
    p.add_argument("--labels",    nargs="+", default=None,
                   help="Legend labels (one per directory)")
    p.add_argument("--time-unit", default="h",
                   choices=["s", "min", "h", "d"])
    p.add_argument("--normalise", action="store_true",
                   help="Normalise quantities by their initial values")
    p.add_argument("--profiles",  action="store_true",
                   help="Overlay 1D field profiles instead of scalars")
    p.add_argument("--snapshot",  default="last",
                   help="Which snapshot to use for profile comparison: 'first', 'last', or step number")
    p.add_argument("--save",      default=None,
                   help="Save figure to this path (omit to display)")
    return p.parse_args()


def main():
    args = parse_args()
    if args.profiles:
        compare_1D_profiles(args.dirs, labels=args.labels,
                            snapshot=args.snapshot, save_path=args.save)
    else:
        compare(args.dirs, labels=args.labels, time_unit=args.time_unit,
                normalise=args.normalise, save_path=args.save)


if __name__ == "__main__":
    main()

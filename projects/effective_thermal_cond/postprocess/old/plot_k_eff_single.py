#!/usr/bin/env python3
"""
Plot effective thermal conductivity (k_xx and k_yy) vs time for a single simulation.

Reads from one directory (cwd by default):
  - SSA_evo.dat          # contains time and step columns
  - k_eff.csv            # contains k_00 (xx), k_11 (yy), and sol_index step column

Usage:
  python plot_k_eff_single.py               # use current directory
  python plot_k_eff_single.py --dir /path/to/sim
  python plot_k_eff_single.py --no-save     # show interactively instead of saving
  python plot_k_eff_single.py --dpi 300
  python plot_k_eff_single.py --normalize   # normalize k_xx and k_yy by their respective max values
  python plot_k_eff_single.py --dir /path/to/simA --dir2 /path/to/simB  # overlay both simulations
"""

import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ---- Publication-quality aesthetics ----
plt.rcParams.update({
    "figure.dpi": 100,
    "savefig.dpi": 600,
    "font.size": 9,                # base size; suitable for single-column figures
    "font.family": "serif",
    "mathtext.fontset": "stix",
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "axes.linewidth": 0.8,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.minor.width": 0.6,
    "ytick.minor.width": 0.6,
    "legend.fontsize": 8,
    "legend.frameon": False,
})

from matplotlib.ticker import AutoMinorLocator, ScalarFormatter

def beautify_axes(ax):
    ax.grid(True, which="major", linestyle="--", alpha=0.35)
    ax.grid(True, which="minor", linestyle=":", alpha=0.20)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    # Scientific formatting without offset box
    for axis in (ax.xaxis, ax.yaxis):
        fmt = ScalarFormatter(useMathText=True)
        fmt.set_powerlimits((-3, 4))
        fmt.set_scientific(True)
        axis.set_major_formatter(fmt)
    return ax

def find_file_or_die(path, name):
    if not os.path.exists(path):
        raise FileNotFoundError(f"Required file not found: {path}\n(looking for {name})")
    return path

def get_col(df, candidates, required=True):
    """Return the first matching column from candidate names (case-sensitive)."""
    for c in candidates:
        if c in df.columns:
            return df[c].values
    if required:
        raise KeyError(f"None of the expected columns found: {candidates}\nAvailable: {list(df.columns)}")
    return None

def load_match_data(sim_dir):
    # --- Load SSA_evo.dat ---
    ssa_path = find_file_or_die(os.path.join(sim_dir, "SSA_evo.dat"), "SSA_evo.dat")
    # Expecting columns like: [ssa_value, ?, time, step, ...]
    ssa = np.loadtxt(ssa_path)
    if ssa.ndim == 1:
        ssa = ssa.reshape(1, -1)
    if ssa.shape[1] < 4:
        raise ValueError(f"SSA_evo.dat has {ssa.shape[1]} columns; expected at least 4.")
    ssa_val  = ssa[:, 0]
    ssa_time = ssa[:, 2]
    ssa_step = ssa[:, 3].astype(int)

    # --- Load k_eff.csv ---
    keff_path = find_file_or_die(os.path.join(sim_dir, "k_eff.csv"), "k_eff.csv")
    df = pd.read_csv(keff_path)

    # Some codebases use different names; try several:
    steps = get_col(df, ["sol_index", "step", "time_step"])
    k_xx  = get_col(df, ["k_00", "k_xx", "kxx", "k11"])  # include common alternates just in case
    k_yy  = get_col(df, ["k_11", "k_yy", "kyy", "k22"])

    # --- Build mapping from step -> time using SSA ---
    step_to_time = {int(s): float(t) for s, t in zip(ssa_step, ssa_time)}
    step_to_ssa  = {int(s): float(v) for s, v in zip(ssa_step, ssa_val)}

    # --- Align times for k_eff rows using step ---
    times = []
    kxx_m = []
    kyy_m = []
    ssa_m = []
    missing = 0
    for i, st in enumerate(steps.astype(int)):
        t = step_to_time.get(st)
        if t is None:
            missing += 1
            continue
        times.append(t)
        kxx_m.append(k_xx[i])
        kyy_m.append(k_yy[i])
        ssa_m.append(step_to_ssa[st])

    if not times:
        raise RuntimeError("No overlapping steps between SSA_evo.dat and k_eff.csv.")

    times = np.array(times)
    kxx_m = np.array(kxx_m)
    kyy_m = np.array(kyy_m)

    # Sort by time in case CSV order differs
    order = np.argsort(times)
    times, kxx_m, kyy_m = times[order], kxx_m[order], kyy_m[order]
    ssa_m = np.array(ssa_m)[order]

    return times, kxx_m, kyy_m, ssa_m, missing

def autoscale_time_axis(times):
    """Return scaled_times and an appropriate label with units."""
    if len(times) == 0:
        return times, "Time (s)"
    tmax = float(np.nanmax(times))
    if tmax >= 2*24*3600:
        return times/86400.0, "Time (days)"
    elif tmax >= 2*3600:
        return times/3600.0, "Time (hours)"
    elif tmax >= 2*60:
        return times/60.0, "Time (minutes)"
    else:
        return times, "Time (s)"

def autoscale_time_axis_two(times_a, times_b):
    """Return (times_a_scaled, times_b_scaled, label) using a common unit."""
    tmax_a = 0.0 if len(times_a) == 0 else float(np.nanmax(times_a))
    tmax_b = 0.0 if len(times_b) == 0 else float(np.nanmax(times_b))
    tmax = max(tmax_a, tmax_b)
    if tmax >= 2*24*3600:
        return times_a/86400.0, times_b/86400.0, "Time (days)"
    elif tmax >= 2*3600:
        return times_a/3600.0, times_b/3600.0, "Time (hours)"
    elif tmax >= 2*60:
        return times_a/60.0, times_b/60.0, "Time (minutes)"
    else:
        return times_a, times_b, "Time (s)"

def main():
    ap = argparse.ArgumentParser(description="Plot k_xx and k_yy vs time for a single simulation directory.")
    ap.add_argument("--dir", default=".", help="Simulation directory containing SSA_evo.dat and k_eff.csv")
    ap.add_argument("--dir2", default=None, help="Optional second simulation directory to overlay")
    ap.add_argument("--save", dest="save", action="store_true", default=True, help="Save figures instead of showing (default)")
    ap.add_argument("--no-save", dest="save", action="store_false", help="Show interactively instead of saving")
    ap.add_argument("--dpi", type=int, default=300, help="Figure DPI when saving")
    ap.add_argument("--normalize", action="store_true", help="Normalize k_xx and k_yy by their respective maximum values")
    ap.add_argument("--plot-ssa", action="store_true", help="Also plot k_xx and k_yy against SSA (from SSA_evo.dat)")
    args = ap.parse_args()

    sim_dir = os.path.abspath(args.dir)
    title = os.path.basename(sim_dir.rstrip("/"))

    times, kxx, kyy, ssa_vals, missing = load_match_data(sim_dir)

    # Optional second dataset
    sim_dir2 = None
    title2 = None
    times2 = kxx2 = kyy2 = ssa_vals2 = None
    missing2 = 0
    if args.dir2:
        sim_dir2 = os.path.abspath(args.dir2)
        title2 = os.path.basename(sim_dir2.rstrip("/"))
        times2, kxx2, kyy2, ssa_vals2, missing2 = load_match_data(sim_dir2)

    # Optional normalization (per-dataset)
    if args.normalize:
        if np.nanmax(kxx) != 0: kxx = kxx / np.nanmax(kxx)
        if np.nanmax(kyy) != 0: kyy = kyy / np.nanmax(kyy)
        if sim_dir2 is not None:
            if np.nanmax(kxx2) != 0: kxx2 = kxx2 / np.nanmax(kxx2)
            if np.nanmax(kyy2) != 0: kyy2 = kyy2 / np.nanmax(kyy2)

    # Choose a common time scaling if overlaying, else scale single set
    if sim_dir2 is not None:
        times_scaled, times2_scaled, time_label = autoscale_time_axis_two(times, times2)
    else:
        times_scaled, time_label = autoscale_time_axis(times)

    if missing:
        print(f"[WARN] {missing} k_eff rows had steps not present in SSA_evo.dat for {title} and were skipped.")
    if sim_dir2 is not None and missing2:
        print(f"[WARN] {missing2} k_eff rows had steps not present in SSA_evo.dat for {title2} and were skipped.")

    # --- Plot 1: both on the same axes (single-column friendly) ---
    fig1, ax1 = plt.subplots(figsize=(3.5, 2.4))  # inches
    ax1.plot(times_scaled, kxx, marker='o', ms=3.0, lw=1.0, label=r"$k_{xx}$")
    ax1.plot(times_scaled, kyy, marker='s', ms=3.0, lw=1.0, label=r"$k_{yy}$")
    if sim_dir2 is not None:
        ax1.plot(times2_scaled, kxx2, marker='o', ms=3.0, lw=1.0, linestyle='--', label=rf"{title2} $k_{{xx}}$")
        ax1.plot(times2_scaled, kyy2, marker='s', ms=3.0, lw=1.0, linestyle='--', label=rf"{title2} $k_{{yy}}$")
    ax1.set_xlabel(time_label)
    ax1.set_ylabel(r"$k$ (W m$^{-1}$ K$^{-1}$)")
    ax1.legend(loc="best", ncols=1 if sim_dir2 is not None else 2, handlelength=1.4, columnspacing=0.8)
    beautify_axes(ax1)
    fig1.tight_layout()

    # --- Plot 2: stacked panels (clear comparison) ---
    fig2, (ax2, ax3) = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(3.5, 3.6))
    ax2.plot(times_scaled, kxx, marker='o', ms=3.0, lw=1.0)
    ax2.set_ylabel(r"$k_{xx}$ (W m$^{-1}$ K$^{-1}$)")
    beautify_axes(ax2)
    if sim_dir2 is not None:
        ax2.plot(times2_scaled, kxx2, marker='o', ms=3.0, lw=1.0, linestyle='--')

    ax3.plot(times_scaled, kyy, marker='s', ms=3.0, lw=1.0)
    ax3.set_xlabel(time_label)
    ax3.set_ylabel(r"$k_{yy}$ (W m$^{-1}$ K$^{-1}$)")
    beautify_axes(ax3)
    if sim_dir2 is not None:
        ax3.plot(times2_scaled, kyy2, marker='s', ms=3.0, lw=1.0, linestyle='--')

    fig2.suptitle("")
    fig2.tight_layout()

    if args.plot_ssa:
        fig_ssa, ax_ssa = plt.subplots(figsize=(3.5, 2.4))
        ax_ssa.plot(ssa_vals, kxx, marker='o', ms=3.0, lw=1.0, label=r"$k_{xx}$")
        ax_ssa.plot(ssa_vals, kyy, marker='s', ms=3.0, lw=1.0, label=r"$k_{yy}$")
        if sim_dir2 is not None:
            ax_ssa.plot(ssa_vals2, kxx2, marker='o', ms=3.0, lw=1.0, linestyle='--', label=rf"{title2} $k_{{xx}}$")
            ax_ssa.plot(ssa_vals2, kyy2, marker='s', ms=3.0, lw=1.0, linestyle='--', label=rf"{title2} $k_{{yy}}$")
        ax_ssa.set_xlabel("SSA")
        ax_ssa.set_ylabel(r"$k$ (W m$^{-1}$ K$^{-1}$)")
        ax_ssa.legend(loc="best", ncols=1 if sim_dir2 is not None else 2, handlelength=1.4, columnspacing=0.8)
        beautify_axes(ax_ssa)
        fig_ssa.tight_layout()

    if args.save:
        suffix = "_overlay" if sim_dir2 is not None else ""
        out1_png = os.path.join(sim_dir, f"k_eff_vs_time_xy{suffix}.png")
        out2_png = os.path.join(sim_dir, f"k_eff_vs_time_stacked{suffix}.png")
        out1_pdf = os.path.join(sim_dir, f"k_eff_vs_time_xy{suffix}.pdf")
        out2_pdf = os.path.join(sim_dir, f"k_eff_vs_time_stacked{suffix}.pdf")
        out_ssa_png = os.path.join(sim_dir, f"k_eff_vs_ssa_xy{suffix}.png")
        out_ssa_pdf = os.path.join(sim_dir, f"k_eff_vs_ssa_xy{suffix}.pdf")
        for path in (out1_png, out2_png):
            plt.gcf().set_facecolor('white')
        fig1.savefig(out1_png, dpi=args.dpi, bbox_inches="tight", transparent=True)
        fig2.savefig(out2_png, dpi=args.dpi, bbox_inches="tight", transparent=True)
        # Vector versions for publication
        fig1.savefig(out1_pdf, bbox_inches="tight", transparent=True)
        fig2.savefig(out2_pdf, bbox_inches="tight", transparent=True)
        if args.plot_ssa:
            fig_ssa.savefig(out_ssa_png, dpi=args.dpi, bbox_inches="tight", transparent=True)
            fig_ssa.savefig(out_ssa_pdf, bbox_inches="tight", transparent=True)
        saved_list = [out1_png, out2_png, out1_pdf, out2_pdf]
        if args.plot_ssa:
            saved_list.extend([out_ssa_png, out_ssa_pdf])
        print("✅ Saved:\n  " + "\n  ".join(saved_list))
        plt.close(fig1); plt.close(fig2)
        if args.plot_ssa:
            plt.close(fig_ssa)
    else:
        plt.show()

if __name__ == "__main__":
    main()
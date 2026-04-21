#!/usr/bin/env python3
"""
plot_relax.py — Interface relaxation analysis

Reads relax_monitor.dat (written at every time step) and produces a two-panel
figure that shows:
  - Top:    Each interface metric normalized by its step-0 value vs. step number
  - Bottom: Step-to-step relative change |Δq/q₀| vs. step number (log scale)

A dashed threshold line marks the convergence criterion (default 0.1 %).
The step at which ALL metrics cross below the threshold is printed to stdout
and annotated on the plot — use that as -nsteps_sed.

Usage:
    python plot_relax.py [--file relax_monitor.dat] [--threshold 1e-3]
                         [--out relax_analysis.png] [--show]

    --file      : path to relax_monitor.dat  (default: relax_monitor.dat)
    --threshold : relative-change threshold for "relaxed" (default: 1e-3)
    --out       : output PNG filename        (default: relax_analysis.png)
    --show      : display the plot interactively in addition to saving
"""

import argparse
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--file",      default="relax_monitor.dat",
                   help="Path to relax_monitor.dat  (default: relax_monitor.dat)")
    p.add_argument("--threshold", type=float, default=1e-3,
                   help="Relative-change threshold for 'relaxed'  (default: 1e-3)")
    p.add_argument("--out",       default="relax_analysis.png",
                   help="Output PNG filename  (default: relax_analysis.png)")
    p.add_argument("--show",      action="store_true",
                   help="Display the figure interactively after saving")
    return p.parse_args()


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_relax_data(path):
    """
    Parse relax_monitor.dat.
    Expected columns (space-separated, comment lines start with #):
        step  t  tot_ice  tot_sed  ice_air_interf  sed_air_interf  ice_sed_interf  tot_trip
    Returns a dict of 1-D numpy arrays keyed by column name.
    """
    if not os.path.isfile(path):
        sys.exit(f"ERROR: file not found: {path}\n"
                 "Run a simulation first, then call this script from the output folder.")

    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data[np.newaxis, :]  # single-row edge case

    names = ["step", "t",
             "tot_ice", "tot_sed",
             "ice_air_interf", "sed_air_interf", "ice_sed_interf",
             "tot_trip"]

    if data.shape[1] < len(names):
        sys.exit(f"ERROR: expected {len(names)} columns, found {data.shape[1]}.\n"
                 "Recompile with the updated assembly.c / monitoring.c and rerun.")

    return {name: data[:, i] for i, name in enumerate(names)}


# ---------------------------------------------------------------------------
# Analysis
# ---------------------------------------------------------------------------

def compute_metrics(d, threshold):
    """
    Returns:
        norm_vals  : dict  metric -> normalized value array  (val / val[0])
        rel_change : dict  metric -> per-step relative change |Δq / q₀|
        converged  : dict  metric -> first step below threshold  (or None)
        relax_step : int   step at which ALL metrics are converged  (or None)
    """
    interface_keys = ["ice_air_interf", "sed_air_interf", "ice_sed_interf", "tot_trip"]

    norm_vals  = {}
    rel_change = {}
    converged  = {}

    for key in interface_keys:
        v   = d[key]
        v0  = v[0] if v[0] != 0.0 else 1.0          # guard against divide-by-zero
        norm_vals[key] = v / v0

        # |Δq| / |q₀|  — first entry is NaN (no previous step)
        delta = np.abs(np.diff(v)) / abs(v0)
        rel_change[key] = np.concatenate([[np.nan], delta])

        # First step where change drops (and stays) below threshold
        below = np.where(delta < threshold)[0]
        converged[key] = int(d["step"][below[0] + 1]) if len(below) > 0 else None

    # Global convergence: the last metric to converge
    conv_steps = [s for s in converged.values() if s is not None]
    relax_step = max(conv_steps) if len(conv_steps) == len(interface_keys) else None

    return norm_vals, rel_change, converged, relax_step


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

COLORS = {
    "ice_air_interf": "#4C9BE8",   # blue
    "sed_air_interf": "#E8834C",   # orange
    "ice_sed_interf": "#4CE87A",   # green
    "tot_trip":       "#C44CE8",   # purple
}

LABELS = {
    "ice_air_interf": r"$\int \phi_i^2\,\phi_a^2\,dV$  (ice–air)",
    "sed_air_interf": r"$\int \phi_s^2\,\phi_a^2\,dV$  (sed–air)",
    "ice_sed_interf": r"$\int \phi_s^2\,\phi_i^2\,dV$  (ice–sed)",
    "tot_trip":       r"$\int \phi_i^2\,\phi_s^2\,\phi_a^2\,dV$  (triple junc.)",
}


def make_figure(d, norm_vals, rel_change, converged, relax_step, threshold, out_path):
    steps = d["step"]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 8),
                                   sharex=True,
                                   gridspec_kw={"hspace": 0.08})

    # --- Top panel: normalized values ---
    for key, vals in norm_vals.items():
        ax1.plot(steps, vals, color=COLORS[key], lw=1.8, label=LABELS[key])

    ax1.axhline(1.0, color="gray", lw=0.8, ls="--")
    ax1.set_ylabel("Interface metric  (normalized to step 0)", fontsize=12)
    ax1.legend(fontsize=9, loc="upper right")
    ax1.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax1.grid(True, which="major", ls=":", alpha=0.4)
    ax1.set_title("Phase-field interface relaxation", fontsize=14, pad=8)

    # --- Bottom panel: step-to-step relative change ---
    for key, vals in rel_change.items():
        ax2.semilogy(steps, vals, color=COLORS[key], lw=1.8, label=LABELS[key])

    ax2.axhline(threshold, color="red", lw=1.2, ls="--",
                label=f"Threshold  {threshold:.0e}")

    if relax_step is not None:
        ax2.axvline(relax_step, color="black", lw=1.5, ls="-.",
                    label=f"Converged at step {relax_step}")
        ax1.axvline(relax_step, color="black", lw=1.5, ls="-.")

    ax2.set_xlabel("Time step", fontsize=12)
    ax2.set_ylabel(r"Step-to-step relative change  $|\Delta q\,/\,q_0|$", fontsize=12)
    ax2.legend(fontsize=9, loc="upper right")
    ax2.grid(True, which="major", ls=":", alpha=0.4)
    ax2.grid(True, which="minor", ls=":", alpha=0.2)
    ax2.set_xlim(left=steps[0])

    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Saved: {out_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    d = load_relax_data(args.file)
    print(f"Loaded {len(d['step'])} steps from '{args.file}'")

    norm_vals, rel_change, converged, relax_step = compute_metrics(d, args.threshold)

    print(f"\nInterface convergence (threshold = {args.threshold:.1e}):")
    for key, step in converged.items():
        label = LABELS[key].split("(")[1].rstrip(")")
        print(f"  {label:20s}: step {step}" if step is not None else
              f"  {label:20s}: never converged within {int(d['step'][-1])} steps")

    if relax_step is not None:
        print(f"\n>>> Recommended -nsteps_sed = {relax_step}  "
              f"(all interfaces converged by this step)")
    else:
        print("\n>>> Not all interfaces converged — run more steps before choosing nsteps_sed")

    make_figure(d, norm_vals, rel_change, converged, relax_step, args.threshold, args.out)

    if args.show:
        matplotlib.use("TkAgg")
        plt.show()


if __name__ == "__main__":
    main()

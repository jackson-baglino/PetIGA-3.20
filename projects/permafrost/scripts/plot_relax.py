#!/usr/bin/env python3
"""
plot_relax.py — Interface relaxation analysis

Reads relax_monitor.dat (written at every time step) and produces figures that show:
  - Normalized interface metrics vs. step number
  - Step-to-step relative change |Δq/q₀| vs. step number (log scale)

A dashed threshold line marks the convergence criterion (default 0.1 %).
The step at which ALL metrics cross below the threshold is printed to stdout
and annotated on the plot — use that as -nsteps_sed.

Usage:
    python plot_relax.py [--file relax_monitor.dat] [--threshold 1e-3]
                         [--out relax_analysis.png] [--show]
                         [--slideshow] [--outdir DIR] [--fmt {png,pdf}]

    --file      : path to relax_monitor.dat  (default: relax_monitor.dat)
    --threshold : relative-change threshold for "relaxed" (default: 1e-3)
    --out       : output filename for the combined figure (default: relax_analysis.png)
    --show      : display the plot interactively in addition to saving
    --slideshow : widescreen 16:9, large fonts, high DPI; saves 5 separate figures
    --outdir    : directory to save all output files (default: same dir as --out)
    --fmt       : output format, png or pdf (default: png)
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
                   help="Output filename for combined figure  (default: relax_analysis.png)")
    p.add_argument("--show",      action="store_true",
                   help="Display the figure interactively after saving")
    p.add_argument("--slideshow", action="store_true",
                   help="Widescreen 16:9 figures with large fonts for presentations")
    p.add_argument("--outdir",    default=None,
                   help="Directory to save all output files (default: same dir as --out)")
    p.add_argument("--fmt",       default="png", choices=["png", "pdf"],
                   help="Output file format  (default: png)")
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
        data = data[np.newaxis, :]

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
        v0  = v[0] if v[0] != 0.0 else 1.0
        norm_vals[key] = v / v0

        delta = np.abs(np.diff(v)) / abs(v0)
        rel_change[key] = np.concatenate([[np.nan], delta])

        below = np.where(delta < threshold)[0]
        converged[key] = int(d["step"][below[0] + 1]) if len(below) > 0 else None

    conv_steps = [s for s in converged.values() if s is not None]
    relax_step = max(conv_steps) if len(conv_steps) == len(interface_keys) else None

    return norm_vals, rel_change, converged, relax_step


# ---------------------------------------------------------------------------
# Style helpers
# ---------------------------------------------------------------------------

COLORS = {
    "ice_air_interf": "#4C9BE8",
    "sed_air_interf": "#E8834C",
    "ice_sed_interf": "#4CE87A",
    "tot_trip":       "#C44CE8",
}

LABELS = {
    "ice_air_interf": r"$\int \phi_i^2\,\phi_a^2\,dV$  (ice–air)",
    "sed_air_interf": r"$\int \phi_s^2\,\phi_a^2\,dV$  (sed–air)",
    "ice_sed_interf": r"$\int \phi_s^2\,\phi_i^2\,dV$  (ice–sed)",
    "tot_trip":       r"$\int \phi_i^2\,\phi_s^2\,\phi_a^2\,dV$  (triple junc.)",
}


def _slide_params():
    return dict(figsize=(13.3, 7.5), dpi=250, lw=2.5,
                label_fs=18, title_fs=20, legend_fs=13, tick_fs=14)


def _default_params():
    return dict(figsize=(9, 8), dpi=150, lw=1.8,
                label_fs=12, title_fs=14, legend_fs=9, tick_fs=10)


def _apply_style(ax, tick_fs):
    ax.tick_params(axis='both', labelsize=tick_fs)
    ax.grid(True, which="major", ls=":", alpha=0.4)


def _add_convergence_annotation(ax, relax_step, slideshow):
    """Add a text box noting the recommended nsteps_sed."""
    if relax_step is None:
        return
    fs = 14 if slideshow else 10
    ax.annotate(
        f"Recommended\nnsteps_sed = {relax_step}",
        xy=(relax_step, ax.get_ylim()[1]),
        xycoords="data",
        xytext=(0.72, 0.95), textcoords="axes fraction",
        fontsize=fs,
        ha="center", va="top",
        bbox=dict(boxstyle="round,pad=0.4", fc="lightyellow", ec="gray", alpha=0.9),
        arrowprops=dict(arrowstyle="->", color="black", lw=1.2),
    )


def _add_time_secondary_axis(ax, steps, times, label_fs):
    """Add a secondary x-axis showing physical time [s] above the top spine."""
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    # Pick ~5 evenly-spaced tick positions (by step index)
    n = len(steps)
    idx = np.linspace(0, n - 1, min(6, n), dtype=int)
    ax2.set_xticks(steps[idx])
    ax2.set_xticklabels([f"{times[i]:.2g}" for i in idx], fontsize=label_fs - 4)
    ax2.set_xlabel("Physical time  [s]", fontsize=label_fs - 2)
    return ax2


# ---------------------------------------------------------------------------
# Figure builders
# ---------------------------------------------------------------------------

def _build_combined(d, norm_vals, rel_change, relax_step, threshold,
                    slideshow, params):
    """Two-panel combined figure (normalized + rate-of-change)."""
    steps = d["step"]
    lw    = params["lw"]
    lfs   = params["label_fs"]
    tfs   = params["title_fs"]
    lefs  = params["legend_fs"]
    tkfs  = params["tick_fs"]

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=params["figsize"],
        sharex=True,
        gridspec_kw={"hspace": 0.08},
        facecolor="white",
        layout="constrained",
    )

    for key, vals in norm_vals.items():
        ax1.plot(steps, vals, color=COLORS[key], lw=lw, label=LABELS[key])
    ax1.axhline(1.0, color="gray", lw=0.8, ls="--")
    if relax_step is not None:
        ax1.axvline(relax_step, color="black", lw=1.5, ls="-.")
    ax1.set_ylabel("Interface metric  (norm. to step 0)", fontsize=lfs)
    ax1.legend(fontsize=lefs, loc="upper right")
    ax1.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax1.set_title("Phase-field interface relaxation — test3 EnclosedGrainPair",
                  fontsize=tfs, pad=8)
    _apply_style(ax1, tkfs)

    for key, vals in rel_change.items():
        ax2.semilogy(steps, vals, color=COLORS[key], lw=lw, label=LABELS[key])
    ax2.axhline(threshold, color="red", lw=1.2, ls="--",
                label=f"Threshold  {threshold:.0e}")
    if relax_step is not None:
        ax2.axvline(relax_step, color="black", lw=1.5, ls="-.",
                    label=f"Converged at step {relax_step}")
    ax2.set_xlabel("Time step", fontsize=lfs)
    ax2.set_ylabel(r"Step-to-step rel. change  $|\Delta q/q_0|$", fontsize=lfs)
    ax2.legend(fontsize=lefs, loc="upper right")
    ax2.grid(True, which="major", ls=":", alpha=0.4)
    ax2.grid(True, which="minor", ls=":", alpha=0.2)
    ax2.set_xlim(left=steps[0])
    _apply_style(ax2, tkfs)

    if slideshow:
        _add_time_secondary_axis(ax1, steps, d["t"], lfs)

    return fig


def _build_normalized(d, norm_vals, relax_step, slideshow, params, step_xlim=None):
    """Single-panel: normalized metrics vs. step."""
    steps = d["step"]
    lw    = params["lw"]
    lfs   = params["label_fs"]
    tfs   = params["title_fs"]
    lefs  = params["legend_fs"]
    tkfs  = params["tick_fs"]

    fig, ax = plt.subplots(figsize=params["figsize"], facecolor="white",
                           layout="constrained")

    for key, vals in norm_vals.items():
        ax.plot(steps, vals, color=COLORS[key], lw=lw, label=LABELS[key])
    ax.axhline(1.0, color="gray", lw=0.8, ls="--")
    if relax_step is not None:
        ax.axvline(relax_step, color="black", lw=1.5, ls="-.",
                   label=f"nsteps_sed = {relax_step}")
        _add_convergence_annotation(ax, relax_step, slideshow)

    if step_xlim is not None:
        ax.set_xlim(*step_xlim)

    ax.set_xlabel("Time step", fontsize=lfs)
    ax.set_ylabel("Interface metric  (norm. to step 0)", fontsize=lfs)
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.legend(fontsize=lefs, loc="upper right")
    ax.set_title("Phase-field interface relaxation — normalized metrics",
                 fontsize=tfs, pad=8)
    _apply_style(ax, tkfs)

    if slideshow:
        _add_time_secondary_axis(ax, steps, d["t"], lfs)

    return fig


def _build_rate(d, rel_change, relax_step, threshold, slideshow, params):
    """Single-panel: step-to-step relative change (log scale)."""
    steps = d["step"]
    lw    = params["lw"]
    lfs   = params["label_fs"]
    tfs   = params["title_fs"]
    lefs  = params["legend_fs"]
    tkfs  = params["tick_fs"]

    fig, ax = plt.subplots(figsize=params["figsize"], facecolor="white",
                           layout="constrained")

    for key, vals in rel_change.items():
        ax.semilogy(steps, vals, color=COLORS[key], lw=lw, label=LABELS[key])
    ax.axhline(threshold, color="red", lw=1.2, ls="--",
               label=f"Threshold  {threshold:.0e}")
    if relax_step is not None:
        ax.axvline(relax_step, color="black", lw=1.5, ls="-.",
                   label=f"Converged at step {relax_step}")

    ax.set_xlabel("Time step", fontsize=lfs)
    ax.set_ylabel(r"Step-to-step rel. change  $|\Delta q / q_0|$", fontsize=lfs)
    ax.legend(fontsize=lefs, loc="upper right")
    ax.set_title("Interface relaxation — rate of change", fontsize=tfs, pad=8)
    ax.grid(True, which="major", ls=":", alpha=0.4)
    ax.grid(True, which="minor", ls=":", alpha=0.2)
    ax.set_xlim(left=steps[0])
    _apply_style(ax, tkfs)

    if slideshow:
        _add_time_secondary_axis(ax, steps, d["t"], lfs)

    return fig


def _build_time_axis(d, norm_vals, relax_step, slideshow, params):
    """Single-panel: normalized metrics with physical time [s] on the x-axis."""
    times = d["t"]
    lw    = params["lw"]
    lfs   = params["label_fs"]
    tfs   = params["title_fs"]
    lefs  = params["legend_fs"]
    tkfs  = params["tick_fs"]

    fig, ax = plt.subplots(figsize=params["figsize"], facecolor="white",
                           layout="constrained")

    for key, vals in norm_vals.items():
        ax.plot(times, vals, color=COLORS[key], lw=lw, label=LABELS[key])
    ax.axhline(1.0, color="gray", lw=0.8, ls="--")

    if relax_step is not None:
        step_arr = d["step"]
        idx = np.searchsorted(step_arr, relax_step)
        if idx < len(times):
            t_relax = times[idx]
            ax.axvline(t_relax, color="black", lw=1.5, ls="-.",
                       label=f"nsteps_sed = {relax_step}  (t = {t_relax:.2g} s)")

    ax.set_xlabel("Physical time  [s]", fontsize=lfs)
    ax.set_ylabel("Interface metric  (norm. to step 0)", fontsize=lfs)
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.legend(fontsize=lefs, loc="upper right")
    ax.set_title("Phase-field interface relaxation — physical time axis",
                 fontsize=tfs, pad=8)
    _apply_style(ax, tkfs)

    return fig


# ---------------------------------------------------------------------------
# Save helper
# ---------------------------------------------------------------------------

def save_fig(fig, outdir, basename, fmt):
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, f"{basename}.{fmt}")
    fig.savefig(path, dpi=fig.get_dpi(), bbox_inches="tight")
    print(f"Saved: {path}")
    plt.close(fig)


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
        if step is not None:
            print(f"  {label:20s}: step {step}")
        else:
            print(f"  {label:20s}: never converged within {int(d['step'][-1])} steps")

    if relax_step is not None:
        print(f"\n>>> Recommended -nsteps_sed = {relax_step}  "
              f"(all interfaces converged by this step)")
    else:
        print("\n>>> Not all interfaces converged — run more steps before choosing nsteps_sed")

    params = _slide_params() if args.slideshow else _default_params()

    # Resolve output directory
    if args.outdir is not None:
        outdir = args.outdir
    else:
        outdir = os.path.dirname(os.path.abspath(args.out)) if args.out else "."

    fmt = args.fmt

    if args.slideshow:
        # --- 1. Combined 2-panel ---
        fig = _build_combined(d, norm_vals, rel_change, relax_step,
                               args.threshold, True, params)
        save_fig(fig, outdir, "relax_combined", fmt)

        # --- 2. Normalized metrics (full run) ---
        fig = _build_normalized(d, norm_vals, relax_step, True, params)
        save_fig(fig, outdir, "relax_normalized", fmt)

        # --- 3. Rate of change ---
        fig = _build_rate(d, rel_change, relax_step, args.threshold, True, params)
        save_fig(fig, outdir, "relax_rate", fmt)

        # --- 4. Early zoom (first 10 % of steps) ---
        n_steps = len(d["step"])
        early_end = max(int(n_steps * 0.10), 5)
        xlim_early = (d["step"][0], d["step"][min(early_end, n_steps - 1)])
        fig = _build_normalized(d, norm_vals, relax_step, True, params,
                                step_xlim=xlim_early)
        # Re-title to indicate zoom
        fig.axes[0].set_title(
            "Phase-field interface relaxation — early transient (first 10 % of steps)",
            fontsize=params["title_fs"], pad=8)
        save_fig(fig, outdir, "relax_early", fmt)

        # --- 5. Physical time x-axis ---
        fig = _build_time_axis(d, norm_vals, relax_step, True, params)
        save_fig(fig, outdir, "relax_time", fmt)

    else:
        # Default: single combined figure at the path specified by --out
        fig = _build_combined(d, norm_vals, rel_change, relax_step,
                               args.threshold, False, params)
        out_path = args.out
        fig.savefig(out_path, dpi=params["dpi"], bbox_inches="tight")
        print(f"Saved: {out_path}")
        plt.close(fig)

    if args.show:
        matplotlib.use("TkAgg")
        plt.show()


if __name__ == "__main__":
    main()

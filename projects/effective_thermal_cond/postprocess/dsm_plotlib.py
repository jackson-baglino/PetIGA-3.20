#!/usr/bin/env python3
"""
plot_heatmaps.py — Build masked heatmaps of k_xx and k_yy vs time across many DSM runs.

USAGE
-----
python postprocess/plot_heatmaps.py PARENT_DIR \
  [--pattern 'REGEX_FOR_SUBFOLDERS'] \
  [--xaxis temperature|porosity] \
  [--no-normalize] \
  [--outdir OUTPUT_DIR]

WHAT IT DOES
------------
- Scans PARENT_DIR for immediate subfolders, optionally filtered by a Python regex (--pattern).
- Each subfolder is loaded as a "run" via dsm_plotlib.collect_runs():
    • reads SSA_evo.dat and k_eff.csv
    • aligns k to time using the shared "step/sol_index" column
    • extracts metadata (porosity, temperature, seed) from metadata.json or folder name
    • optionally normalizes (default: k_xx, k_yy and SSA divided by their first values)
- Creates a common time grid over the overlap of all runs, interpolates each run onto it,
  and produces masked heatmaps (no extrapolation) for k_xx and k_yy.
- Saves figures into OUTPUT_DIR (default: PARENT_DIR/output_plots).

OUTPUTS
-------
- OUTPUT_DIR/heatmap_kyy_vs_time_by_{xaxis}.png
- OUTPUT_DIR/heatmap_kxx_vs_time_by_{xaxis}.png  (if k_xx present in data)
- Also prints a short summary of loaded runs and where files are written.

NOTES
-----
- Use --xaxis porosity when you want vertical series sorted by φ rather than by temperature.
- The --pattern must be quoted in the shell to avoid expansion by your shell.
"""

from __future__ import annotations

import argparse
import os
from typing import List, Tuple

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# Local utilities
from dsm_plotlib import (
    set_default_rc,
    beautify_axes,
    autoscale_time,
    collect_runs,
    Run,
)


def build_common_time(runs: List[Run], Nt: int = 200) -> Tuple[np.ndarray, np.ndarray, str]:
    """Return (common_time, common_time_scaled, unit) over the *overlap* of all runs."""
    if not runs:
        return np.array([]), np.array([]), "seconds"
    starts = np.array([r.time[0] for r in runs if r.time.size > 0], dtype=float)
    ends   = np.array([r.time[-1] for r in runs if r.time.size > 0], dtype=float)
    if starts.size == 0 or ends.size == 0:
        return np.array([]), np.array([]), "seconds"

    t_min = float(np.nanmax(starts))  # latest start
    t_max = float(np.nanmin(ends))    # earliest end
    if not np.isfinite(t_min) or not np.isfinite(t_max) or t_max <= t_min:
        return np.array([]), np.array([]), "seconds"

    common_time = np.linspace(t_min, t_max, Nt)
    common_scaled, unit = autoscale_time(common_time)
    return common_time, common_scaled, unit


def interpolate_runs_onto_time(
    runs: List[Run], common_time: np.ndarray
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Interpolate k_yy (required) and k_xx (optional) onto common_time.
    Returns (kyy_map, kxx_map, mask) where each is shape [Nruns, Nt] masked array,
    and mask is True for valid (non-NaN) interpolated points.
    """
    Nr = len(runs)
    Nt = common_time.size
    kyy = np.full((Nr, Nt), np.nan, dtype=float)
    kxx = np.full((Nr, Nt), np.nan, dtype=float)

    for i, r in enumerate(runs):
        if r.time.size < 2:
            continue
        # Interpolate with NaN outside range (we only evaluate within overlap, so should be fine)
        fyy = interp1d(r.time, r.k_yy, bounds_error=False, fill_value=np.nan)
        kyy[i, :] = fyy(common_time)
        if r.k_xx is not None and r.k_xx.size == r.time.size:
            fxx = interp1d(r.time, r.k_xx, bounds_error=False, fill_value=np.nan)
            kxx[i, :] = fxx(common_time)

    # Masks: valid where finite
    mask_yy = np.isfinite(kyy)
    mask_xx = np.isfinite(kxx)

    kyy_ma = np.ma.array(kyy, mask=~mask_yy)
    kxx_ma = np.ma.array(kxx, mask=~mask_xx)
    return kyy_ma, kxx_ma, mask_yy


def main():
    ap = argparse.ArgumentParser(description="Plot masked heatmaps of k_xx and k_yy vs time across runs.")
    ap.add_argument("parent_dir", help="Directory containing run subfolders")
    ap.add_argument("--pattern", default=None, help="Python regex for subfolder names (quote in the shell)")
    ap.add_argument("--xaxis", choices=["temperature", "porosity"], default="temperature",
                    help="Variable used on the x-axis of the heatmap")
    ap.add_argument("--no-normalize", action="store_true", help="Disable normalization of k and SSA")
    ap.add_argument("--outdir", default=None, help="Output directory (default: PARENT_DIR/output_plots)")
    args = ap.parse_args()

    parent = os.path.abspath(os.path.expanduser(args.parent_dir))
    normalize = not args.no_normalize
    outdir = os.path.abspath(os.path.expanduser(args.outdir)) if args.outdir else os.path.join(parent, "output_plots")
    os.makedirs(outdir, exist_ok=True)

    set_default_rc()

    runs = collect_runs(parent, pattern=args.pattern, normalize=normalize)
    if not runs:
        print(f"No runs found in {parent}. Use --pattern to filter, or ensure folders contain SSA_evo.dat and k_eff.csv.")
        return

    # Sort runs by selected x-axis var
    if args.xaxis == "porosity":
        # use None-safe sort key; missing porosity to the end
        runs.sort(key=lambda r: (r.porosity is None, r.porosity))
        x_vals = np.array([r.porosity if r.porosity is not None else np.nan for r in runs], dtype=float)
        x_label = r"Porosity $\phi$"
    else:
        runs.sort(key=lambda r: (r.temperature is None, r.temperature))
        x_vals = np.array([r.temperature if r.temperature is not None else np.nan for r in runs], dtype=float)
        x_label = "Temperature (°C)"

    # Build common time
    common_time, common_time_scaled, unit = build_common_time(runs, Nt=240)
    if common_time.size == 0:
        print("⚠️  Not enough overlap in time across simulations to build heatmaps.")
        return

    # Interpolate onto the common grid
    kyy_map, kxx_map, _ = interpolate_runs_onto_time(runs, common_time)

    # --- Plot k_yy heatmap ---
    fig, ax = plt.subplots(figsize=(11.5, 4.0))
    # We want x along the horizontal axis. Our arrays are [Nruns, Nt]; pcolormesh expects Z with shape (Ny, Nx)
    # So we plot with Z.T where X=x_vals and Y=common_time_scaled
    c = ax.pcolormesh(x_vals, common_time_scaled, kyy_map.T, shading="auto", cmap="plasma")
    cb_label = r"$k_{yy}/k_{yy,0}$" if normalize else r"$k_{yy}$ (W/m·K)"
    fig.colorbar(c, ax=ax, label=cb_label)
    ax.set_xlabel(x_label)
    ax.set_ylabel(f"Time ({unit})")
    ax.set_title(r"$k_{yy}$ vs Time")
    beautify_axes(ax)
    fig.tight_layout()
    out_kyy = os.path.join(outdir, f"heatmap_kyy_vs_time_by_{'phi' if args.xaxis=='porosity' else 'temp'}.png")
    fig.savefig(out_kyy, dpi=300, bbox_inches="tight", transparent=True)
    plt.close(fig)

    # --- Plot k_xx heatmap (if any k_xx present) ---
    if np.isfinite(np.ma.masked_invalid(kxx_map).compressed()).any():
        fig, ax = plt.subplots(figsize=(11.5, 4.0))
        c = ax.pcolormesh(x_vals, common_time_scaled, kxx_map.T, shading="auto", cmap="plasma")
        cb_label = r"$k_{xx}/k_{xx,0}$" if normalize else r"$k_{xx}$ (W/m·K)"
        fig.colorbar(c, ax=ax, label=cb_label)
        ax.set_xlabel(x_label)
        ax.set_ylabel(f"Time ({unit})")
        ax.set_title(r"$k_{xx}$ vs Time")
        beautify_axes(ax)
        fig.tight_layout()
        out_kxx = os.path.join(outdir, f"heatmap_kxx_vs_time_by_{'phi' if args.xaxis=='porosity' else 'temp'}.png")
        fig.savefig(out_kxx, dpi=300, bbox_inches="tight", transparent=True)
        plt.close(fig)
    else:
        out_kxx = None

    # Summary
    print(f"\n✅ Loaded {len(runs)} runs from: {parent}")
    filled = lambda arr: int(np.isfinite(arr).sum())
    print(f"   • x-axis values ({'phi' if args.xaxis=='porosity' else 'temp'}): {filled(x_vals)}/{x_vals.size} finite")
    print(f"   • common time grid: Nt={common_time.size}, unit={unit}")
    print("\n🖼️  Saved:")
    print(f"   • {out_kyy}")
    if out_kxx:
        print(f"   • {out_kxx}")


if __name__ == "__main__":
    main()
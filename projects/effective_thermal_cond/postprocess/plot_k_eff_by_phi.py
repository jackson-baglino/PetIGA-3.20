#!/usr/bin/env python3
"""
Make color-plots (heatmaps) of k_eff for xx and yy across POROSITY and TIME
for a set of simulation result folders in a single parent directory.

Each subdirectory should contain:
  - k_eff.csv (columns: k_00/k_xx, k_11/k_yy, sol_index/step)
  - SSA_evo.dat (cols: SSA, ?, time, step)  [used to map step -> time]
  - metadata.json (optional; may contain porosity fields)

Porosity (phi) is taken from the folder name if possible, else from metadata.json.

Outputs:
  <outdir>/heatmap_kxx_vs_time_by_phi.png
  <outdir>/heatmap_kyy_vs_time_by_phi.png   (if data present)

Usage:
  python postprocess/plot_k_eff_by_phi.py <PARENT_DIR> \
    [--pattern 'REGEX'] [--outdir PATH] [--normalize true|false] [--dpi 300] [--cmap plasma] [--theme light|dark|both]

Examples:
  # All subfolders:
  python postprocess/plot_k_eff_by_phi.py ~/SimulationResults/effective_thermal_cond/const_temp_varying_phi

  # Only match DSM_phi... folders (quote to avoid shell globbing):
  python postprocess/plot_k_eff_by_phi.py ~/.../const_temp_varying_phi \
    --pattern '^DSM_phi[0-9]+\\.[0-9]+_Lx[0-9]+mm_Ly[0-9]+mm_seed[0-9]+_Tm-?[0-9]+_hum[0-9]+_tf[0-9]+d$'

Notes:
- If SSA_evo.dat is missing or malformed, the script falls back to a synthetic
  monotonic "time" based on the k_eff row index (plots will still render).
- Normalization (default true) divides each k-series by its first value in that run.
- Theme: use --theme light (white background, dark text), --theme dark (dark background, light text), or --theme both to save both versions with _light/_dark suffixes.
"""

from __future__ import annotations

import argparse
import json
import os
import re
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# --- Stable discrete color mapping for porosity values ---
def build_phi_color_lookup(phis: np.ndarray):
    """
    Return a function color_for(phi) that assigns a stable, discrete color
    to each unique porosity (rounded to 3 decimals). This avoids shifts in
    color when the set of runs changes.
    """
    phis = np.asarray(phis, dtype=float)
    if phis.size == 0:
        def _fallback(_):
            return plt.get_cmap("viridis")(0.5)
        return _fallback
    # Round to stabilize keys, then unique + sorted for deterministic ordering
    keys = np.unique(np.round(phis, 3))
    # Build a long discrete palette (tab20 + tab20b + tab20c)
    palette = list(plt.cm.tab20.colors) + list(plt.cm.tab20b.colors) + list(plt.cm.tab20c.colors)
    # If still too short, fall back to viridis samples
    if len(palette) < len(keys):
        extra = plt.get_cmap("viridis")(np.linspace(0.05, 0.95, len(keys) - len(palette)))
        palette.extend([tuple(c) for c in extra])
    mapping = {k: palette[i % len(palette)] for i, k in enumerate(keys)}
    def color_for(phi_value: float):
        k = float(np.round(phi_value, 3))
        return mapping.get(k, plt.get_cmap("viridis")(0.5))
    return color_for



def set_theme(theme: str):
    """Apply a light or dark Matplotlib theme."""
    import matplotlib as mpl
    plt.rcParams.update(plt.rcParamsDefault)
    if theme == "dark":
        plt.rcParams.update({
            'figure.facecolor': '#000000',
            'axes.facecolor': '#000000',
            'savefig.facecolor': '#000000',
            'text.color': 'white',
            'axes.labelcolor': 'white',
            'axes.edgecolor': 'white',
            'xtick.color': 'white',
            'ytick.color': 'white',
            'grid.color': '0.35',
        })
    else:  # light
        plt.rcParams.update({
            'figure.facecolor': 'white',
            'axes.facecolor': 'white',
            'savefig.facecolor': 'white',
            'text.color': 'black',
            'axes.labelcolor': 'black',
            'axes.edgecolor': 'black',
            'xtick.color': 'black',
            'ytick.color': 'black',
            'grid.color': '0.7',
        })
    plt.rcParams.update({
        'font.family': 'DejaVu Sans',
        'font.size': 14,
        'axes.linewidth': 1.2,
        'xtick.direction': 'in', 'ytick.direction': 'in',
        'xtick.major.size': 6,   'ytick.major.size': 6,
        'xtick.major.width': 1.2,'ytick.major.width': 1.2,
        'legend.fontsize': 12,
    })


def bool_arg(x: str) -> bool:
    return str(x).lower() in ("1", "true", "yes", "y", "on")


def discover_folders(parent_dir: str, pattern: Optional[str]) -> List[str]:
    names = []
    for n in sorted(os.listdir(parent_dir)):
        full = os.path.join(parent_dir, n)
        if os.path.isdir(full):
            names.append(n)
    if not pattern:
        return names
    rx = re.compile(pattern)
    return [n for n in names if rx.match(n)]


def load_json_if_exists(p: str) -> dict:
    try:
        with open(p, "r") as f:
            return json.load(f)
    except Exception:
        return {}


def coerce_float(x):
    try:
        return float(x)
    except Exception:
        return None


def extract_phi_from_name(name: str) -> Optional[float]:
    m = re.search(r"phi(\d+(?:\.\d+)?)", name)
    return coerce_float(m.group(1)) if m else None


def extract_phi_from_meta(meta: dict) -> Optional[float]:
    # try a few likely keys
    for k in ("porosity", "porosity_target", "phi", "phi0"):
        if k in meta:
            v = meta[k]
            if isinstance(v, (int, float)):
                return float(v)
            if isinstance(v, str):
                try:
                    return float(v.strip())
                except Exception:
                    pass
    # sometimes nested under structure/…
    try:
        s = meta.get("structure", {})
        for k in ("porosity_target", "porosity"):
            if k in s:
                return float(s[k])
    except Exception:
        pass
    return None


def load_run(folder: str, normalize: bool) -> Optional[dict]:
    """
    Return dict with fields:
      name, path, phi, time (1D), kxx (1D), kyy (1D)
    or None if not loadable.
    """
    name = os.path.basename(folder)
    keff = os.path.join(folder, "k_eff.csv")
    ssa  = os.path.join(folder, "SSA_evo.dat")

    if not os.path.exists(keff):
        return None

    # Map step -> time and SSA from SSA_evo.dat
    step_to_time = None
    step_to_ssa = None
    if os.path.exists(ssa):
        try:
            arr = np.loadtxt(ssa)
            if arr.ndim == 1:
                arr = arr[None, :]
            # SSA_evo.dat (cols: SSA, ?, time, step)
            ssa_val = arr[:, 0]
            ssa_time = arr[:, 2]
            ssa_step = arr[:, 3].astype(int)
            step_to_time = {int(s): float(t) for s, t in zip(ssa_step, ssa_time)}
            step_to_ssa  = {int(s): float(v) for s, v in zip(ssa_step, ssa_val)}
        except Exception:
            step_to_time = None
            step_to_ssa = None

    # Load k_eff.csv
    df = pd.read_csv(keff)
    col_kxx = "k_00" if "k_00" in df.columns else ("k_xx" if "k_xx" in df.columns else None)
    col_kyy = "k_11" if "k_11" in df.columns else ("k_yy" if "k_yy" in df.columns else None)
    col_step = "sol_index" if "sol_index" in df.columns else ("step" if "step" in df.columns else None)
    if not (col_kxx and col_kyy and col_step):
        # missing required columns
        return None

    kxx = df[col_kxx].to_numpy(dtype=float)
    kyy = df[col_kyy].to_numpy(dtype=float)
    steps = df[col_step].to_numpy(dtype=int)

    # SSA aligned to k_eff steps (if available)
    ssa_for_steps = None
    if step_to_ssa is not None:
        ssa_for_steps = np.array([step_to_ssa.get(int(s), np.nan) for s in steps], dtype=float)

    # Build time vector
    if step_to_time:
        time = np.array([step_to_time.get(int(s), np.nan) for s in steps], dtype=float)
        # trim leading/trailing NaNs consistently
        valid = np.isfinite(time) & np.isfinite(kxx) & np.isfinite(kyy)
        if valid.sum() >= 2:
            i0 = int(np.argmax(valid))
            i1 = len(valid) - int(np.argmax(valid[::-1]))
            time = time[i0:i1]
            kxx  = kxx[i0:i1]
            kyy  = kyy[i0:i1]
            if ssa_for_steps is not None:
                ssa_for_steps = ssa_for_steps[i0:i1]
        else:
            # fall back to synthetic increasing "time"
            time = np.arange(len(steps), dtype=float)
    else:
        time = np.arange(len(steps), dtype=float)

    # Normalize (optional)
    if normalize:
        if kxx.size:
            kxx = kxx / kxx[0]
        if kyy.size:
            kyy = kyy / kyy[0]

    # Porosity
    phi = extract_phi_from_name(name)
    if phi is None:
        meta = {}
        for cand in ("metadata.json", "Metadata.json", "meta.json", "config.json"):
            p = os.path.join(folder, cand)
            if os.path.exists(p):
                meta = load_json_if_exists(p)
                break
        phi = extract_phi_from_meta(meta)

    if phi is None:
        return None

    return dict(name=name, path=folder, phi=float(phi), time=time, kxx=kxx, kyy=kyy, ssa=ssa_for_steps, steps=steps)


def autoscale_time(arr: np.ndarray) -> Tuple[np.ndarray, str]:
    unit = "seconds"
    scaled = arr
    if arr.size:
        tmax = float(np.nanmax(arr))
        if tmax >= 2 * 24 * 3600:
            scaled, unit = arr / 86400.0, "days"
        elif tmax >= 2 * 3600:
            scaled, unit = arr / 3600.0, "hours"
        elif tmax >= 120:
            scaled, unit = arr / 60.0, "minutes"
    return scaled, unit


def to_edges_from_centers(x: np.ndarray) -> np.ndarray:
    """Convert 1D array of centers to edges for pcolormesh."""
    x = np.asarray(x, dtype=float)
    if x.size < 2:
        dx = 1.0
        return np.array([x[0] - dx/2, x[0] + dx/2])
    mids = 0.5 * (x[1:] + x[:-1])
    edges = np.concatenate(([x[0] - (mids[0] - x[0])], mids, [x[-1] + (x[-1] - mids[-1])]))
    return edges


def main():
    ap = argparse.ArgumentParser(description="Plot k_eff heatmaps vs porosity (x) and time (y).")
    ap.add_argument("parent_dir", help="Folder containing run subdirectories")
    ap.add_argument("--pattern", help="Python regex to filter subfolder names", default=None)
    ap.add_argument("--outdir", help="Where to save plots (default: <parent_dir>/output_plots)", default=None)
    ap.add_argument("--normalize", type=bool_arg, default=True, help="Normalize k by first value (default: true)")
    ap.add_argument("--dpi", type=int, default=300)
    ap.add_argument("--cmap", type=str, default="plasma")
    ap.add_argument("--theme", type=str, default="light", choices=["light", "dark", "both"],
                    help="Color theme: 'light' (default), 'dark', or 'both' to export both.")
    args = ap.parse_args()

    parent_dir = os.path.abspath(os.path.expanduser(args.parent_dir))
    outdir = os.path.abspath(os.path.expanduser(args.outdir)) if args.outdir else os.path.join(parent_dir, "output_plots")
    os.makedirs(outdir, exist_ok=True)

    # Determine which themes to render
    themes = [args.theme] if args.theme in ("light", "dark") else ["light", "dark"]

    folders = discover_folders(parent_dir, args.pattern)
    if not folders:
        print(f"[WARN] No subfolders found under {parent_dir}. Check --pattern.")
        return

    runs = []
    for n in folders:
        r = load_run(os.path.join(parent_dir, n), normalize=args.normalize)
        if r is not None:
            runs.append(r)

    if len(runs) < 2:
        print(f"[WARN] Found only {len(runs)} usable run(s). Need ≥2 to make a φ×time heatmap.")
        return

    # Sort runs by porosity
    runs.sort(key=lambda d: d["phi"])
    phis = np.array([d["phi"] for d in runs], dtype=float)

    # Build stable color lookup for porosity
    color_for_phi = build_phi_color_lookup(phis)

    # Build common time grid across overlap (to avoid extrapolation)
    t_mins = []
    t_maxs = []
    for d in runs:
        t = d["time"]
        if t.size:
            t_mins.append(np.nanmin(t))
            t_maxs.append(np.nanmax(t))
    if not t_mins or not t_maxs:
        print("[WARN] Could not determine time extents.")
        return
    t0 = max(t_mins)
    t1 = min(t_maxs)
    if not np.isfinite(t0) or not np.isfinite(t1) or t1 <= t0:
        print("[WARN] Not enough common time overlap across runs.")
        return

    Nt = 250
    common_time = np.linspace(t0, t1, Nt)

    # Interpolate each run's k onto the common time
    kxx_rows = []
    kyy_rows = []
    for d in runs:
        t = d["time"]
        kxx = d["kxx"]
        kyy = d["kyy"]

        # guard against duplicate times
        ord_idx = np.argsort(t)
        t_sorted = t[ord_idx]
        kxx_sorted = kxx[ord_idx]
        kyy_sorted = kyy[ord_idx]

        # numpy interp requires finite values
        valid_xx = np.isfinite(t_sorted) & np.isfinite(kxx_sorted)
        valid_yy = np.isfinite(t_sorted) & np.isfinite(kyy_sorted)

        if valid_xx.sum() >= 2:
            kxxi = np.interp(common_time, t_sorted[valid_xx], kxx_sorted[valid_xx], left=np.nan, right=np.nan)
        else:
            kxxi = np.full_like(common_time, np.nan)

        if valid_yy.sum() >= 2:
            kyyi = np.interp(common_time, t_sorted[valid_yy], kyy_sorted[valid_yy], left=np.nan, right=np.nan)
        else:
            kyyi = np.full_like(common_time, np.nan)

        kxx_rows.append(kxxi)
        kyy_rows.append(kyyi)

    kxx_map = np.array(kxx_rows)   # shape: (Nphi, Nt)
    kyy_map = np.array(kyy_rows)

    # Time scaling for nicer labels
    t_scaled, t_unit = autoscale_time(common_time)

    # Build edges for pcolormesh (non-uniform φ is OK)
    phi_edges = to_edges_from_centers(phis)
    t_edges = to_edges_from_centers(t_scaled)

    def plot_heat(map2d: np.ndarray, label: str, fname: str, theme: str):
        # Apply theme-specific rcParams
        set_theme(theme)

        fig, ax = plt.subplots(figsize=(10.0, 4.4))
        # pcolormesh expects Z.shape == (len(y)-1, len(x)-1) in (Y, X) order; our map is (Nphi, Nt).
        # So transpose to (Nt, Nphi) to align Y (time) x X (phi).
        c = ax.pcolormesh(phi_edges, t_edges, map2d.T, shading="auto", cmap=args.cmap)
        cb = fig.colorbar(c, ax=ax, label=label)
        ax.set_xlabel(r"Porosity $\phi$")
        ax.set_ylabel(f"Time ({t_unit})")
        ax.set_title(label + r" vs $\phi$ and time")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(True, alpha=0.35)
        fig.tight_layout()
        base, ext = os.path.splitext(fname)
        suffix = "" if len(themes) == 1 else f"_{theme}"
        out = os.path.join(outdir, f"{base}{suffix}{ext}")
        fig.savefig(out, dpi=args.dpi, bbox_inches="tight", transparent=True)
        plt.close(fig)
        print("[OK]", out)

    # Labels depend on normalization
    if args.normalize:
        for theme in themes:
            plot_heat(kxx_map, r"$k_{xx}/k_{xx,0}$", "heatmap_kxx_vs_time_by_phi.png", theme)
            plot_heat(kyy_map, r"$k_{yy}/k_{yy,0}$", "heatmap_kyy_vs_time_by_phi.png", theme)
    else:
        for theme in themes:
            plot_heat(kxx_map, r"$k_{xx}$ (W m$^{-1}$ K$^{-1}$)", "heatmap_kxx_vs_time_by_phi.png", theme)
            plot_heat(kyy_map, r"$k_{yy}$ (W m$^{-1}$ K$^{-1}$)", "heatmap_kyy_vs_time_by_phi.png", theme)

    # Quick summary
    print(f"\nDone. Runs used: {len(runs)}")
    for d in runs:
        print(f"  • {d['name']}  (phi={d['phi']:.2f}, Nt={len(d['time'])})")

    # ---------- Additional curve plots (conference-ready) ----------

    def _global_time_scale(runs_list):
        # Concatenate times to choose a consistent unit
        all_t = np.concatenate([d["time"] for d in runs_list if d["time"].size > 0]) if runs_list else np.array([])
        scaled, unit = autoscale_time(all_t if all_t.size else np.array([0.0]))
        # Determine factor to convert each run's time to the chosen unit
        if unit == "days":
            factor = 1.0 / 86400.0
        elif unit == "hours":
            factor = 1.0 / 3600.0
        elif unit == "minutes":
            factor = 1.0 / 60.0
        else:
            factor = 1.0
            unit = "seconds"
        return (lambda arr: arr * factor), unit

    time_scale, unit_str = _global_time_scale(runs)

    def plot_time_curves(runs_list, ylabel, key, fname, theme):
        set_theme(theme)
        fig, ax = plt.subplots(figsize=(10.0, 4.4))
        # Color by porosity for consistency (discrete, stable)
        for d in runs_list:
            col = color_for_phi(d["phi"])
            t_plot = time_scale(d["time"])
            y = d[key]
            ax.plot(
                t_plot, y, lw=2.0, alpha=0.95, label=fr"$\phi={d['phi']:.2f}$", color=col,
                marker='o', ms=3, markevery=max(1, len(t_plot)//18)
            )
        ax.set_xlabel(f"Time ({unit_str})")
        ax.set_ylabel(ylabel)
        ax.set_title(f"{ylabel} vs time")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(True, alpha=0.35)
        n = len(runs_list)
        if n > 6:
            # Figure-level legend at right
            handles, labels = ax.get_legend_handles_labels()
            fig.legend(handles, labels, loc="center right", bbox_to_anchor=(1.02, 0.5), title="Porosity")
            fig.subplots_adjust(right=0.84)
            fig.tight_layout(rect=[0.0, 0.0, 0.84, 1.0])
        else:
            ax.legend()
            fig.tight_layout()
        base, ext = os.path.splitext(fname)
        suffix = "" if len(themes) == 1 else f"_{theme}"
        out = os.path.join(outdir, f"{base}{suffix}{ext}")
        fig.savefig(out, dpi=args.dpi, bbox_inches="tight", transparent=True)
        plt.close(fig)
        print("[OK]", out)

    def plot_norm_vs_ssa_curves(runs_list, ylabel, key, fname, theme):
        # Only include runs with valid SSA
        valid_runs = []
        for d in runs_list:
            ssa = d.get("ssa", None)
            if ssa is None:
                continue
            # Need at least two finite points and finite start
            finite = np.isfinite(ssa) & np.isfinite(d[key])
            if finite.sum() >= 2:
                ssa = ssa[finite]
                kval = d[key][finite]
                # Normalize
                ssa0 = ssa[0]
                k0 = kval[0]
                if np.isfinite(ssa0) and np.isfinite(k0) and k0 != 0:
                    ssa_n = ssa / ssa0
                    k_n = kval / k0
                    valid_runs.append((d["phi"], ssa_n, k_n))

        if not valid_runs:
            print("[WARN] No runs with usable SSA for normalized SSA plots; skipping", fname)
            return

        set_theme(theme)
        fig, ax = plt.subplots(figsize=(10.0, 4.4))
        for phi, ssa_n, k_n in valid_runs:
            col = color_for_phi(phi)
            ax.plot(
                ssa_n, k_n, lw=2.0, alpha=0.95, label=fr"$\phi={phi:.2f}$", color=col,
                marker='o', ms=3, markevery=max(1, len(ssa_n)//18)
            )
        ax.set_xlabel(r"SSA/SSA$_0$")
        ax.set_ylabel(ylabel)
        ax.set_title(f"{ylabel} vs SSA/SSA$_0$")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(True, alpha=0.35)
        n = len(valid_runs)
        if n > 6:
            handles, labels = ax.get_legend_handles_labels()
            fig.legend(handles, labels, loc="center right", bbox_to_anchor=(1.02, 0.5), title="Porosity")
            fig.subplots_adjust(right=0.84)
            fig.tight_layout(rect=[0.0, 0.0, 0.84, 1.0])
        else:
            ax.legend()
            fig.tight_layout()
        base, ext = os.path.splitext(fname)
        suffix = "" if len(themes) == 1 else f"_{theme}"
        out = os.path.join(outdir, f"{base}{suffix}{ext}")
        fig.savefig(out, dpi=args.dpi, bbox_inches="tight", transparent=True)
        plt.close(fig)
        print("[OK]", out)

    # Generate the new curve plots per theme
    for theme in themes:
        # Absolute (not normalized)
        if args.normalize:
            # If heatmaps were normalized, still show *absolute* time curves using current arrays
            plot_time_curves(runs, r"$k_{xx}/k_{xx,0}$" if args.normalize else r"$k_{xx}$ (W m$^{-1}$ K$^{-1}$)",
                             key="kxx", fname="curves_kxx_vs_time.png", theme=theme)
            plot_time_curves(runs, r"$k_{yy}/k_{yy,0}$" if args.normalize else r"$k_{yy}$ (W m$^{-1}$ K$^{-1}$)",
                             key="kyy", fname="curves_kyy_vs_time.png", theme=theme)
        else:
            plot_time_curves(runs, r"$k_{xx}$ (W m$^{-1}$ K$^{-1}$)",
                             key="kxx", fname="curves_kxx_vs_time.png", theme=theme)
            plot_time_curves(runs, r"$k_{yy}$ (W m$^{-1}$ K$^{-1}$)",
                             key="kyy", fname="curves_kyy_vs_time.png", theme=theme)

        # Normalized vs SSA/SSA0
        plot_norm_vs_ssa_curves(runs, r"$k_{xx}/k_{xx,0}$", key="kxx", fname="curves_kxx_norm_vs_SSA.png", theme=theme)
        plot_norm_vs_ssa_curves(runs, r"$k_{yy}/k_{yy,0}$", key="kyy", fname="curves_kyy_norm_vs_SSA.png", theme=theme)


if __name__ == "__main__":
    main()
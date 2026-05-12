#!/usr/bin/env python3
"""
plot_single_run_k_tensor.py

Generate per-simulation figures:
  1) k_xx vs time
  2) k_yy vs time
  3) k_xx vs SSA
  4) k_yy vs SSA
  5) k_xx & k_yy vs time (same y-axis)
  6) k_xx & k_yy vs SSA (same y-axis)

Features
- Light & dark themed outputs (both by default; choose via --theme)
- Optional normalization of k-components and SSA by each series' first value
- Robust to 2D (k_00, k_11) or 3D (k_00, k_11, k_22) CSVs
- Aligns k rows to SSA using 'sol_index' (or 'step') as the join key
- Auto-scales time axis to seconds/minutes/hours/days

Usage
  python plot_single_run_k_tensor.py /path/to/one/sim/folder
  python plot_single_run_k_tensor.py /path/to/folder --theme dark
  python plot_single_run_k_tensor.py /path/to/folder --no-norm-k --no-norm-ssa
"""

import os
import argparse
import glob
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ======= Flags (also exposed via CLI) =======
NORMALIZE_K_DEFAULT = True     # normalize k_ij by its first value
NORMALIZE_SSA_DEFAULT = True   # normalize SSA by its first value

# ======= Minimal styling helpers =======

def set_theme(theme: str):
    """Apply a light or dark theme."""
    if theme not in ("light", "dark"):
        theme = "light"

    # Reset first
    plt.rcParams.update(plt.rcParamsDefault)

    base = {
        'font.family': 'Arial',
        'font.size': 18,
        'axes.labelsize': 18,
        'axes.titlesize': 18,
        'legend.fontsize': 16,
        'xtick.labelsize': 16,
        'ytick.labelsize': 16,
        'axes.linewidth': 1.2,
        'xtick.direction': 'in',
        'ytick.direction': 'in',
        'xtick.major.width': 1.2,
        'ytick.major.width': 1.2,
        'figure.dpi': 100,
    }
    plt.rcParams.update(base)

    if theme == "light":
        plt.rcParams.update({
            'figure.facecolor': 'white',
            'axes.facecolor': 'white',
            'axes.edgecolor': 'black',
            'axes.labelcolor': 'black',
            'xtick.color': 'black',
            'ytick.color': 'black',
            'text.color': 'black',
            'savefig.facecolor': 'white',
            'savefig.edgecolor': 'white',
        })
        # default color cycle is fine
    else:
        # dark
        plt.rcParams.update({
            'figure.facecolor': 'black',
            'axes.facecolor': 'black',
            'axes.edgecolor': 'white',
            'axes.labelcolor': 'white',
            'xtick.color': 'white',
            'ytick.color': 'white',
            'text.color': 'white',
            'savefig.facecolor': 'black',
            'savefig.edgecolor': 'black',
        })
        # High-contrast cycle
        plt.rcParams['axes.prop_cycle'] = plt.cycler(color=[
            '#8ecae6','#ffd166','#ef476f','#06d6a0','#adb5bd','#f4a261','#e76f51'
        ])

def beautify(ax):
    ax.grid(False)
    # spines already in theme colors
    return ax

# ======= IO / Alignment =======

def autoscale_time(arr):
    unit = 'seconds'
    scaled = arr
    if arr.size > 0:
        tmax = float(np.nanmax(arr))
        if tmax >= 2 * 24 * 3600:
            scaled, unit = arr / 86400.0, 'days'
        elif tmax >= 2 * 3600:
            scaled, unit = arr / 3600.0, 'hours'
        elif tmax >= 120:
            scaled, unit = arr / 60.0, 'minutes'
    return scaled, unit

def load_ssa(folder):
    """Locate and load SSA time series.
       Prefer `SSA_evo.dat`, but accept common variants via globbing.
       Expected columns (order will be read as in SSA_evo.dat spec):
         [0]=SSA, [2]=time(s), [3]=step (int)
    """
    # Preferred filename
    pref = os.path.join(folder, 'SSA_evo.dat')
    candidates = []
    if os.path.exists(pref):
        candidates.append(pref)
    else:
        # Try a few permissive patterns
        patterns = [
            'SSA*_evo*.dat', 'SSA*evo*.dat', '*SSA*evo*.dat',
            'SSA*.dat', '*ssa*.dat'
        ]
        for pat in patterns:
            for p in glob.glob(os.path.join(folder, pat)):
                candidates.append(p)

    for path in candidates:
        try:
            data = np.loadtxt(path)
            if data.ndim != 2 or data.shape[1] < 4:
                print(f"[load_ssa] '{os.path.basename(path)}' has unexpected shape {data.shape}, skipping")
                continue
            ssa      = data[:, 0]
            time_s   = data[:, 2]
            step_int = data[:, 3].astype(int)
            print(f"[load_ssa] using '{os.path.basename(path)}' with {len(step_int)} rows")
            return dict(ssa=ssa, time=time_s, step=step_int)
        except Exception as e:
            print(f"[load_ssa] failed to read {path}: {e}")
            continue

    return None

def _first_present(cols, names):
    for n in names:
        if n in cols: return n
    return None

def load_k(folder):
    """Locate and load k_eff CSV.
       Accept `k_eff.csv` or variants (e.g., 'keff.csv', 'k-effective.csv').
       Tries columns in this order for the join key and components.
    """
    pref = os.path.join(folder, 'k_eff.csv')
    candidates = []
    if os.path.exists(pref):
        candidates.append(pref)
    else:
        for pat in ['k*eff*.csv', 'keff*.csv', '*k*.csv']:
            for p in glob.glob(os.path.join(folder, pat)):
                candidates.append(p)

    for path in candidates:
        try:
            df = pd.read_csv(path)
        except Exception as e:
            print(f"[load_k] failed to read {path}: {e}")
            continue

        cols = df.columns
        step_col = _first_present(cols, ['sol_index', 'step', 'index'])
        if step_col is None:
            print(f"[load_k] '{os.path.basename(path)}' has no step/index column, skipping")
            continue

        # Direct names first
        k11 = _first_present(cols, ['k_11', 'k11', 'k_yy'])
        k22 = _first_present(cols, ['k_22', 'k22', 'k_zz'])

        # Robust 2D fallback: if k_00 & k_11 exist but no k_22, map them to two series
        if 'k_00' in cols and 'k_11' in cols and k22 is None and k11 is None:
            k11 = 'k_00'
            k22 = 'k_11'
        elif 'k_00' in cols and 'k_11' in cols and k22 is None and k11 is not None:
            # if k11 already chosen (e.g., 'k_11'), use k_00 as the second series
            if k11 == 'k_11':
                k22 = 'k_00'
            else:
                k22 = 'k_11'

        # If still missing both, try to guess from any k_* columns
        if k11 is None and k22 is None:
            kcols = [c for c in cols if c.startswith('k_')]
            if len(kcols) >= 2:
                k11, k22 = kcols[0], kcols[1]
            elif len(kcols) == 1:
                k11 = kcols[0]

        data = dict(
            step=df[step_col].astype(int).values,
            k11=df[k11].values if k11 else None,
            k22=df[k22].values if k22 else None,
        )
        print(f"[load_k] using '{os.path.basename(path)}' | step='{step_col}', k11='{k11}', k22='{k22}'")
        return data

    return None

def align_by_step(ssa_dict, k_dict, norm_k=True, norm_ssa=True):
    """Return aligned arrays: time, ssa, k11, k22, plus time unit label."""
    step_to_t = {int(s): t for s, t in zip(ssa_dict['step'], ssa_dict['time'])}
    step_to_s = {int(s): v for s, v in zip(ssa_dict['step'], ssa_dict['ssa'])}

    times, ssa, k11, k22 = [], [], [], []
    for i, st in enumerate(k_dict['step']):
        if st in step_to_t:
            times.append(step_to_t[st])
            ssa.append(step_to_s.get(st, np.nan))
            k11.append(k_dict['k11'][i] if k_dict['k11'] is not None else np.nan)
            k22.append(k_dict['k22'][i] if k_dict['k22'] is not None else np.nan)

    times = np.asarray(times, float)
    ssa   = np.asarray(ssa, float)
    k11   = np.asarray(k11, float)
    k22   = np.asarray(k22, float)

    # Normalize (per-series) if requested and have at least one value
    if norm_k and k11.size: k11 = k11 / k11[0]
    if norm_k and k22.size and np.isfinite(k22[0]): k22 = k22 / k22[0]
    if norm_ssa and ssa.size: ssa = ssa / ssa[0]

    t_scaled, unit = autoscale_time(times)
    print(f"[align] aligned points: n={times.size}, with SSA={np.isfinite(ssa).sum()}, k11={np.isfinite(k11).sum()}, k22={np.isfinite(k22).sum()}")
    return times, t_scaled, unit, ssa, k11, k22

# ======= Plotting =======

def savefig(fig, folder, stem, theme, dpi):
    out_dir = os.path.join(folder, "output_plots")
    os.makedirs(out_dir, exist_ok=True)
    name = f"{stem}_{theme}.png"
    fig.savefig(os.path.join(out_dir, name), dpi=dpi, bbox_inches='tight', transparent=True)
    plt.close(fig)

def plot_series_vs_time(folder, t_scaled, unit, y, ylabel, title, theme, color='C0', marker='o', dpi=600):
    set_theme(theme)
    fig, ax = plt.subplots(figsize=(11, 4))
    ax.plot(t_scaled, y, lw=1.8, marker=marker)
    ax.set_xlabel(f"Time ({unit})")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    beautify(ax)
    fig.tight_layout()
    return fig

def plot_series_vs_ssa(folder, ssa, y, xlabel, ylabel, title, theme, color='C0', marker='o', dpi=600):
    set_theme(theme)
    fig, ax = plt.subplots(figsize=(11, 4))
    ax.plot(ssa, y, lw=1.8, marker=marker)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    beautify(ax)
    fig.tight_layout()
    return fig

def plot_dual_same_y_vs_time(folder, t_scaled, unit, y1, y2, labels, theme, dpi=600):
    set_theme(theme)
    fig, ax = plt.subplots(figsize=(11.5, 4.2))
    ax.plot(t_scaled, y1, lw=2.0, marker='o', label=labels[0])
    ax.plot(t_scaled, y2, lw=2.0, marker='s', label=labels[1])
    ax.set_xlabel(f"Time ({unit})")
    ax.set_ylabel(labels[2])  # common ylabel
    ax.set_title(labels[3])   # title
    beautify(ax)
    ax.legend()
    fig.tight_layout()
    return fig

def plot_dual_same_y_vs_ssa(folder, ssa, y1, y2, labels, theme, dpi=600):
    set_theme(theme)
    fig, ax = plt.subplots(figsize=(11.5, 4.2))
    ax.plot(ssa, y1, lw=2.0, marker='o', label=labels[0])
    ax.plot(ssa, y2, lw=2.0, marker='s', label=labels[1])
    ax.set_xlabel(labels[2])  # x label (SSA or SSA/SSA0)
    ax.set_ylabel(labels[3])  # y label (k/k0 or W/m·K)
    ax.set_title(labels[4])   # title
    beautify(ax)
    ax.legend()
    fig.tight_layout()
    return fig

# ======= Main =======

def main():
    ap = argparse.ArgumentParser(description="Per-simulation k_xx, k_yy plotting (time & SSA), light/dark themes")
    ap.add_argument("folder", help="Path to a single simulation result folder")
    ap.add_argument("--theme", choices=["light","dark","both"], default="both")
    ap.add_argument("--no-norm-k", action="store_true", help="Disable normalization of k-values")
    ap.add_argument("--no-norm-ssa", action="store_true", help="Disable normalization of SSA")
    ap.add_argument("--dpi", type=int, default=600)
    args = ap.parse_args()

    norm_k = NORMALIZE_K_DEFAULT and not args.no_norm_k
    norm_s = NORMALIZE_SSA_DEFAULT and not args.no_norm_ssa

    folder = os.path.abspath(args.folder)
    ssa_dict = load_ssa(folder)
    k_dict   = load_k(folder)

    if ssa_dict is None or k_dict is None:
        present = ', '.join(sorted(os.listdir(folder)))
        missing = []
        if ssa_dict is None: missing.append("SSA_evo.dat (or variant)")
        if k_dict is None:   missing.append("k_eff.csv (or variant)")
        print(f"[error] In '{folder}' — could not find: {', '.join(missing)}")
        print(f"[error] Files present: {present}")
        raise SystemExit("Missing required input files. See log above for details.")

    times, t_scaled, unit, ssa, k11, k22 = align_by_step(ssa_dict, k_dict, norm_k, norm_s)

    # Build labels based on normalization flags
    yk = r'$k/k_0$' if norm_k else r'$k$ (W/m·K)'
    xssa = r'SSA/SSA$_0$' if norm_s else 'SSA'

    # Helpful title (from metadata if present)
    title_extra = ""
    meta_path = None
    for name in ("metadata.json","Metadata.json","meta.json","config.json"):
        p = os.path.join(folder, name)
        if os.path.exists(p): meta_path = p; break
    if meta_path:
        try:
            with open(meta_path,'r') as f: meta = json.load(f)
            # Try a few fields for compact context
            T = None
            for k in ("temperature_C","temperature","temp","T"):
                if isinstance(meta, dict):
                    # safe dive into nested dicts if needed would be overkill; keep simple
                    v = meta.get(k)
                    if v is not None: T = v; break
            if T is not None:
                title_extra = f" (T={T})"
        except Exception:
            pass

    themes = [args.theme] if args.theme in ("light","dark") else ["light","dark"]

    # 1) k_11 vs time
    if np.isfinite(k11).any():
        for th in themes:
            lab = r'$k_{xx}/k_{xx,0}$' if norm_k else r'$k_{xx}$ (W/m·K)'
            fig = plot_series_vs_time(folder, t_scaled, unit, k11, lab, f'$k_{{xx}}$ vs Time{title_extra}', th, dpi=args.dpi)
            savefig(fig, folder, 'k11_vs_time', th, args.dpi)

    # 2) k_22 vs time
    if np.isfinite(k22).any():
        for th in themes:
            lab = r'$k_{yy}/k_{yy,0}$' if norm_k else r'$k_{yy}$ (W/m·K)'
            fig = plot_series_vs_time(folder, t_scaled, unit, k22, lab, f'$k_{{yy}}$ vs Time{title_extra}', th, dpi=args.dpi)
            savefig(fig, folder, 'k22_vs_time', th, args.dpi)

    # 3) k_11 vs SSA
    valid = np.isfinite(ssa) & np.isfinite(k11)
    if np.count_nonzero(valid):
        for th in themes:
            ylab = r'$k_{xx}/k_{xx,0}$' if norm_k else r'$k_{xx}$ (W/m·K)'
            fig = plot_series_vs_ssa(folder, ssa[valid], k11[valid], xssa, ylab, f'$k_{{xx}}$ vs SSA{title_extra}', th, dpi=args.dpi)
            savefig(fig, folder, 'k11_vs_ssa', th, args.dpi)

    # 4) k_22 vs SSA
    valid = np.isfinite(ssa) & np.isfinite(k22)
    if np.count_nonzero(valid):
        for th in themes:
            ylab = r'$k_{yy}/k_{yy,0}$' if norm_k else r'$k_{yy}$ (W/m·K)'
            fig = plot_series_vs_ssa(folder, ssa[valid], k22[valid], xssa, ylab, f'$k_{{yy}}$ vs SSA{title_extra}', th, dpi=args.dpi)
            savefig(fig, folder, 'k22_vs_ssa', th, args.dpi)

    # 5) k_11 and k_22 vs time (same axes) — plot what's available
    has11 = np.isfinite(k11).any()
    has22 = np.isfinite(k22).any()
    if has11 or has22:
        valid = np.isfinite(t_scaled)
        # Build series lists conditionally
        t_for_plot = t_scaled[valid]
        y1 = k11[valid] if has11 else None
        y2 = k22[valid] if has22 else None
        for th in themes:
            set_theme(th)
            fig, ax = plt.subplots(figsize=(11.5, 4.2))
            if has11:
                ax.plot(t_for_plot, y1, lw=2.0, marker='o', label=r'$k_{xx}$')
            if has22:
                ax.plot(t_for_plot, y2, lw=2.0, marker='s', label=r'$k_{yy}$')
            ax.set_xlabel(f"Time ({unit})")
            ax.set_ylabel(yk)
            ax.set_title(f'$k_{{xx}}$ & $k_{{yy}}$ vs Time{title_extra}')
            beautify(ax)
            ax.legend()
            fig.tight_layout()
            savefig(fig, folder, 'k11_k22_vs_time', th, args.dpi)

    # 6) k_11 and k_22 vs SSA (same axes) — plot what's available
    has11_ssa = (np.isfinite(ssa) & np.isfinite(k11)).any()
    has22_ssa = (np.isfinite(ssa) & np.isfinite(k22)).any()
    if has11_ssa or has22_ssa:
        valid = np.isfinite(ssa)
        x_for_plot = ssa[valid]
        y1 = k11[valid] if has11_ssa else None
        y2 = k22[valid] if has22_ssa else None
        for th in themes:
            set_theme(th)
            fig, ax = plt.subplots(figsize=(11.5, 4.2))
            if has11_ssa:
                ax.plot(x_for_plot, y1, lw=2.0, marker='o', label=r'$k_{xx}$')
            if has22_ssa:
                ax.plot(x_for_plot, y2, lw=2.0, marker='s', label=r'$k_{yy}$')
            ax.set_xlabel(xssa)
            ax.set_ylabel(yk)
            ax.set_title(f'$k_{{xx}}$ & $k_{{yy}}$ vs SSA{title_extra}')
            beautify(ax)
            ax.legend()
            fig.tight_layout()
            savefig(fig, folder, 'k11_k22_vs_ssa', th, args.dpi)

    print(f"Saved plots to: {os.path.join(folder, 'output_plots')}")
    print("Done.")

if __name__ == "__main__":
    main()
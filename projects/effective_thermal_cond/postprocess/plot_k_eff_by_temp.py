#!/usr/bin/env python3
"""
Plot k_eff (k_xx and k_yy) vs time as LINE PLOTS across multiple simulation runs
in a parent directory. Each subfolder corresponds (typically) to a different temperature.

This version reads BOTH temperature (°C) and porosity from metadata.json when available.

Data read per run:
  - k_eff.csv
      columns: prefers ['k_00','k_11'] (k_xx, k_yy). Accepts fallbacks like k_xx/k_yy/k11/k22.
      time: uses 'time' or 't' if present; else 'sol_index' mapped via SSA_evo.dat; else index.
  - SSA_evo.dat (optional)
      for mapping sol_index -> time when k_eff.csv lacks time,
      and for plotting normalized conductivity vs. normalized SSA.
      expected columns: [SSA, (optional), time, step]
  - metadata.json (optional but preferred)
      temperature: prefers dsm_run.parameters.temperature_C; otherwise any key containing 'temp'
      porosity: prefers structure.porosity_target, else structure.porosity_achieved

Output:
  - Creates <parent_dir>/output_plots/ (or <parent_dir>/log_output_plots/ when --log_plots is enabled)
  - Saves 2 figures (x2 if --theme=both):
      lines_kxx_vs_time_by_temp_[light|dark].png
      lines_kyy_vs_time_by_temp_[light|dark].png
      lines_kxx_vs_ssa_by_temp_[light|dark].png
      lines_kyy_vs_ssa_by_temp_[light|dark].png

Usage:
  python postprocess/plot_k_eff_by_temp.py PARENT_DIR
     [--pattern REGEX]        # filter subfolder basenames with regex
     [--normalize]            # divide each series by its first value
     [--theme {light,dark,both}]  # default: both
     [--dpi 300]

Example:
  python postprocess/plot_k_eff_by_temp.py \
      ~/SimulationResults/effective_thermal_cond/const_porosity_varying_temp \
      --pattern '^DSM_.*$' --normalize --theme both
"""

from __future__ import annotations

import os, re, sys, json, argparse
from typing import Dict, Any, List, Optional, Tuple
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ------------------------ Global plot styling ------------------------ #
# Target figure size (inches)
FIG_W_IN = 6.0
FIG_H_IN = 3.5

# Font sizes (tuned for 6"x3.5" figures)
TITLE_FONTSIZE = 16
AXIS_FONTSIZE = 18
TICK_FONTSIZE = 10
LEGEND_FONTSIZE = 10

# A small negative-time pad so axes don't look "cut off" at t=0

TIME_AXIS_PAD_FRACTION = 0.02

# ------------------------ Axis label wrapping helpers ------------------------ #
# Keep axis labels compact so they fit within the figure bounds
MAX_XLABEL_CHARS = 34
MAX_YLABEL_CHARS = 26

def _wrap_label(text: str, max_chars: int) -> str:
    """Soft-wrap a label by inserting a newline at the last space before max_chars."""
    if text is None:
        return ""
    s = str(text)
    if len(s) <= max_chars or "\n" in s:
        return s
    # Try to break on the last space before max_chars
    cut = s.rfind(" ", 0, max_chars + 1)
    if cut == -1:
        return s
    return s[:cut].rstrip() + "\n" + s[cut + 1 :].lstrip()

def set_compact_axis_labels(ax: plt.Axes, xlabel: str, ylabel: str):
    """Set labels with wrapping so they stay within the figure bounds."""
    xl = _wrap_label(xlabel, MAX_XLABEL_CHARS)
    yl = _wrap_label(ylabel, MAX_YLABEL_CHARS)
    ax.set_xlabel(xl)
    ax.set_ylabel(yl)
    # Enable wrapping for backends that respect it
    try:
        ax.xaxis.label.set_wrap(True)
        ax.yaxis.label.set_wrap(True)
    except Exception:
        pass

# ------------------------ Small helpers ------------------------ #

def autoscale_time(arr: np.ndarray) -> Tuple[np.ndarray, str]:
    """Scale seconds to minutes/hours/days for nicer axes; returns (scaled, unit)."""
    if arr.size == 0 or not np.isfinite(arr).any():
        return arr, "seconds"
    tmax = float(np.nanmax(arr))
    if tmax >= 2 * 24 * 3600:
        return arr / 86400.0, "days"
    if tmax >= 2 * 3600:
        return arr / 3600.0, "hours"
    if tmax >= 120:
        return arr / 60.0, "minutes"
    return arr, "seconds"

# ------------------------ Normalization & formatting helpers ------------------------ #

def _is_normalized_series(arr: np.ndarray) -> bool:
    """Return True if arr is non-empty, finite at index 0, and abs(arr[0] - 1.0) < 1e-6."""
    if arr is None or not hasattr(arr, "size") or arr.size == 0:
        return False
    if not np.isfinite(arr[0]):
        return False
    return abs(float(arr[0]) - 1.0) < 1e-6

def _fmt_temp_phi(tempC: Optional[float], por: Optional[float]) -> str:
    """Return string like 'T=-20.0 °C, φ=0.24' using unicode phi. Only include finite values."""
    bits = []
    if tempC is not None and np.isfinite(tempC):
        bits.append(f"T={tempC:.1f} °C")
    if por is not None and np.isfinite(por):
        bits.append(f"φ={por:.2f}")
    return ", ".join(bits)

def set_theme(theme: str, use_cmocean_thermal: bool = False, n_colors: int = 10):
    """Apply a consistent Matplotlib theme.

    Requirements:
      - No grid
      - No background (transparent)
      - Title 28pt, axes labels 24pt
      - Optional cmocean thermal color scheme
    """
    plt.rcParams.update(plt.rcParamsDefault)

    # Transparent backgrounds (no plot/figure background)
    plt.rcParams.update({
        'figure.facecolor': 'none',
        'axes.facecolor': 'none',
        'savefig.facecolor': 'none',
        'savefig.transparent': True,
    })

    # Foreground colors (light vs dark) for text/spines/ticks
    if theme == "dark":
        fg = 'white'
        plt.rcParams.update({'text.color': fg, 'axes.labelcolor': fg, 'axes.edgecolor': fg, 'xtick.color': fg, 'ytick.color': fg})
    else:
        fg = 'black'
        plt.rcParams.update({'text.color': fg, 'axes.labelcolor': fg, 'axes.edgecolor': fg, 'xtick.color': fg, 'ytick.color': fg})

    # Typography + axes
    plt.rcParams.update({
        'font.family': 'DejaVu Sans',
        'axes.titlesize': TITLE_FONTSIZE,
        'axes.labelsize': AXIS_FONTSIZE,
        'xtick.labelsize': TICK_FONTSIZE,
        'ytick.labelsize': TICK_FONTSIZE,
        'legend.fontsize': LEGEND_FONTSIZE,
        'axes.linewidth': 1.2,
        'axes.grid': False,
        'xtick.direction': 'in',
        'ytick.direction': 'in',
        'xtick.major.size': 6,
        'ytick.major.size': 6,
        'xtick.major.width': 1.2,
        'ytick.major.width': 1.2,
    })

    # Note: we intentionally do NOT set a color cycle here when using cmocean.
    # We instead map each run's temperature to a color spanning the full colormap
    # range (coldest -> warmest) inside the plotting functions.
    if use_cmocean_thermal:
        try:
            import cmocean  # noqa: F401
        except Exception as e:
            print(f"[WARN] Could not import cmocean ({e}); falling back to Matplotlib defaults.")
# ------------------------ Temperature colormap helpers ------------------------ #

def get_temp_colormap(use_cmocean_thermal: bool):
    """Return a colormap for temperature-based coloring."""
    if use_cmocean_thermal:
        try:
            import cmocean
            return cmocean.cm.thermal
        except Exception:
            # Fallback that exists in Matplotlib
            return plt.get_cmap("plasma")
    return None


def temp_to_color(tempC: float, tmin: float, tmax: float, cmap) -> Optional[Tuple[float, float, float, float]]:
    """Map temperature to a color spanning the full colormap range.

    coldest (tmin) -> cmap(0.0)
    warmest (tmax) -> cmap(1.0)
    """
    if cmap is None or tempC is None:
        return None
    if not np.isfinite(tempC) or not np.isfinite(tmin) or not np.isfinite(tmax):
        return None
    if tmax == tmin:
        u = 0.5
    else:
        u = (float(tempC) - float(tmin)) / (float(tmax) - float(tmin))
    u = float(np.clip(u, 0.0, 1.0))
    return cmap(u)

# ------------------------ φ summarization helper ------------------------ #
def summarize_phi_for_title(runs: List[Dict[str, Any]]) -> Optional[str]:
    """
    Return a short text like 'φ=0.24' for use in a figure title if porosity is
    effectively constant across runs. If porosity is absent, returns None.
    If porosity varies noticeably across runs, returns 'φ varies'.
    """
    phis = []
    for d in runs:
        p = d.get("porosity")
        if p is not None and np.isfinite(p):
            phis.append(float(p))
    if not phis:
        return None
    if (max(phis) - min(phis)) < 5e-4:  # ~constant within tolerance
        return f"φ={phis[0]:.2f}"
    return "φ varies"

def parse_temp_from_name(name: str) -> Optional[float]:
    """Extract temperature (°C) from folder name like ..._Tm-25_... or ..._Tm-25C_..."""
    m = re.search(r'_Tm(-?\d+(?:\.\d+)?)C?(_|$)', name)
    if m:
        try:
            return float(m.group(1))
        except Exception:
            return None
    return None

def read_json(path: str) -> Optional[Dict[str, Any]]:
    try:
        with open(path, 'r') as f:
            return json.load(f)
    except Exception:
        return None

def _walk_items(obj, prefix=""):
    if isinstance(obj, dict):
        for k, v in obj.items():
            p = f"{prefix}.{k}" if prefix else k
            if isinstance(v, (dict, list)):
                yield from _walk_items(v, p)
            else:
                yield (p, k, v)
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            p = f"{prefix}[{i}]"
            if isinstance(v, (dict, list)):
                yield from _walk_items(v, p)
            else:
                yield (p, str(i), v)

def _coerce_float(v) -> Optional[float]:
    try:
        return float(v)
    except Exception:
        return None

def find_temp_from_meta(meta: Dict[str, Any]) -> Optional[float]:
    """Prefer dsm_run.parameters.temperature_C; else any key containing 'temp' parsed as float."""
    # Priority exact path if present
    try:
        p = meta.get("dsm_run", {}).get("parameters", {}).get("temperature_C", None)
        if p is not None:
            val = _coerce_float(p)
            if val is not None:
                return val
    except Exception:
        pass
    # Common alternative exact keys
    for k in ("temperature_C", "tempC", "T_C", "TdegC"):
        if k in meta:
            val = _coerce_float(meta[k])
            if val is not None:
                return val
    # Any key containing temp
    for _, k, v in _walk_items(meta):
        if 'temp' in k.lower():
            val = _coerce_float(v)
            if val is not None:
                return val
    return None

def find_porosity_from_meta(meta: Dict[str, Any]) -> Optional[float]:
    """Prefer structure.porosity_target; else structure.porosity_achieved; else any key containing 'poros' or 'phi'."""
    try:
        s = meta.get("structure", {})
        for key in ("porosity_target", "porosity_achieved"):
            if key in s:
                val = _coerce_float(s[key])
                if val is not None:
                    return val
    except Exception:
        pass
    # Any key containing 'poros' or 'phi'
    for _, k, v in _walk_items(meta):
        kl = k.lower()
        if "poros" in kl or kl == "phi" or kl.endswith(".phi"):
            val = _coerce_float(v)
            if val is not None:
                return val
    return None

def read_ssa_map(ssa_path: str) -> Tuple[Optional[np.ndarray], Optional[np.ndarray], Optional[np.ndarray]]:
    """Read SSA_evo.dat; returns (steps, times, ssa_vals) or (None,None,None).
    Expected columns: [SSA, (optional), time, step]
    """
    try:
        data = np.loadtxt(ssa_path)
        if data.ndim != 2 or data.shape[1] < 4:
            return None, None, None
        ssa_vals = data[:, 0]
        times = data[:, 2]
        steps = data[:, 3].astype(int)
        return steps, times, ssa_vals
    except Exception:
        return None, None, None

# ------------------------ IO per run ------------------------ #

def load_run(folder: str, normalize: bool) -> Optional[Dict[str, Any]]:
    """
    Load one run folder and return:
      dict(name, tempC, porosity, time, kxx, kyy, ssa)
    """
    name = os.path.basename(folder.rstrip(os.sep))

    # ---- metadata.json: temperature & porosity ----
    meta = None
    for mname in ("metadata.json", "Metadata.json", "meta.json"):
        mp = os.path.join(folder, mname)
        if os.path.exists(mp):
            meta = read_json(mp)
            if meta:
                break

    tempC = None
    porosity = None
    if meta:
        tempC = find_temp_from_meta(meta)
        porosity = find_porosity_from_meta(meta)

    # Fallback temperature from folder name if metadata missing
    if tempC is None:
        tempC = parse_temp_from_name(name)
    if tempC is None:
        print(f"[WARN] Could not determine temperature for '{name}', skipping.")
        return None

    # ---- k_eff.csv ----
    kcsv = os.path.join(folder, "k_eff.csv")
    if not os.path.exists(kcsv):
        print(f"[WARN] Missing k_eff.csv in '{name}', skipping.")
        return None

    try:
        df = pd.read_csv(kcsv)
    except Exception as e:
        print(f"[WARN] Failed to read k_eff.csv in '{name}': {e}")
        return None

    def pick_col(cands):
        for c in cands:
            if c in df.columns:
                return df[c].to_numpy()
        return None

    kxx = pick_col(["k_00", "k_xx", "k11", "k_11_xx"])
    kyy = pick_col(["k_11", "k_yy", "k22", "k_22_yy"])
    if kxx is None or kyy is None:
        print(f"[WARN] k columns not found in '{name}', columns present: {list(df.columns)}")
        return None

    # ---- Time: time or sol_index->SSA map or index; SSA alignment ----
    steps_k = pick_col(["sol_index", "step", "iter", "n"])
    time_k = pick_col(["time", "t"])
    ssa_aligned = None
    ssa_path = os.path.join(folder, "SSA_evo.dat")
    ssa_steps, ssa_times, ssa_vals = read_ssa_map(ssa_path)
    time = None
    if time_k is None and steps_k is not None:
        if ssa_steps is not None and ssa_times is not None:
            step2time = {int(s): float(t) for s, t in zip(ssa_steps, ssa_times)}
            time = np.array([step2time.get(int(s), np.nan) for s in steps_k], dtype=float)
        else:
            print(f"[INFO] No SSA_evo.dat mapping for '{name}', using step index as time (s).")
            time = np.arange(len(kxx), dtype=float)
    elif time_k is not None:
        time = np.asarray(time_k, dtype=float)
    else:
        print(f"[INFO] No time/step columns in k_eff.csv for '{name}', using index as time (s).")
        time = np.arange(len(kxx), dtype=float)

    # SSA alignment
    if steps_k is not None and ssa_steps is not None and ssa_vals is not None:
        step2ssa = {int(s): float(v) for s, v in zip(ssa_steps, ssa_vals)}
        ssa_aligned = np.array([step2ssa.get(int(s), np.nan) for s in steps_k], dtype=float)
    elif time is not None and ssa_times is not None and ssa_vals is not None:
        # Interpolate SSA to time
        ssa_times_sorted = np.asarray(ssa_times)
        ssa_vals_sorted = np.asarray(ssa_vals)
        order_ssa = np.argsort(ssa_times_sorted)
        ssa_times_sorted = ssa_times_sorted[order_ssa]
        ssa_vals_sorted = ssa_vals_sorted[order_ssa]
        ssa_aligned = np.interp(time, ssa_times_sorted, ssa_vals_sorted, left=np.nan, right=np.nan)

    time = np.asarray(time, dtype=float)
    kxx = np.asarray(kxx, dtype=float)
    kyy = np.asarray(kyy, dtype=float)

    valid = np.isfinite(time) & np.isfinite(kxx) & np.isfinite(kyy)
    if not np.any(valid):
        print(f"[WARN] No finite data in '{name}', skipping.")
        return None
    time, kxx, kyy = time[valid], kxx[valid], kyy[valid]
    # Also filter ssa_aligned if present
    if ssa_aligned is not None:
        ssa_aligned = np.asarray(ssa_aligned, dtype=float)
        if ssa_aligned.shape[0] == valid.shape[0]:
            ssa_aligned = ssa_aligned[valid]
        else:
            # fallback: ignore mismatch
            ssa_aligned = None

    # Sort by time to avoid zig-zags
    order = np.argsort(time)
    time, kxx, kyy = time[order], kxx[order], kyy[order]
    if ssa_aligned is not None:
        ssa_aligned = ssa_aligned[order]

    if normalize:
        if kxx.size > 0 and kxx[0] != 0:
            kxx = kxx / kxx[0]
        if kyy.size > 0 and kyy[0] != 0:
            kyy = kyy / kyy[0]

    return dict(
        name=name,
        tempC=float(tempC),
        porosity=(None if porosity is None else float(porosity)),
        time=time,
        kxx=kxx,
        kyy=kyy,
        ssa=ssa_aligned,
    )
# ------------------------ Plotting ------------------------ #
def plot_lines_k_vs_ssa(runs: List[Dict[str, Any]],
                        which: str,
                        y_label_norm: str,
                        outdir: str,
                        base_name: str,
                        theme: str,
                        dpi: int,
                        use_cmocean_thermal: bool = False,
                        log_plots: bool = False):
    """
    Plot normalized k (k/k0) vs normalized SSA (SSA/SSA0) for each run, lines colored by temperature.
    Requires SSA_evo.dat to be present (ssa in run dict). Runs lacking SSA are skipped.
    """
    def do_plot(theme_name: str, suffix: str):
        set_theme(theme_name, use_cmocean_thermal=use_cmocean_thermal)
        fig, ax = plt.subplots(figsize=(FIG_W_IN, FIG_H_IN))
        runs_sorted = sorted(runs, key=lambda d: d["tempC"])
        temps = [float(d["tempC"]) for d in runs_sorted if d.get("tempC") is not None and np.isfinite(d.get("tempC"))]
        tmin = float(np.nanmin(temps)) if temps else np.nan
        tmax = float(np.nanmax(temps)) if temps else np.nan
        cmap = get_temp_colormap(use_cmocean_thermal)
        any_series = False
        if log_plots:
            ax.set_xscale('log')
            ax.set_yscale('log')
        for d in runs_sorted:
            ssa = d.get("ssa")
            y = d.get(which)
            if ssa is None or y is None:
                continue
            ssa = np.asarray(ssa, dtype=float)
            y = np.asarray(y, dtype=float)
            valid = np.isfinite(ssa) & np.isfinite(y)
            if log_plots:
                valid = valid & (ssa > 0) & (y > 0)
            if np.count_nonzero(valid) < 2:
                continue
            ssa = ssa[valid]
            y = y[valid]
            # Normalize both
            if ssa[0] == 0 or y[0] == 0:
                continue
            ssa_n = ssa / ssa[0]
            y_n = y / y[0]
            if log_plots:
                mpos = (ssa_n > 0) & (y_n > 0)
                ssa_n = ssa_n[mpos]
                y_n = y_n[mpos]
                if ssa_n.size < 2:
                    continue
            lab = f"{d['tempC']:.1f} °C"
            color = temp_to_color(d["tempC"], tmin, tmax, cmap)
            ax.plot(ssa_n, y_n, lw=2.0, label=lab, color=color)
            any_series = True
        set_compact_axis_labels(
            ax,
            "Normalized Specific Surface Area",
            "Normalized Effective Thermal Conductivity",
        )
        # Title: two lines (no component in title)
        title = "Normalized Effective Thermal Conductivity\nvs. Normalized Specific Surface Area"
        ax.set_title(title)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(False)
        if any_series:
            if len(runs) > 6:
                # Many entries: keep legend outside
                fig.legend(loc="center right", bbox_to_anchor=(1.02, 0.5))
                fig.subplots_adjust(right=0.84)
                fig.tight_layout(rect=[0.0, 0.0, 0.84, 1.0])
            else:
                # Use Matplotlib's "best" placement and a slightly smaller legend font
                ax.legend(loc="best", fontsize=max(8, LEGEND_FONTSIZE - 1))
                fig.tight_layout()
        out_path = os.path.join(outdir, base_name.replace('.png', f'_{suffix}.png'))
        fig.savefig(out_path, dpi=dpi, bbox_inches='tight', transparent=True)
        plt.close(fig)
        print("[OK]", out_path)

    if theme in ("light", "both"):
        do_plot("light", "light")
    if theme in ("dark", "both"):
        do_plot("dark", "dark")

# ------------------------ Discovery & scaling ------------------------ #

def discover_folders(parent: str, pattern: Optional[str]) -> List[str]:
    """List immediate subfolders (optionally filter by regex pattern).
    Skips utility/hidden folders such as 'output_plots' and names starting with '.'
    """
    try:
        entries = sorted(os.listdir(parent))
    except FileNotFoundError:
        return []

    # Exclude hidden folders and our plot output folder by default
    def _is_valid_dir(name: str) -> bool:
        if name.startswith('.'):  # hidden
            return False
        if name in {"output_plots", "plots", "figures"}:
            return False
        return True

    subs = [os.path.join(parent, f) for f in entries
            if _is_valid_dir(f) and os.path.isdir(os.path.join(parent, f))]

    if pattern:
        rx = re.compile(pattern)
        subs = [p for p in subs if rx.match(os.path.basename(p))]

    return subs

def choose_global_time_unit(runs: List[Dict[str, Any]]) -> str:
    """Pick a single time unit across all runs for consistent x-axes."""
    all_t = np.concatenate([r["time"] for r in runs if isinstance(r.get("time"), np.ndarray) and r["time"].size > 0])
    if all_t.size == 0:
        return "seconds"
    _, unit = autoscale_time(all_t)
    return unit

def scale_to_unit(arr: np.ndarray, unit: str) -> np.ndarray:
    if unit == "days": return arr / 86400.0
    if unit == "hours": return arr / 3600.0
    if unit == "minutes": return arr / 60.0
    return arr

# ------------------------ Plotting ------------------------ #

def plot_lines_vs_time(runs: List[Dict[str, Any]],
                       which: str,
                       y_label: str,
                       outdir: str,
                       base_name: str,
                       theme: str,
                       dpi: int,
                       use_cmocean_thermal: bool = False,
                       log_plots: bool = False):
    """
    Make a line plot of 'which' in {'kxx','kyy'} vs time, with one line per temperature.
    Legend includes temperature only.
    """
    unit = choose_global_time_unit(runs)

    def do_plot(theme_name: str, suffix: str):
        set_theme(theme_name, use_cmocean_thermal=use_cmocean_thermal)
        fig, ax = plt.subplots(figsize=(FIG_W_IN, FIG_H_IN))

        # Sort by temperature for legend order
        runs_sorted = sorted(runs, key=lambda d: d["tempC"])
        temps = [float(d["tempC"]) for d in runs_sorted if d.get("tempC") is not None and np.isfinite(d.get("tempC"))]
        tmin = float(np.nanmin(temps)) if temps else np.nan
        tmax = float(np.nanmax(temps)) if temps else np.nan
        cmap = get_temp_colormap(use_cmocean_thermal)

        # Decide normalization for y-axis label and title (use first run with data)
        normalized = False
        if runs_sorted and which in runs_sorted[0]:
            arr = runs_sorted[0][which]
            if isinstance(arr, np.ndarray) and _is_normalized_series(arr):
                normalized = True

        # Y-axis label, compact with line break
        if normalized:
            ylab = "Normalized\nEffective Thermal Conductivity"
        else:
            ylab = "Effective\nThermal Conductivity"

        # Title: two lines (no component in title)
        norm_prefix = "Normalized " if normalized else ""
        title = f"{norm_prefix}Effective Thermal Conductivity\nvs. Time"

        if log_plots:
            ax.set_xscale('log')
            ax.set_yscale('log')

        for d in runs_sorted:
            t = scale_to_unit(d["time"], unit)
            y = d[which]
            if t.size == 0 or y.size == 0:
                continue
            lab = f"{d['tempC']:.1f} °C"
            color = temp_to_color(d["tempC"], tmin, tmax, cmap)
            if log_plots:
                m = np.isfinite(t) & np.isfinite(y) & (t > 0) & (y > 0)
                t_plot = t[m]
                y_plot = y[m]
            else:
                t_plot = t
                y_plot = y
            if t_plot.size < 2:
                continue
            ax.plot(t_plot, y_plot, lw=2.0, label=lab, color=color)

        set_compact_axis_labels(ax, f"Time ({unit})", ylab)
        ax.set_title(title)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(False)
        # Extend x-axis slightly before 0 so it doesn't look cut off (skip for log-x)
        if not log_plots:
            if t.size > 0 and np.isfinite(t).any():
                xmin, xmax = np.nanmin(t), np.nanmax(t)
                pad = TIME_AXIS_PAD_FRACTION * (xmax - xmin if xmax > xmin else 1.0)
                ax.set_xlim(xmin - pad, xmax)

        # Place legend
        if len(runs) > 6:
            # Many entries: keep legend outside
            fig.legend(loc="center right", bbox_to_anchor=(1.02, 0.5))
            fig.subplots_adjust(right=0.84)
            fig.tight_layout(rect=[0.0, 0.0, 0.84, 1.0])
        else:
            # Use Matplotlib's "best" placement and a slightly smaller legend font
            ax.legend(loc="best", fontsize=max(8, LEGEND_FONTSIZE - 1))
            fig.tight_layout()

        out_path = os.path.join(outdir, base_name.replace(".png", f"_{suffix}.png"))
        fig.savefig(out_path, dpi=dpi, bbox_inches="tight", transparent=True)
        plt.close(fig)
        print("[OK]", out_path)

    if theme in ("light", "both"):
        do_plot("light", "light")
    if theme in ("dark", "both"):
        do_plot("dark", "dark")


# ------------------------ Per-run plotting helper ------------------------ #
def plot_individual_runs(runs: List[Dict[str, Any]],
                         outdir: str,
                         theme: str,
                         dpi: int,
                         use_cmocean_thermal: bool = False,
                         log_plots: bool = False):
    """Save per-run plots (overlaying kxx and kyy vs time) into output_plots/<run_name>/.

    For each run:
      - Creates output_plots/<run_name>/
      - Saves k_eff_vs_time_[light|dark].png (overlaying kxx and kyy)
    """

    def _safe_name(s: str) -> str:
        # Keep folder names stable and filesystem-friendly
        return re.sub(r"[^A-Za-z0-9._-]+", "_", s).strip("_") or "run"

    def _plot_both(run: Dict[str, Any], unit: str, theme_name: str, out_path: str):
        # Always use default color cycle for individual overlays, do not use cmocean thermal
        set_theme(theme_name, use_cmocean_thermal=False)
        fig, ax = plt.subplots(figsize=(FIG_W_IN, FIG_H_IN))
        t = scale_to_unit(run["time"], unit)
        kxx = np.asarray(run.get("kxx", np.array([])), dtype=float)
        kyy = np.asarray(run.get("kyy", np.array([])), dtype=float)
        if log_plots:
            ax.set_xscale('log')
            ax.set_yscale('log')
        # Plot both kxx and kyy, filtering nonpositive for log if needed
        if log_plots:
            mxx = np.isfinite(t) & np.isfinite(kxx) & (t > 0) & (kxx > 0)
            myy = np.isfinite(t) & np.isfinite(kyy) & (t > 0) & (kyy > 0)
            ax.plot(t[mxx], kxx[mxx], lw=2.2, label=r"$k_{xx}$")
            ax.plot(t[myy], kyy[myy], lw=2.2, label=r"$k_{yy}$")
        else:
            ax.plot(t, kxx, lw=2.2, label=r"$k_{xx}$")
            ax.plot(t, kyy, lw=2.2, label=r"$k_{yy}$")
        # Decide y-label based on normalization (use kxx as heuristic)
        if _is_normalized_series(kxx):
            ylab = "Normalized\nEffective Thermal Conductivity"
            main_title = "Normalized Effective Thermal Conductivity\nvs. Time"
        else:
            ylab = "Effective\nThermal Conductivity"
            main_title = "Effective Thermal Conductivity\nvs. Time"
        set_compact_axis_labels(ax, f"Time ({unit})", ylab)
        # Subtitle: run metadata (temp, phi)
        tempC = run.get("tempC", None)
        por = run.get("porosity", None)
        subtitle = _fmt_temp_phi(tempC, por)
        if subtitle:
            ax.set_title(f"{main_title}\n{subtitle}")
        else:
            ax.set_title(main_title)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(False)
        # Extend x-axis slightly before 0 so it doesn't look cut off (skip for log-x)
        if not log_plots:
            if t.size > 0 and np.isfinite(t).any():
                xmin, xmax = np.nanmin(t), np.nanmax(t)
                pad = TIME_AXIS_PAD_FRACTION * (xmax - xmin if xmax > xmin else 1.0)
                ax.set_xlim(xmin - pad, xmax)
        ax.legend(loc="best", fontsize=max(8, LEGEND_FONTSIZE - 1))
        fig.tight_layout()
        fig.savefig(out_path, dpi=dpi, bbox_inches='tight', transparent=True)
        plt.close(fig)

    def _plot_both_vs_ssa(run: Dict[str, Any], theme_name: str, out_path: str):
        """Plot normalized kxx and kyy vs normalized SSA for a single run."""
        # Always use default color cycle for individual overlays, do not use cmocean thermal
        set_theme(theme_name, use_cmocean_thermal=False)
        fig, ax = plt.subplots(figsize=(FIG_W_IN, FIG_H_IN))

        ssa = run.get("ssa")
        if ssa is None:
            plt.close(fig)
            return

        ssa = np.asarray(ssa, dtype=float)
        kxx = np.asarray(run.get("kxx", np.array([])), dtype=float)
        kyy = np.asarray(run.get("kyy", np.array([])), dtype=float)

        valid = np.isfinite(ssa) & np.isfinite(kxx) & np.isfinite(kyy)
        if np.count_nonzero(valid) < 2:
            plt.close(fig)
            return

        ssa = ssa[valid]
        kxx = kxx[valid]
        kyy = kyy[valid]

        # Normalize both axes (matches multi-temperature SSA plots)
        if ssa.size == 0 or ssa[0] == 0 or kxx.size == 0 or kxx[0] == 0 or kyy.size == 0 or kyy[0] == 0:
            plt.close(fig)
            return

        ssa_n = ssa / ssa[0]
        kxx_n = kxx / kxx[0]
        kyy_n = kyy / kyy[0]

        if log_plots:
            ax.set_xscale('log')
            ax.set_yscale('log')
            mxx = np.isfinite(ssa_n) & np.isfinite(kxx_n) & (ssa_n > 0) & (kxx_n > 0)
            myy = np.isfinite(ssa_n) & np.isfinite(kyy_n) & (ssa_n > 0) & (kyy_n > 0)
            if np.count_nonzero(mxx) >= 2:
                ax.plot(ssa_n[mxx], kxx_n[mxx], lw=2.2, label=r"$k_{xx}$")
            if np.count_nonzero(myy) >= 2:
                ax.plot(ssa_n[myy], kyy_n[myy], lw=2.2, label=r"$k_{yy}$")
        else:
            ax.plot(ssa_n, kxx_n, lw=2.2, label=r"$k_{xx}$")
            ax.plot(ssa_n, kyy_n, lw=2.2, label=r"$k_{yy}$")

        set_compact_axis_labels(
            ax,
            "Normalized Specific Surface Area",
            "Normalized Effective Thermal Conductivity",
        )

        main_title = "Normalized Effective Thermal Conductivity\nvs. Normalized Specific Surface Area"
        tempC = run.get("tempC", None)
        por = run.get("porosity", None)
        subtitle = _fmt_temp_phi(tempC, por)
        if subtitle:
            ax.set_title(f"{main_title}\n{subtitle}")
        else:
            ax.set_title(main_title)

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(False)
        ax.legend(loc="best", fontsize=max(8, LEGEND_FONTSIZE - 1))
        fig.tight_layout()
        fig.savefig(out_path, dpi=dpi, bbox_inches='tight', transparent=True)
        plt.close(fig)

    for run in runs:
        run_dir = os.path.join(outdir, _safe_name(run.get("name", "run")))
        os.makedirs(run_dir, exist_ok=True)
        _, unit = autoscale_time(run["time"])

        def _emit(theme_name: str, suffix: str):
            out_path = os.path.join(run_dir, f"k_eff_vs_time_{suffix}.png")
            _plot_both(run, unit, theme_name, out_path)
            print("[OK]", out_path)

            # Also save normalized k vs normalized SSA, if SSA is available
            out_path_ssa = os.path.join(run_dir, f"k_eff_vs_ssa_norm_{suffix}.png")
            _plot_both_vs_ssa(run, theme_name, out_path_ssa)
            if os.path.exists(out_path_ssa):
                print("[OK]", out_path_ssa)

        if theme in ("light", "both"):
            _emit("light", "light")
        if theme in ("dark", "both"):
            _emit("dark", "dark")

# ------------------------ Main ------------------------ #

def main():
    ap = argparse.ArgumentParser(description="Line plots of k_eff vs time by temperature (reads temp & porosity from metadata.json when available).")
    ap.add_argument("parent_dir", help="Directory containing run subfolders.")
    ap.add_argument("--pattern", default=None,
                    help="Regex to filter subfolder names. Match is applied to the folder basename.")
    ap.add_argument("--normalize", action="store_true",
                    help="Normalize k_xx and k_yy by their initial values.")
    ap.add_argument("--theme", choices=["light", "dark", "both"], default="both",
                    help="Color theme for saved plots.")
    ap.add_argument("--dpi", type=int, default=300, help="Output DPI.")
    ap.add_argument("--use_thermal_cmocean", action="store_true",
                    help="Use cmocean 'thermal' colormap as the line color cycle (requires cmocean).")
    ap.add_argument("--log_plots", action="store_true",
                    help="Save log-scale versions of all plots (log y-axis). Output goes to log_output_plots.")
    args = ap.parse_args()

    parent = os.path.abspath(os.path.expanduser(args.parent_dir))
    if not os.path.isdir(parent):
        print(f"[ERROR] Parent directory not found: {parent}")
        sys.exit(2)

    subfolders = discover_folders(parent, args.pattern)
    if not subfolders:
        print(f"[ERROR] No subfolders found in {parent}. Use --pattern to filter if needed.")
        sys.exit(1)

    runs: List[Dict[str, Any]] = []
    for sub in subfolders:
        r = load_run(sub, normalize=args.normalize)
        if r is not None:
            runs.append(r)

    if not runs:
        print(f"[ERROR] No valid runs in {parent}. Ensure each folder has k_eff.csv and (optionally) SSA_evo.dat & metadata.json.")
        sys.exit(1)

    outdir_name = "log_output_plots" if args.log_plots else "output_plots"
    outdir = os.path.join(parent, outdir_name)
    os.makedirs(outdir, exist_ok=True)

    plot_lines_vs_time(runs, "kxx", "", outdir, "lines_kxx_vs_time_by_temp.png", args.theme, args.dpi,
                       use_cmocean_thermal=args.use_thermal_cmocean, log_plots=args.log_plots)
    plot_lines_vs_time(runs, "kyy", "", outdir, "lines_kyy_vs_time_by_temp.png", args.theme, args.dpi,
                       use_cmocean_thermal=args.use_thermal_cmocean, log_plots=args.log_plots)

    # Per-run plots (saved into output_plots/<run_name>/)
    plot_individual_runs(runs, outdir, args.theme, args.dpi,
                         use_cmocean_thermal=args.use_thermal_cmocean, log_plots=args.log_plots)

    # Also plot normalized k vs normalized SSA (always normalized)
    ylab_xx_norm = r"$k_{xx}/k_{xx,0}$"
    ylab_yy_norm = r"$k_{yy}/k_{yy,0}$"
    plot_lines_k_vs_ssa(runs, "kxx", ylab_xx_norm, outdir, "lines_kxx_vs_ssa_by_temp.png", args.theme, args.dpi,
                        use_cmocean_thermal=args.use_thermal_cmocean, log_plots=args.log_plots)
    plot_lines_k_vs_ssa(runs, "kyy", ylab_yy_norm, outdir, "lines_kyy_vs_ssa_by_temp.png", args.theme, args.dpi,
                        use_cmocean_thermal=args.use_thermal_cmocean, log_plots=args.log_plots)

    print("\n✅ Done.")

if __name__ == "__main__":
    main()
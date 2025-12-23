#!/usr/bin/env python3
"""
Produce conference-ready curve plots of effective thermal conductivity vs time and vs (normalized) SSA
for varying porosity runs. This script loads simulation result folders from a parent directory,
extracts k_eff and SSA evolution, and generates high-quality plots suitable for publication or conference use.

Each subdirectory should contain:
  - k_eff.csv (columns: k_00/k_xx, k_11/k_yy, sol_index/step)
  - SSA_evo.dat (cols: SSA, ?, time, step)  [used to map step -> time]
  - metadata.json (optional; may contain porosity fields)

Porosity (phi) is taken from the folder name if possible, else from metadata.json.

Outputs:
  <outdir>/curves_kxx_vs_time[_normalized][_loglog][_light|_dark].png
  <outdir>/curves_kyy_vs_time[_normalized][_loglog][_light|_dark].png
  <outdir>/curves_kxx_norm_vs_SSA_normalized[_loglog][_light|_dark].png
  <outdir>/curves_kyy_norm_vs_SSA_normalized[_loglog][_light|_dark].png
  <outdir>/individual_plots/<run>/k_eff_vs_time[_normalized][_loglog][_light|_dark].png
  <outdir>/individual_plots/<run>/k_eff_vs_norm_SSA_normalized[_loglog][_light|_dark].png

Usage:
  python postprocess/plot_k_eff_by_phi.py <PARENT_DIR> \
    [--pattern 'REGEX'] [--outdir PATH] [--normalize true|false] [--dpi 300] [--theme light|dark|both]

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

# ----------------- Style constants (conference-ready) -----------------
FIG_W_IN = 6.0
FIG_H_IN = 3.5
TITLE_FS = 18
LABEL_FS = 16
TICK_FS  = 14
LEGEND_FS = 13
LINE_W = 2.5
SAVE_DPI = 300

# ----------------- Label wrapping helpers -----------------
def _wrap_label(text: str) -> str:
    """
    Wrap known long axis labels onto two lines using \n.
    """
    mapping = {
        # k_eff labels
        "Normalized Effective Thermal Conductivity": "Normalized Effective\nThermal Conductivity",
        "Effective Thermal Conductivity": "Effective\nThermal Conductivity",
        "Effective Thermal Conductivity (W m$^{-1}$ K$^{-1}$)": "Effective Thermal Conductivity\n(W m$^{-1}$ K$^{-1}$)",

        # SSA labels
        "Normalized Specific Surface Area": "Normalized Specific\nSurface Area",
        "Specific Surface Area": "Specific\nSurface Area",
        "Specific Surface Area (m$^{-1}$)": "Specific Surface Area\n(m$^{-1}$)",
    }
    return mapping.get(text, text)


def _apply_label_wrap(ax):
    ax.title.set_wrap(True)
    ax.xaxis.label.set_wrap(True)
    ax.yaxis.label.set_wrap(True)

# ----------------- Axis-label fitting helpers -----------------
def _soft_wrap_at_space(text: str) -> str:
    """If `text` has no newline, insert one near the middle at a space."""
    if "\n" in text:
        return text
    s = text.strip()
    if len(s) < 12:
        return text
    # Find a space closest to the midpoint
    mid = len(s) // 2
    left = s.rfind(" ", 0, mid)
    right = s.find(" ", mid)
    if left == -1 and right == -1:
        return text
    if left == -1:
        cut = right
    elif right == -1:
        cut = left
    else:
        cut = left if (mid - left) <= (right - mid) else right
    if cut <= 0 or cut >= len(s) - 1:
        return text
    return s[:cut].rstrip() + "\n" + s[cut + 1 :].lstrip()


def _ensure_axis_labels_fit(fig, ax, *, max_passes: int = 6):
    """
    Ensure axis labels/titles fit within axes/figure by soft-wrapping at spaces.
    Only wraps text, does not adjust subplot paddings.
    """
    def exceeds(ax, txtobj):
        """Check if txtobj's window extent exceeds the axes in width or height."""
        if not txtobj or not hasattr(txtobj, "get_window_extent"):
            return False
        renderer = fig.canvas.get_renderer()
        try:
            bb = txtobj.get_window_extent(renderer=renderer)
            axbb = ax.get_window_extent(renderer=renderer)
            # Allow a few px of margin
            w_exceeds = bb.width > axbb.width - 4
            h_exceeds = bb.height > axbb.height - 4
            return w_exceeds or h_exceeds
        except Exception:
            return False

    def exceeds_fig(fig, txtobj):
        if not txtobj or not hasattr(txtobj, "get_window_extent"):
            return False
        renderer = fig.canvas.get_renderer()
        try:
            bb = txtobj.get_window_extent(renderer=renderer)
            figbb = fig.bbox
            w_exceeds = bb.width > figbb.width - 4
            h_exceeds = bb.height > figbb.height - 4
            return w_exceeds or h_exceeds
        except Exception:
            return False

    for _ in range(max_passes):
        changed = False
        fig.canvas.draw()
        # X label
        xlabel_obj = ax.xaxis.label
        xlab = xlabel_obj.get_text()
        if exceeds(ax, xlabel_obj) and ("\n" not in xlab):
            newlab = _soft_wrap_at_space(xlab)
            if newlab != xlab:
                ax.set_xlabel(newlab)
                changed = True
        # Y label
        ylabel_obj = ax.yaxis.label
        ylab = ylabel_obj.get_text()
        if exceeds(ax, ylabel_obj) and ("\n" not in ylab):
            newlab = _soft_wrap_at_space(ylab)
            if newlab != ylab:
                ax.set_ylabel(newlab)
                changed = True
        # Subtitle (axes title)
        title_obj = ax.title
        tit = title_obj.get_text()
        if exceeds(ax, title_obj) and ("\n" not in tit):
            newtit = _soft_wrap_at_space(tit)
            if newtit != tit:
                ax.set_title(newtit)
                changed = True
        # Suptitle (figure)
        suptitle_obj = getattr(fig, "_suptitle", None)
        if suptitle_obj is not None:
            stext = suptitle_obj.get_text()
            if exceeds_fig(fig, suptitle_obj) and ("\n" not in stext):
                newstext = _soft_wrap_at_space(stext)
                if newstext != stext:
                    fig.suptitle(newstext, fontsize=suptitle_obj.get_fontsize(), y=suptitle_obj.get_position()[1])
                    changed = True
        if not changed:
            break


# ----------------- Title and layout helpers -----------------
def _set_titles(fig, ax, title: str, subtitle: str | None = None):
    """Set a consistent (multi-line safe) title style.

    We use fig.suptitle for the main title and (optionally) ax.set_title for a subtitle,
    matching the style used in the temperature-sweep script.
    """
    # Clear any existing titles
    ax.set_title("")
    # Main title
    st = fig.suptitle(title, fontsize=TITLE_FS, y=0.98)
    try:
        st.set_wrap(True)
    except Exception:
        pass
    # Optional subtitle
    if subtitle:
        t2 = ax.set_title(subtitle, fontsize=max(10, LABEL_FS - 2), pad=4)
        try:
            t2.set_wrap(True)
        except Exception:
            pass


def _finalize_layout(fig, ax, *, top: float = 0.90, max_passes: int = 8):
    fig.tight_layout(rect=[0.0, 0.0, 1.0, top])

#
# ----------------- Legend placement helper -----------------
def _place_legend_no_overlap(fig, ax, *, prefer: str = "best", title: str | None = None, legend_mode: int = 1):
    """
    Place legend so it does not overlap plotted curves by expanding y-limits downward if needed.
    Never moves legend outside axes or changes subplot paddings.
    """
    if legend_mode == 2:
        # Do not draw legend on plot
        return
    font_main = max(9, LEGEND_FS - 2)
    font_small = max(8, LEGEND_FS - 4)
    for attempt in range(8):
        leg = ax.legend(loc=prefer, fontsize=font_main, frameon=True, title=title)
        if leg is None:
            return
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        leg_bb = leg.get_window_extent(renderer=renderer)
        # Get all line artists (curves)
        lines = ax.get_lines()
        overlap = False
        for line in lines:
            try:
                line_bb = line.get_window_extent(renderer=renderer)
                if leg_bb.overlaps(line_bb):
                    overlap = True
                    break
            except Exception:
                continue
        if not overlap:
            return
        # If overlap, expand y-limits downward (decrease lower limit)
        y0, y1 = ax.get_ylim()
        if ax.get_yscale() == "log":
            # Downward expansion by ~8% in log space
            # Keep y0 > 0 (required for log scale)
            y0_pos = max(float(y0), np.finfo(float).tiny)
            dy = np.log10(float(y1)) - np.log10(y0_pos)
            # Move the lower bound down by 8% of the current log-range
            y0_new = 10 ** (np.log10(y0_pos) - dy * 0.08)
            ax.set_ylim(y0_new, y1)
        else:
            # Downward expansion by 8% of the current range
            rng = float(y1) - float(y0)
            if rng <= 0:
                rng = abs(float(y1)) if float(y1) != 0 else 1.0
            y0_new = float(y0) - 0.08 * rng
            ax.set_ylim(y0_new, y1)
        leg.remove()
    # Fallback: smaller font and "best" location
    leg = ax.legend(loc="best", fontsize=font_small, frameon=True, title=title)
    fig.canvas.draw()


# --- Save legend as a separate PNG ---
def _save_legend_only(fig, ax, out_path, title=None):
    """Save the legend (handles/labels) as a standalone PNG file (transparent), no axes.

    NOTE: We intentionally avoid passing a tuple to `bbox_inches` (matplotlib expects a
    Bbox-like object or the strings 'tight'/'standard'). We instead draw the legend on a
    tiny figure and save with `bbox_inches='tight'`.
    """
    handles, labels = ax.get_legend_handles_labels()
    if not handles or not labels:
        return

    import matplotlib.pyplot as plt

    # Make a small figure that only contains the legend.
    tmp_fig = plt.figure(figsize=(2.0, 2.0))
    tmp_ax = tmp_fig.add_axes([0, 0, 1, 1])
    tmp_ax.axis("off")

    leg = tmp_ax.legend(
        handles,
        labels,
        fontsize=LEGEND_FS,
        frameon=True,
        title=title,
        loc="center",
    )

    # Draw once so legend has a computed size, then save tightly around it.
    tmp_fig.canvas.draw()

    tmp_fig.savefig(
        out_path,
        dpi=SAVE_DPI,
        bbox_inches="tight",
        pad_inches=0.02,
        transparent=True,
    )
    plt.close(tmp_fig)

# --- Stable discrete color mapping for porosity values ---
def build_phi_color_lookup(phis: np.ndarray):
    """
    Return a function color_for(phi) that assigns a stable, discrete color
    to each unique porosity (rounded to 3 decimals), with fixed mapping for
    common φ values, else fallback to discrete palette.
    """
    phis = np.asarray(phis, dtype=float)
    if phis.size == 0:
        def _fallback(_):
            return plt.get_cmap("viridis")(0.5)
        return _fallback
    # Fixed mapping for common porosities (rounded to 2 decimals)
    fixed_colors = {
        0.24: "#34B11E",
        0.26: "#9080A0",
        0.28: "#002060",
        0.30: "#FE5500",
    }
    # Round to stabilize keys, then unique + sorted for deterministic ordering
    keys = np.unique(np.round(phis, 3))
    # Build a long discrete palette (tab20 + tab20b + tab20c)
    palette = list(plt.cm.tab20.colors) + list(plt.cm.tab20b.colors) + list(plt.cm.tab20c.colors)
    if len(palette) < len(keys):
        extra = plt.get_cmap("viridis")(np.linspace(0.05, 0.95, len(keys) - len(palette)))
        palette.extend([tuple(c) for c in extra])
    mapping = {k: palette[i % len(palette)] for i, k in enumerate(keys)}
    def color_for(phi_value: float):
        k2 = float(np.round(phi_value, 2))
        if k2 in fixed_colors:
            return fixed_colors[k2]
        k3 = float(np.round(phi_value, 3))
        return mapping.get(k3, plt.get_cmap("viridis")(0.5))
    return color_for



def set_theme(theme: str):
    """Apply a light or dark Matplotlib theme (conference-ready, transparent, small fonts, no grid)."""
    import matplotlib as mpl
    plt.rcParams.update(plt.rcParamsDefault)
    # Always use transparent background for figure, axes, and savefig
    bg_none = 'none'
    if theme == "dark":
        plt.rcParams.update({
            'figure.facecolor': bg_none,
            'axes.facecolor': bg_none,
            'savefig.facecolor': bg_none,
            'text.color': 'white',
            'axes.labelcolor': 'white',
            'axes.edgecolor': 'white',
            'xtick.color': 'white',
            'ytick.color': 'white',
            'axes.grid': False,
            'grid.color': '0.35',
        })
    else:  # light
        plt.rcParams.update({
            'figure.facecolor': bg_none,
            'axes.facecolor': bg_none,
            'savefig.facecolor': bg_none,
            'text.color': 'black',
            'axes.labelcolor': 'black',
            'axes.edgecolor': 'black',
            'xtick.color': 'black',
            'ytick.color': 'black',
            'axes.grid': False,
            'grid.color': '0.7',
        })
    plt.rcParams.update({
        'font.family': 'DejaVu Sans',
        'axes.titlesize': TITLE_FS,
        'axes.labelsize': LABEL_FS,
        'xtick.labelsize': TICK_FS,
        'ytick.labelsize': TICK_FS,
        'legend.fontsize': LEGEND_FS,
        'axes.linewidth': 1.2,
        'xtick.direction': 'in', 'ytick.direction': 'in',
        'xtick.major.size': 6,   'ytick.major.size': 6,
        'xtick.major.width': 1.2,'ytick.major.width': 1.2,
    })
def extract_temp_from_name(name: str) -> Optional[float]:
    """
    Extract temperature from folder name patterns like 'Tm-20', 'Tm-05', 'T-05.0', 'T-5.0'.
    Returns float in Celsius, or None if not found.
    """
    # Patterns: Tm-20, Tm-05, T-05.0, T-5.0, T-20, etc.
    m = re.search(r"Tm-?(\d+)", name)
    if m:
        return -float(m.group(1))
    m = re.search(r"T-(-?\d+(?:\.\d+)?)", name)
    if m:
        return float(m.group(1))
    return None


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





def main():
    ap = argparse.ArgumentParser(description="Produce conference-ready curve plots of effective thermal conductivity vs time and vs (normalized) SSA for varying porosity runs.")
    ap.add_argument("parent_dir", help="Folder containing run subdirectories")
    ap.add_argument("--pattern", help="Python regex to filter subfolder names", default=None)
    ap.add_argument("--outdir", help="Where to save plots (default: <parent_dir>/output_plots)", default=None)
    ap.add_argument("--normalize", type=bool_arg, default=True, help="Normalize k by first value (default: true)")
    ap.add_argument("--dpi", type=int, default=300)
    ap.add_argument("--theme", type=str, default="light", choices=["light", "dark", "both"],
                    help="Color theme: 'light' (default), 'dark', or 'both' to export both.")
    ap.add_argument("--loglog", action="store_true", help="Also save log-log versions of the curve plots to <outdir>/../log_output_plots (or <parent_dir>/log_output_plots).")
    ap.add_argument("--legend_mode", type=int, choices=[1, 2], default=1,
                    help="Legend rendering mode: 1=draw on plot (default), 2=save legend as separate PNG")
    args = ap.parse_args()

    parent_dir = os.path.abspath(os.path.expanduser(args.parent_dir))
    outdir = os.path.abspath(os.path.expanduser(args.outdir)) if args.outdir else os.path.join(parent_dir, "output_plots")
    os.makedirs(outdir, exist_ok=True)

    # Set up log-log output directory if needed
    if not args.outdir:
        log_outdir = os.path.join(os.path.dirname(outdir), "log_output_plots")
    else:
        log_outdir = os.path.join(os.path.dirname(outdir), "log_output_plots")
    if args.loglog:
        os.makedirs(log_outdir, exist_ok=True)

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
        print(f"[WARN] Found only {len(runs)} usable run(s). Need ≥2 for meaningful plots.")
        return

    # Sort runs by porosity
    runs.sort(key=lambda d: d["phi"])
    phis = np.array([d["phi"] for d in runs], dtype=float)

    # Build stable color lookup for porosity
    color_for_phi = build_phi_color_lookup(phis)

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

    def plot_time_curves(runs_list, ylabel, key, fname, theme, out_base_dir, use_loglog=False):
        set_theme(theme)
        fig, ax = plt.subplots(figsize=(FIG_W_IN, FIG_H_IN))
        # Color by porosity for consistency (discrete, stable)
        for d in runs_list:
            col = color_for_phi(d["phi"])
            t_plot = time_scale(d["time"])
            y = d[key]
            if use_loglog:
                mask = np.isfinite(t_plot) & np.isfinite(y) & (t_plot > 0) & (y > 0)
                if np.count_nonzero(mask) < 2:
                    continue
                t_plot = t_plot[mask]
                y = y[mask]
            # No markers for combined curves
            ax.plot(
                t_plot, y, lw=LINE_W, alpha=0.98, label=fr"$\phi={d['phi']:.2f}$", color=col
            )
        ax.set_xlabel(f"Time ({unit_str})", fontsize=LABEL_FS)
        ax.set_ylabel(_wrap_label(ylabel), fontsize=LABEL_FS)
        # Title and legend location
        ylab_for_title = _wrap_label(ylabel).replace('\n', ' ')
        if "time" in fname:
            _set_titles(fig, ax, title=f"{ylab_for_title} vs. Time")
            legend_loc = 'lower right'
        else:
            _set_titles(fig, ax, title=f"{ylab_for_title} vs. Specific Surface Area")
            legend_loc = 'upper right'
        if use_loglog:
            ax.set_xscale('log')
            ax.set_yscale('log')
        _apply_label_wrap(ax)
        _ensure_axis_labels_fit(fig, ax)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        # Legend handling
        if args.legend_mode == 1:
            _place_legend_no_overlap(fig, ax, prefer=legend_loc, title="Porosity", legend_mode=args.legend_mode)
        else:
            # Remove any legend, save legend separately after plot is saved
            pass
        _finalize_layout(fig, ax, top=0.88)
        base, ext = os.path.splitext(fname)
        suffix = "" if len(themes) == 1 else f"_{theme}"
        loglog_suffix = "_loglog" if use_loglog else ""
        # Add _normalized to base name if normalized plot
        normalized = args.normalize or ("Normalized" in ylabel)
        if normalized and "_normalized" not in base:
            base += "_normalized"
        out = os.path.join(out_base_dir, f"{base}{loglog_suffix}{suffix}{ext}")
        fig.savefig(out, dpi=SAVE_DPI, bbox_inches="tight", transparent=True)
        # Save legend separately if legend_mode==2
        if args.legend_mode == 2:
            legend_out = os.path.join(out_base_dir, f"{base}{loglog_suffix}{suffix}_legend.png")
            _save_legend_only(fig, ax, legend_out, title="Porosity")
        plt.close(fig)
        print("[OK]", out)

    def plot_norm_vs_ssa_curves(runs_list, ylabel, key, fname, theme, out_base_dir, use_loglog=False):
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
                    if use_loglog:
                        mask = (ssa_n > 0) & (k_n > 0)
                        if np.count_nonzero(mask) < 2:
                            continue
                        ssa_n = ssa_n[mask]
                        k_n = k_n[mask]
                    valid_runs.append((d["phi"], ssa_n, k_n))

        if not valid_runs:
            print("[WARN] No runs with usable SSA for normalized SSA plots; skipping", fname)
            return

        set_theme(theme)
        fig, ax = plt.subplots(figsize=(FIG_W_IN, FIG_H_IN))
        for phi, ssa_n, k_n in valid_runs:
            col = color_for_phi(phi)
            # No markers for combined curves
            ax.plot(
                ssa_n, k_n, lw=LINE_W, alpha=0.98, label=fr"$\phi={phi:.2f}$", color=col
            )
        # Axis labels and title
        ax.set_xlabel(_wrap_label("Normalized Specific Surface Area"), fontsize=LABEL_FS)
        ax.set_ylabel(_wrap_label(ylabel), fontsize=LABEL_FS)
        ylab_for_title = _wrap_label(ylabel).replace('\n', ' ')
        _set_titles(fig, ax, title=f"{ylab_for_title} vs. Normalized Specific Surface Area")
        if use_loglog:
            ax.set_xscale('log')
            ax.set_yscale('log')
        _apply_label_wrap(ax)
        _ensure_axis_labels_fit(fig, ax)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        legend_loc = 'upper right'
        if args.legend_mode == 1:
            _place_legend_no_overlap(fig, ax, prefer=legend_loc, title="Porosity", legend_mode=args.legend_mode)
        else:
            pass
        _finalize_layout(fig, ax, top=0.88)
        base, ext = os.path.splitext(fname)
        suffix = "" if len(themes) == 1 else f"_{theme}"
        loglog_suffix = "_loglog" if use_loglog else ""
        # Always append _normalized for these plots
        if not base.endswith("_normalized"):
            base += "_normalized"
        out = os.path.join(out_base_dir, f"{base}{loglog_suffix}{suffix}{ext}")
        fig.savefig(out, dpi=SAVE_DPI, bbox_inches="tight", transparent=True)
        if args.legend_mode == 2:
            legend_out = os.path.join(out_base_dir, f"{base}{loglog_suffix}{suffix}_legend.png")
            _save_legend_only(fig, ax, legend_out, title="Porosity")
        plt.close(fig)
        print("[OK]", out)

    # --------- Individual-run (per-run) conference-ready plots ---------
    import shutil
    from matplotlib.lines import Line2D
    def plot_individual_run(run, theme, normalize, out_base_dir, use_loglog=False):
        # Make output dir
        run_dir = os.path.join(out_base_dir, "individual_plots", run["name"])
        os.makedirs(run_dir, exist_ok=True)
        set_theme(theme)
        color = color_for_phi(run["phi"])
        # --- k_xx and k_yy vs time (on same axes) ---
        fig, ax = plt.subplots(figsize=(FIG_W_IN, FIG_H_IN))
        t_plot = time_scale(run["time"])
        kxx = run["kxx"]
        kyy = run["kyy"]
        # Line styles for kxx and kyy
        if use_loglog:
            mask = np.isfinite(t_plot) & np.isfinite(kxx) & np.isfinite(kyy) & (t_plot > 0) & (kxx > 0) & (kyy > 0)
            if np.count_nonzero(mask) < 2:
                plt.close(fig)
                return
            t_plot2 = t_plot[mask]
            kxx2 = kxx[mask]
            kyy2 = kyy[mask]
        else:
            t_plot2 = t_plot
            kxx2 = kxx
            kyy2 = kyy
        ax.plot(t_plot2, kxx2, lw=LINE_W, color=color, label=r"$k_{xx}$", linestyle='-')
        ax.plot(t_plot2, kyy2, lw=LINE_W, color=color, label=r"$k_{yy}$", linestyle='--')
        # Labels
        ax.set_xlabel(f"Time ({unit_str})", fontsize=LABEL_FS)
        if normalize:
            ylab = "Normalized Effective Thermal Conductivity"
        else:
            ylab = "Effective Thermal Conductivity"
        ax.set_ylabel(_wrap_label(ylab), fontsize=LABEL_FS)
        # Title/subtitle (consistent with combined plots)
        main_title = f"{_wrap_label(ylab).replace(chr(10),' ')} vs. Time"

        temp = extract_temp_from_name(run["name"])
        subtitle = ""
        if temp is not None:
            subtitle += f"T={temp:.0f}°C"
        if run.get("phi", None) is not None:
            if subtitle:
                subtitle += ", "
            subtitle += f"φ={run['phi']:.2f}"

        _set_titles(fig, ax, title=main_title, subtitle=subtitle if subtitle else None)
        if use_loglog:
            ax.set_xscale('log')
            ax.set_yscale('log')
        _apply_label_wrap(ax)
        _ensure_axis_labels_fit(fig, ax)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        if args.legend_mode == 1:
            _place_legend_no_overlap(fig, ax, prefer='lower right', legend_mode=args.legend_mode)
        else:
            pass
        _finalize_layout(fig, ax, top=0.86)
        loglog_suffix = "_loglog" if use_loglog else ""
        # Add _normalized to filename if normalized
        norm_suffix = "_normalized" if normalize else ""
        out1 = os.path.join(run_dir, f"k_eff_vs_time{norm_suffix}{loglog_suffix}{'_'+theme if len(themes)>1 else ''}.png")
        fig.savefig(out1, dpi=SAVE_DPI, bbox_inches="tight", transparent=True)
        if args.legend_mode == 2:
            legend_out = os.path.join(run_dir, f"k_eff_vs_time{norm_suffix}{loglog_suffix}{'_'+theme if len(themes)>1 else ''}_legend.png")
            _save_legend_only(fig, ax, legend_out, title=r"$k_{xx}$/$k_{yy}$")
        plt.close(fig)

        # --- k_xx and k_yy vs normalized SSA (if available) ---
        ssa = run.get("ssa", None)
        if ssa is not None and np.isfinite(ssa).sum() >= 2:
            finite = np.isfinite(ssa) & np.isfinite(run["kxx"]) & np.isfinite(run["kyy"])
            if finite.sum() >= 2:
                ssa_n = ssa[finite] / ssa[finite][0] if ssa[finite][0] != 0 else ssa[finite]
                kxx_n = run["kxx"][finite]
                kyy_n = run["kyy"][finite]
                if use_loglog:
                    mask = (ssa_n > 0) & (kxx_n > 0) & (kyy_n > 0)
                    if np.count_nonzero(mask) < 2:
                        return
                    ssa_n2 = ssa_n[mask]
                    kxx_n2 = kxx_n[mask]
                    kyy_n2 = kyy_n[mask]
                else:
                    ssa_n2 = ssa_n
                    kxx_n2 = kxx_n
                    kyy_n2 = kyy_n
                fig2, ax2 = plt.subplots(figsize=(FIG_W_IN, FIG_H_IN))
                ax2.plot(ssa_n2, kxx_n2, lw=LINE_W, color=color, label=r"$k_{xx}$", linestyle='-')
                ax2.plot(ssa_n2, kyy_n2, lw=LINE_W, color=color, label=r"$k_{yy}$", linestyle='--')
                ax2.set_xlabel(_wrap_label("Normalized Specific Surface Area"), fontsize=LABEL_FS)
                if normalize:
                    ylab2 = "Normalized Effective Thermal Conductivity"
                else:
                    ylab2 = "Effective Thermal Conductivity"
                ax2.set_ylabel(_wrap_label(ylab2), fontsize=LABEL_FS)
                main_title2 = f"{_wrap_label(ylab2).replace(chr(10),' ')} vs. Normalized Specific Surface Area"
                _set_titles(fig2, ax2, title=main_title2, subtitle=subtitle if subtitle else None)
                if use_loglog:
                    ax2.set_xscale('log')
                    ax2.set_yscale('log')
                _apply_label_wrap(ax2)
                _ensure_axis_labels_fit(fig2, ax2)
                ax2.spines["top"].set_visible(False)
                ax2.spines["right"].set_visible(False)
                if args.legend_mode == 1:
                    _place_legend_no_overlap(fig2, ax2, prefer='upper right', legend_mode=args.legend_mode)
                else:
                    pass
                _finalize_layout(fig2, ax2, top=0.86)
                loglog_suffix = "_loglog" if use_loglog else ""
                # Always append _normalized to filename for these plots
                out2 = os.path.join(run_dir, f"k_eff_vs_norm_SSA_normalized{loglog_suffix}{'_'+theme if len(themes)>1 else ''}.png")
                fig2.savefig(out2, dpi=SAVE_DPI, bbox_inches="tight", transparent=True)
                if args.legend_mode == 2:
                    legend_out = os.path.join(run_dir, f"k_eff_vs_norm_SSA_normalized{loglog_suffix}{'_'+theme if len(themes)>1 else ''}_legend.png")
                    _save_legend_only(fig2, ax2, legend_out, title=r"$k_{xx}$/$k_{yy}$")
                plt.close(fig2)

    # Generate the new curve plots per theme
    for theme in themes:
        # Y labels and titles
        if args.normalize:
            ylab_kxx = "Normalized Effective Thermal Conductivity"
            ylab_kyy = "Normalized Effective Thermal Conductivity"
        else:
            ylab_kxx = "Effective Thermal Conductivity (W m$^{-1}$ K$^{-1}$)"
            ylab_kyy = "Effective Thermal Conductivity (W m$^{-1}$ K$^{-1}$)"

        # Combined curves (all runs, colored by porosity)
        plot_time_curves(runs, ylab_kxx, key="kxx", fname="curves_kxx_vs_time.png", theme=theme, out_base_dir=outdir, use_loglog=False)
        plot_time_curves(runs, ylab_kyy, key="kyy", fname="curves_kyy_vs_time.png", theme=theme, out_base_dir=outdir, use_loglog=False)
        plot_norm_vs_ssa_curves(runs, ylab_kxx, key="kxx", fname="curves_kxx_norm_vs_SSA.png", theme=theme, out_base_dir=outdir, use_loglog=False)
        plot_norm_vs_ssa_curves(runs, ylab_kyy, key="kyy", fname="curves_kyy_norm_vs_SSA.png", theme=theme, out_base_dir=outdir, use_loglog=False)

        # Individual-run plots (kxx and kyy on same axes)
        for run in runs:
            plot_individual_run(run, theme, args.normalize, out_base_dir=outdir, use_loglog=False)

    # Log-log versions if requested
    if args.loglog:
        for theme in themes:
            if args.normalize:
                ylab_kxx = "Normalized Effective Thermal Conductivity"
                ylab_kyy = "Normalized Effective Thermal Conductivity"
            else:
                ylab_kxx = "Effective Thermal Conductivity (W m$^{-1}$ K$^{-1}$)"
                ylab_kyy = "Effective Thermal Conductivity (W m$^{-1}$ K$^{-1}$)"
            plot_time_curves(runs, ylab_kxx, key="kxx", fname="curves_kxx_vs_time.png", theme=theme, out_base_dir=log_outdir, use_loglog=True)
            plot_time_curves(runs, ylab_kyy, key="kyy", fname="curves_kyy_vs_time.png", theme=theme, out_base_dir=log_outdir, use_loglog=True)
            plot_norm_vs_ssa_curves(runs, ylab_kxx, key="kxx", fname="curves_kxx_norm_vs_SSA.png", theme=theme, out_base_dir=log_outdir, use_loglog=True)
            plot_norm_vs_ssa_curves(runs, ylab_kyy, key="kyy", fname="curves_kyy_norm_vs_SSA.png", theme=theme, out_base_dir=log_outdir, use_loglog=True)
            for run in runs:
                plot_individual_run(run, theme, args.normalize, out_base_dir=log_outdir, use_loglog=True)


if __name__ == "__main__":
    main()
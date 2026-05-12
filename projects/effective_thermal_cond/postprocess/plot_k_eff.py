#!/usr/bin/env python3
"""
Plot effective thermal conductivity for a SINGLE simulation run.

Expected files in the run directory:
  - k_eff.csv
      Columns:
        time or t (optional)
        k_00 / k_xx
        k_11 / k_yy
        sol_index or step (optional)

  - SSA_evo.dat (optional)
      Columns:
        SSA, (unused), time, step

  - metadata.json (optional)
      Used for temperature and porosity annotation.

Outputs (written inside the run directory):
  - single_run_plots_<format>/
      - k_xx_vs_time_<theme>.<format>
      - k_yy_vs_time_<theme>.<format>
      - log_k_xx_vs_time_<theme>.<format>          (if --log_plots)
      - log_k_yy_vs_time_<theme>.<format>          (if --log_plots)
      - k_xx_vs_ssa_norm_<theme>.<format>          (if SSA_evo.dat exists)
      - k_yy_vs_ssa_norm_<theme>.<format>          (if SSA_evo.dat exists)
      - log_k_xx_vs_ssa_norm_<theme>.<format>      (if --log_plots and SSA exists)
      - log_k_yy_vs_ssa_norm_<theme>.<format>      (if --log_plots and SSA exists)

Usage:
  python postprocess/plot_k_eff.py <run_dir> [options]

Required positional arguments:
  run_dir
      Path to a single simulation run directory containing k_eff.csv

Options:
  --theme {light,dark,both}
      Color theme for plots. Default: both

  --dpi <int>
      DPI for saved figures (only meaningful for raster outputs like PNG).
      Default: 300

  --format {svg,pdf,png}
      Output file format.
        - svg (default): vector, editable in Inkscape
        - pdf           : vector, journal-friendly
        - png           : raster

  --log_plots
      Generate log-log versions of the plots.

  --no_legend
      Disable legends on plots.

  --no_labels
      Remove axis labels and titles (keeps ticks and tick labels).

  --mark_time_points
      Draw marker points along k vs time curves at approximately evenly spaced times.
      The same marker points are also marked on the k vs SSA plots (if SSA is available).

  --n_markers <int>
      Number of marker points to draw when --mark_time_points is enabled.
      Default: 4

  --fit_powerlaw
      Fit a power-law curve of the form y = C * t^a to the k vs time data (log-log least squares)
      and overlay the best-fit curve on the time plots. Also prints the fitted exponent a.
      (Requires positive time and positive y values.)

  --dry_run
      Run the script and report what would be generated, but do not write any output files.

Examples:
  # Vector SVG output (recommended for Inkscape editing)
  python postprocess/plot_k_eff.py /path/to/run --theme both --format svg

  # PDF output with log-log plots
  python postprocess/plot_k_eff.py /path/to/run --theme light --format pdf --log_plots

  # Minimal figure styling for panel assembly
  python postprocess/plot_k_eff.py /path/to/run --no_labels --no_legend --format svg

  # Mark the 4 time points used for microstructure snapshots
  python postprocess/plot_k_eff.py /path/to/run --mark_time_points --format svg

  # Use 6 marker points and dry-run mode (no files written)
  python postprocess/plot_k_eff.py /path/to/run --mark_time_points --n_markers 6 --dry_run
"""

import os, json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse, sys
from typing import Optional

# ------------------- Styling ------------------- #
FIG_W, FIG_H = 6.0, 3.5
TITLE_FONTSIZE = 16
AXIS_FONTSIZE = 18
TICK_FONTSIZE = 10
LEGEND_FONTSIZE = 10
TIME_AXIS_PAD_FRACTION = 0.02

# Side-specific padding fractions (used when you want extra room for insets)
# Linear: additive padding as a fraction of data range
# Log: multiplicative padding in decades as a fraction of log-range
Y_PAD_FRAC_LINEAR = 0.45
Y_PAD_FRAC_LOG    = 0.45

# Minimum padding on ALL sides so plots never look clipped (applies to x and y)
PAD_MIN_FRAC = 0.03

MAX_XLABEL_CHARS = 34
MAX_YLABEL_CHARS = 26

def set_theme(theme):
    plt.rcParams.update(plt.rcParamsDefault)
    fg = "black" if theme == "light" else "white"
    plt.rcParams.update({
        "figure.facecolor": "none",
        "axes.facecolor": "none",
        "savefig.facecolor": "none",
        "savefig.transparent": True,
        # Keep text as text in vector outputs (do NOT convert glyphs to paths)
        "svg.fonttype": "none",
        # Embed TrueType fonts in PDF/PS so text remains editable in Inkscape
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "text.color": fg,
        "axes.labelcolor": fg,
        "axes.edgecolor": fg,
        "xtick.color": fg,
        "ytick.color": fg,
        "font.family": "DejaVu Sans",
        "axes.titlesize": TITLE_FONTSIZE,
        "axes.labelsize": AXIS_FONTSIZE,
        "xtick.labelsize": TICK_FONTSIZE,
        "ytick.labelsize": TICK_FONTSIZE,
        "legend.fontsize": LEGEND_FONTSIZE,
        "axes.linewidth": 1.2,
        "axes.grid": False,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.major.size": 6,
        "ytick.major.size": 6,
        "xtick.major.width": 1.2,
        "ytick.major.width": 1.2,
        "xtick.minor.size": 3,
        "ytick.minor.size": 3,
        "xtick.minor.width": 1.2,
        "ytick.minor.width": 1.2,
    })

def _wrap_label(text, max_chars):
    if len(text) <= max_chars:
        return text
    words = text.split()
    lines = []
    current_line = ""
    for w in words:
        if len(current_line) + len(w) + 1 <= max_chars:
            if current_line:
                current_line += " "
            current_line += w
        else:
            lines.append(current_line)
            current_line = w
    if current_line:
        lines.append(current_line)
    return "\n".join(lines)

def set_compact_axis_labels(ax, xlabel, ylabel):
    xlabel_wrapped = _wrap_label(xlabel, MAX_XLABEL_CHARS)
    ylabel_wrapped = _wrap_label(ylabel, MAX_YLABEL_CHARS)
    ax.set_xlabel(xlabel_wrapped)
    ax.set_ylabel(ylabel_wrapped)
    # Try to enable wrapping if supported
    for label in (ax.xaxis.label, ax.yaxis.label):
        try:
            label.set_wrap(True)
        except Exception:
            pass

def autoscale_time(arr):
    arr = np.asarray(arr, float)
    if arr.size == 0:
        return arr, "s"
    max_t = np.nanmax(arr)
    if max_t < 60:
        return arr, "s"
    elif max_t < 3600:
        return arr / 60, "min"
    elif max_t < 86400:
        return arr / 3600, "hr"
    else:
        return arr / 86400, "d"

def scale_to_unit(arr, unit):
    arr = np.asarray(arr, float)
    if unit == "s":
        return arr
    elif unit == "min":
        return arr / 60
    elif unit == "hr":
        return arr / 3600
    elif unit == "d":
        return arr / 86400
    else:
        return arr

# ------------------- Helpers ------------------- #
def read_metadata(run_dir):
    meta_path = os.path.join(run_dir, "metadata.json")
    if not os.path.exists(meta_path):
        return None, None
    with open(meta_path) as f:
        meta = json.load(f)

    temp = meta.get("dsm_run", {}).get("parameters", {}).get("temperature_C")
    por  = meta.get("structure", {}).get("porosity_target")
    return temp, por

def read_ssa(run_dir):
    """Read SSA_evo.dat.

    We support multiple column conventions seen in your workflows.

    Primary convention (your current runs):
      - Column 0: SSA
      - Column 2 (3rd column): physical time (seconds)
      - Column 3 (4th column): step index

    Fallback behavior (robust):
      - Column 0 is SSA.
      - If there are >= 2 columns and the primary convention is not available,
        we try to infer which column is time and which is step:
          * time: the column that looks most monotone increasing (not necessarily integer)
          * step: the column that looks most integer-like.

    Returns:
      (steps, times, ssa_vals) where steps may be None if it cannot be inferred.
    """
    ssa_path = os.path.join(run_dir, "SSA_evo.dat")
    if not os.path.exists(ssa_path):
        return None, None, None

    try:
        data = np.loadtxt(ssa_path)
    except Exception:
        return None, None, None

    if data.ndim != 2 or data.shape[1] < 2:
        return None, None, None

    ssa_vals = np.asarray(data[:, 0], float)

    # --- Primary convention: [SSA, (unused), time, step] ---
    if data.shape[1] >= 4:
        times = np.asarray(data[:, 2], float)  # 3rd column is physical time
        try:
            steps = np.asarray(np.round(data[:, 3]), int)  # 4th column is step index
        except Exception:
            steps = None
        return steps, times, ssa_vals

    # --- Fallback / inference mode ---
    # If we do not have 4 columns, fall back to inference.
    times = None
    steps = None

    # Candidates exclude SSA column 0
    cand_cols = list(range(1, data.shape[1]))
    if not cand_cols:
        return None, None, None

    # Score monotonicity for time-like behavior
    def _mono_score(col: np.ndarray) -> float:
        finite = np.isfinite(col)
        if np.count_nonzero(finite) < max(2, int(0.5 * col.size)):
            return -np.inf
        c = col[finite]
        if c.size < 2:
            return -np.inf
        dif = np.diff(c)
        # fraction of nondecreasing steps
        frac_nondec = float(np.sum(dif >= 0)) / float(dif.size)
        # prefer columns that actually change
        span = float(np.nanmax(c) - np.nanmin(c))
        span_score = 0.0 if not np.isfinite(span) else np.tanh(span)
        return frac_nondec + 0.05 * span_score

    # Score integer-likeness for step-like behavior
    def _int_score(col: np.ndarray) -> float:
        finite = np.isfinite(col)
        if np.count_nonzero(finite) < max(2, int(0.5 * col.size)):
            return -np.inf
        c = col[finite]
        frac_part = np.abs(c - np.round(c))
        return -float(np.nanmedian(frac_part))

    # Pick time column
    best_time_col = None
    best_time_score = -np.inf
    for c in cand_cols:
        score = _mono_score(np.asarray(data[:, c], float))
        if score > best_time_score:
            best_time_score = score
            best_time_col = c

    if best_time_col is None:
        return None, None, None

    times = np.asarray(data[:, best_time_col], float)

    # Pick step column among remaining (if any)
    remaining = [c for c in cand_cols if c != best_time_col]
    best_step_col = None
    best_step_score = -np.inf
    for c in remaining:
        score = _int_score(np.asarray(data[:, c], float))
        if score > best_step_score:
            best_step_score = score
            best_step_col = c

    if best_step_col is not None:
        try:
            steps = np.asarray(np.round(data[:, best_step_col]), int)
        except Exception:
            steps = None

    return steps, times, ssa_vals

def print_usage_and_exit(msg):
    print(f"ERROR: {msg}\n")
    print("Required CLI arguments:")
    print("  run_dir            Path to the simulation run directory containing k_eff.csv\n")
    print("Optional files expected in run_dir:")
    print("  SSA_evo.dat")
    print("  metadata.json\n")
    print("Usage example:")
    print("  python plot_k_eff.py <run_dir> [--theme both] [--dpi 300] [--log_plots] [--no_legend] [--format svg] [--mark_time_points]")
    sys.exit(1)

# ------------------- Argument parsing ------------------- #
parser = argparse.ArgumentParser(
    description="Plot effective thermal conductivity for a single simulation run.\n"
                "Required files in the run directory:\n"
                "  - k_eff.csv\n"
                "Optional files:\n"
                "  - SSA_evo.dat\n"
                "  - metadata.json",
    add_help=True
)
parser.add_argument("run_dir", nargs='?', help="Path to the simulation run directory containing k_eff.csv")
parser.add_argument("--theme", choices=["light","dark","both"], default="both", help="Color theme for plots")
parser.add_argument("--dpi", type=int, default=300, help="DPI for saved figures")
parser.add_argument("--log_plots", action="store_true", help="Generate log-log plots")
parser.add_argument("--no_legend", action="store_true", help="Disable legends on plots")
parser.add_argument("--format", choices=["svg", "pdf", "png"], default="svg",
                    help="Output figure format. Use 'svg' (default) or 'pdf' for vector graphics; 'png' for raster.")
parser.add_argument("--mark_time_points", action="store_true",
                    help="Draw marker points at evenly spaced times along k vs time curves.")
parser.add_argument("--n_markers", type=int, default=4,
                    help="Number of marker points to draw when --mark_time_points is enabled (default: 4).")
parser.add_argument("--fit_powerlaw", action="store_true",
                    help="Fit y = C * t^a to k vs time (log-log least squares) and overlay the best-fit curve.")
parser.add_argument("--dry_run", action="store_true",
                    help="Do not write any output files; just report what would be generated.")
parser.add_argument("--no_labels", action="store_true",
                    help="Remove axis labels and titles (keep ticks and tick numbers).")
args = parser.parse_args()

run_dir = args.run_dir
theme_arg = args.theme
dpi = args.dpi
log_plots = args.log_plots
no_legend = args.no_legend
out_format = args.format

mark_time_points = args.mark_time_points
no_labels = args.no_labels
n_markers = int(args.n_markers)
dry_run = bool(args.dry_run)
fit_powerlaw = bool(args.fit_powerlaw)
# ------------------- Power-law fit helper ------------------- #
def fit_powerlaw_exponent(t: np.ndarray, y: np.ndarray):
    """Fit y = C * x^a via log-log least squares.

    Returns (a, C, r2). Requires t>0 and y>0.
    """
    t = np.asarray(t, float)
    y = np.asarray(y, float)
    mask = np.isfinite(t) & np.isfinite(y) & (t > 0) & (y > 0)
    if np.count_nonzero(mask) < 3:
        return None, None, None

    lt = np.log(t[mask])
    ly = np.log(y[mask])

    # Linear fit: ly = a*lt + b
    a, b = np.polyfit(lt, ly, 1)
    C = float(np.exp(b))

    # R^2 on log-space fit
    ly_hat = a * lt + b
    ss_res = float(np.sum((ly - ly_hat) ** 2))
    ss_tot = float(np.sum((ly - float(np.mean(ly))) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan

    return float(a), C, r2

# --- New helper: fit inverse power-law y = C / x^a ---
def fit_inverse_powerlaw_exponent(x: np.ndarray, y: np.ndarray):
    """Fit y = C / x^a via log-log least squares.

    Equivalent to: log(y) = log(C) - a*log(x)

    Returns (a, C, r2). Requires x>0 and y>0.
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    mask = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    if np.count_nonzero(mask) < 3:
        return None, None, None

    lx = np.log(x[mask])
    ly = np.log(y[mask])

    # Linear fit: ly = m*lx + b  where m = -a
    m, b = np.polyfit(lx, ly, 1)
    a = float(-m)
    C = float(np.exp(b))

    # R^2 on log-space fit
    ly_hat = m * lx + b
    ss_res = float(np.sum((ly - ly_hat) ** 2))
    ss_tot = float(np.sum((ly - float(np.mean(ly))) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan

    return a, C, r2

if n_markers < 1:
    print_usage_and_exit("--n_markers must be >= 1")

if run_dir is None:
    print_usage_and_exit("Missing required argument: run_dir")

if not os.path.isdir(run_dir):
    print_usage_and_exit(f"Directory not found: {run_dir}")

# Define output directory INSIDE the run directory so plots stay with the results
outdir = os.path.join(os.path.abspath(run_dir), f"single_run_plots_{out_format}")
if not dry_run:
    os.makedirs(outdir, exist_ok=True)

k_eff_path = os.path.join(run_dir, "k_eff.csv")
if not os.path.isfile(k_eff_path):
    print_usage_and_exit(f"Missing required file k_eff.csv in {run_dir}")

# ------------------- Load data ------------------- #
df = pd.read_csv(k_eff_path)

def pick(colnames):
    for c in colnames:
        if c in df.columns:
            return df[c].to_numpy()
    return None

kxx = pick(["k_00", "k_xx", "k11"])
kyy = pick(["k_11", "k_yy", "k22"])
time = pick(["time", "t"])
steps_k = pick(["sol_index", "step", "iter", "n"])

# NOTE: We prefer physical time from SSA_evo.dat (3rd column) when available.
# If k_eff.csv has time, we will still override it with SSA time when we can align.
using_time_fallback_index = False
if time is None:
    # Temporary placeholder; we will try to replace with SSA-aligned physical time below.
    using_time_fallback_index = True
    time = np.arange(len(kxx))

time = np.asarray(time, float)
kxx  = np.asarray(kxx, float)
kyy  = np.asarray(kyy, float)

tempC, porosity = read_metadata(run_dir)
ssa_steps, ssa_times, ssa_vals = read_ssa(run_dir)

# -----------------------------------------------------------------------------
# Prefer SSA physical time (last column of SSA_evo.dat) for time axis whenever possible.
# We do NOT normalize time; normalization applies only to SSA and k.
# -----------------------------------------------------------------------------
if ssa_times is not None and np.asarray(ssa_times).size >= 2:
    # Try to align SSA time to k_eff sampling using step indices (preferred)
    time_aligned_from_ssa = None

    if steps_k is not None and ssa_steps is not None:
        try:
            step2time = {int(s): float(t) for s, t in zip(ssa_steps, ssa_times)}
            time_aligned_from_ssa = np.array([step2time.get(int(s), np.nan) for s in steps_k], dtype=float)
        except Exception:
            time_aligned_from_ssa = None

    # If we cannot step-align, fall back to interpolating SSA time onto k_eff index/time ordering
    # (useful when k_eff.csv has no reliable time column).
    if time_aligned_from_ssa is None and using_time_fallback_index:
        # Map evenly by sample index (not ideal, but better than pretending index is time).
        try:
            # Sort SSA by time to ensure monotone x for interp
            o = np.argsort(np.asarray(ssa_times, float))
            ssa_t_sorted = np.asarray(ssa_times, float)[o]
            ssa_idx = np.linspace(0.0, 1.0, ssa_t_sorted.size)
            k_idx = np.linspace(0.0, 1.0, len(time))
            time_aligned_from_ssa = np.interp(k_idx, ssa_idx, ssa_t_sorted, left=np.nan, right=np.nan)
        except Exception:
            time_aligned_from_ssa = None

    if time_aligned_from_ssa is not None:
        tssa = np.asarray(time_aligned_from_ssa, float)
        finite = np.isfinite(tssa)
        uniq = np.unique(tssa[finite]) if np.count_nonzero(finite) else np.array([])
        if np.count_nonzero(finite) >= 2 and uniq.size >= 2:
            # Override time with SSA-aligned physical time (seconds)
            time = tssa
        else:
            print("[WARN] SSA_evo.dat present, but SSA time alignment was degenerate (all-NaN or constant); using time from k_eff.csv (or index fallback).")
    else:
        print("[WARN] SSA_evo.dat present, but could not align SSA time to k_eff sampling; using time from k_eff.csv (or index fallback).")

# Normalize SSA & k for SSA plot
# IMPORTANT: SSA_evo.dat is often sampled at many more steps than k_eff.csv.
# We align SSA to the k_eff sampling using the same logic as plot_k_eff_by_temp.py:
#   1) If k_eff.csv has sol_index/step and SSA_evo.dat has steps -> map step->SSA
#   2) Else if both have time -> interpolate SSA(time) onto k_eff time
#   3) Else -> skip SSA plot
ssa_n = kxx_n = kyy_n = None
ssa_use = kxx_use = kyy_use = time_use = None
idx_use = None

# Coerce arrays
if steps_k is not None:
    steps_k = np.asarray(steps_k, int)

# Build a working mask for finite k and time
valid_k = np.isfinite(kxx) & np.isfinite(kyy) & np.isfinite(time)
if np.count_nonzero(valid_k) >= 2:
    # Keep only valid entries for alignment/plotting
    time = time[valid_k]
    kxx = kxx[valid_k]
    kyy = kyy[valid_k]
    if steps_k is not None and steps_k.shape[0] == valid_k.shape[0]:
        steps_k = steps_k[valid_k]
    elif steps_k is not None:
        # Mismatch: disable step alignment
        steps_k = None

    # Sort by time (matches multi-temp script behavior)
    order = np.argsort(time)
    time = time[order]
    kxx = kxx[order]
    kyy = kyy[order]
    if steps_k is not None:
        steps_k = steps_k[order]

    ssa_aligned = None

    # (1) Step-based alignment (preferred)
    if steps_k is not None and ssa_steps is not None and ssa_vals is not None:
        try:
            step2ssa = {int(s): float(v) for s, v in zip(ssa_steps, ssa_vals)}
            ssa_aligned = np.array([step2ssa.get(int(s), np.nan) for s in steps_k], dtype=float)
        except Exception:
            ssa_aligned = None

    # (2) Time-based interpolation alignment
    if ssa_aligned is None and ssa_times is not None and ssa_vals is not None:
        ssa_times = np.asarray(ssa_times, float)
        ssa_vals = np.asarray(ssa_vals, float)
        # Sort SSA by time for interpolation safety
        o = np.argsort(ssa_times)
        ssa_times_s = ssa_times[o]
        ssa_vals_s = ssa_vals[o]
        ssa_aligned = np.interp(time, ssa_times_s, ssa_vals_s, left=np.nan, right=np.nan)

    # Normalize (guard against NaN/zero)
    # Also keep aligned (unnormalized) arrays so we can place the SAME time-markers on SSA plots.
    ssa_use = kxx_use = kyy_use = time_use = None
    idx_use = None

    if ssa_aligned is not None:
        ssa_aligned = np.asarray(ssa_aligned, float)
        valid_ssa = np.isfinite(ssa_aligned) & np.isfinite(kxx) & np.isfinite(kyy) & np.isfinite(time)
        if np.count_nonzero(valid_ssa) >= 2:
            ssa_use = ssa_aligned[valid_ssa]
            kxx_use = kxx[valid_ssa]
            kyy_use = kyy[valid_ssa]
            time_use = time[valid_ssa]
            # Map: SSA-aligned index -> index into the (sorted/filtered) k/time arrays
            idx_use = np.nonzero(valid_ssa)[0]

            if ssa_use[0] != 0 and kxx_use[0] != 0 and kyy_use[0] != 0:
                ssa_n = ssa_use / ssa_use[0]
                kxx_n = kxx_use / kxx_use[0]
                kyy_n = kyy_use / kyy_use[0]


# ------------------- Helper for marking time points ------------------- #
def choose_even_time_indices(t: np.ndarray, n: int = 4) -> np.ndarray:
    """Return indices of points whose times are approximately evenly spaced over [tmin, tmax].

    Includes endpoints when possible. Works with nonuniform sampling.
    Assumes t is sorted ascending.
    """
    t = np.asarray(t, float)
    if t.size == 0:
        return np.array([], dtype=int)
    if t.size <= n:
        return np.arange(t.size, dtype=int)

    tmin = float(np.nanmin(t))
    tmax = float(np.nanmax(t))
    if not np.isfinite(tmin) or not np.isfinite(tmax) or tmax == tmin:
        # Degenerate time axis: fall back to evenly spaced indices
        return np.unique(np.round(np.linspace(0, t.size - 1, n)).astype(int))

    targets = np.linspace(tmin, tmax, n)
    idx = []
    for tt in targets:
        j = int(np.searchsorted(t, tt, side="left"))
        if j <= 0:
            idx.append(0)
        elif j >= t.size:
            idx.append(t.size - 1)
        else:
            # pick closer of j-1 and j
            if abs(t[j] - tt) < abs(t[j - 1] - tt): 
                idx.append(j)
            else:
                idx.append(j - 1)

    # Ensure uniqueness (can happen if sampling is coarse)
    idx = np.unique(np.array(idx, dtype=int))
    return idx


# Compute marker indices ONCE and expose mapping for SSA plots as well
def compute_time_marker_indices(time_scaled_arr: np.ndarray, log_plots_flag: bool, n: int = 4) -> tuple[np.ndarray, np.ndarray]:
    """Return (idx_global, t_marked).

    idx_global are indices into the current (sorted/filtered) time/k arrays.
    t_marked are the corresponding time values (in scaled units).

    We pick points evenly spaced in time, and for log-x plots we require t > 0.
    """
    t = np.asarray(time_scaled_arr, float)
    if t.size == 0:
        return np.array([], dtype=int), np.array([], dtype=float)

    if log_plots_flag:
        mask = np.isfinite(t) & (t > 0)
    else:
        mask = np.isfinite(t)

    if np.count_nonzero(mask) < 2:
        return np.array([], dtype=int), np.array([], dtype=float)

    t_valid = t[mask]
    idx_in_valid = choose_even_time_indices(t_valid, n=n)
    idx_in_valid = np.asarray(idx_in_valid, dtype=int)

    idx_global_all = np.nonzero(mask)[0]
    idx_global = idx_global_all[idx_in_valid]
    return idx_global, t[idx_global]

# ----------- Marker reporting helper -----------
def print_time_marker_report(idx_mark_global: np.ndarray,
                             time_seconds_sorted: np.ndarray,
                             time_scaled_sorted: np.ndarray,
                             unit: str,
                             steps_sorted: Optional[np.ndarray],
                             *,
                             context: str = "time-plots",
                             ssa_steps_all: Optional[np.ndarray] = None,
                             ssa_times_all: Optional[np.ndarray] = None):
    """Print a small report of which points were chosen for marker placement.

    idx_mark_global indexes into the *sorted/filtered* time/k arrays.

    `context` is a label so it's obvious whether the report is for k-vs-time or k-vs-SSA.

    If SSA steps/times are provided, we also sanity-check that the `step` values map to the
    expected SSA time values (useful to catch step mismatches between files).
    """
    idx_mark_global = np.asarray(idx_mark_global, dtype=int)
    if idx_mark_global.size == 0:
        print(f"[mark_time_points:{context}] No valid marker indices found (time axis may be degenerate).")
        return

    print(f"\n[mark_time_points:{context}] Selected marker points (evenly spaced in time):")
    print(f"  (n_markers = {idx_mark_global.size})")
    header = "  idx   time[s]        time[{}]".format(unit)
    if steps_sorted is not None:
        header += "        step"
    print(header)

    marked_steps = []
    for idx in idx_mark_global:
        tsec = float(time_seconds_sorted[idx])
        tsc  = float(time_scaled_sorted[idx])
        if steps_sorted is not None:
            st = int(steps_sorted[idx])
            marked_steps.append(st)
            print(f"  {idx:4d}  {tsec: .6e}  {tsc: .6e}  {st:d}")
        else:
            print(f"  {idx:4d}  {tsec: .6e}  {tsc: .6e}")

    # Optional: sanity-check step->time mapping using SSA_evo.dat
    if steps_sorted is not None and ssa_steps_all is not None and ssa_times_all is not None:
        try:
            step2time = {int(s): float(t) for s, t in zip(ssa_steps_all, ssa_times_all)}
            missing = [s for s in marked_steps if s not in step2time]
            if missing:
                print(f"[WARN] Some marker steps are missing from SSA_evo.dat step list: {missing}")
        except Exception:
            pass

# ------------------- XY Padding Helpers ------------------- #
def _linear_pad_amount(vmin: float, vmax: float, frac: float) -> float:
    """Additive pad amount for linear axes."""
    rng = vmax - vmin
    if not np.isfinite(rng) or rng <= 0:
        # Degenerate range: use a scale based on magnitude or 1.0
        scale = max(abs(vmin), abs(vmax), 1.0)
        return frac * scale
    return frac * rng

def apply_xy_padding(ax, x, y, logx: bool, logy: bool,
                     xpad_min_frac: float,
                     ypad_min_frac: float,
                     xpad_left_frac: float = None,
                     xpad_right_frac: float = None,
                     ypad_bottom_frac: float = None,
                     ypad_top_frac: float = None):
    """
    Apply padding on all four sides.
    - Always applies a minimum padding (xpad_min_frac / ypad_min_frac).
    - If a side-specific padding is provided and is larger than the minimum, it overrides.
    - For log axes, padding is multiplicative in decades based on the log-range.
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)

    # -------- X limits --------
    if x.size > 0 and np.isfinite(x).any():
        if logx:
            xv = x[np.isfinite(x) & (x > 0)]
            if xv.size >= 2:
                xmin, xmax = float(np.min(xv)), float(np.max(xv))
                lxmin, lxmax = np.log10(xmin), np.log10(xmax)
                lrange = lxmax - lxmin
                if not np.isfinite(lrange) or lrange <= 0:
                    lrange = 1.0
                # minimum pad (both sides)
                lp_min = xpad_min_frac * lrange
                lp_left = lp_min if xpad_left_frac is None else max(lp_min, xpad_left_frac * lrange)
                lp_right = lp_min if xpad_right_frac is None else max(lp_min, xpad_right_frac * lrange)
                ax.set_xlim(10**(lxmin - lp_left), 10**(lxmax + lp_right))
        else:
            xv = x[np.isfinite(x)]
            if xv.size >= 2:
                xmin, xmax = float(np.min(xv)), float(np.max(xv))
                pmin = _linear_pad_amount(xmin, xmax, xpad_min_frac)
                pL = pmin if xpad_left_frac is None else max(pmin, _linear_pad_amount(xmin, xmax, xpad_left_frac))
                pR = pmin if xpad_right_frac is None else max(pmin, _linear_pad_amount(xmin, xmax, xpad_right_frac))
                ax.set_xlim(xmin - pL, xmax + pR)

    # -------- Y limits --------
    if y.size > 0 and np.isfinite(y).any():
        if logy:
            yv = y[np.isfinite(y) & (y > 0)]
            if yv.size >= 2:
                ymin, ymax = float(np.min(yv)), float(np.max(yv))
                lymin, lymax = np.log10(ymin), np.log10(ymax)
                lrange = lymax - lymin
                if not np.isfinite(lrange) or lrange <= 0:
                    lrange = 1.0
                lp_min = ypad_min_frac * lrange
                lp_bot = lp_min if ypad_bottom_frac is None else max(lp_min, ypad_bottom_frac * lrange)
                lp_top = lp_min if ypad_top_frac is None else max(lp_min, ypad_top_frac * lrange)
                ax.set_ylim(10**(lymin - lp_bot), 10**(lymax + lp_top))
        else:
            yv = y[np.isfinite(y)]
            if yv.size >= 2:
                ymin, ymax = float(np.min(yv)), float(np.max(yv))
                pmin = _linear_pad_amount(ymin, ymax, ypad_min_frac)
                pB = pmin if ypad_bottom_frac is None else max(pmin, _linear_pad_amount(ymin, ymax, ypad_bottom_frac))
                pT = pmin if ypad_top_frac is None else max(pmin, _linear_pad_amount(ymin, ymax, ypad_top_frac))
                ax.set_ylim(ymin - pB, ymax + pT)


# ------------------- Save helper (for dry_run) ------------------- #
def maybe_savefig(fig, path: str, dpi_val: int, fmt: str, dry_run_flag: bool) -> None:
    """Save a figure unless dry_run is enabled."""
    if dry_run_flag:
        print(f"[dry_run] Would save: {path}")
        return
    fig.savefig(path, dpi=dpi_val, transparent=True, format=fmt)

# ------------------- Plots ------------------- #

themes = []
if theme_arg == "both":
    themes = ["light", "dark"]
else:
    themes = [theme_arg]

for theme in themes:
    set_theme(theme)

    # ---- kxx vs time ----
    # For time plots, force physical time to be shown in days (we do NOT normalize time).
    unit = "d"
    time_scaled = scale_to_unit(time, unit)

    # For time plots, normalize K_eff by its initial value: K_eff(t) / K_eff0
    def _safe_norm_by_first(arr: np.ndarray) -> np.ndarray:
        arr = np.asarray(arr, float)
        if arr.size == 0:
            return arr
        k0 = float(arr[0])
        if (not np.isfinite(k0)) or k0 == 0.0:
            return arr
        return arr / k0

    kxx_time = _safe_norm_by_first(kxx)
    kyy_time = _safe_norm_by_first(kyy)
    # Pre-compute the marker indices ONCE so time plots and SSA plots share the same time points
    idx_mark_global, _t_marked = (np.array([], dtype=int), np.array([], dtype=float))
    if mark_time_points:
        idx_mark_global, _t_marked = compute_time_marker_indices(time_scaled, log_plots, n=n_markers)
    if mark_time_points:
        # Print which points were selected (for the k-vs-time plots). Uses sorted/filtered arrays.
        steps_sorted = steps_k if (steps_k is not None and steps_k.shape[0] == time.shape[0]) else None
        if steps_k is not None and steps_sorted is None:
            print(f"[WARN] steps_k present but length mismatch vs time array (steps={steps_k.shape[0]} time={time.shape[0]}). Step column will be ignored for reporting.")

        print_time_marker_report(
            idx_mark_global,
            time,            # seconds (after SSA alignment override)
            time_scaled,
            unit,
            steps_sorted,
            context="time-plots",
            ssa_steps_all=ssa_steps,
            ssa_times_all=ssa_times,
        )

    def _title_with_meta(base: str) -> str:
        t = base
        if tempC is not None or porosity is not None:
            subtitle = []
            if tempC is not None:
                subtitle.append(f"T={tempC:.1f} °C")
            if porosity is not None:
                subtitle.append(f"φ={porosity:.2f}")
            t += "\n" + ", ".join(subtitle)
        return t

    def _plot_single_vs_time(y, label_tex: str, outname: str):
        # Determine if normalized (heuristic: y[0] close to 1.0)
        is_normalized = False
        if y.size > 0:
            is_normalized = np.isclose(y[0], 1.0, atol=1e-3)

        if is_normalized:
            ylab = "Normalized Effective\nThermal Conductivity ($K_{eff}/K_{eff,0}$)"
            title_main = "Normalized Effective Thermal Conductivity\nvs. Time"
        else:
            ylab = "Effective\nThermal Conductivity"
            title_main = "Effective Thermal Conductivity\nvs. Time"

        fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))

        if log_plots:
            mask = (time_scaled > 0) & (y > 0)
            t_plot = time_scaled[mask]
            y_plot = y[mask]
            ax.set_xscale("log")
            ax.set_yscale("log")
        else:
            t_plot = time_scaled
            y_plot = y

        line = ax.plot(t_plot, y_plot, lw=2.2, label=label_tex)[0]
        # Optional power-law fit overlay: y = C * t^a
        if fit_powerlaw:
            a_fit, C_fit, r2_fit = fit_powerlaw_exponent(t_plot, y_plot)
            if a_fit is not None and C_fit is not None:
                # Build a smooth fit curve over the plotted time range
                t_fit = np.asarray(t_plot, float)
                t_fit = t_fit[np.isfinite(t_fit) & (t_fit > 0)]
                if t_fit.size >= 2:
                    t_fit = np.logspace(np.log10(np.min(t_fit)), np.log10(np.max(t_fit)), 200)
                    y_fit = C_fit * (t_fit ** a_fit)
                    ax.plot(t_fit, y_fit, lw=1.8, ls="--",
                            label=fr"fit: $t^{{{a_fit:.3f}}}$ (R$^2$={r2_fit:.3f})")
                    print(f"[fit_powerlaw] {outname}: a={a_fit:.6f}, C={C_fit:.6e}, R^2={r2_fit:.6f}")
            else:
                print(f"[fit_powerlaw] {outname}: not enough positive (t,y) points to fit.")
        line_color = line.get_color()
        if mark_time_points and idx_mark_global.size > 0:
            # Convert global indices (into time_scaled/y arrays) to indices in the plotted arrays
            if log_plots:
                # t_plot/y_plot are masked versions; reconstruct mapping
                mask_plot = (time_scaled > 0) & np.isfinite(time_scaled) & np.isfinite(y) & (y > 0)
                idx_global_all = np.nonzero(mask_plot)[0]
                # Keep only marker indices that survived the mask
                keep = np.isin(idx_global_all, idx_mark_global)
                idx_in_plot = np.nonzero(keep)[0]
            else:
                # Linear: plotted arrays match the full arrays
                idx_in_plot = idx_mark_global

            if idx_in_plot.size > 0:
                ax.plot(
                    t_plot[idx_in_plot], y_plot[idx_in_plot],
                    linestyle="None",
                    marker="o",
                    markersize=5.5,
                    markerfacecolor=line_color,
                    markeredgecolor=line_color,
                    markeredgewidth=0.0,
                    zorder=5,
                )

        # Apply padding on all four sides:
        # - Always minimum padding (PAD_MIN_FRAC)
        # - Extra room for microstructure snapshots on the y-bottom for time plots
        if log_plots:
            apply_xy_padding(
                ax, t_plot, y_plot,
                logx=True, logy=True,
                xpad_min_frac=PAD_MIN_FRAC,
                ypad_min_frac=PAD_MIN_FRAC,
                xpad_left_frac=TIME_AXIS_PAD_FRACTION,
                xpad_right_frac=TIME_AXIS_PAD_FRACTION,
                ypad_bottom_frac=PAD_MIN_FRAC,
                ypad_top_frac=Y_PAD_FRAC_LOG
            )
        else:
            apply_xy_padding(
                ax, t_plot, y_plot,
                logx=False, logy=False,
                xpad_min_frac=PAD_MIN_FRAC,
                ypad_min_frac=PAD_MIN_FRAC,
                xpad_left_frac=TIME_AXIS_PAD_FRACTION,
                xpad_right_frac=TIME_AXIS_PAD_FRACTION,
                ypad_bottom_frac=Y_PAD_FRAC_LINEAR,
                ypad_top_frac=PAD_MIN_FRAC
            )

        # Axis labels and title logic
        if not no_labels:
            set_compact_axis_labels(ax, "Time (d)", ylab)
            ax.set_title(_title_with_meta(title_main))
        else:
            ax.set_xlabel("")
            ax.set_ylabel("")
            ax.set_title("")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(False)

        if not no_legend:
            ax.legend(loc="best", fontsize=max(8, LEGEND_FONTSIZE - 1))

        fig.tight_layout()
        maybe_savefig(fig, os.path.join(outdir, outname), dpi, out_format, dry_run)
        plt.close(fig)

    fname_kxx_time = f"k_xx_vs_time_{theme}.{out_format}"
    fname_kyy_time = f"k_yy_vs_time_{theme}.{out_format}"
    if log_plots:
        fname_kxx_time = f"log_{fname_kxx_time}"
        fname_kyy_time = f"log_{fname_kyy_time}"

    _plot_single_vs_time(kxx_time, r"$k_{xx}$", fname_kxx_time)
    _plot_single_vs_time(kyy_time, r"$k_{yy}$", fname_kyy_time)

    # ---- k vs SSA ----
    if ssa_n is not None and kxx_n is not None and kyy_n is not None:
        # Filter data for log plots if needed
        # Also determine marker indices in SSA-plot space that correspond to idx_mark_global
        if log_plots:
            mask = (ssa_n > 0) & (kxx_n > 0) & (kyy_n > 0)
            ssa_plot = ssa_n[mask]
            kxx_plot = kxx_n[mask]
            kyy_plot = kyy_n[mask]
        else:
            mask = np.ones_like(ssa_n, dtype=bool)
            ssa_plot = ssa_n
            kxx_plot = kxx_n
            kyy_plot = kyy_n

        # Determine which SSA points correspond to the SAME global time-marker indices
        idx_mark_ssa_plot = np.array([], dtype=int)
        if mark_time_points and idx_mark_global.size > 0 and idx_use is not None:
            # idx_use maps SSA-aligned indices -> indices into time/k arrays
            idx_mark_ssa = np.nonzero(np.isin(idx_use, idx_mark_global))[0]
            # Apply the SSA plotting mask (e.g., positivity for log plots)
            if idx_mark_ssa.size > 0:
                idx_in_masked = np.nonzero(mask)[0]  # positions kept in ssa_plot
                keep = np.isin(idx_in_masked, idx_mark_ssa)
                idx_mark_ssa_plot = np.nonzero(keep)[0]

        def _plot_single_vs_ssa(y, label_tex: str, outname: str):
            fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
            line = ax.plot(ssa_plot, y, lw=2.2, label=label_tex)[0]
            line_color = line.get_color()

            # Optional inverse power-law fit overlay on SSA plots: y = C / (SSA)^a
            if fit_powerlaw:
                a_fit, C_fit, r2_fit = fit_inverse_powerlaw_exponent(ssa_plot, np.asarray(y, float))
                if a_fit is not None and C_fit is not None:
                    x_fit = np.asarray(ssa_plot, float)
                    x_fit = x_fit[np.isfinite(x_fit) & (x_fit > 0)]
                    if x_fit.size >= 2:
                        # Smooth fit curve over SSA range
                        x_fit = np.logspace(np.log10(np.min(x_fit)), np.log10(np.max(x_fit)), 200)
                        y_fit = C_fit / (x_fit ** a_fit)
                        ax.plot(
                            x_fit, y_fit,
                            lw=1.8, ls="--",
                            label=fr"fit: $\mathrm{{SSA}}^{{-{a_fit:.3f}}}$ (R$^2$={r2_fit:.3f})"
                        )
                        print(f"[fit_powerlaw:ssa] {outname}: a={a_fit:.6f}, C={C_fit:.6e}, R^2={r2_fit:.6f}  (fit: y=C/SSA^a)")
                else:
                    print(f"[fit_powerlaw:ssa] {outname}: not enough positive (SSA,y) points to fit.")

            # Add the SAME time-series marker points (mapped onto SSA space)
            if mark_time_points and idx_mark_ssa_plot.size > 0:
                ax.plot(
                    ssa_plot[idx_mark_ssa_plot], np.asarray(y, float)[idx_mark_ssa_plot],
                    linestyle="None",
                    marker="o",
                    markersize=5.5,
                    markerfacecolor=line_color,
                    markeredgecolor=line_color,
                    markeredgewidth=0.0,
                    zorder=5,
                )

            # Apply padding on all four sides:
            # SSA plots: you wanted headroom above the curve (especially for insets),
            # while still keeping a minimum padding everywhere.
            y_plot = np.asarray(y, float)
            if log_plots:
                apply_xy_padding(
                    ax, ssa_plot, y_plot,
                    logx=True, logy=True,
                    xpad_min_frac=PAD_MIN_FRAC,
                    ypad_min_frac=PAD_MIN_FRAC,
                    xpad_left_frac=PAD_MIN_FRAC,
                    xpad_right_frac=PAD_MIN_FRAC,
                    ypad_bottom_frac=PAD_MIN_FRAC,
                    ypad_top_frac=Y_PAD_FRAC_LOG
                )
            else:
                apply_xy_padding(
                    ax, ssa_plot, y_plot,
                    logx=False, logy=False,
                    xpad_min_frac=PAD_MIN_FRAC,
                    ypad_min_frac=PAD_MIN_FRAC,
                    xpad_left_frac=PAD_MIN_FRAC,
                    xpad_right_frac=PAD_MIN_FRAC,
                    ypad_bottom_frac=PAD_MIN_FRAC,
                    ypad_top_frac=Y_PAD_FRAC_LINEAR
                )

            # Axis labels and title logic
            if not no_labels:
                set_compact_axis_labels(ax, "Normalized Specific Surface Area", "Normalized Effective Thermal Conductivity")
                title_main = "Normalized Effective Thermal Conductivity\nvs. Normalized Specific Surface Area"
                ax.set_title(_title_with_meta(title_main))
            else:
                ax.set_xlabel("")
                ax.set_ylabel("")
                ax.set_title("")

            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.grid(False)

            if log_plots:
                ax.set_xscale("log")
                ax.set_yscale("log")

            if not no_legend:
                ax.legend(loc="best", fontsize=max(8, LEGEND_FONTSIZE - 1))

            fig.tight_layout()
            maybe_savefig(fig, os.path.join(outdir, outname), dpi, out_format, dry_run)
            plt.close(fig)

        fname_kxx_ssa = f"k_xx_vs_ssa_norm_{theme}.{out_format}"
        fname_kyy_ssa = f"k_yy_vs_ssa_norm_{theme}.{out_format}"
        if log_plots:
            fname_kxx_ssa = f"log_{fname_kxx_ssa}"
            fname_kyy_ssa = f"log_{fname_kyy_ssa}"

        _plot_single_vs_ssa(kxx_plot, r"$k_{xx}$", fname_kxx_ssa)
        _plot_single_vs_ssa(kyy_plot, r"$k_{yy}$", fname_kyy_ssa)

if dry_run:
    print(f"✅ Dry run complete. No files were written. (Would have used output dir: {outdir})")
else:
    print(f"✅ Single-run plots ({out_format}) generated in: {outdir}")
#!/usr/bin/env python3
"""
Single-run plotting for k_eff results.

Reads:
  - <RUN_DIR>/k_eff.csv
  - <RUN_DIR>/SSA_evo.dat (optional)
  - <RUN_DIR>/metadata.json (optional)

Outputs:
  - <RUN_DIR>/output_plots_single_run/
      kxx_vs_time_[light|dark].png
      kyy_vs_time_[light|dark].png
      k_eff_overlay_vs_time_[light|dark].png
    If SSA_evo.dat present:
      kxx_vs_ssa_norm_[light|dark].png
      kyy_vs_ssa_norm_[light|dark].png
      k_eff_overlay_vs_ssa_norm_[light|dark].png
"""

from __future__ import annotations

import os, re, json
from typing import Dict, Any, Optional, Tuple
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ======================== USER SETTINGS ========================

# Point this to ONE simulation output folder (the one containing k_eff.csv)
RUN_DIR = "/path/to/one/run_folder"

# Plot options
THEME = "both"          # "light", "dark", or "both"
DPI = 300
NORMALIZE = False       # normalize kxx/kyy by first value for time plots
LOG_PLOTS = False       # log-log versions for SSA plots + time plots (skips nonpositive)
USE_CMOCEAN_THERMAL = False  # only affects temperature color in multi-run; kept for compatibility

# ======================== GLOBAL PLOT STYLING ========================

FIG_W_IN = 6.0
FIG_H_IN = 3.5

TITLE_FONTSIZE = 16
AXIS_FONTSIZE = 18
TICK_FONTSIZE = 10
LEGEND_FONTSIZE = 10
TIME_AXIS_PAD_FRACTION = 0.02

MAX_XLABEL_CHARS = 34
MAX_YLABEL_CHARS = 26

def _wrap_label(text: str, max_chars: int) -> str:
    if text is None:
        return ""
    s = str(text)
    if len(s) <= max_chars or "\n" in s:
        return s
    cut = s.rfind(" ", 0, max_chars + 1)
    if cut == -1:
        return s
    return s[:cut].rstrip() + "\n" + s[cut + 1 :].lstrip()

def set_compact_axis_labels(ax: plt.Axes, xlabel: str, ylabel: str):
    ax.set_xlabel(_wrap_label(xlabel, MAX_XLABEL_CHARS))
    ax.set_ylabel(_wrap_label(ylabel, MAX_YLABEL_CHARS))
    try:
        ax.xaxis.label.set_wrap(True)
        ax.yaxis.label.set_wrap(True)
    except Exception:
        pass

def set_theme(theme: str):
    plt.rcParams.update(plt.rcParamsDefault)
    plt.rcParams.update({
        "figure.facecolor": "none",
        "axes.facecolor": "none",
        "savefig.facecolor": "none",
        "savefig.transparent": True,
    })

    if theme == "dark":
        fg = "white"
        plt.rcParams.update({
            "text.color": fg,
            "axes.labelcolor": fg,
            "axes.edgecolor": fg,
            "xtick.color": fg,
            "ytick.color": fg,
        })
    else:
        fg = "black"
        plt.rcParams.update({
            "text.color": fg,
            "axes.labelcolor": fg,
            "axes.edgecolor": fg,
            "xtick.color": fg,
            "ytick.color": fg,
        })

    plt.rcParams.update({
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
    })

def autoscale_time(arr: np.ndarray) -> Tuple[np.ndarray, str]:
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

def _is_normalized_series(arr: np.ndarray) -> bool:
    if arr is None or not hasattr(arr, "size") or arr.size == 0:
        return False
    if not np.isfinite(arr[0]):
        return False
    return abs(float(arr[0]) - 1.0) < 1e-6

def _coerce_float(v) -> Optional[float]:
    try:
        return float(v)
    except Exception:
        return None

def read_json(path: str) -> Optional[Dict[str, Any]]:
    try:
        with open(path, "r") as f:
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

def find_temp_from_meta(meta: Dict[str, Any]) -> Optional[float]:
    try:
        p = meta.get("dsm_run", {}).get("parameters", {}).get("temperature_C", None)
        if p is not None:
            val = _coerce_float(p)
            if val is not None:
                return val
    except Exception:
        pass
    for k in ("temperature_C", "tempC", "T_C", "TdegC"):
        if k in meta:
            val = _coerce_float(meta[k])
            if val is not None:
                return val
    for _, k, v in _walk_items(meta):
        if "temp" in k.lower():
            val = _coerce_float(v)
            if val is not None:
                return val
    return None

def find_porosity_from_meta(meta: Dict[str, Any]) -> Optional[float]:
    try:
        s = meta.get("structure", {})
        for key in ("porosity_target", "porosity_achieved"):
            if key in s:
                val = _coerce_float(s[key])
                if val is not None:
                    return val
    except Exception:
        pass
    for _, k, v in _walk_items(meta):
        kl = k.lower()
        if "poros" in kl or kl == "phi" or kl.endswith(".phi"):
            val = _coerce_float(v)
            if val is not None:
                return val
    return None

def parse_temp_from_name(name: str) -> Optional[float]:
    m = re.search(r"_Tm(-?\d+(?:\.\d+)?)C?(_|$)", name)
    if m:
        try:
            return float(m.group(1))
        except Exception:
            return None
    return None

def read_ssa_map(ssa_path: str) -> Tuple[Optional[np.ndarray], Optional[np.ndarray], Optional[np.ndarray]]:
    """Expected columns: [SSA, (optional), time, step]"""
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

def _fmt_temp_phi(tempC: Optional[float], por: Optional[float]) -> str:
    bits = []
    if tempC is not None and np.isfinite(tempC):
        bits.append(f"T={tempC:.1f} °C")
    if por is not None and np.isfinite(por):
        bits.append(f"φ={por:.2f}")
    return ", ".join(bits)

def load_single_run(run_dir: str, normalize: bool) -> Dict[str, Any]:
    run_dir = os.path.abspath(os.path.expanduser(run_dir))
    name = os.path.basename(run_dir.rstrip(os.sep))

    # metadata
    meta = None
    for mname in ("metadata.json", "Metadata.json", "meta.json"):
        mp = os.path.join(run_dir, mname)
        if os.path.exists(mp):
            meta = read_json(mp)
            if meta:
                break

    tempC = find_temp_from_meta(meta) if meta else None
    porosity = find_porosity_from_meta(meta) if meta else None

    if tempC is None:
        tempC = parse_temp_from_name(name)

    # k_eff.csv
    kcsv = os.path.join(run_dir, "k_eff.csv")
    if not os.path.exists(kcsv):
        raise FileNotFoundError(f"Missing k_eff.csv in {run_dir}")

    df = pd.read_csv(kcsv)

    def pick_col(cands):
        for c in cands:
            if c in df.columns:
                return df[c].to_numpy()
        return None

    kxx = pick_col(["k_00", "k_xx", "k11", "k_11_xx"])
    kyy = pick_col(["k_11", "k_yy", "k22", "k_22_yy"])
    if kxx is None or kyy is None:
        raise ValueError(f"k columns not found. Columns present: {list(df.columns)}")

    steps_k = pick_col(["sol_index", "step", "iter", "n"])
    time_k = pick_col(["time", "t"])

    # SSA mapping if needed
    ssa_aligned = None
    ssa_path = os.path.join(run_dir, "SSA_evo.dat")
    ssa_steps, ssa_times, ssa_vals = read_ssa_map(ssa_path)

    if time_k is None and steps_k is not None:
        if ssa_steps is not None and ssa_times is not None:
            step2time = {int(s): float(t) for s, t in zip(ssa_steps, ssa_times)}
            time = np.array([step2time.get(int(s), np.nan) for s in steps_k], dtype=float)
        else:
            time = np.arange(len(kxx), dtype=float)
    elif time_k is not None:
        time = np.asarray(time_k, dtype=float)
    else:
        time = np.arange(len(kxx), dtype=float)

    # SSA alignment to k points
    if steps_k is not None and ssa_steps is not None and ssa_vals is not None:
        step2ssa = {int(s): float(v) for s, v in zip(ssa_steps, ssa_vals)}
        ssa_aligned = np.array([step2ssa.get(int(s), np.nan) for s in steps_k], dtype=float)
    elif time is not None and ssa_times is not None and ssa_vals is not None:
        order = np.argsort(ssa_times)
        ssa_times_sorted = np.asarray(ssa_times)[order]
        ssa_vals_sorted = np.asarray(ssa_vals)[order]
        ssa_aligned = np.interp(time, ssa_times_sorted, ssa_vals_sorted, left=np.nan, right=np.nan)

    time = np.asarray(time, dtype=float)
    kxx = np.asarray(kxx, dtype=float)
    kyy = np.asarray(kyy, dtype=float)

    valid = np.isfinite(time) & np.isfinite(kxx) & np.isfinite(kyy)
    time, kxx, kyy = time[valid], kxx[valid], kyy[valid]
    if ssa_aligned is not None:
        ssa_aligned = np.asarray(ssa_aligned, dtype=float)
        if ssa_aligned.shape[0] == valid.shape[0]:
            ssa_aligned = ssa_aligned[valid]
        else:
            ssa_aligned = None

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
        run_dir=run_dir,
        tempC=(None if tempC is None else float(tempC)),
        porosity=(None if porosity is None else float(porosity)),
        time=time,
        kxx=kxx,
        kyy=kyy,
        ssa=ssa_aligned,
    )

def _save(fig, out_path: str, dpi: int):
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print("[OK]", out_path)

def plot_vs_time(run: Dict[str, Any], outdir: str, theme: str, dpi: int, which: str):
    def emit(theme_name: str, suffix: str):
        set_theme(theme_name)
        fig, ax = plt.subplots(figsize=(FIG_W_IN, FIG_H_IN))

        t_scaled, unit = autoscale_time(run["time"])
        y = run[which]

        if LOG_PLOTS:
            ax.set_xscale("log")
            ax.set_yscale("log")
            m = np.isfinite(t_scaled) & np.isfinite(y) & (t_scaled > 0) & (y > 0)
            t_plot = t_scaled[m]
            y_plot = y[m]
        else:
            t_plot = t_scaled
            y_plot = y

        ax.plot(t_plot, y_plot, lw=2.2)

        normalized = _is_normalized_series(y)
        ylab = "Normalized\nEffective Thermal Conductivity" if normalized else "Effective\nThermal Conductivity"
        title_main = ("Normalized Effective Thermal Conductivity\nvs. Time" if normalized
                      else "Effective Thermal Conductivity\nvs. Time")
        subtitle = _fmt_temp_phi(run.get("tempC"), run.get("porosity"))
        ax.set_title(f"{title_main}\n{subtitle}" if subtitle else title_main)

        set_compact_axis_labels(ax, f"Time ({unit})", ylab)

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(False)

        if not LOG_PLOTS and t_plot.size > 0 and np.isfinite(t_plot).any():
            xmin, xmax = float(np.nanmin(t_plot)), float(np.nanmax(t_plot))
            pad = TIME_AXIS_PAD_FRACTION * (xmax - xmin if xmax > xmin else 1.0)
            ax.set_xlim(xmin - pad, xmax)

        out_path = os.path.join(outdir, f"{which}_vs_time_{suffix}.png")
        _save(fig, out_path, dpi)

    if theme in ("light", "both"):
        emit("light", "light")
    if theme in ("dark", "both"):
        emit("dark", "dark")

def plot_overlay_vs_time(run: Dict[str, Any], outdir: str, theme: str, dpi: int):
    def emit(theme_name: str, suffix: str):
        set_theme(theme_name)
        fig, ax = plt.subplots(figsize=(FIG_W_IN, FIG_H_IN))

        t_scaled, unit = autoscale_time(run["time"])
        kxx = np.asarray(run["kxx"], dtype=float)
        kyy = np.asarray(run["kyy"], dtype=float)

        if LOG_PLOTS:
            ax.set_xscale("log")
            ax.set_yscale("log")
            mxx = np.isfinite(t_scaled) & np.isfinite(kxx) & (t_scaled > 0) & (kxx > 0)
            myy = np.isfinite(t_scaled) & np.isfinite(kyy) & (t_scaled > 0) & (kyy > 0)
            ax.plot(t_scaled[mxx], kxx[mxx], lw=2.2, label=r"$k_{xx}$")
            ax.plot(t_scaled[myy], kyy[myy], lw=2.2, label=r"$k_{yy}$")
        else:
            ax.plot(t_scaled, kxx, lw=2.2, label=r"$k_{xx}$")
            ax.plot(t_scaled, kyy, lw=2.2, label=r"$k_{yy}$")

        normalized = _is_normalized_series(kxx)
        ylab = "Normalized\nEffective Thermal Conductivity" if normalized else "Effective\nThermal Conductivity"
        title_main = ("Normalized Effective Thermal Conductivity\nvs. Time" if normalized
                      else "Effective Thermal Conductivity\nvs. Time")
        subtitle = _fmt_temp_phi(run.get("tempC"), run.get("porosity"))
        ax.set_title(f"{title_main}\n{subtitle}" if subtitle else title_main)

        set_compact_axis_labels(ax, f"Time ({unit})", ylab)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(False)
        ax.legend(loc="best", fontsize=max(8, LEGEND_FONTSIZE - 1))

        if not LOG_PLOTS and t_scaled.size > 0 and np.isfinite(t_scaled).any():
            xmin, xmax = float(np.nanmin(t_scaled)), float(np.nanmax(t_scaled))
            pad = TIME_AXIS_PAD_FRACTION * (xmax - xmin if xmax > xmin else 1.0)
            ax.set_xlim(xmin - pad, xmax)

        out_path = os.path.join(outdir, f"k_eff_overlay_vs_time_{suffix}.png")
        _save(fig, out_path, dpi)

    if theme in ("light", "both"):
        emit("light", "light")
    if theme in ("dark", "both"):
        emit("dark", "dark")

def plot_vs_ssa_norm(run: Dict[str, Any], outdir: str, theme: str, dpi: int, which: str):
    ssa = run.get("ssa")
    if ssa is None:
        return

    def emit(theme_name: str, suffix: str):
        set_theme(theme_name)
        fig, ax = plt.subplots(figsize=(FIG_W_IN, FIG_H_IN))

        ssa = np.asarray(run["ssa"], dtype=float)
        y = np.asarray(run[which], dtype=float)

        valid = np.isfinite(ssa) & np.isfinite(y)
        if LOG_PLOTS:
            valid = valid & (ssa > 0) & (y > 0)
        if np.count_nonzero(valid) < 2:
            plt.close(fig)
            return

        ssa = ssa[valid]
        y = y[valid]
        if ssa[0] == 0 or y[0] == 0:
            plt.close(fig)
            return

        ssa_n = ssa / ssa[0]
        y_n = y / y[0]

        if LOG_PLOTS:
            ax.set_xscale("log")
            ax.set_yscale("log")
            m = (ssa_n > 0) & (y_n > 0)
            ssa_n = ssa_n[m]
            y_n = y_n[m]
            if ssa_n.size < 2:
                plt.close(fig)
                return

        ax.plot(ssa_n, y_n, lw=2.2)

        set_compact_axis_labels(
            ax,
            "Normalized Specific Surface Area",
            "Normalized Effective Thermal Conductivity",
        )

        main_title = "Normalized Effective Thermal Conductivity\nvs. Normalized Specific Surface Area"
        subtitle = _fmt_temp_phi(run.get("tempC"), run.get("porosity"))
        ax.set_title(f"{main_title}\n{subtitle}" if subtitle else main_title)

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(False)

        out_path = os.path.join(outdir, f"{which}_vs_ssa_norm_{suffix}.png")
        _save(fig, out_path, dpi)

    if theme in ("light", "both"):
        emit("light", "light")
    if theme in ("dark", "both"):
        emit("dark", "dark")

def plot_overlay_vs_ssa_norm(run: Dict[str, Any], outdir: str, theme: str, dpi: int):
    if run.get("ssa") is None:
        return

    def emit(theme_name: str, suffix: str):
        set_theme(theme_name)
        fig, ax = plt.subplots(figsize=(FIG_W_IN, FIG_H_IN))

        ssa = np.asarray(run["ssa"], dtype=float)
        kxx = np.asarray(run["kxx"], dtype=float)
        kyy = np.asarray(run["kyy"], dtype=float)

        valid = np.isfinite(ssa) & np.isfinite(kxx) & np.isfinite(kyy)
        if np.count_nonzero(valid) < 2:
            plt.close(fig)
            return

        ssa = ssa[valid]
        kxx = kxx[valid]
        kyy = kyy[valid]

        if ssa[0] == 0 or kxx[0] == 0 or kyy[0] == 0:
            plt.close(fig)
            return

        ssa_n = ssa / ssa[0]
        kxx_n = kxx / kxx[0]
        kyy_n = kyy / kyy[0]

        if LOG_PLOTS:
            ax.set_xscale("log")
            ax.set_yscale("log")
            mxx = (ssa_n > 0) & (kxx_n > 0)
            myy = (ssa_n > 0) & (kyy_n > 0)
            ax.plot(ssa_n[mxx], kxx_n[mxx], lw=2.2, label=r"$k_{xx}$")
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
        subtitle = _fmt_temp_phi(run.get("tempC"), run.get("porosity"))
        ax.set_title(f"{main_title}\n{subtitle}" if subtitle else main_title)

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(False)
        ax.legend(loc="best", fontsize=max(8, LEGEND_FONTSIZE - 1))

        out_path = os.path.join(outdir, f"k_eff_overlay_vs_ssa_norm_{suffix}.png")
        _save(fig, out_path, dpi)

    if theme in ("light", "both"):
        emit("light", "light")
    if theme in ("dark", "both"):
        emit("dark", "dark")

def main():
    run = load_single_run(RUN_DIR, normalize=NORMALIZE)

    outdir = os.path.join(run["run_dir"], "output_plots_single_run")
    os.makedirs(outdir, exist_ok=True)

    # time plots
    plot_vs_time(run, outdir, THEME, DPI, "kxx")
    plot_vs_time(run, outdir, THEME, DPI, "kyy")
    plot_overlay_vs_time(run, outdir, THEME, DPI)

    # SSA plots (normalized)
    plot_vs_ssa_norm(run, outdir, THEME, DPI, "kxx")
    plot_vs_ssa_norm(run, outdir, THEME, DPI, "kyy")
    plot_overlay_vs_ssa_norm(run, outdir, THEME, DPI)

    print("\n✅ Done.")

if __name__ == "__main__":
    main()
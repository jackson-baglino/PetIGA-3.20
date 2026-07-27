#!/usr/bin/env python3
"""
analyze_interface.py  —  Quantify ice-air interface geometry over time.

Works for both 2D and 1D runs.  For each solution snapshot it computes:

  1D metrics (from 1D sol_*.dat):
    - Interface position(s)  [m]  — φ_i = 0.5 crossing(s)
    - Interface width  [m]       — distance between φ_i = 0.1 and φ_i = 0.9 crossings
    - Total ice volume fraction
    - Peak vapor supersaturation  max((ρ_v − ρ_vs)/ρ_vs)

  2D metrics (from SSA_evo.dat, available for any dimension):
    - Interface density  Σ/ε  (proxy for specific surface area)
    - Total ice volume
    - Triple-junction density  (ice·air·sediment overlap)
    - Ice-air interface area  = (Σ/ε) × ε

Results are saved to a CSV and optionally plotted.

Usage
-----
  # 1D analysis
  python analyze_interface.py --dir /path/to/1D/run --dim 1

  # 2D/3D scalar analysis (from SSA_evo.dat)
  python analyze_interface.py --dir /path/to/2D/run --dim 2

  # Save CSV + figure
  python analyze_interface.py --dir . --dim 1 --save-csv metrics.csv --save-fig metrics.png
"""

import argparse
import os
import sys
import glob

import numpy as np
import matplotlib.pyplot as plt

try:
    from igakit.io import PetIGA
except ImportError:
    PetIGA = None   # only needed for 1D field analysis


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def rho_vs(T_C: np.ndarray) -> np.ndarray:
    """Saturation vapor density over ice [kg/m³]."""
    return 3.25e-3 * np.exp(-6150.0 / (T_C + 273.15))


def _crossings(x: np.ndarray, y: np.ndarray, level: float = 0.5) -> np.ndarray:
    """Linear-interpolation crossings of y=level in array y(x)."""
    out = []
    for i in range(len(y) - 1):
        y0, y1 = y[i], y[i + 1]
        if (y0 - level) * (y1 - level) < 0:
            t = (level - y0) / (y1 - y0)
            out.append(x[i] + t * (x[i + 1] - x[i]))
    return np.array(out)


def _step_number(path: str) -> int:
    base = os.path.splitext(os.path.basename(path))[0]
    digits = base.lstrip("abcdefghijklmnopqrstuvwxyz_")
    return int(digits) if digits else 0


# ---------------------------------------------------------------------------
# 1D analysis
# ---------------------------------------------------------------------------

def analyze_1D(run_dir: str, eps: float = None,
               iga_file: str = "igasol.dat") -> dict:
    """
    Read all sol_*.dat files from a 1D run and return a dict of time-series
    numpy arrays.

    eps : interface width parameter (used to compute normalised interface
          density if provided).  If None, raw widths are returned.
    """
    if PetIGA is None:
        sys.exit("ERROR: igakit required for 1D analysis. pip install igakit")

    iga_path = os.path.join(run_dir, iga_file)
    if not os.path.isfile(iga_path):
        sys.exit(f"IGA geometry file not found: {iga_path}")

    nrb  = PetIGA().read(iga_path)
    ctrl = nrb.control
    x    = ctrl[:, 0]   # physical x-coords [m]

    sol_files = sorted(glob.glob(os.path.join(run_dir, "sol_*.dat")),
                       key=_step_number)
    if not sol_files:
        sys.exit(f"No sol_*.dat files found in '{run_dir}'")

    # Load time labels from SSA_evo if available
    times = _load_times(run_dir)

    records = []
    for sf in sol_files:
        try:
            sol = PetIGA().read_vec(sf, nrb)
        except Exception as e:
            print(f"  WARNING: skipping {sf}: {e}")
            continue
        if sol.ndim < 2:
            continue

        step = _step_number(sf)
        t    = times.get(step, np.nan)

        ice  = sol[:, 0]
        tem  = sol[:, 1]
        rhov = sol[:, 2]
        sed  = np.clip(sol[:, 3], 0.0, 1.0) if sol.shape[1] > 3 else np.zeros_like(ice)

        # Volume fractions (trapezoidal)
        Lx        = x[-1] - x[0]
        vol_frac  = np.trapz(ice, x) / Lx
        sed_frac  = np.trapz(sed, x) / Lx

        # Interface positions (φ_i = 0.5 crossings)
        c05  = _crossings(x, ice, 0.5)
        c01  = _crossings(x, ice, 0.1)
        c09  = _crossings(x, ice, 0.9)

        # Width: distance between outermost 0.1 and 0.9 crossings
        if len(c09) >= 1 and len(c01) >= 1:
            width_left  = abs(c01[0]  - c09[0])
            width_right = abs(c09[-1] - c01[-1])
            width_mean  = np.nanmean([width_left, width_right])
        else:
            width_left = width_right = width_mean = np.nan

        # Interface count (number of 0.5 crossings)
        n_intf = len(c05)

        # Vapor supersaturation: peak and mean over air region (air = 1 - ice - sed)
        air       = np.clip(1.0 - ice - sed, 0.0, 1.0)
        rvs       = rho_vs(tem)
        supersat  = np.where(rvs > 0, (rhov - rvs) / rvs, 0.0)
        ss_peak   = float(np.max(supersat))
        ss_mean   = float(np.trapz(supersat * air, x) / max(np.trapz(air, x), 1e-30))

        # Temperature range
        T_min = float(np.min(tem))
        T_max = float(np.max(tem))

        records.append({
            "step":            step,
            "t_s":             t,
            "vol_frac":        vol_frac,
            "sed_frac":        sed_frac,
            "n_interfaces":    n_intf,
            "x_left_mm":       c05[0]  * 1e3 if len(c05) >= 1 else np.nan,
            "x_right_mm":      c05[-1] * 1e3 if len(c05) >= 2 else np.nan,
            "width_left_um":   width_left  * 1e6,
            "width_right_um":  width_right * 1e6,
            "width_mean_um":   width_mean  * 1e6,
            "ss_peak":         ss_peak,
            "ss_mean_air":     ss_mean,
            "T_min_C":         T_min,
            "T_max_C":         T_max,
        })

    if not records:
        sys.exit("No valid solution files could be processed.")

    # Convert list-of-dicts to dict-of-arrays
    keys = list(records[0].keys())
    return {k: np.array([r[k] for r in records]) for k in keys}


def _load_times(run_dir: str) -> dict:
    """Return {step: t_s} from SSA_evo.dat if available."""
    path = os.path.join(run_dir, "SSA_evo.dat")
    times = {}
    if not os.path.isfile(path):
        return times
    try:
        data = np.genfromtxt(path, dtype=float, comments="#", invalid_raise=False)
        if data.ndim == 1:
            data = data[np.newaxis, :]
        if data.shape[1] >= 4:
            for row in data:
                times[int(row[3])] = float(row[2])
    except Exception:
        pass
    return times


# ---------------------------------------------------------------------------
# 2D / generic scalar analysis (from SSA_evo.dat)
# ---------------------------------------------------------------------------

def analyze_scalars(run_dir: str, eps: float = 9.3295e-7) -> dict:
    """
    Parse SSA_evo.dat and compute derived metrics.
    """
    path = os.path.join(run_dir, "SSA_evo.dat")
    if not os.path.isfile(path):
        sys.exit(f"SSA_evo.dat not found in '{run_dir}'")

    data = np.genfromtxt(path, dtype=float, comments="#", invalid_raise=False)
    if data.ndim == 1:
        data = data[np.newaxis, :]
    data = data[~np.isnan(data).any(axis=1)]

    return {
        "step":          data[:, 3].astype(int),
        "t_s":           data[:, 2],
        "t_h":           data[:, 2] / 3600.0,
        "intf_density":  data[:, 0],         # Σ/ε
        "intf_area":     data[:, 0] * eps,    # Σ  [m or m²]
        "ice_vol":       data[:, 1],
        "d_ice_vol":     data[:, 1] - data[0, 1],   # change from initial
    }


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_1D_metrics(metrics: dict, save_path: str = None, title: str = ""):
    t = metrics["t_s"]
    t_h = t / 3600.0
    valid = ~np.isnan(t_h)
    t_h = t_h[valid]

    fig, axes = plt.subplots(4, 2, figsize=(12, 13))

    def _plot(ax, y, ylabel, color="tab:blue"):
        yv = y[valid]
        ax.plot(t_h, yv, "o-", color=color, ms=4, lw=2)
        ax.set_ylabel(ylabel, fontsize=11)
        ax.grid(True, alpha=0.3)
        ax.tick_params(labelsize=10)

    _plot(axes[0, 0], metrics["vol_frac"],       r"Mean ice fraction $\langle\phi_i\rangle$",       "tab:blue")
    _plot(axes[0, 1], metrics["sed_frac"],        r"Mean sediment fraction $\langle\phi_s\rangle$",  "tab:brown")
    _plot(axes[1, 0], metrics["n_interfaces"],   "Number of interfaces",                             "tab:orange")
    _plot(axes[1, 1], metrics["ss_peak"],        "Peak supersaturation",                             "tab:red")
    _plot(axes[2, 0], metrics["x_left_mm"],      "Left interface pos.  [mm]",                        "tab:green")
    _plot(axes[2, 1], metrics["x_right_mm"],     "Right interface pos. [mm]",                        "tab:red")
    _plot(axes[3, 0], metrics["width_mean_um"],  "Mean interface width  [μm]",                       "tab:purple")
    _plot(axes[3, 1], metrics["T_min_C"],        r"$T_{\min}$  [°C]",                               "tab:cyan")

    for row in axes:
        row[-1].set_xlabel("Time  [h]", fontsize=11)
        row[0].set_xlabel("Time  [h]", fontsize=11)

    fig.suptitle(f"1D interface metrics{' — ' + title if title else ''}",
                 fontsize=13, y=1.01)
    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Figure saved to: {save_path}")
    else:
        plt.show()


def plot_2D_metrics(metrics: dict, save_path: str = None, title: str = ""):
    t_h = metrics["t_h"]

    fig, axes = plt.subplots(2, 2, figsize=(11, 8))

    def _plot(ax, y, ylabel, color):
        ax.plot(t_h, y, "o-", color=color, ms=4, lw=2)
        ax.set_ylabel(ylabel, fontsize=11)
        ax.set_xlabel("Time  [h]", fontsize=11)
        ax.grid(True, alpha=0.3)

    _plot(axes[0, 0], metrics["ice_vol"],     r"$\int \phi_i \, dV$",       "tab:blue")
    _plot(axes[0, 1], metrics["d_ice_vol"],   r"$\Delta \int \phi_i \, dV$","tab:red")
    _plot(axes[1, 0], metrics["intf_density"],r"Interface density $\Sigma/\varepsilon$","tab:green")
    _plot(axes[1, 1], metrics["intf_area"],   r"Interface area $\Sigma$  [m or m²]",   "tab:purple")

    axes[0, 1].axhline(0, color="gray", ls="--", lw=1)

    fig.suptitle(f"Scalar metrics{' — ' + title if title else ''}",
                 fontsize=13, y=1.01)
    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Figure saved to: {save_path}")
    else:
        plt.show()


# ---------------------------------------------------------------------------
# CSV export
# ---------------------------------------------------------------------------

def save_csv(metrics: dict, csv_path: str):
    keys = list(metrics.keys())
    n    = len(metrics[keys[0]])
    header = ",".join(keys)
    rows   = np.column_stack([metrics[k] for k in keys])
    np.savetxt(csv_path, rows, delimiter=",", header=header, comments="")
    print(f"CSV saved to: {csv_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Quantify ice-air interface geometry and scalar evolution."
    )
    p.add_argument("--dir",      default=".",
                   help="Output directory containing sol_*.dat and/or SSA_evo.dat")
    p.add_argument("--dim",      type=int, default=2, choices=[1, 2, 3],
                   help="Simulation dimension (1 uses sol_*.dat; 2/3 use SSA_evo.dat)")
    p.add_argument("--eps",      type=float, default=9.3295e-7,
                   help="Interface width ε [m] (used for interface area calculation)")
    p.add_argument("--iga",      default="igasol.dat",
                   help="IGA geometry file (1D only)")
    p.add_argument("--save-csv", default=None,
                   help="Save metrics to this CSV file")
    p.add_argument("--save-fig", default=None,
                   help="Save figure to this path (omit to display)")
    p.add_argument("--title",    default="",
                   help="Figure title suffix")
    return p.parse_args()


def main():
    args = parse_args()

    if args.dim == 1:
        metrics = analyze_1D(args.dir, eps=args.eps, iga_file=args.iga)
        if args.save_csv:
            save_csv(metrics, args.save_csv)
        plot_1D_metrics(metrics, save_path=args.save_fig, title=args.title)
    else:
        metrics = analyze_scalars(args.dir, eps=args.eps)
        if args.save_csv:
            save_csv(metrics, args.save_csv)
        plot_2D_metrics(metrics, save_path=args.save_fig, title=args.title)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
plot1D_profiles.py  —  Visualize 1D permafrost simulation field profiles.

Reads sol_*.dat files from a PetIGA 1D run and plots spatial profiles of:
  - Ice phase field  φ_i(x)
  - Temperature      T(x)
  - Vapor density    ρ_v(x)

Multiple time snapshots are overlaid in a single figure per field.

Usage
-----
  # Plot all snapshots in current directory
  python plot1D_profiles.py

  # Specify output directory
  python plot1D_profiles.py --dir /path/to/output

  # Limit to first 5 snapshots; save figure instead of showing
  python plot1D_profiles.py --dir . --max-steps 5 --save profiles.png
"""

import argparse
import glob
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# ---------------------------------------------------------------------------
# Try importing igakit; give a clear error if missing
# ---------------------------------------------------------------------------
try:
    from igakit.io import PetIGA
except ImportError:
    sys.exit(
        "ERROR: igakit is not installed in the active Python environment.\n"
        "Install it with:  pip install igakit"
    )


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_geometry(iga_path: str):
    if not os.path.isfile(iga_path):
        sys.exit(f"ERROR: IGA geometry file not found: {iga_path}")
    return PetIGA().read(iga_path)


def load_solution(sol_path: str, nrb):
    return PetIGA().read_vec(sol_path, nrb)


def get_x_coords(nrb) -> np.ndarray:
    """Return physical x-coordinates of control points (1D uniform mesh)."""
    ctrl = nrb.control  # shape (n, 2) for 1D: [x, w]
    return ctrl[:, 0]


def step_number(path: str) -> int:
    """Extract integer step number from sol_NNNNN.dat filename."""
    base = os.path.splitext(os.path.basename(path))[0]  # 'sol_00010'
    digits = base.lstrip("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_")
    return int(digits) if digits else 0


# ---------------------------------------------------------------------------
# Saturation vapor density (Clausius-Clapeyron, matches material_properties.c)
# ---------------------------------------------------------------------------

def rho_vs(T_C: np.ndarray) -> np.ndarray:
    """
    Saturation vapor density over ice [kg/m³] as a function of temperature.
    Uses the same formula as RhoVS_I in material_properties.c:
        ρ_vs = A * exp(B / (T_C + 273.15))
    with A = 3.25e-3 kg/m³, B = -6150 K  (ice-vapor Clausius-Clapeyron).
    """
    T_K = T_C + 273.15
    return 3.25e-3 * np.exp(-6150.0 / T_K)


# ---------------------------------------------------------------------------
# Main plotting routine
# ---------------------------------------------------------------------------

def plot_profiles(run_dir: str, max_steps: int = None, save_path: str = None,
                  iga_file: str = "igasol.dat"):
    """
    Read all solution files in run_dir and produce a 3-panel profile plot.
    """
    iga_path = os.path.join(run_dir, iga_file)
    nrb = load_geometry(iga_path)
    x   = get_x_coords(nrb)  # (n,)

    # Collect and sort solution files
    pattern = os.path.join(run_dir, "sol_*.dat")
    sol_files = sorted(glob.glob(pattern), key=step_number)
    if not sol_files:
        sys.exit(f"No sol_*.dat files found in '{run_dir}'")

    if max_steps is not None:
        sol_files = sol_files[:max_steps]

    n_steps  = len(sol_files)
    colormap = cm.viridis(np.linspace(0, 1, n_steps))

    print(f"Found {n_steps} solution file(s) in '{run_dir}'")

    # ------------------------------------------------------------------
    # Read SSA_evo.dat for time labels (optional)
    # ------------------------------------------------------------------
    times = {}  # step -> time (s)
    ssa_path = os.path.join(run_dir, "SSA_evo.dat")
    if os.path.isfile(ssa_path):
        try:
            ssa = np.loadtxt(ssa_path)
            if ssa.ndim == 1:
                ssa = ssa[np.newaxis, :]
            # columns: sub_interf/eps  tot_ice  t  step
            for row in ssa:
                if row.shape[0] >= 4:
                    times[int(row[3])] = row[2]
        except Exception:
            pass  # time labels are optional

    # ------------------------------------------------------------------
    # Load all solutions
    # ------------------------------------------------------------------
    ice_snaps  = []
    tem_snaps  = []
    rhov_snaps = []
    step_labels = []

    for sf in sol_files:
        try:
            sol = load_solution(sf, nrb)  # shape (nx, 3) for 1D
        except Exception as e:
            print(f"  WARNING: Could not read {sf}: {e}")
            continue

        if sol.ndim == 1:
            print(f"  WARNING: {sf} returned 1D array — skipping.")
            continue

        ice_snaps.append(sol[:, 0])
        tem_snaps.append(sol[:, 1])
        rhov_snaps.append(sol[:, 2])

        step = step_number(sf)
        t    = times.get(step)
        if t is not None:
            label = f"step {step}  t={_fmt_time(t)}"
        else:
            label = f"step {step}"
        step_labels.append(label)

    if not ice_snaps:
        sys.exit("No valid solution files could be read.")

    n_valid = len(ice_snaps)

    # ------------------------------------------------------------------
    # Figure: 3 rows × 1 column
    # ------------------------------------------------------------------
    fig, axes = plt.subplots(3, 1, figsize=(9, 10), sharex=True)

    x_mm = x * 1e3  # convert m → mm for display

    for i in range(n_valid):
        c = cm.viridis(i / max(n_valid - 1, 1))
        lw = 1.5 if i not in (0, n_valid - 1) else 2.0

        axes[0].plot(x_mm, ice_snaps[i],  color=c, lw=lw, label=step_labels[i])
        axes[1].plot(x_mm, tem_snaps[i],  color=c, lw=lw)
        axes[2].plot(x_mm, rhov_snaps[i], color=c, lw=lw)

    # Also plot saturation vapor density at the mean temperature
    if tem_snaps:
        T_mean = np.mean(tem_snaps[-1])
        rhoVS_mean = rho_vs(T_mean)
        axes[2].axhline(rhoVS_mean, color="k", ls="--", lw=1.0,
                        label=rf"$\rho_{{vs}}(T_{{mean}}={T_mean:.1f}°C)$")

    # Formatting
    axes[0].set_ylabel(r"Ice phase  $\phi_i$",      fontsize=13)
    axes[1].set_ylabel("Temperature  [°C]",          fontsize=13)
    axes[2].set_ylabel(r"Vapor density  [kg m$^{-3}$]", fontsize=13)
    axes[2].set_xlabel("x  [mm]",                    fontsize=13)

    axes[0].set_ylim(-0.05, 1.1)
    axes[0].set_title(f"1D Permafrost — field profiles\n(dir: {run_dir})", fontsize=13)
    axes[0].legend(fontsize=7, loc="upper right", ncol=2)
    axes[2].legend(fontsize=9, loc="upper right")

    for ax in axes:
        ax.grid(True, alpha=0.3)
        ax.tick_params(labelsize=11)

    # Colorbar as time axis (only if > 1 snapshot)
    if n_valid > 1:
        sm = plt.cm.ScalarMappable(cmap="viridis",
                                    norm=plt.Normalize(0, n_valid - 1))
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=axes, orientation="vertical",
                            fraction=0.02, pad=0.02)
        cbar.set_label("Snapshot index (0 = earliest)", fontsize=10)

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Figure saved to: {save_path}")
    else:
        plt.show()


# ---------------------------------------------------------------------------
# Auxiliary: derived quantities plot (interface position + ice volume)
# ---------------------------------------------------------------------------

def plot_derived(run_dir: str, save_path: str = None,
                 iga_file: str = "igasol.dat"):
    """
    Plot derived 1D scalar quantities from all snapshots:
      - Total ice volume fraction vs. snapshot index
      - Interface position(s) vs. snapshot index (defined as φ_i = 0.5 contour)
      - Ice slab width vs. snapshot index
    """
    iga_path = os.path.join(run_dir, iga_file)
    nrb = load_geometry(iga_path)
    x   = get_x_coords(nrb)
    Lx  = x[-1] - x[0]

    pattern  = os.path.join(run_dir, "sol_*.dat")
    sol_files = sorted(glob.glob(pattern), key=step_number)
    if not sol_files:
        sys.exit(f"No sol_*.dat files found in '{run_dir}'")

    steps         = []
    ice_volumes   = []
    left_edges    = []
    right_edges   = []
    slab_widths   = []

    for sf in sol_files:
        try:
            sol = load_solution(sf, nrb)
        except Exception:
            continue
        if sol.ndim < 2:
            continue

        step = step_number(sf)
        ice  = sol[:, 0]

        # Total ice volume fraction (trapz integration, normalised by Lx)
        vol = np.trapz(ice, x) / Lx
        ice_volumes.append(vol)
        steps.append(step)

        # Interface positions: interpolate φ_i = 0.5 crossings
        crossings = _find_crossings(x, ice, level=0.5)
        if len(crossings) >= 2:
            left_edges.append(crossings[0] * 1e3)
            right_edges.append(crossings[-1] * 1e3)
            slab_widths.append((crossings[-1] - crossings[0]) * 1e3)
        elif len(crossings) == 1:
            left_edges.append(crossings[0] * 1e3)
            right_edges.append(np.nan)
            slab_widths.append(np.nan)
        else:
            left_edges.append(np.nan)
            right_edges.append(np.nan)
            slab_widths.append(np.nan)

    if not steps:
        sys.exit("No valid solutions found.")

    steps       = np.array(steps)
    ice_volumes = np.array(ice_volumes)
    left_edges  = np.array(left_edges)
    right_edges = np.array(right_edges)
    slab_widths = np.array(slab_widths)

    fig, axes = plt.subplots(3, 1, figsize=(8, 9), sharex=True)

    axes[0].plot(steps, ice_volumes, "o-", color="#1f77b4", ms=4)
    axes[0].set_ylabel("Mean ice fraction  $\\langle\\phi_i\\rangle$", fontsize=13)
    axes[0].axhline(ice_volumes[0], color="gray", ls="--", lw=1, label="initial")
    axes[0].legend(fontsize=10)

    axes[1].plot(steps, left_edges,  "s-", color="#2ca02c", ms=4, label="left edge")
    axes[1].plot(steps, right_edges, "D-", color="#d62728", ms=4, label="right edge")
    axes[1].set_ylabel("Interface position  [mm]", fontsize=13)
    axes[1].legend(fontsize=10)

    axes[2].plot(steps, slab_widths, "^-", color="#9467bd", ms=4)
    axes[2].set_ylabel("Slab width  [mm]", fontsize=13)
    axes[2].set_xlabel("Step index", fontsize=13)

    axes[0].set_title(f"1D Permafrost — derived quantities\n(dir: {run_dir})", fontsize=13)

    for ax in axes:
        ax.grid(True, alpha=0.3)
        ax.tick_params(labelsize=11)

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Figure saved to: {save_path}")
    else:
        plt.show()


# ---------------------------------------------------------------------------
# Small utilities
# ---------------------------------------------------------------------------

def _find_crossings(x: np.ndarray, y: np.ndarray, level: float = 0.5) -> np.ndarray:
    """Return x positions where y crosses `level` (linear interpolation)."""
    crossings = []
    for i in range(len(y) - 1):
        y0, y1 = y[i], y[i + 1]
        if (y0 - level) * (y1 - level) < 0:
            # linear interpolation
            t = (level - y0) / (y1 - y0)
            crossings.append(x[i] + t * (x[i + 1] - x[i]))
    return np.array(crossings)


def _fmt_time(t_sec: float) -> str:
    """Human-readable time string."""
    if t_sec < 60:
        return f"{t_sec:.1f} s"
    elif t_sec < 3600:
        return f"{t_sec / 60:.1f} min"
    else:
        return f"{t_sec / 3600:.2f} h"


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description="Plot 1D permafrost field profiles.")
    p.add_argument("--dir",        default=".",     help="Output directory with sol_*.dat files")
    p.add_argument("--iga",        default="igasol.dat", help="IGA geometry file")
    p.add_argument("--max-steps",  type=int, default=None, help="Maximum snapshots to load")
    p.add_argument("--save",       default=None, help="Save figure to this path (omit to display)")
    p.add_argument("--derived",    action="store_true",
                   help="Plot derived quantities (ice volume, interface positions) instead")
    return p.parse_args()


def main():
    args = parse_args()
    if args.derived:
        plot_derived(args.dir, save_path=args.save, iga_file=args.iga)
    else:
        plot_profiles(args.dir, max_steps=args.max_steps,
                      save_path=args.save, iga_file=args.iga)


if __name__ == "__main__":
    main()

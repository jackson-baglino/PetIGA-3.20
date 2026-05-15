#!/usr/bin/env python3
"""
plot_mass.py  —  Plot phase mass vs. time for a permafrost simulation.

Reads sol_*.dat solution snapshots and computes the integrated mass of each
phase over the simulation domain at every output step:

  mass_ice(t)   = rho_ice × ∫ φ_i(x,t) dV
  mass_sed(t)   = rho_sed × ∫ φ_s(x,t) dV
  mass_vap(t)   =           ∫ ρ_v(x,t) · (1 − φ_i − φ_s) dV
  mass_total(t) = mass_ice + mass_sed + mass_vap

A vertical dashed line is drawn at t = t_sed_freeze (the moment the model
switches from the full three-phase formulation to the frozen-sediment two-phase
formulation). t_sed_freeze is read from the .opts files in the run folder.

Unit notes
----------
  1D run: masses are *per unit cross-sectional area*  [kg/m²]
  2D run: masses are *per unit depth*                 [kg/m]
  3D run: masses are in kg

Physical constants (rho_ice, rho_sed) must match src/permafrost2.c.

Usage
-----
  # From inside a run folder:
  python postprocess/plot_mass.py

  # Explicit path:
  python postprocess/plot_mass.py --dir /path/to/run --save /path/to/mass.png

  # Choose time unit:
  python postprocess/plot_mass.py --time-unit h
"""

import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    from igakit.io import PetIGA
except ImportError:
    sys.exit(
        "ERROR: igakit is not installed.\n"
        "Install with:  pip install igakit"
    )

# ---------------------------------------------------------------------------
# Physical constants — must match src/permafrost2.c
# ---------------------------------------------------------------------------
RHO_ICE = 919.0    # kg/m³  — density of ice
RHO_SED = 7753.0   # kg/m³  — density of sediment (metal placeholder; update if needed)

# ---------------------------------------------------------------------------
# Plot colours (consistent with existing postprocess scripts)
# ---------------------------------------------------------------------------
COLOR_TOTAL = "black"
COLOR_ICE   = "#2166ac"   # steel blue
COLOR_SED   = "#8c510a"   # warm brown
COLOR_VAP   = "#7b3294"   # purple

# ---------------------------------------------------------------------------
# Time-unit helpers
# ---------------------------------------------------------------------------
TIME_SCALES = {
    "s":   (1.0,     "Time  [s]"),
    "min": (60.0,    "Time  [min]"),
    "h":   (3600.0,  "Time  [h]"),
    "d":   (86400.0, "Time  [days]"),
}


def auto_time_unit(t_max_sec: float) -> str:
    """Choose a sensible time unit for the x-axis."""
    if t_max_sec <= 600:
        return "s"
    if t_max_sec <= 7200:
        return "min"
    if t_max_sec <= 3 * 86400:
        return "h"
    return "d"


# ---------------------------------------------------------------------------
# Opts-file helpers
# ---------------------------------------------------------------------------

def find_opts_files(run_dir: str):
    """Return opts files in run_dir: [universal.opts, sim-specific.opts].
    PETSc processes left-to-right so later entries override earlier ones."""
    result = []
    u = os.path.join(run_dir, "universal.opts")
    if os.path.isfile(u):
        result.append(u)
    for f in sorted(os.listdir(run_dir)):
        if f.endswith(".opts") and f != "universal.opts":
            result.append(os.path.join(run_dir, f))
    return result


def parse_opts_float(opts_files, key: str, default=None):
    """Return the last value of -key found across opts_files (sim file wins)."""
    value = default
    for path in opts_files:
        try:
            with open(path) as fh:
                for line in fh:
                    tokens = line.split()
                    if len(tokens) >= 2 and tokens[0] == key:
                        value = float(tokens[1])
        except (OSError, ValueError):
            pass
    return value


# ---------------------------------------------------------------------------
# SSA_evo.dat loader (duplicated from plot_scalars.py for standalone use)
# ---------------------------------------------------------------------------

def load_ssa(path: str) -> np.ndarray:
    """Load SSA_evo.dat → (N, 5) array [sub_interf/eps, tot_ice, t, step, dt]."""
    if not os.path.isfile(path):
        return None
    try:
        data = np.genfromtxt(path, dtype=float, comments="#", invalid_raise=False)
    except Exception:
        return None
    if data.ndim == 1:
        data = data[np.newaxis, :]
    if data.shape[1] < 4:
        return None
    if data.shape[1] == 4:
        data = np.hstack([data, np.full((len(data), 1), np.nan)])
    mask = ~np.isnan(data[:, :4]).any(axis=1)
    data = data[mask]
    if len(data) == 0:
        return None
    # Deduplicate by step number
    steps = data[:, 3].astype(int)
    _, last_idx = np.unique(steps[::-1], return_index=True)
    keep = np.sort(len(steps) - 1 - last_idx)
    return data[keep]


# ---------------------------------------------------------------------------
# Core computation
# ---------------------------------------------------------------------------

def compute_masses(run_dir: str):
    """
    Read all sol_*.dat files in run_dir and return:
      times      [s]       — 1-D array, length N
      mass_ice   [kg/m^k]  — 1-D array, length N   (k = 3-dim)
      mass_sed   [kg/m^k]
      mass_vap   [kg/m^k]
      dim        int        — spatial dimension of the run
    """
    geo_file = os.path.join(run_dir, "igasol.dat")
    if not os.path.isfile(geo_file):
        sys.exit(f"ERROR: IGA geometry file not found: {geo_file}")

    sol_files = sorted(
        os.path.join(run_dir, f)
        for f in os.listdir(run_dir)
        if f.startswith("sol_") and f.endswith(".dat")
    )
    if not sol_files:
        sys.exit(f"ERROR: No sol_*.dat files found in {run_dir}")

    # ── IGA geometry ──────────────────────────────────────────────────────
    nrb = PetIGA().read(geo_file)
    dim = nrb.dim

    # Physical domain extent and cell volume
    # nrb.points has shape (*shape, sdim) where sdim >= dim
    ctrl = nrb.points[..., :dim]               # physical coordinates
    flat = ctrl.reshape(-1, dim)
    extents = flat.max(axis=0) - flat.min(axis=0)
    V_domain = float(np.prod(extents))

    shape = nrb.points.shape[:-1]              # (Nx,) or (Nx, Ny) or (Nx, Ny, Nz)
    N_total = int(np.prod(shape))
    dV = V_domain / N_total                    # approximate cell volume per ctrl point

    # ── Time axis from SSA_evo.dat ────────────────────────────────────────
    ssa = load_ssa(os.path.join(run_dir, "SSA_evo.dat"))

    n_sol = len(sol_files)
    if ssa is not None:
        ssa_times = ssa[:, 2]
        n_ssa = len(ssa_times)
        if n_sol == n_ssa + 1:
            # sol_00000 is the IC at t=0; SSA rows correspond to sol_00001 onward
            times = np.concatenate([[0.0], ssa_times])
        elif n_sol == n_ssa:
            times = ssa_times
        else:
            # Mismatch: use what we can
            n = min(n_sol, n_ssa)
            times = ssa_times[:n]
            sol_files = sol_files[:n]
    else:
        # No SSA file: assume uniform dt, reconstruct from opts
        times = np.arange(n_sol, dtype=float)  # step index as fallback

    # ── Integrate each snapshot ───────────────────────────────────────────
    n = len(times)
    mass_ice = np.zeros(n)
    mass_sed = np.zeros(n)
    mass_vap = np.zeros(n)

    for k, sf in enumerate(sol_files[:n]):
        sol = PetIGA().read_vec(sf, nrb)
        sol = sol.reshape(-1, 4)               # ensure (N_total, 4) even for 2D

        phi_i = sol[:, 0].clip(0.0, 1.0)
        phi_s = sol[:, 3].clip(0.0, 1.0)
        rho_v = sol[:, 2].clip(0.0)
        phi_a = (1.0 - phi_i - phi_s).clip(0.0, 1.0)

        mass_ice[k] = RHO_ICE * np.sum(phi_i) * dV
        mass_sed[k] = RHO_SED * np.sum(phi_s) * dV
        mass_vap[k] = np.sum(rho_v * phi_a) * dV

    return times, mass_ice, mass_sed, mass_vap, dim


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_mass(run_dir: str, time_unit: str = None, save_path: str = None):

    times, mass_ice, mass_sed, mass_vap, dim = compute_masses(run_dir)
    mass_total = mass_ice + mass_sed + mass_vap

    # Auto time unit if not given
    t_max = times[-1] if len(times) > 0 else 1.0
    unit = time_unit or auto_time_unit(t_max)
    scale, xlabel = TIME_SCALES.get(unit, (1.0, "Time  [s]"))
    t_plot = times / scale

    # Parse t_sed_freeze
    opts_files = find_opts_files(run_dir)
    t_freeze = parse_opts_float(opts_files, "-t_sed_freeze")

    # Mass units label
    unit_suffix = {1: "kg m⁻²", 2: "kg m⁻¹", 3: "kg"}.get(dim, "kg")

    # ── Figure ────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(8, 5))

    ax.plot(t_plot, mass_total, color=COLOR_TOTAL, lw=2.2, zorder=4,
            label="Total")
    ax.plot(t_plot, mass_ice,   color=COLOR_ICE,   lw=1.6, zorder=3,
            label=r"Ice  ($\rho_i \int \phi_i \, dV$)")
    ax.plot(t_plot, mass_sed,   color=COLOR_SED,   lw=1.6, zorder=3,
            label=r"Sediment  ($\rho_s \int \phi_s \, dV$)")
    ax.plot(t_plot, mass_vap,   color=COLOR_VAP,   lw=1.6, zorder=3,
            label=r"Vapor  ($\int \rho_v \, \phi_a \, dV$)")

    # Dashed semi-opaque reference lines at the initial (t=0) mass values.
    # These make it easy to see whether each phase is gaining or losing mass.
    ax.axhline(mass_total[0], color=COLOR_TOTAL, lw=1.0, ls="--", alpha=0.45, zorder=2)
    ax.axhline(mass_ice[0],   color=COLOR_ICE,   lw=1.0, ls="--", alpha=0.45, zorder=2)
    ax.axhline(mass_sed[0],   color=COLOR_SED,   lw=1.0, ls="--", alpha=0.45, zorder=2)
    ax.axhline(mass_vap[0],   color=COLOR_VAP,   lw=1.0, ls="--", alpha=0.45, zorder=2)

    # Vertical dashed line at t_sed_freeze
    if t_freeze is not None:
        xf = t_freeze / scale
        if t_plot[0] <= xf <= t_plot[-1] * 1.05:
            ax.axvline(x=xf, color="gray", ls="--", lw=1.2, zorder=2,
                       label=rf"$t_{{\rm freeze}}$ = {t_freeze:.0f} s")
            # Small annotation above the line
            ylim = ax.get_ylim()
            ax.text(xf, ylim[1], r"$t_{\rm freeze}$",
                    ha="center", va="bottom", fontsize=9, color="gray")

    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(f"Mass  [{unit_suffix}]", fontsize=12)
    ax.set_title("Phase mass vs. time", fontsize=13)
    ax.legend(fontsize=10, loc="best")
    ax.grid(True, alpha=0.3)
    ax.tick_params(labelsize=10)
    fig.tight_layout()

    out = save_path or os.path.join(run_dir, "mass.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Figure saved to: {out}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Plot integrated phase mass vs. time (1D/2D/3D)."
    )
    p.add_argument("--dir",       default=".",
                   help="Run folder containing igasol.dat and sol_*.dat (default: .)")
    p.add_argument("--save",      default=None,
                   help="Output PNG path (default: <dir>/mass.png)")
    p.add_argument("--time-unit", default=None,
                   choices=["s", "min", "h", "d"],
                   help="Time unit for x-axis (default: auto)")
    return p.parse_args()


def main():
    args = parse_args()
    run_dir = os.path.abspath(args.dir)
    plot_mass(run_dir, time_unit=args.time_unit, save_path=args.save)


if __name__ == "__main__":
    main()

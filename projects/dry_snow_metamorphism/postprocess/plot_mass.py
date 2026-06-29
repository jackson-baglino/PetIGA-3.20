#!/usr/bin/env python3
"""
plot_mass.py  —  Plot phase mass vs. time for a dry snow metamorphism simulation.

Reads sol_*.dat solution snapshots and computes the integrated mass of each
phase over the simulation domain at every output step:

  mass_ice(t)   = rho_ice × ∫ φ_i(x,t) dV
  mass_vap(t)   =           ∫ ρ_v(x,t) · (1 − φ_i) dV
  mass_total(t) = mass_ice + mass_vap

Unit notes
----------
  1D run: masses are *per unit cross-sectional area*  [kg/m²]
  2D run: masses are *per unit depth*                 [kg/m]
  3D run: masses are in kg

Physical constant (rho_ice) must match src/dry_snow_metamorphism.c.

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

# ---------------------------------------------------------------------------
# Plot colours (consistent with existing postprocess scripts)
# ---------------------------------------------------------------------------
COLOR_TOTAL = "black"
COLOR_ICE   = "#2166ac"   # steel blue
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
# SSA_evo.dat loader (duplicated from plot_scalars.py for standalone use)
# ---------------------------------------------------------------------------

def load_ssa(path: str) -> np.ndarray:
    """Load SSA_evo.dat → (N, >=5) array.

    Columns: [sub_interf/eps, tot_ice, t, step, dt, tot_air, tot_rhov, tot_mass].
    The last three columns (tot_air, tot_rhov, tot_mass) are written directly
    from monitoring.c's IGA-quadrature integrals and are only present in runs
    produced after that change; older runs have just the first 5 columns.
    """
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
    Return:
      times      [s]       — 1-D array, length N
      mass_ice   [kg/m^k]  — 1-D array, length N   (k = 3-dim)
      mass_vap   [kg/m^k]
      dim        int        — spatial dimension of the run

    Preferred path: SSA_evo.dat written by monitoring.c carries tot_ice,
    tot_air, and tot_rhov computed via proper IGA-quadrature integration
    (matching outp.txt's TOTAL_MASS exactly) for *every* time step — so we
    use those directly instead of recomputing from the sparsely-written
    sol_*.dat snapshots. This avoids both (a) a step/time mismatch when
    sol_*.dat is written less often than every step, and (b) a flawed
    uniform-dV approximation that overcounts non-rectangular geometries.

    Fallback (older runs without the extra SSA_evo.dat columns): recompute
    from sol_*.dat with a uniform-dV approximation, matching each snapshot
    to its SSA_evo.dat row by step number (parsed from the filename).
    """
    ssa = load_ssa(os.path.join(run_dir, "SSA_evo.dat"))

    if ssa is not None and ssa.shape[1] >= 8:
        times    = ssa[:, 2]
        tot_ice  = ssa[:, 1]
        tot_rhov = ssa[:, 6]
        mass_ice = RHO_ICE * tot_ice
        mass_vap = tot_rhov

        # dim: read from igasol.dat if available, else assume 2D.
        geo_file = os.path.join(run_dir, "igasol.dat")
        dim = 2
        if os.path.isfile(geo_file):
            dim = PetIGA().read(geo_file).dim

        return times, mass_ice, mass_vap, dim

    # ── Fallback: integrate sol_*.dat with uniform-dV approximation ────────
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

    # ── Time axis: match each sol_NNNNN.dat to its SSA_evo.dat row by step ──
    if ssa is not None:
        ssa_steps = ssa[:, 3].astype(int)
        ssa_times = ssa[:, 2]
        step_to_time = dict(zip(ssa_steps, ssa_times))

        times = []
        matched_files = []
        for sf in sol_files:
            base = os.path.basename(sf)
            digits = "".join(ch for ch in base if ch.isdigit())
            if not digits:
                continue
            step = int(digits)
            if step == 0 and step not in step_to_time:
                t = 0.0
            elif step in step_to_time:
                t = step_to_time[step]
            else:
                continue  # no matching SSA row for this snapshot's step
            times.append(t)
            matched_files.append(sf)
        times = np.array(times, dtype=float)
        sol_files = matched_files
    else:
        # No SSA file: assume uniform dt, reconstruct from opts
        times = np.arange(len(sol_files), dtype=float)  # step index as fallback

    # ── Integrate each snapshot ───────────────────────────────────────────
    n = len(times)
    mass_ice = np.zeros(n)
    mass_vap = np.zeros(n)

    for k, sf in enumerate(sol_files[:n]):
        sol = PetIGA().read_vec(sf, nrb)
        sol = sol.reshape(-1, 3)               # ensure (N_total, 3) even for 2D

        phi_i = sol[:, 0].clip(0.0, 1.0)
        rho_v = sol[:, 2].clip(0.0)
        phi_a = (1.0 - phi_i).clip(0.0, 1.0)

        mass_ice[k] = RHO_ICE * np.sum(phi_i) * dV
        mass_vap[k] = np.sum(rho_v * phi_a) * dV

    return times, mass_ice, mass_vap, dim


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def _pct_change(mass: np.ndarray) -> float:
    """Final-vs-initial percent change of a mass time series.

    Robust to mass[0] == 0 (returns NaN in that case so the formatter shows
    "n/a" rather than dividing by zero).
    """
    if mass[0] == 0.0:
        return float("nan")
    return (mass[-1] - mass[0]) / mass[0] * 100.0


def _fmt_pct(pct: float) -> str:
    """Format a percent change like '+0.12%' / '-0.45e-5%' / 'n/a' (NaN)."""
    if not np.isfinite(pct):
        return "n/a"
    if abs(pct) < 1e-3:
        return f"{pct:+.3e} %"
    return f"{pct:+.4f} %"


def _per_phase_plot(t_plot, mass, color, label_math, xlabel, unit_suffix,
                    title, save_path):
    """One self-contained plot for a single phase with its initial-value reference.

    Legend now reports the final-vs-initial percent change so the magnitude of
    drift is visible without having to read the y-axis values."""
    fig, ax = plt.subplots(figsize=(8, 5))
    pct = _pct_change(mass)
    ax.plot(t_plot, mass, color=color, lw=1.8, zorder=3,
            label=f"{label_math}   Δ = {_fmt_pct(pct)}")
    ax.axhline(mass[0], color=color, lw=1.0, ls="--", alpha=0.45, zorder=2,
               label=f"initial = {mass[0]:.3e}")
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(f"Mass  [{unit_suffix}]", fontsize=12)
    ax.set_title(title, fontsize=13)
    ax.legend(fontsize=10, loc="best")
    ax.grid(True, alpha=0.3)
    ax.tick_params(labelsize=10)
    fig.tight_layout()
    fig.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def _change_loglog_plot(times_sec, masses, labels, colors, save_path):
    """Log-log plot of |mass(t) - mass(0)| / mass(0) vs time (in seconds).

    Each entry in ``masses`` is a 1-D array (same length as ``times_sec``).
    The y-axis is the *absolute relative change* so positive- and negative-going
    drifts both appear on a log scale. Legend reports the *signed* final
    percent change so the direction is preserved."""
    fig, ax = plt.subplots(figsize=(8, 5))

    # Drop t=0 since log(0) is undefined; keep only strictly positive times.
    pos = times_sec > 0
    if not np.any(pos):
        plt.close(fig)
        return
    t = times_sec[pos]

    plotted_any = False
    for mass, label, color in zip(masses, labels, colors):
        if mass[0] == 0.0:
            continue   # nothing to normalize against
        rel = np.abs(mass[pos] - mass[0]) / abs(mass[0])
        # Hide exact zeros (which become -inf in log) by clipping a floor.
        rel = np.where(rel > 0, rel, np.nan)
        pct = _pct_change(mass)
        ax.loglog(t, rel, color=color, lw=1.6, zorder=3,
                  label=f"{label}   Δ = {_fmt_pct(pct)}")
        plotted_any = True

    if not plotted_any:
        plt.close(fig)
        return

    ax.set_xlabel("Time  [s]", fontsize=12)
    ax.set_ylabel(r"$|m(t) - m(0)| / |m(0)|$", fontsize=12)
    ax.set_title("Relative change in mass vs. time (log-log)", fontsize=13)
    ax.legend(fontsize=10, loc="best")
    ax.grid(True, which="both", alpha=0.3)
    ax.tick_params(labelsize=10)
    fig.tight_layout()
    fig.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_mass(run_dir: str, time_unit: str = None, save_path: str = None,
              per_phase_dir: str = None):
    """Generate the combined mass plot and (by default) per-phase subplots.

    Per-phase plots go into <run_dir>/mass_plots/ unless `per_phase_dir` is
    given. Each phase (total, ice, sed, vapor) gets its own PNG on its own
    y-axis so subtle drifts that are hidden by the dominant-phase scale on
    the combined plot become visible.
    """

    times, mass_ice, mass_vap, dim = compute_masses(run_dir)
    mass_total = mass_ice + mass_vap

    # Auto time unit if not given
    t_max = times[-1] if len(times) > 0 else 1.0
    unit = time_unit or auto_time_unit(t_max)
    scale, xlabel = TIME_SCALES.get(unit, (1.0, "Time  [s]"))
    t_plot = times / scale

    # Mass units label
    unit_suffix = {1: "kg m⁻²", 2: "kg m⁻¹", 3: "kg"}.get(dim, "kg")

    # ── Combined figure (unchanged behaviour) ─────────────────────────────
    fig, ax = plt.subplots(figsize=(8, 5))

    pct_total = _pct_change(mass_total)
    pct_ice   = _pct_change(mass_ice)
    pct_vap   = _pct_change(mass_vap)

    ax.plot(t_plot, mass_total, color=COLOR_TOTAL, lw=2.2, zorder=4,
            label=f"Total   Δ = {_fmt_pct(pct_total)}")
    ax.plot(t_plot, mass_ice,   color=COLOR_ICE,   lw=1.6, zorder=3,
            label=fr"Ice  ($\rho_i \int \phi_i \, dV$)   Δ = {_fmt_pct(pct_ice)}")
    ax.plot(t_plot, mass_vap,   color=COLOR_VAP,   lw=1.6, zorder=3,
            label=fr"Vapor  ($\int \rho_v \, \phi_a \, dV$)   Δ = {_fmt_pct(pct_vap)}")

    # Dashed semi-opaque reference lines at the initial (t=0) mass values.
    # These make it easy to see whether each phase is gaining or losing mass.
    ax.axhline(mass_total[0], color=COLOR_TOTAL, lw=1.0, ls="--", alpha=0.45, zorder=2)
    ax.axhline(mass_ice[0],   color=COLOR_ICE,   lw=1.0, ls="--", alpha=0.45, zorder=2)
    ax.axhline(mass_vap[0],   color=COLOR_VAP,   lw=1.0, ls="--", alpha=0.45, zorder=2)

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

    # ── Per-phase figures ─────────────────────────────────────────────────
    pp_dir = per_phase_dir or os.path.join(run_dir, "mass_plots")
    os.makedirs(pp_dir, exist_ok=True)

    _per_phase_plot(t_plot, mass_total, COLOR_TOTAL, "Total",
                    xlabel, unit_suffix, "Total mass vs. time",
                    os.path.join(pp_dir, "total.png"))
    _per_phase_plot(t_plot, mass_ice, COLOR_ICE,
                    r"Ice  ($\rho_i \int \phi_i \, dV$)",
                    xlabel, unit_suffix, "Ice mass vs. time",
                    os.path.join(pp_dir, "ice.png"))
    _per_phase_plot(t_plot, mass_vap, COLOR_VAP,
                    r"Vapor  ($\int \rho_v \, \phi_a \, dV$)",
                    xlabel, unit_suffix, "Vapor mass vs. time",
                    os.path.join(pp_dir, "vapor.png"))

    print(f"Per-phase plots saved to: {pp_dir}/")

    # ── Log-log relative-change plot ──────────────────────────────────────
    _change_loglog_plot(
        times_sec=times,
        masses=[mass_total, mass_ice, mass_vap],
        labels=["Total",
                r"Ice  ($\rho_i \int \phi_i \, dV$)",
                r"Vapor  ($\int \rho_v \, \phi_a \, dV$)"],
        colors=[COLOR_TOTAL, COLOR_ICE, COLOR_VAP],
        save_path=os.path.join(pp_dir, "change_loglog.png"),
    )
    print(f"Log-log change plot saved to: {pp_dir}/change_loglog.png")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Plot integrated phase mass vs. time (1D/2D/3D)."
    )
    p.add_argument("--dir",            default=".",
                   help="Run folder containing igasol.dat and sol_*.dat (default: .)")
    p.add_argument("--save",           default=None,
                   help="Output PNG path for the combined plot (default: <dir>/mass.png)")
    p.add_argument("--per-phase-dir",  default=None,
                   help="Directory for per-phase plots (default: <dir>/mass_plots/)")
    p.add_argument("--time-unit",      default=None,
                   choices=["s", "min", "h", "d"],
                   help="Time unit for x-axis (default: auto)")
    return p.parse_args()


def main():
    args = parse_args()
    run_dir = os.path.abspath(args.dir)
    plot_mass(run_dir,
              time_unit=args.time_unit,
              save_path=args.save,
              per_phase_dir=args.per_phase_dir)


if __name__ == "__main__":
    main()

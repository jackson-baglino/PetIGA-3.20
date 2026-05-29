#!/usr/bin/env python3
"""
plot_energy.py  —  Plot the phase-field free energy of the system vs. time.

Reads sol_*.dat snapshots and evaluates the Allen-Cahn free-energy functional
that the phase-field dynamics minimise:

    F(t) = F_bulk(t) + F_grad(t)

    F_bulk = ∫ [ Σ_k η_k φ_k²(1−φ_k)²  +  Λ φ_i² φ_s² φ_a² ] dΩ      (multi-well + triple junction)
    F_grad = ∫ [ (ε²/2) Σ_k |∇φ_k|² ] dΩ                            (diffuse-interface gradient energy)

with k ∈ {ice, sed, air}, φ_a = 1 − φ_i − φ_s.  The η_k combination energies
and Λ, ε are read from the staged .opts files (falling back to the
src/permafrost2.c defaults).

Why this is useful as a diagnostic
----------------------------------
  * F is a Lyapunov functional for the Allen-Cahn part: under pure AC
    relaxation it decreases monotonically and plateaus at equilibrium.
    A flat F over the last steps = the system has stalled (nothing evolving).
  * The bulk/gradient split is diagnostic on its own: if F_grad rises while
    F_bulk is flat, the interface is spreading (becoming more diffuse) without
    the phases changing — which is exactly the "neck looks too diffuse"
    symptom.  If F_grad falls, interfaces are sharpening/shrinking.

Unit notes (same convention as plot_mass.py)
--------------------------------------------
  η_k, Λ have units J/m²; ε has units m. With dV in m^dim the reported energy
  is per-unit-out-of-plane-measure (1D: J/m², 2D: J/m, 3D: J). Absolute scale
  depends on the model's normalisation; the *trend* is what's diagnostic.

Usage
-----
  python postprocess/plot_energy.py --dir /path/to/run [--save energy.png] [--time-unit h]
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
    sys.exit("ERROR: igakit is not installed.  Install with:  pip install igakit")

# ---------------------------------------------------------------------------
# Defaults — must match src/permafrost2.c (overridden by .opts when present)
# ---------------------------------------------------------------------------
GAMMA_IV = 0.109    # ice–vapor      surface energy [J/m²]
GAMMA_IM = 0.057    # ice–sediment   surface energy [J/m²]
GAMMA_MV = 0.080    # sediment–vapor surface energy [J/m²]
LAMBDA_DEFAULT = 1.0e4
EPS_DEFAULT = 7.12e-7

COLOR_TOTAL = "black"
COLOR_BULK  = "#b2182b"   # red
COLOR_GRAD  = "#2166ac"   # blue

TIME_SCALES = {
    "s":   (1.0,     "Time  [s]"),
    "min": (60.0,    "Time  [min]"),
    "h":   (3600.0,  "Time  [h]"),
    "d":   (86400.0, "Time  [days]"),
}


def auto_time_unit(t_max_sec: float) -> str:
    if t_max_sec <= 600:
        return "s"
    if t_max_sec <= 7200:
        return "min"
    if t_max_sec <= 3 * 86400:
        return "h"
    return "d"


# ---------------------------------------------------------------------------
# Opts-file helpers (mirror plot_mass.py so behaviour is consistent)
# ---------------------------------------------------------------------------

def find_opts_files(run_dir: str):
    result = []
    for base in ("solver.opts", "universal.opts"):
        p = os.path.join(run_dir, base)
        if os.path.isfile(p):
            result.append(p)
            break
    for f in sorted(os.listdir(run_dir)):
        if f.endswith(".opts") and f not in ("solver.opts", "universal.opts"):
            result.append(os.path.join(run_dir, f))
    return result


def parse_opts_float(opts_files, key: str, default=None):
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


def load_ssa(path: str):
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
    mask = ~np.isnan(data[:, :4]).any(axis=1)
    data = data[mask]
    if len(data) == 0:
        return None
    steps = data[:, 3].astype(int)
    _, last_idx = np.unique(steps[::-1], return_index=True)
    keep = np.sort(len(steps) - 1 - last_idx)
    return data[keep]


# ---------------------------------------------------------------------------
# Core computation
# ---------------------------------------------------------------------------

def compute_energy(run_dir: str):
    """Return (times, F_total, F_bulk, F_grad, dim)."""
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

    opts = find_opts_files(run_dir)
    g_iv = parse_opts_float(opts, "-gamma_iv", GAMMA_IV)
    g_im = parse_opts_float(opts, "-gamma_im", GAMMA_IM)
    g_mv = parse_opts_float(opts, "-gamma_mv", GAMMA_MV)
    eta_i = g_iv + g_im - g_mv
    eta_s = g_mv + g_im - g_iv
    eta_a = g_iv + g_mv - g_im
    Lambda = parse_opts_float(opts, "-Lambda", LAMBDA_DEFAULT)
    eps = parse_opts_float(opts, "-eps", EPS_DEFAULT)

    nrb = PetIGA().read(geo_file)
    dim = nrb.dim

    coords = nrb.points[..., :dim]
    shape = nrb.points.shape[:-1]              # (Nx,) / (Nx,Ny) / (Nx,Ny,Nz)
    flat = coords.reshape(-1, dim)
    extents = flat.max(axis=0) - flat.min(axis=0)
    V_domain = float(np.prod(extents))
    N_total = int(np.prod(shape))
    dV = V_domain / N_total

    # Per-axis spacing for gradients (uniform-grid approximation on the
    # control net, consistent with plot_mass.py's grid treatment).
    spacing = [extents[d] / max(shape[d] - 1, 1) for d in range(dim)]

    ssa = load_ssa(os.path.join(run_dir, "SSA_evo.dat"))
    n_sol = len(sol_files)
    if ssa is not None:
        ssa_times = ssa[:, 2]
        if n_sol == len(ssa_times) + 1:
            times = np.concatenate([[0.0], ssa_times])
        elif n_sol == len(ssa_times):
            times = ssa_times
        else:
            n = min(n_sol, len(ssa_times))
            times = ssa_times[:n]
            sol_files = sol_files[:n]
    else:
        times = np.arange(n_sol, dtype=float)

    n = len(times)
    F_bulk = np.zeros(n)
    F_grad = np.zeros(n)

    for k, sf in enumerate(sol_files[:n]):
        sol = PetIGA().read_vec(sf, nrb).reshape(-1, 4)
        phi_i = sol[:, 0].clip(0.0, 1.0)
        phi_s = sol[:, 3].clip(0.0, 1.0)
        phi_a = (1.0 - phi_i - phi_s).clip(0.0, 1.0)

        # --- bulk multi-well + triple-junction energy density ---
        well = (eta_i * phi_i**2 * (1.0 - phi_i)**2
                + eta_s * phi_s**2 * (1.0 - phi_s)**2
                + eta_a * phi_a**2 * (1.0 - phi_a)**2
                + Lambda * phi_i**2 * phi_s**2 * phi_a**2)
        F_bulk[k] = float(np.sum(well) * dV)

        # --- gradient (diffuse-interface) energy density ---
        grad_sq = np.zeros(N_total)
        for fld in (phi_i, phi_s, phi_a):
            g = np.gradient(fld.reshape(shape), *spacing)
            if dim == 1:
                g = [g]
            for gc in g:
                grad_sq += (gc.reshape(-1))**2
        F_grad[k] = float(0.5 * eps**2 * np.sum(grad_sq) * dV)

    return times, F_bulk + F_grad, F_bulk, F_grad, dim


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def _unit_suffix(dim):
    return {1: "J/m²", 2: "J/m", 3: "J"}.get(dim, "a.u.")


def plot_energy(run_dir, save_path, time_unit):
    times, F_tot, F_bulk, F_grad, dim = compute_energy(run_dir)
    if len(times) == 0:
        sys.exit("ERROR: no snapshots to plot.")

    unit = time_unit or auto_time_unit(times[-1] if times[-1] > 0 else 1.0)
    scale, xlabel = TIME_SCALES[unit]
    t_plot = times / scale
    t_freeze = parse_opts_float(find_opts_files(run_dir), "-t_sed_freeze", None)
    usuf = _unit_suffix(dim)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    # Total + components
    ax1.plot(t_plot, F_tot,  color=COLOR_TOTAL, lw=2.0, label="F total")
    ax1.plot(t_plot, F_bulk, color=COLOR_BULK,  lw=1.5, ls="-",  label="F bulk (wells + triple jn)")
    ax1.plot(t_plot, F_grad, color=COLOR_GRAD,  lw=1.5, ls="-",  label="F gradient (interface)")
    if F_tot[0] != 0:
        d_tot = (F_tot[-1] - F_tot[0]) / abs(F_tot[0]) * 100.0
        ax1.set_title(f"Phase-field free energy   (ΔF = {d_tot:+.3f} %)", fontsize=12)
    else:
        ax1.set_title("Phase-field free energy", fontsize=12)
    ax1.set_xlabel(xlabel); ax1.set_ylabel(f"Free energy  [{usuf}]")
    ax1.grid(alpha=0.3)

    # Components normalised to their initial value (trend, regardless of scale)
    if F_bulk[0] != 0:
        ax2.plot(t_plot, F_bulk / F_bulk[0], color=COLOR_BULK, lw=1.8, label="F bulk / F bulk(0)")
    if F_grad[0] != 0:
        ax2.plot(t_plot, F_grad / F_grad[0], color=COLOR_GRAD, lw=1.8, label="F grad / F grad(0)")
    ax2.axhline(1.0, color="gray", lw=0.8, ls=":")
    ax2.set_title("Components, normalised to t=0", fontsize=12)
    ax2.set_xlabel(xlabel); ax2.set_ylabel("relative energy")
    ax2.grid(alpha=0.3)

    for ax in (ax1, ax2):
        if t_freeze is not None:
            xf = t_freeze / scale
            if t_plot[0] <= xf <= t_plot[-1] * 1.05:
                ax.axvline(xf, color="green", ls="--", lw=1.0, alpha=0.6)
        ax.legend(fontsize=9)

    fig.tight_layout()
    fig.savefig(save_path, dpi=150)
    print(f"  Wrote {save_path}")
    print(f"  F_total: {F_tot[0]:.4e} -> {F_tot[-1]:.4e} [{usuf}]   "
          f"(bulk {F_bulk[-1]/F_bulk[0]*100 if F_bulk[0] else float('nan'):.1f}% , "
          f"grad {F_grad[-1]/F_grad[0]*100 if F_grad[0] else float('nan'):.1f}% of t=0)")


def main():
    ap = argparse.ArgumentParser(description="Plot phase-field free energy vs time.")
    ap.add_argument("--dir", default=".", help="run directory (default: cwd)")
    ap.add_argument("--save", default=None, help="output PNG (default: <dir>/energy.png)")
    ap.add_argument("--time-unit", choices=list(TIME_SCALES.keys()), default=None)
    args = ap.parse_args()

    run_dir = os.path.abspath(args.dir)
    save_path = args.save or os.path.join(run_dir, "energy.png")
    plot_energy(run_dir, save_path, args.time_unit)


if __name__ == "__main__":
    main()

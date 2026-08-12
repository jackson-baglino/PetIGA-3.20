#!/usr/bin/env python3
"""
plot_energy.py  —  Plot the phase-field free energy of the system vs. time.

Reads sol_*.dat snapshots and evaluates the Allen-Cahn free-energy functional
that the phase-field dynamics minimise:

    F(t) = F_bulk(t) + F_grad(t)

The number of fields per node is detected from the snapshots, and the
functional follows the solver that wrote them.

Two-phase solvers (dof=3, Field = {ice, tem, rhov} — the current
lunar_regolith_DSM and enceladus_DSM mains):

    F_bulk = √2 γ_iv ∫ (3/ε) φ²(1−φ)²  dΩ
    F_grad = √2 γ_iv ∫ (3ε/2) |∇φ|²    dΩ

This is exactly the Lyapunov functional of the ice residual in
src/assembly.c (R_ice = φ_t + 3M[ε ∇N·∇φ + f′(φ)N/ε], f = φ²(1−φ)²), scaled
so that an equilibrated flat interface carries γ_iv per unit area.

Legacy 3-phase solver (dof=4, with a fourth sediment field):

    F_bulk = ∫ [ Σ_k η_k φ_k²(1−φ_k)²  +  Λ φ_i² φ_s² φ_a² ] dΩ      (multi-well + triple junction)
    F_grad = ∫ [ (ε²/2) Σ_k |∇φ_k|² ] dΩ                            (diffuse-interface gradient energy)

with k ∈ {ice, sed, air}, φ_a = 1 − φ_i − φ_s.  The η_k combination energies
and Λ, ε are read from the staged .opts files (falling back to the
src/<project>_main.c defaults).

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
from pplib import auto_time_unit, load_ssa

try:
    from igakit.io import PetIGA
except ImportError:
    sys.exit("ERROR: igakit is not installed.  Install with:  pip install igakit")

# ---------------------------------------------------------------------------
# Defaults — must match src/<project>_main.c (overridden by .opts when present)
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
    N_total = int(np.prod(shape))

    # --- metric of the control net -----------------------------------------
    # The net is NOT a uniform Cartesian grid: open knot vectors bunch control
    # points near the boundaries, and mapped geometries (e.g. the wedge_band
    # family, whose channel widens from 5e-5 m to 2e-4 m along x) make the
    # physical spacing a function of position.  Treating it as uniform gets
    # both the cell volume and |grad phi| wrong -- on the wedge that is a 1.6x
    # error in the domain measure alone.
    #
    # So build the Jacobian of the index -> physical map, J[..., d, k] =
    # dx_d/dxi_k, using unit spacing in index space.  |det J| is then the cell
    # volume and J^-1 converts index-space gradients into physical ones.
    def _idx_grad(arr):
        """d(arr)/d(index) along each axis -> list of dim arrays."""
        g = np.gradient(arr)
        return [g] if dim == 1 else list(g)

    J = np.empty(shape + (dim, dim))
    for d in range(dim):
        for k, gk in enumerate(_idx_grad(coords[..., d])):
            J[..., d, k] = gk

    cellV = np.abs(np.linalg.det(J)).reshape(-1)     # per-node volume [m^dim]
    Jinv = np.linalg.inv(J)                          # Jinv[..., k, d] = dxi_k/dx_d
    V_domain = float(cellV.sum())

    def phys_grad_sq(fld):
        """|grad fld|^2 in physical space, flattened."""
        dfd = _idx_grad(fld.reshape(shape))
        out = np.zeros(shape)
        for d in range(dim):
            gd = np.zeros(shape)
            for k in range(dim):
                gd += dfd[k] * Jinv[..., k, d]
            out += gd**2
        return out.reshape(-1)

    # Time axis: match each sol_NNNNN.dat to its SSA_evo.dat row BY STEP NUMBER,
    # exactly as plot_mass.py does.  sol_*.dat is written every -outp N steps
    # while SSA_evo.dat carries a row per step, so the two are not index-aligned
    # -- pairing them positionally (or truncating to the shorter one) silently
    # compresses the time axis, e.g. reporting 6.3 days for a 90-day run.
    ssa = load_ssa(os.path.join(run_dir, "SSA_evo.dat"))
    if ssa is not None:
        step_to_time = dict(zip(ssa[:, 3].astype(int), ssa[:, 2]))
        times = []
        matched_files = []
        for sf in sol_files:
            digits = "".join(ch for ch in os.path.basename(sf) if ch.isdigit())
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
        times = np.arange(len(sol_files), dtype=float)

    if not sol_files:
        sys.exit("ERROR: no sol_*.dat snapshot could be matched to an SSA_evo.dat step.")

    # Fields per node.  The current solvers carry Field = {ice, tem, rhov}
    # (dof=3, two phases: ice and air); the legacy 3-phase solver carried a
    # fourth sediment field (dof=4).  Detect it from the first snapshot rather
    # than hard-coding, so the same script serves both.
    n_first = PetIGA().read_vec(sol_files[0], nrb).size
    dof, rem = divmod(n_first, N_total)
    if rem != 0 or dof not in (3, 4):
        sys.exit(
            f"ERROR: {os.path.basename(sol_files[0])} holds {n_first} values on a "
            f"{N_total}-node mesh ({'x'.join(str(s) for s in shape)}); expected an "
            f"integer 3 or 4 fields per node, got {n_first / N_total:g}. "
            "Geometry and solution files are from different runs?"
        )
    three_phase = (dof == 4)

    # Two-phase calibration.  src/assembly.c integrates
    #     R_ice = phi_t + 3 M [ eps grad_N.grad_phi + f1(phi) N / eps ]
    # i.e. phi_t = -M dF/dphi with the Lyapunov functional
    #     F = C * int 3 [ g(phi)/eps + (eps/2)|grad phi|^2 ] dOmega,
    # where g is the antiderivative of the code's f1.  NOTE that f1 there is
    # phi(1-phi)(1-2phi), which is HALF of d/dphi[phi^2 (1-phi)^2] -- the
    # comment at assembly.c:10 drops that factor.  The well the solver
    # actually descends is g = (1/2) phi^2 (1-phi)^2, which gives the
    # equilibrium profile phi = 1/(1+exp(-x/eps)) (decay length exactly eps,
    # as documented in CLAUDE.md) and equipartition F_bulk == F_grad.
    #
    # That bracket integrates to 1/2 across an equilibrated flat profile, so
    # C = 2*gamma_iv makes F equal exactly gamma_iv per unit interface area
    # — giving the same J/m (2D) / J (3D)
    # units as _unit_suffix() reports.
    C2 = 2.0 * g_iv

    n = len(times)
    F_bulk = np.zeros(n)
    F_grad = np.zeros(n)

    for k, sf in enumerate(sol_files[:n]):
        sol = PetIGA().read_vec(sf, nrb).reshape(-1, dof)
        phi_i = sol[:, 0].clip(0.0, 1.0)

        if three_phase:
            phi_s = sol[:, 3].clip(0.0, 1.0)
            phi_a = (1.0 - phi_i - phi_s).clip(0.0, 1.0)

            # --- bulk multi-well + triple-junction energy density ---
            well = (eta_i * phi_i**2 * (1.0 - phi_i)**2
                    + eta_s * phi_s**2 * (1.0 - phi_s)**2
                    + eta_a * phi_a**2 * (1.0 - phi_a)**2
                    + Lambda * phi_i**2 * phi_s**2 * phi_a**2)
            grad_fields = (phi_i, phi_s, phi_a)
            grad_pref = 0.5 * eps**2
        else:
            # --- two-phase (ice/air) double well, see C2 above ---
            well = C2 * (1.5 / eps) * phi_i**2 * (1.0 - phi_i)**2
            grad_fields = (phi_i,)
            grad_pref = C2 * 1.5 * eps

        F_bulk[k] = float(np.sum(well * cellV))

        # --- gradient (diffuse-interface) energy density ---
        grad_sq = np.zeros(N_total)
        for fld in grad_fields:
            grad_sq += phys_grad_sq(fld)
        F_grad[k] = float(grad_pref * np.sum(grad_sq * cellV))

    return times, F_bulk + F_grad, F_bulk, F_grad, dim, three_phase


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def _unit_suffix(dim):
    return {1: "J/m²", 2: "J/m", 3: "J"}.get(dim, "a.u.")


def plot_energy(run_dir, save_path, time_unit):
    times, F_tot, F_bulk, F_grad, dim, three_phase = compute_energy(run_dir)
    if len(times) == 0:
        sys.exit("ERROR: no snapshots to plot.")
    bulk_label = ("F bulk (wells + triple jn)" if three_phase
                  else "F bulk (double well)")

    unit = time_unit or auto_time_unit(times[-1] if times[-1] > 0 else 1.0)
    scale, xlabel = TIME_SCALES[unit]
    t_plot = times / scale
    t_freeze = parse_opts_float(find_opts_files(run_dir), "-t_sed_freeze", None)
    usuf = _unit_suffix(dim)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    # Total + components
    ax1.plot(t_plot, F_tot,  color=COLOR_TOTAL, lw=2.0, label="F total")
    ax1.plot(t_plot, F_bulk, color=COLOR_BULK,  lw=1.5, ls="-",  label=bulk_label)
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

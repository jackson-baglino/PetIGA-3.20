#!/usr/bin/env python3
"""
plot2D_snapshot.py  —  2D field snapshots without ParaView/VTK.

Reads a PetIGA 2D solution file and produces a multi-panel matplotlib figure
showing all four fields evaluated on a regular grid:
  - Ice phase field  φ_i
  - Sediment phase   φ_s  (from soil.dat)
  - Temperature      T
  - Vapor density    ρ_v
  - Air phase        φ_a = 1 − φ_i − φ_s  (derived)
  - Vapor super-saturation  (ρ_v − ρ_vs(T)) / ρ_vs(T)

Also optionally plots 1D cross-sectional profiles along a horizontal or
vertical cut through the domain.

Usage
-----
  # Plot step 10 from current directory (uses sol_00010.dat)
  python plot2D_snapshot.py --step 10

  # Specify directory and grid resolution
  python plot2D_snapshot.py --dir /path/to/output --step 50 --nx 200 --ny 200

  # Also plot x- and y-cuts through the domain midpoint
  python plot2D_snapshot.py --step 10 --cuts

  # Save figure
  python plot2D_snapshot.py --step 10 --save snap_010.png
"""

import argparse
import os
import sys
import glob

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1 import make_axes_locatable

try:
    from igakit.io import PetIGA
    from igakit.nurbs import NURBS
except ImportError:
    sys.exit("ERROR: igakit is required.  Install with:  pip install igakit")


# ---------------------------------------------------------------------------
# Saturation vapor density (matches material_properties.c: RhoVS_I)
# ---------------------------------------------------------------------------

def rho_vs(T_C: np.ndarray) -> np.ndarray:
    """Saturation vapor density over ice [kg/m³] at temperature T_C [°C]."""
    T_K = T_C + 273.15
    return 3.25e-3 * np.exp(-6150.0 / T_K)


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def _find_sol_file(run_dir: str, step: int) -> str:
    """Return path to sol_NNNNN.dat for the given step number."""
    pattern = os.path.join(run_dir, f"sol_{step:05d}.dat")
    if os.path.isfile(pattern):
        return pattern
    # Fallback: any file ending with the right number
    candidates = sorted(glob.glob(os.path.join(run_dir, "sol_*.dat")))
    if not candidates:
        sys.exit(f"No sol_*.dat files found in '{run_dir}'")
    if step < 0 or step >= len(candidates):
        return candidates[-1]
    return candidates[step]


def _load_2d(run_dir: str, step: int, nx_eval: int, ny_eval: int,
             iga_file: str = "igasol.dat"):
    """
    Read the 2D solution at a given step and evaluate on an (nx_eval × ny_eval)
    regular grid.  Also loads the sediment field if available.

    Returns
    -------
    X, Y       : (ny_eval, nx_eval) meshgrid arrays in metres
    ice        : (ny_eval, nx_eval)
    tem        : (ny_eval, nx_eval)
    rhov       : (ny_eval, nx_eval)
    sed        : (ny_eval, nx_eval)  or zeros if soil.dat unavailable
    Lx, Ly     : domain extents [m]
    """
    iga_path = os.path.join(run_dir, iga_file)
    if not os.path.isfile(iga_path):
        sys.exit(f"IGA geometry file not found: {iga_path}")

    nrb = PetIGA().read(iga_path)
    sol_path = _find_sol_file(run_dir, step)
    print(f"  Reading solution: {sol_path}")
    sol_ctrl = PetIGA().read_vec(sol_path, nrb)  # (nx, ny, 3)

    # Evaluate on a regular grid ------------------------------------------
    u_vals = np.linspace(0, 1, nx_eval)
    v_vals = np.linspace(0, 1, ny_eval)
    UV     = np.meshgrid(u_vals, v_vals, indexing="ij")  # both (nx_eval, ny_eval)
    params = np.stack([UV[0].ravel(), UV[1].ravel()], axis=-1)  # (n, 2)

    # Evaluate physical coordinates and solution fields
    xyz  = nrb.evaluate(*[u_vals, v_vals])   # (nx_eval, ny_eval, ndim+1) or similar
    # igakit NURBS.evaluate returns (nx, ny, physical_dim+1) for 2D:
    # last axis is homogeneous weight; first ndim are x, y
    Xe = xyz[..., 0]   # (nx_eval, ny_eval)
    Ye = xyz[..., 1]

    # Evaluate solution fields
    ice_e  = _eval_field(nrb, sol_ctrl, u_vals, v_vals, comp=0)
    tem_e  = _eval_field(nrb, sol_ctrl, u_vals, v_vals, comp=1)
    rhov_e = _eval_field(nrb, sol_ctrl, u_vals, v_vals, comp=2)

    Lx = float(Xe.max() - Xe.min())
    Ly = float(Ye.max() - Ye.min())

    # Sediment ----------------------------------------------------------------
    sed_e = np.zeros_like(ice_e)
    soil_dat  = os.path.join(run_dir, "soil.dat")
    igasoil   = os.path.join(run_dir, "igasoil.dat")
    if os.path.isfile(soil_dat) and os.path.isfile(igasoil):
        nrb_s   = PetIGA().read(igasoil)
        sed_ctrl = PetIGA().read_vec(soil_dat, nrb_s)
        sed_e   = _eval_field(nrb_s, sed_ctrl, u_vals, v_vals, comp=0)
        print("  Sediment field loaded.")
    else:
        print("  soil.dat / igasoil.dat not found — sediment set to zero.")

    # Transpose to (ny, nx) for imshow convention (row = y, col = x)
    return (Xe.T, Ye.T,
            ice_e.T, tem_e.T, rhov_e.T, sed_e.T,
            Lx, Ly)


def _eval_field(nrb: NURBS, ctrl: np.ndarray,
                u_vals: np.ndarray, v_vals: np.ndarray, comp: int) -> np.ndarray:
    """
    Evaluate a single field component on a (len(u_vals), len(v_vals)) grid.

    ctrl is shape (nx_ctrl, ny_ctrl, ncomp).  We use scipy's B-spline
    interpolation over the control point grid as an approximation.
    For a pure isogeometric field this is exact at control points and
    a good approximation elsewhere.
    """
    from scipy.interpolate import RegularGridInterpolator

    nx_c, ny_c = ctrl.shape[:2]
    if ctrl.ndim == 2:
        # single-component field
        vals = ctrl
    else:
        vals = ctrl[..., comp]  # (nx_c, ny_c)

    # Parametric coords of control points (uniform knot vector → uniform in [0,1])
    u_c = np.linspace(0, 1, nx_c)
    v_c = np.linspace(0, 1, ny_c)

    interp = RegularGridInterpolator((u_c, v_c), vals,
                                     method="linear", bounds_error=False,
                                     fill_value=None)
    UU, VV = np.meshgrid(u_vals, v_vals, indexing="ij")
    pts    = np.stack([UU.ravel(), VV.ravel()], axis=-1)
    result = interp(pts).reshape(len(u_vals), len(v_vals))
    return result


# ---------------------------------------------------------------------------
# Main plot
# ---------------------------------------------------------------------------

def plot_snapshot(run_dir: str, step: int, nx_eval: int = 200, ny_eval: int = 200,
                  save_path: str = None, do_cuts: bool = False,
                  iga_file: str = "igasol.dat"):
    """
    6-panel figure: ice, sediment, air, temperature, vapor density,
    vapor supersaturation.
    """
    X, Y, ice, tem, rhov, sed, Lx, Ly = _load_2d(
        run_dir, step, nx_eval, ny_eval, iga_file
    )

    # Derived fields
    air   = np.clip(1.0 - ice - sed, 0.0, 1.0)
    rvs   = rho_vs(tem)
    supersat = np.where(rvs > 0, (rhov - rvs) / rvs, 0.0)

    # Convert to mm for display
    Xmm = X * 1e3
    Ymm = Y * 1e3

    fields = [
        (ice,       r"$\phi_i$  (ice)",             "Blues",        (0, 1)),
        (sed,       r"$\phi_s$  (sediment)",         "YlOrBr",       (0, 1)),
        (air,       r"$\phi_a$  (air)",              "Greens",       (0, 1)),
        (tem,       r"$T$  [°C]",                    "RdBu_r",       None),
        (rhov,      r"$\rho_v$  [kg m$^{-3}$]",     "viridis",      None),
        (supersat,  r"$(ρ_v - ρ_{vs})/ρ_{vs}$",     "coolwarm",     None),
    ]

    fig, axes = plt.subplots(2, 3, figsize=(14, 9))
    axes = axes.ravel()

    for ax, (field, label, cmap, clim) in zip(axes, fields):
        if clim:
            im = ax.pcolormesh(Xmm, Ymm, field, cmap=cmap,
                               vmin=clim[0], vmax=clim[1], shading="auto")
        else:
            im = ax.pcolormesh(Xmm, Ymm, field, cmap=cmap, shading="auto")

        divider = make_axes_locatable(ax)
        cax     = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(im, cax=cax)

        ax.set_title(label, fontsize=12)
        ax.set_xlabel("x  [mm]", fontsize=10)
        ax.set_ylabel("y  [mm]", fontsize=10)
        ax.set_aspect("equal")
        ax.tick_params(labelsize=9)

    run_label = os.path.basename(run_dir.rstrip("/")) or run_dir
    fig.suptitle(f"2D snapshot — step {step} — {run_label}", fontsize=13, y=1.01)
    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Figure saved to: {save_path}")
    else:
        plt.show()

    if do_cuts:
        _plot_cuts(Xmm, Ymm, ice, tem, rhov, sed, step, run_dir, save_path)


def _plot_cuts(Xmm, Ymm, ice, tem, rhov, sed, step, run_dir, save_path):
    """Horizontal (y=Ly/2) and vertical (x=Lx/2) cross-section profiles."""
    mid_row = Ymm.shape[0] // 2
    mid_col = Xmm.shape[1] // 2

    x_cut = Xmm[mid_row, :]
    y_cut = Ymm[:, mid_col]

    air = np.clip(1.0 - ice - sed, 0.0, 1.0)

    fig, axes = plt.subplots(2, 3, figsize=(14, 7))

    def _cut_row(ax, field, title):
        ax.plot(x_cut, field[mid_row, :], lw=2)
        ax.set_xlabel("x  [mm]", fontsize=10)
        ax.set_title(f"{title}  (y = {Ymm[mid_row, 0]:.2f} mm)", fontsize=11)
        ax.grid(True, alpha=0.3)

    def _cut_col(ax, field, title):
        ax.plot(field[:, mid_col], y_cut, lw=2)
        ax.set_ylabel("y  [mm]", fontsize=10)
        ax.set_title(f"{title}  (x = {Xmm[0, mid_col]:.2f} mm)", fontsize=11)
        ax.grid(True, alpha=0.3)

    _cut_row(axes[0, 0], ice,  r"$\phi_i$")
    _cut_row(axes[0, 1], tem,  r"$T$ [°C]")
    _cut_row(axes[0, 2], rhov, r"$\rho_v$")

    _cut_col(axes[1, 0], ice,  r"$\phi_i$")
    _cut_col(axes[1, 1], tem,  r"$T$ [°C]")
    _cut_col(axes[1, 2], rhov, r"$\rho_v$")

    fig.suptitle(f"Cross-section profiles — step {step}", fontsize=13)
    plt.tight_layout()

    if save_path:
        cut_path = save_path.replace(".png", "_cuts.png").replace(".pdf", "_cuts.pdf")
        fig.savefig(cut_path, dpi=150, bbox_inches="tight")
        print(f"Cuts figure saved to: {cut_path}")
    else:
        plt.show()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Plot 2D field snapshot from PetIGA output."
    )
    p.add_argument("--dir",   default=".",         help="Output directory")
    p.add_argument("--step",  type=int, default=0, help="Step number to plot (0 = first)")
    p.add_argument("--nx",    type=int, default=200, help="Evaluation grid resolution in x")
    p.add_argument("--ny",    type=int, default=200, help="Evaluation grid resolution in y")
    p.add_argument("--iga",   default="igasol.dat", help="IGA geometry file")
    p.add_argument("--cuts",  action="store_true",  help="Also plot cross-section profiles")
    p.add_argument("--save",  default=None,         help="Save figure to this path")
    return p.parse_args()


def main():
    args = parse_args()
    plot_snapshot(
        run_dir=args.dir, step=args.step,
        nx_eval=args.nx, ny_eval=args.ny,
        save_path=args.save, do_cuts=args.cuts,
        iga_file=args.iga,
    )


if __name__ == "__main__":
    main()

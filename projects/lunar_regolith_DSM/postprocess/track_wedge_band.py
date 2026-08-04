#!/usr/bin/env python3
"""track_wedge_band.py — does the ice band migrate along the wedge?

The wedge run (geometry/2D_wedge_band.opts, zero temperature gradient) asks one
question: does pore taper alone redistribute ice? The answer is one number per
snapshot — the ice-weighted centroid's distance from the wedge apex, r_c(t):

    FALLING r_c  -> the band migrates toward the NARROW end, as predicted from
                    the meniscus curvatures (inner face concave, kappa=-1/r1;
                    outer convex, kappa=+1/r2; so the inner face is the vapour
                    sink).
    FLAT r_c     -> no migration. Before calling it a null, check the total ice
                    is also flat: a band that is simply evaporating uniformly
                    holds its centroid while losing mass.
    RISING r_c   -> migration toward the WIDE end, contradicting the curvature
                    argument. Worth chasing rather than reporting.

AREA WEIGHTING. The sampling grid is uniform in PARAMETRIC (u,v), but the wedge
mesh maps it to physical cells whose area varies 4x between throat and mouth
(Jacobian min/max = 0.25), so parametric sums are not proportional to mass.
Every sum below is therefore weighted by the local Jacobian, finite-differenced
from the sampled physical coordinates.

This matters most for the ICE MASS row: an unweighted sum is not proportional
to mass on this mesh, so the conservation gate would be reading the wrong
quantity. For r_c it is a smaller correction — measured against a synthetic
analytic band, an unweighted centroid lands 1.3% low (biased toward the narrow
end). Because the geometry is fixed, that is a constant offset rather than a
spurious trend, so it could not by itself manufacture a migration signal — but
it would misstate how far the band actually sits from the apex.

Validated against the analytic annular sector on the production mesh: band area
and area-weighted mean radius both recovered to ~1e-3 and ~2e-5 relative.

Usage:
    python3 postprocess/track_wedge_band.py --dir <run folder> \\
        --save-csv wedge_band.csv --save-fig wedge_band.png --no-show

The apex is read from the staged geometry .opts in the run folder (run_lunar.sh
copies it there); override with --apex-x/--apex-y.
"""

import argparse
import glob
import os
import re
import sys

import numpy as np
from igakit.io import PetIGA


def _dense_uv(nrb, n_per_elem):
    """Sample n_per_elem points per element in each direction (matches
    postprocess/track_grain_mass.py)."""
    out = []
    for d in range(2):
        b = nrb.breaks(d)
        pts = [np.linspace(b[i], b[i + 1], n_per_elem, endpoint=False)
               for i in range(len(b) - 1)]
        out.append(np.unique(np.concatenate(pts + [[b[-1]]])))
    return out[0], out[1]


def _find_apex(run_dir):
    """Read -wedge_apex_x/-wedge_apex_y from whichever staged .opts has them."""
    for path in sorted(glob.glob(os.path.join(run_dir, "*.opts"))):
        txt = open(path).read()
        mx = re.search(r"^-wedge_apex_x\s+(\S+)", txt, re.M)
        my = re.search(r"^-wedge_apex_y\s+(\S+)", txt, re.M)
        if mx and my:
            return float(mx.group(1)), float(my.group(1)), os.path.basename(path)
    return None, None, None


def _cell_areas(X, Y, u, v):
    """Jacobian |d(x,y)/d(u,v)| * du * dv on the sampled grid — the physical
    area each sample point stands for."""
    du = np.gradient(u)
    dv = np.gradient(v)
    xu = np.gradient(X, u, axis=0); yu = np.gradient(Y, u, axis=0)
    xv = np.gradient(X, v, axis=1); yv = np.gradient(Y, v, axis=1)
    J = np.abs(xu * yv - xv * yu)
    return J * du[:, None] * dv[None, :]


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dir", required=True,
                   help="Run output directory (contains igasol.dat, sol_*.dat)")
    p.add_argument("--apex-x", type=float, default=None, help="Override wedge apex x [m]")
    p.add_argument("--apex-y", type=float, default=None, help="Override wedge apex y [m]")
    p.add_argument("--n-per-elem", type=int, default=3,
                   help="Dense-sample points per element (default 3)")
    p.add_argument("--save-csv", default=None, help="Write the table to this CSV path")
    p.add_argument("--save-fig", default=None, help="Save the plot to this path")
    p.add_argument("--no-show", action="store_true", help="Don't open a plot window")
    args = p.parse_args()

    igasol = os.path.join(args.dir, "igasol.dat")
    if not os.path.isfile(igasol):
        sys.exit(f"Not found: {igasol}")
    nrb = PetIGA().read(igasol)
    u, v = _dense_uv(nrb, args.n_per_elem)

    apex_x, apex_y = args.apex_x, args.apex_y
    if apex_x is None or apex_y is None:
        ax, ay, src = _find_apex(args.dir)
        if ax is None:
            sys.exit("No -wedge_apex_x/-wedge_apex_y in any staged .opts; "
                     "pass --apex-x/--apex-y explicitly.")
        apex_x = ax if apex_x is None else apex_x
        apex_y = ay if apex_y is None else apex_y
        print(f"apex ({apex_x:.4e}, {apex_y:.4e}) from {src}")

    sol_files = sorted(glob.glob(os.path.join(args.dir, "sol_*.dat")))
    if not sol_files:
        sys.exit(f"No sol_*.dat files found in {args.dir}")

    # SSA_evo.dat row i matches sol_i.dat when -outp 1 (see track_grain_mass.py).
    ssa_path = os.path.join(args.dir, "SSA_evo.dat")
    times = None
    if os.path.isfile(ssa_path):
        ssa = np.loadtxt(ssa_path)
        if ssa.ndim == 1:
            ssa = ssa.reshape(1, -1)
        if len(ssa) == len(sol_files):
            times = ssa[:, 2]
        else:
            print(f"Warning: SSA_evo.dat has {len(ssa)} rows but {len(sol_files)} "
                  f"sol_*.dat files -- falling back to step index for the time axis.",
                  file=sys.stderr)

    # Geometry is fixed, so coordinates, areas and radii are computed once.
    C0, _ = nrb(u, v, fields=PetIGA().read_vec(sol_files[0], nrb))
    X, Y = C0[..., 0], C0[..., 1]
    dA = _cell_areas(X, Y, u, v)
    rho = np.hypot(X - apex_x, Y - apex_y)

    rc, xc, mass = [], [], []
    for f in sol_files:
        _, F = nrb(u, v, fields=PetIGA().read_vec(f, nrb))
        w = F[..., 0] * dA                      # ice fraction x physical area
        m = float(w.sum())
        mass.append(m)
        rc.append(float((rho * w).sum() / m) if m > 0 else np.nan)
        xc.append(float((X * w).sum() / m) if m > 0 else np.nan)

    rc, xc, mass = np.array(rc), np.array(xc), np.array(mass)
    steps = np.arange(len(sol_files))
    if times is None:
        times = steps.astype(float)

    print(f"{'step':>6} {'time [s]':>12} {'r_c [m]':>12} {'x_c [m]':>12} {'ice area [m^2]':>15}")
    for s, t, r, x, m in zip(steps, times, rc, xc, mass):
        print(f"{s:6d} {t:12.4e} {r:12.6e} {x:12.6e} {m:15.6e}")

    if len(rc) > 1:
        dr = rc[-1] - rc[0]
        dm = (mass[-1] - mass[0]) / mass[0] * 100.0
        print(f"\nr_c change:   {dr:+.4e} m  ({dr / rc[0] * 100:+.3f}%)")
        print(f"x_c change:   {xc[-1] - xc[0]:+.4e} m")
        print(f"ice change:   {dm:+.4f}%   (should be ~0; everything is Neumann)")
        if dr < 0:
            print("\n-> r_c FELL: the band migrated toward the NARROW end, as predicted.")
        elif dr > 0:
            print("\n-> r_c ROSE: migration toward the WIDE end. This contradicts the "
                  "curvature argument -- check the contact angle in the first snapshot "
                  "and the sign of the geometry before reporting it.")
        if abs(dr) < 1e-3 * rc[0] and abs(dm) > 1.0:
            print("   NOTE: r_c is flat but ice mass moved %.2f%% -- the band is losing "
                  "mass rather than translating." % dm)

    if args.save_csv:
        with open(args.save_csv, "w") as fh:
            fh.write("step,time,r_c,x_c,ice_area\n")
            for s, t, r, x, m in zip(steps, times, rc, xc, mass):
                fh.write(f"{s},{t},{r},{x},{m}\n")
        print(f"\nWrote {args.save_csv}")

    if args.save_fig or not args.no_show:
        import matplotlib
        if args.no_show:
            matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, (a1, a2) = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
        td = times / 86400.0
        a1.plot(td, rc * 1e6, "o-", ms=3, color="#1f6fd0")
        a1.set_ylabel("centroid distance from apex  r$_c$ [µm]")
        a1.set_title("Wedge band migration — falling r$_c$ = toward the narrow end")
        a1.grid(alpha=0.3)
        a2.plot(td, mass / mass[0] * 100.0, "o-", ms=3, color="#c05a1f")
        a2.set_ylabel("ice area [% of initial]")
        a2.set_xlabel("time [days]")
        a2.grid(alpha=0.3)
        fig.tight_layout()
        if args.save_fig:
            fig.savefig(args.save_fig, dpi=140)
            print(f"Wrote {args.save_fig}")
        if not args.no_show:
            plt.show()


if __name__ == "__main__":
    main()

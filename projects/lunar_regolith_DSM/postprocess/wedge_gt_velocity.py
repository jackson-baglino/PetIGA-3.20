#!/usr/bin/env python3
"""Gibbs-Thomson check on the wedge-band menisci at the pore centreline.

Does the solver actually obey the interface condition its sharp-interface
asymptotics claim?

    sigma = beta_sub0 * v_n + d0_sub0 * chi
        =>  v_n = ( sigma * rho_vs(T)/rho_vs(T0) - d0_sub0 * chi ) / beta_sub0

This script answers that for the two ice-air menisci a wedge-band run puts on
the horizontal centreline y = Ly/2, by measuring v_n directly (rate of change
of the phi=0.5 crossing) and comparing it against the prediction above.

WHY THE CENTRELINE.  The band IC is an annulus about the VIRTUAL APEX where
the two channel walls would meet, and for the standard wedge geometry that
apex sits at y = Ly/2.  So on the centreline the radial direction IS the x
direction, the interface normal is exactly +/- x_hat, and v_n collapses to a
1-D displacement rate.  No contour differentiation, no normal reconstruction,
no curvature stencils -- which is exactly what makes this comparison clean.
Every other point on the arc would need all three.

CONVENTIONS (from src/lunar_main.c:822-833 and src/assembly.c:166-180)

  * phi = 1 in ice, 0 in vapour.  Outward normal n = -grad(phi)/|grad(phi)|,
    so v_n > 0 means GROWTH (deposition) on both faces.
  * chi is the SUM of principal curvatures, positive where the ice is convex
    into the vapour.  In 2-D planar that is +1/R.  On the centreline the left
    (inner) meniscus is concave, chi = -1/r_L; the right (outer) one is
    convex, chi = +1/r_R; r = x_interface - apex_x.
  * d0_sub0 and beta_sub0 are the UNSCALED option values.  The run header
    also prints "beta_sub (SCALED = beta0*rho_vs/rho_ice = beta_HK)" -- that
    is the Hertz-Knudsen coefficient, NOT the one that pairs with a
    sigma-normalised driving force.  Using it is wrong by rho_ice/rho_vs,
    about 1e6.

SIGMA AT THE INTERFACE.  sigma is read by fitting a line to sigma(x) on the
VAPOUR side, over [lo, hi] * eps out from the crossing, and evaluating that
fit AT the crossing.  This matters more than it looks: sigma and d0*chi agree
to ~10%, so v_n is their small residual and a 1% error in sigma is a ~10%
error in v_n.  Sampling sigma at a fixed standoff without extrapolating (the
supersat_probes.py convention) moves the prediction by ~20%.  The direct
phi=0.5 sample is carried in the CSV as a cross-check; on well-resolved data
the two agree to better than 0.1%.

THE RIPPLE IN THE MEASURED CURVE IS NOT PHYSICS.  A meniscus that travels a
few micrometres over the whole run crosses only a handful of cells, and the
sub-cell interpolated crossing carries a small periodic bias as it does -- the
same discretised-extremum bias analyze_interface.py warns about, differentiated.
It shows up as a few-percent oscillation in v_n at the period of one cell
crossing.  --vn-window widens the sliding fit against it; the summary uses
medians, which are insensitive to it.  Do not read it as a growth cycle.

Usage:
    python3 postprocess/wedge_gt_velocity.py --dir <run> --save <run>/plots/x.png
"""

import argparse
import os
import sys

import numpy as np

import matplotlib
matplotlib.use("Agg")          # headless: must precede the pyplot import
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pplib import (                                    # noqa: E402
    auto_time_unit,
    in_time_unit,
    opt_float,
    read_opts,
    read_vts,
    rho_vs,
    step_of,
    step_times,
    supersaturation,
)

CSV_NAME = "wedge_gt_velocity.csv"


# ---------------------------------------------------------------------------
# Options, with the precedence the run scripts actually use
# ---------------------------------------------------------------------------
def _parse_opts_file(path):
    out = {}
    with open(path, errors="replace") as fh:
        for line in fh:
            line = line.split("#", 1)[0].strip()
            if not line.startswith("-"):
                continue
            parts = line.split(None, 1)
            out[parts[0]] = parts[1].strip() if len(parts) > 1 else ""
    return out


def read_opts_ordered(run_dir):
    """read_opts(), but resolving flags set in more than one .opts correctly.

    pplib.read_opts merges the staged *.opts in sorted-name order and says so:
    that puts solver.opts LAST, so its `-flag_BC_rhovfix 0` beats an experiment
    file that set 1.  run_lunar.sh passes them solver -> geometry -> experiment
    (lines 74-76), so the experiment file is the one that wins.  The run folder
    is named "<geometry>__<experiment>", which is enough to recover that order.

    Only the vapour-BC flags are actually set twice; every physics flag this
    script reads lives in exactly one file, so read_opts would do for those.
    """
    if not os.path.isdir(run_dir):
        return {}
    files = sorted(f for f in os.listdir(run_dir) if f.endswith(".opts"))
    geom, _, exp = os.path.basename(os.path.normpath(run_dir)).partition("__")

    def rank(fn):
        stem = fn[: -len(".opts")]
        if stem == "solver":
            return 0
        if stem == exp:
            return 2
        return 1                      # geometry, or anything unrecognised

    opts = {}
    for fn in sorted(files, key=lambda f: (rank(f), f)):
        try:
            opts.update(_parse_opts_file(os.path.join(run_dir, fn)))
        except OSError:
            continue
    return opts


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------
def crossings(x, y, level=0.5):
    """Linear-interpolation crossings of y=level in y(x).  Never a grid node."""
    s = np.flatnonzero((y[:-1] - level) * (y[1:] - level) < 0.0)
    if s.size == 0:
        return np.empty(0)
    f = (level - y[s]) / (y[s + 1] - y[s])
    return x[s] + f * (x[s + 1] - x[s])


def centreline_weights(Y, y0):
    """(j, w) such that (1-w)*A[j,i] + w*A[j+1,i] is A on the line y = y0.

    The wedge patch is curvilinear: a column of constant index i sits at
    constant x (verified exactly), but a ROW of constant j does not sit at
    constant y and no row lands on Ly/2.  Y[:, i] is monotone, so the line is
    recovered column by column.  The mesh never moves, so this is computed once.
    """
    ny, nx = Y.shape
    j = np.empty(nx, dtype=np.intp)
    w = np.empty(nx)
    for i in range(nx):
        col = Y[:, i]
        k = min(max(int(np.searchsorted(col, y0)) - 1, 0), ny - 2)
        j[i] = k
        w[i] = (y0 - col[k]) / (col[k + 1] - col[k])
    return j, w


def on_centreline(A, j, w):
    i = np.arange(A.shape[1])
    return (1.0 - w) * A[j, i] + w * A[j + 1, i]


def contour_points(x, Y, phi, level=0.5):
    """phi=level crossings row by row, as physical (x, y) pairs.

    Returns (xl, yl, xr, yr) for the rows that have exactly two crossings --
    i.e. the left and right menisci sampled along their whole arcs.
    """
    d = phi - level
    hit = (d[:, :-1] * d[:, 1:]) < 0.0
    jj, ii = np.nonzero(hit)
    if jj.size == 0:
        return (np.empty(0),) * 4
    f = (level - phi[jj, ii]) / (phi[jj, ii + 1] - phi[jj, ii])
    px = x[ii] + f * (x[ii + 1] - x[ii])
    py = Y[jj, ii] + f * (Y[jj, ii + 1] - Y[jj, ii])

    n_per_row = np.bincount(jj, minlength=phi.shape[0])
    keep = n_per_row[jj] == 2
    if not keep.any():
        return (np.empty(0),) * 4
    jj, px, py = jj[keep], px[keep], py[keep]
    order = np.lexsort((px, jj))            # per row: left crossing first
    px, py = px[order], py[order]
    return px[0::2], py[0::2], px[1::2], py[1::2]


def circle_radius(px, py):
    """Least-squares circle through (px, py).  Returns (xc, yc, r) or None.

    Linear formulation 2*a*x + 2*b*y + c = x^2 + y^2, so no iteration and no
    initial guess.  Used only to VALIDATE the analytic apex radius -- second
    differences of the contour would be the noise-amplifying alternative that
    neck_width.py:18-23 rejects for exactly this geometry.
    """
    if px.size < 5:
        return None
    A = np.column_stack((2.0 * px, 2.0 * py, np.ones(px.size)))
    b = px ** 2 + py ** 2
    try:
        sol, *_ = np.linalg.lstsq(A, b, rcond=None)
    except np.linalg.LinAlgError:
        return None
    xc, yc = sol[0], sol[1]
    disc = sol[2] + xc ** 2 + yc ** 2
    if not np.isfinite(disc) or disc <= 0.0:
        return None
    return xc, yc, float(np.sqrt(disc))


def sigma_at_interface(x, sigma, x_iface, sign_out, eps, win):
    """sigma extrapolated to the interface from the vapour side.

    `sign_out` is +1 when vapour lies at larger x, -1 when it lies at smaller x.
    Fits sigma(x) over [win[0], win[1]] * eps out from x_iface and returns the
    fit's value at x_iface, so the diffuse band itself contributes nothing.
    """
    lo, hi = sorted((x_iface + sign_out * win[0] * eps,
                     x_iface + sign_out * win[1] * eps))
    m = (x >= lo) & (x <= hi)
    if m.sum() < 3:
        return float("nan")
    return float(np.polyval(np.polyfit(x[m] - x_iface, sigma[m], 1), 0.0))


def sliding_slope(t, y, win):
    """d y/d t from a local least-squares line over `win` samples.

    Adjacent-sample differencing of an interface position resolved to a small
    fraction of dx is dominated by interpolation jitter (~5% here); a windowed
    fit is what track_wedge_band.py uses for the same reason.
    """
    n = t.size
    out = np.full(n, np.nan)
    h = max(win // 2, 1)
    for k in range(n):
        a, b = max(0, k - h), min(n, k + h + 1)
        if b - a < 3:
            continue
        out[k] = np.polyfit(t[a:b] - t[k], y[a:b], 1)[0]
    return out


# ---------------------------------------------------------------------------
# Core
# ---------------------------------------------------------------------------
def analyse(run_dir, source="vtkOut", stride=1, sigma_win=(4.0, 10.0),
            curv_frac=0.2, vn_win=21, progress=True):
    opts = read_opts_ordered(run_dir)

    for dead in ("-mob_sub", "-alph_sub"):
        if opt_float(opts, dead) not in (None, 0.0):
            raise SystemExit(
                f"{dead} was overridden ({opts[dead]}).  mob_sub and alph_sub are\n"
                "derived from d0_sub0/beta_sub0 in lunar_main.c:822-833; setting either\n"
                "directly breaks that chain, so -d0_sub0/-beta_sub0 no longer describe\n"
                "this run's kinetics and the Gibbs-Thomson prediction is meaningless."
            )

    apex_x = opt_float(opts, "-wedge_apex_x")
    apex_y = opt_float(opts, "-wedge_apex_y")
    if apex_x is None or apex_y is None:
        raise SystemExit(
            "no -wedge_apex_x/-wedge_apex_y in this run's .opts -- not a wedge run.\n"
            "This diagnostic only applies to the apex-centred band geometry."
        )

    eps = opt_float(opts, "-eps")
    Ly = opt_float(opts, "-Ly")
    d0 = opt_float(opts, "-d0_sub0")
    beta = opt_float(opts, "-beta_sub0")
    T0 = opt_float(opts, "-temp")
    rho_air = opt_float(opts, "-rho_air", 1.341)
    missing = [k for k, v in (("-eps", eps), ("-Ly", Ly), ("-d0_sub0", d0),
                              ("-beta_sub0", beta), ("-temp", T0)) if v is None]
    if missing:
        raise SystemExit("missing required options in the run folder: "
                         + ", ".join(missing))

    y0 = 0.5 * Ly
    if abs(apex_y - y0) > 1e-12:
        print(f"  ⚠️  apex_y = {apex_y:.6e} is not Ly/2 = {y0:.6e}; on the centreline "
              "the interface normal is then only approximately +/- x_hat.",
              file=sys.stderr)

    vdir = os.path.join(run_dir, source)
    files = sorted(f for f in os.listdir(vdir)
                   if f.startswith("solV_") and f.endswith(".vts")) if os.path.isdir(vdir) else []
    if not files:
        raise SystemExit(f"no solV_*.vts under {vdir}.  Generate them first:\n"
                         f"    python3 postprocess/plot_fields.py --dir {run_dir}")
    files = files[::max(stride, 1)]

    times = step_times(run_dir)
    rho_vs_T0 = rho_vs(T0, rho_air)

    j = w = None
    rows = []
    for n, fn in enumerate(files):
        path = os.path.join(vdir, fn)
        fields, X, Y = read_vts(path)
        if j is None:
            x = X[0, :]
            if not np.allclose(X, x[None, :], atol=0.0, rtol=1e-12):
                raise SystemExit("mesh columns are not at constant x; the centreline "
                                 "extraction in this script assumes they are.")
            j, w = centreline_weights(Y, y0)

        phi = fields["IcePhase"]
        phi_c = on_centreline(phi, j, w)
        T_c = on_centreline(fields["Temperature"], j, w)
        if "Supersaturation" in fields:
            sig_c = on_centreline(fields["Supersaturation"], j, w)
        else:
            sig_c = supersaturation(on_centreline(fields["VaporDensity"], j, w),
                                    T_c, rho_air)

        c = crossings(x, phi_c, 0.5)
        step = step_of(fn)
        if c.size != 2:
            print(f"  ⚠️  step {step}: {c.size} centreline crossings, expected 2 "
                  "(band broken or ice off the centreline) -- snapshot skipped.",
                  file=sys.stderr)
            continue

        xl, xr = float(c[0]), float(c[1])
        rl, rr = xl - apex_x, xr - apex_x
        # Outward normal points AWAY from the ice: -x on the inner face, +x on
        # the outer one.  chi follows: inner meniscus concave, outer convex.
        chi_l, chi_r = -1.0 / rl, +1.0 / rr

        cxl, cyl, cxr, cyr = contour_points(x, Y, phi, 0.5)
        chi_l_fit = chi_r_fit = float("nan")
        for px, py, r_ana, sgn, slot in ((cxl, cyl, rl, -1.0, "l"),
                                         (cxr, cyr, rr, +1.0, "r")):
            if px.size == 0:
                continue
            m = np.abs(py - y0) < curv_frac * r_ana
            fit = circle_radius(px[m], py[m])
            if fit is None:
                continue
            if slot == "l":
                chi_l_fit = -1.0 / fit[2]
            else:
                chi_r_fit = +1.0 / fit[2]

        sig_l = sigma_at_interface(x, sig_c, xl, -1.0, eps, sigma_win)
        sig_r = sigma_at_interface(x, sig_c, xr, +1.0, eps, sigma_win)
        sig_l_dir = float(np.interp(xl, x, sig_c))
        sig_r_dir = float(np.interp(xr, x, sig_c))

        # sigma is normalised by the LOCAL rho_vs(T) but the solver's d0/beta
        # carry rho_vs(T0); the ratio is the bridge.  It is 1 for isothermal runs.
        fl = rho_vs(float(np.interp(xl, x, T_c)), rho_air) / rho_vs_T0
        fr = rho_vs(float(np.interp(xr, x, T_c)), rho_air) / rho_vs_T0

        rows.append(dict(
            step=step, time=times.get(step, float("nan")),
            r_left=rl, r_right=rr,
            chi_left=chi_l, chi_right=chi_r,
            chi_left_fit=chi_l_fit, chi_right_fit=chi_r_fit,
            sigma_left=sig_l, sigma_right=sig_r,
            sigma_left_direct=sig_l_dir, sigma_right_direct=sig_r_dir,
            vn_left_pred=(sig_l * fl - d0 * chi_l) / beta,
            vn_right_pred=(sig_r * fr - d0 * chi_r) / beta,
            vn_left_pred_direct=(sig_l_dir * fl - d0 * chi_l) / beta,
            vn_right_pred_direct=(sig_r_dir * fr - d0 * chi_r) / beta,
            vn_left_pred_chifit=(sig_l * fl - d0 * chi_l_fit) / beta,
            vn_right_pred_chifit=(sig_r * fr - d0 * chi_r_fit) / beta,
        ))
        if progress and (n % 100 == 0 or n == len(files) - 1):
            print(f"  {n + 1}/{len(files)} snapshots", flush=True)

    if len(rows) < 3:
        raise SystemExit("fewer than 3 usable snapshots -- nothing to differentiate.")

    data = {k: np.array([r[k] for r in rows], dtype=float) for k in rows[0]}
    if not np.isfinite(data["time"]).all():
        raise SystemExit("outp.txt did not supply a time for every snapshot.  Note "
                         "snapshot index is NOT the row index in SSA_evo.dat.")

    t = data["time"]
    # Both faces recede toward the apex, so a shrinking radius is GROWTH on the
    # inner face and RECESSION on the outer one.  Sign both positive-for-growth.
    data["vn_left_meas"] = -sliding_slope(t, data["r_left"], vn_win)
    data["vn_right_meas"] = +sliding_slope(t, data["r_right"], vn_win)

    meta = dict(run_dir=run_dir, source=source, eps=eps, Ly=Ly, d0=d0, beta=beta,
                T0=T0, rho_air=rho_air, apex_x=apex_x, apex_y=apex_y,
                n_snapshots=len(rows), sigma_win=sigma_win, curv_frac=curv_frac,
                vn_win=vn_win,
                rhovfix=(opts.get("-flag_BC_rhovfix", "0"),
                         opts.get("-rhovfix_lo", "-"), opts.get("-rhovfix_hi", "-")))
    return data, meta


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------
COLUMNS = ["step", "time", "r_left", "r_right", "chi_left", "chi_right",
           "chi_left_fit", "chi_right_fit", "sigma_left", "sigma_right",
           "sigma_left_direct", "sigma_right_direct",
           "vn_left_meas", "vn_right_meas", "vn_left_pred", "vn_right_pred"]


def write_csv(data, path):
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    arr = np.column_stack([data[c] for c in COLUMNS])
    np.savetxt(path, arr, delimiter=",", header=",".join(COLUMNS), comments="",
               fmt="%.10e")
    print(f"csv  -> {path}")


def _median(a):
    a = a[np.isfinite(a)]
    return float(np.median(a)) if a.size else float("nan")


def summarise(data, meta, skip):
    keep = slice(skip, None)
    d0, beta = meta["d0"], meta["beta"]
    flag, lo, hi = meta["rhovfix"]

    print()
    print("=" * 74)
    print("  Gibbs-Thomson interface velocity at y = Ly/2")
    print("=" * 74)
    print(f"  run           : {os.path.basename(os.path.normpath(meta['run_dir']))}")
    print(f"  source        : {meta['source']}/  ({meta['n_snapshots']} snapshots, "
          f"first {skip} dropped from the summary)")
    print(f"  apex          : ({meta['apex_x']:.4e}, {meta['apex_y']:.4e}) m")
    print(f"  eps           : {meta['eps']:.4e} m      T0 = {meta['T0']:.2f} C")
    print(f"  d0_sub0       : {d0:.4e} m       beta_sub0 = {beta:.4e} s/m  (UNSCALED)")
    print(f"  vapour BC     : -flag_BC_rhovfix {flag}, rhovfix_lo {lo}, rhovfix_hi {hi}")
    print("-" * 74)
    print(f"  {'':14s}{'LEFT (inner, concave)':>26s}{'RIGHT (outer, convex)':>26s}")

    def row(label, l, r, fmt="{:+.4e}"):
        print(f"  {label:14s}{fmt.format(l):>26s}{fmt.format(r):>26s}")

    row("sigma", _median(data["sigma_left"][keep]), _median(data["sigma_right"][keep]))
    row("d0*chi", d0 * _median(data["chi_left"][keep]),
        d0 * _median(data["chi_right"][keep]))
    row("residual", _median((data["sigma_left"] - d0 * data["chi_left"])[keep]),
        _median((data["sigma_right"] - d0 * data["chi_right"])[keep]))
    print("-" * 74)
    row("v_n measured", _median(data["vn_left_meas"][keep]),
        _median(data["vn_right_meas"][keep]))
    row("v_n predicted", _median(data["vn_left_pred"][keep]),
        _median(data["vn_right_pred"][keep]))

    ratio_l = _median((data["vn_left_meas"] / data["vn_left_pred"])[keep])
    ratio_r = _median((data["vn_right_meas"] / data["vn_right_pred"])[keep])
    row("meas / pred", ratio_l, ratio_r, "{:+.4f}")

    beta_eff_l = _median(((data["sigma_left"] - d0 * data["chi_left"])
                          / data["vn_left_meas"])[keep])
    beta_eff_r = _median(((data["sigma_right"] - d0 * data["chi_right"])
                          / data["vn_right_meas"])[keep])
    row("beta_eff [s/m]", beta_eff_l, beta_eff_r)
    row("beta_eff/beta0", beta_eff_l / beta, beta_eff_r / beta, "{:+.4f}")

    print("-" * 74)
    print("  Sensitivity of v_n,pred -- the residual above is a small difference of")
    print("  two nearly equal terms, so read the ratio with these in view:")

    def spread(pred, alt):
        v = np.abs((alt - pred) / pred)[keep]
        return 100.0 * _median(v)

    print(f"    sigma extrapolated vs. sampled at phi=0.5 : "
          f"{spread(data['vn_left_pred'], data['vn_left_pred_direct']):6.2f} % (L)   "
          f"{spread(data['vn_right_pred'], data['vn_right_pred_direct']):6.2f} % (R)")
    print(f"    chi analytic (apex) vs. circle-fit        : "
          f"{spread(data['vn_left_pred'], data['vn_left_pred_chifit']):6.2f} % (L)   "
          f"{spread(data['vn_right_pred'], data['vn_right_pred_chifit']):6.2f} % (R)")
    print("=" * 74)
    print()


def make_figure(data, meta, skip, save, title=None):
    t_all = data["time"]
    unit = auto_time_unit(float(np.nanmax(t_all)))
    t = in_time_unit(t_all, unit)
    d0 = meta["d0"]
    k = slice(skip, None)

    fig, ax = plt.subplots(3, 1, figsize=(9.0, 10.5), sharex=True)
    cl, cr = "#1f77b4", "#d62728"

    ax[0].axhline(0.0, color="0.7", lw=0.8)
    for c, side, lab in ((cl, "left", "left (inner, concave)"),
                         (cr, "right", "right (outer, convex)")):
        ax[0].plot(t[k], data[f"vn_{side}_meas"][k], "o", ms=2.5, color=c, alpha=0.5,
                   label=f"{lab} — measured")
        ax[0].plot(t[k], data[f"vn_{side}_pred"][k], "-", lw=1.8, color=c,
                   label=f"{lab} — Gibbs-Thomson")
    ax[0].set_ylabel(r"$v_n$  [m/s]   (+ = growth)")
    ax[0].legend(fontsize=8, ncol=2, loc="best")
    ax[0].set_title("Measured vs. predicted normal velocity", fontsize=10, loc="left")

    ax[1].axhline(1.0, color="0.4", lw=1.0, ls="--", label="perfect agreement")
    for c, side, lab in ((cl, "left", "left"), (cr, "right", "right")):
        ax[1].plot(t[k], (data[f"vn_{side}_meas"] / data[f"vn_{side}_pred"])[k],
                   "-", lw=1.4, color=c, label=lab)
    ax[1].set_ylabel(r"$v_{n,\mathrm{meas}} / v_{n,\mathrm{pred}}$")
    ax[1].set_ylim(0.0, 1.5)
    ax[1].legend(fontsize=8, loc="best")
    ax[1].set_title("Fraction of the Gibbs-Thomson rate the solver recovers",
                    fontsize=10, loc="left")

    ax[2].axhline(0.0, color="0.7", lw=0.8)
    for c, side, lab in ((cl, "left", "left"), (cr, "right", "right")):
        ax[2].plot(t[k], data[f"sigma_{side}"][k], "-", lw=1.5, color=c,
                   label=rf"$\sigma$ {lab}")
        ax[2].plot(t[k], d0 * data[f"chi_{side}"][k], "--", lw=1.5, color=c,
                   label=rf"$d_0\chi$ {lab}")
    ax[2].set_ylabel("driving-force terms  [-]")
    ax[2].set_xlabel(f"time [{unit}]")
    ax[2].legend(fontsize=8, ncol=2, loc="best")
    ax[2].set_title(r"$v_n$ is the residual $\sigma - d_0\chi$ of these curves",
                    fontsize=10, loc="left")

    for a in ax:
        a.grid(alpha=0.25)
    flag, lo, hi = meta["rhovfix"]
    head = title or os.path.basename(os.path.normpath(meta["run_dir"]))
    fig.suptitle(f"{head}\nrho_v BC: lo {lo}, hi {hi}  (flag {flag})   "
                 f"d0={meta['d0']:.3e} m, beta={meta['beta']:.3e} s/m",
                 fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    os.makedirs(os.path.dirname(os.path.abspath(save)) or ".", exist_ok=True)
    fig.savefig(save, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"plot -> {save}")


# ---------------------------------------------------------------------------
def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dir", default=".", help="run folder (default: cwd)")
    p.add_argument("--save", default=None,
                   help="figure path (default: <run>/plots/wedge_gt_velocity.png)")
    p.add_argument("--save-csv", default=None,
                   help=f"csv path (default: <run>/{CSV_NAME})")
    p.add_argument("--source", default="vtkOut", choices=("vtkOut", "vtkOut_highres"),
                   help="snapshot directory (default: vtkOut)")
    p.add_argument("--stride", type=int, default=1, help="use every Nth snapshot")
    p.add_argument("--skip", type=int, default=5,
                   help="snapshots dropped from the summary and figure; the IC is not "
                        "at its equilibrium profile (default: 5)")
    p.add_argument("--sigma-window", type=float, nargs=2, default=(4.0, 10.0),
                   metavar=("LO", "HI"),
                   help="vapour-side fit window for sigma, in multiples of eps "
                        "(default: 4 10)")
    p.add_argument("--curv-window", type=float, default=0.2,
                   help="half-height of the circle-fit window as a fraction of the "
                        "meniscus radius (default: 0.2)")
    p.add_argument("--vn-window", type=int, default=21,
                   help="snapshots in the sliding fit for measured v_n (default: 21)")
    p.add_argument("--title", default=None, help="figure title override")
    p.add_argument("--quiet", action="store_true", help="suppress progress lines")
    a = p.parse_args(argv)

    run = os.path.abspath(a.dir)
    data, meta = analyse(run, source=a.source, stride=a.stride,
                         sigma_win=tuple(a.sigma_window), curv_frac=a.curv_window,
                         vn_win=a.vn_window, progress=not a.quiet)

    write_csv(data, a.save_csv or os.path.join(run, CSV_NAME))
    make_figure(data, meta, a.skip,
                a.save or os.path.join(run, "plots", "wedge_gt_velocity.png"),
                a.title)
    summarise(data, meta, a.skip)
    return 0


if __name__ == "__main__":
    sys.exit(main())

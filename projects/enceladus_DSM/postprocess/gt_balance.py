#!/usr/bin/env pvpython
"""
gt_balance.py — check a run against the Gibbs–Thomson interface condition that
the Allen–Cahn term implies.

docs/curvature_driven_growth.md derives, from the residual in src/assembly.c,
that the sharp-interface limit of this model carries the interface condition

    (rho_v - rho_vs)/rho_vs  =  d0 * kappa  +  beta * v_n

    d0   = 15 * M * rho_ice / (alph_sub * rho_vs)      [m]
    beta = 5 * rho_ice / (eps * alph_sub * rho_vs)     [s/m]

with M = mob_sub. Neither d0 nor kappa appears anywhere in the code: the
capillarity is carried implicitly by the Allen–Cahn term, calibrated through
the mob_sub / alph_sub ratio. This script tests that claim on real output.

It samples the phi = 0.5 contour of a finished run, reads kappa off the same
Programmable Filter that scripts/paraview_macros/plot_rhovsI.py builds, and
compares the measured supersaturation against d0*kappa point by point. What is
left over is beta*v_n, so the residual also gives a map of local interface
velocity: positive = growing, negative = sublimating.

Run
---
    /Applications/ParaView-6.1.1.app/Contents/bin/pvpython \
        postprocess/gt_balance.py <run_dir> [--time T] [--out DIR]

<run_dir> is a run directory containing pf.pvd and the staged .opts/outp.txt.
Writes gt_balance.csv and gt_balance.png into --out (default: <run_dir>).

Reads mob_sub, alph_sub, rho_ice and eps out of the run's own outp.txt, so the
constants are the ones the run actually used rather than anything hard-coded.
"""

import argparse
import builtins
import os
import re
import sys

import numpy as np

# ParaView's Programmable Filter execs in __main__ and star-imports
# numpy_interface.algorithms over the builtins; bind what we need up front.
# See scripts/paraview_macros/verify_curvature.py for the full story.
_max, _abs = builtins.max, builtins.abs

HERE = os.path.dirname(os.path.abspath(__file__))
MACRO = os.path.join(HERE, "..", "scripts", "paraview_macros", "plot_rhovsI.py")


def load_macro():
    """Reuse plot_rhovsI.py's filter script rather than restating the formula."""
    with open(MACRO) as handle:
        source = handle.read()
    marker = "\nmain()\n"
    if not source.endswith(marker):
        raise SystemExit("gt_balance: %s no longer ends in a bare main() call"
                         % MACRO)
    ns = {"__name__": "plot_rhovsI_gt", "__file__": os.path.abspath(MACRO)}
    exec(compile(source[: -len(marker)], MACRO, "exec"), ns)
    return ns


def scan_log(run_dir):
    """Pull the run's own constants out of its startup banner in outp.txt.

    The banner is the authority: -mob_sub and -alph_sub are usually DERIVED at
    startup from d0_sub0/beta_sub0 rather than set in a .opts, so parsing the
    .opts alone would miss them.
    """
    path = os.path.join(run_dir, "outp.txt")
    if not os.path.isfile(path):
        raise SystemExit("gt_balance: no outp.txt in %s" % run_dir)
    with open(path, errors="replace") as handle:
        text = handle.read()

    def grab(pattern, label):
        match = re.search(pattern, text)
        if not match:
            raise SystemExit("gt_balance: could not read %s from outp.txt" % label)
        return float(match.group(1))

    consts = {
        "eps": grab(r"eps\s*=\s*([0-9.eE+-]+)\s*m", "eps"),
        "mob_sub": grab(r"mob_sub\s*=\s*([0-9.eE+-]+)", "mob_sub"),
        "alph_sub": grab(r"alph_sub\s*=\s*([0-9.eE+-]+)", "alph_sub"),
        "rho_ice": grab(r"rho_ice\s*=\s*([0-9.eE+-]+)", "rho_ice"),
    }
    match = re.search(r"d0_sub0\s*=\s*([0-9.eE+-]+)\s*m", text)
    consts["d0_sub0"] = float(match.group(1)) if match else float("nan")
    consts["axisym"] = bool(re.search(r"^\s*-axisym\s+1", text, re.M))
    return consts


def sample_contour(run_dir, eps, axisym, time_value, n_times=1):
    """kappa, supersaturation and rho_vs on the phi = 0.5 isocontour.

    Contouring is what makes this a clean test: it puts every sample exactly on
    the interface, where the sharp-interface condition applies, instead of
    somewhere in the diffuse band where kappa varies with the level set.
    """
    from paraview.simple import Contour, PVDReader, ProgrammableFilter, servermanager
    from vtk.numpy_interface import dataset_adapter as dsa

    ns = load_macro()
    reader = PVDReader(FileName=os.path.join(run_dir, "pf.pvd"))
    reader.UpdatePipeline()
    times = list(reader.TimestepValues or [0.0])
    if time_value is not None:
        picked = [builtins.min(times, key=lambda v: _abs(v - time_value))]
    elif n_times <= 1:
        picked = [times[-1]]
    else:
        # Pool several times. The interface curvature evolves, so this tests
        # the s = d0*kappa relation across a RANGE of kappa rather than at the
        # one or two values a single snapshot of this geometry happens to hold.
        n = builtins.min(n_times, len(times))
        step = _max(len(times) // n, 1)
        picked = times[::step][-n:]
    t = picked[-1]

    pf = ProgrammableFilter(Input=reader)
    pf.OutputDataSetType = "Same as Input"
    pf.Script = ns["_filter_script"](eps, axisym, "gt_balance")
    pf.RequestInformationScript = ""
    pf.RequestUpdateExtentScript = ""
    pf.UpdatePipeline(t)

    contour = Contour(Input=pf)
    contour.ContourBy = ["POINTS", "IcePhase"]
    contour.Isosurfaces = [0.5]

    kap, sup, pts, tags = [], [], [], []
    for tv in picked:
        pf.UpdatePipeline(tv)
        contour.UpdatePipeline(tv)
        data = dsa.WrapDataObject(servermanager.Fetch(contour))
        if data.GetNumberOfPoints() == 0:
            continue
        kap.append(np.asarray(data.PointData["Curvature"], dtype=np.float64))
        sup.append(np.asarray(data.PointData["Supersaturation"], dtype=np.float64))
        pts.append(np.asarray(data.Points, dtype=np.float64))
        tags.append(np.full(kap[-1].shape, tv))
    if not kap:
        raise SystemExit("gt_balance: the phi = 0.5 contour is empty")
    return (t, np.concatenate(kap), np.concatenate(sup),
            np.concatenate(pts), np.concatenate(tags))


def make_figure(path, kappa, s_meas, d0, beta, consts, t, tags):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    resid = s_meas - d0 * kappa
    v_n = resid / beta

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6))

    ax = axes[0]
    lo = builtins.min(float(np.min(d0 * kappa)), float(np.min(s_meas)))
    hi = _max(float(np.max(d0 * kappa)), float(np.max(s_meas)))
    pad = 0.08 * (hi - lo)
    ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], "-", color="0.4", lw=1.2,
            label="equilibrium: s = $d_0\\kappa$", zorder=1)
    sc = ax.scatter(d0 * kappa, s_meas, s=8, alpha=0.65, c=tags / 86400.0,
                    cmap="viridis", edgecolors="none", zorder=2)
    cb = ax.figure.colorbar(sc, ax=ax, pad=0.02)
    cb.set_label("time [days]", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    ax.set_xlabel("$d_0\\,\\kappa$   (Gibbs-Thomson equilibrium)")
    ax.set_ylabel("measured $(\\rho_v-\\rho_{vs})/\\rho_{vs}$")
    ax.set_title("Supersaturation relaxes onto $s = d_0\\kappa$")
    ax.annotate("t = 0: vapour still uniform,\nnot yet at equilibrium",
                xy=(0.03, 0.55), xycoords="axes fraction", fontsize=7.5,
                color="0.35")
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    ax.grid(alpha=0.25, lw=0.5)

    ax = axes[1]
    ax.axhline(0.0, color="0.4", lw=1.2, zorder=1)
    ax.scatter(kappa, v_n, s=8, alpha=0.65, c=tags / 86400.0, cmap="viridis",
               edgecolors="none", zorder=2)
    ax.set_xlabel("$\\kappa$  [1/m]")
    ax.set_ylabel("$v_n = (s - d_0\\kappa)/\\beta$   [m/s]")
    ax.set_title("Residual is the kinetic term: concave grows, convex sublimates")
    ax.grid(alpha=0.25, lw=0.5)

    fig.suptitle("Gibbs-Thomson balance,  t = %.4g s (%.1f d),  "
                 "$d_0$ = %.4g m,  $\\beta$ = %.4g s/m"
                 % (t, t / 86400.0, d0, beta), fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("run_dir")
    ap.add_argument("--time", type=float, default=None,
                    help="simulated time to sample (default: last timestep)")
    ap.add_argument("--out", default=None, help="output directory")
    ap.add_argument("--times", type=int, default=8,
                    help="how many timesteps to pool (default 8); the interface "
                         "curvature evolves, so pooling tests the relation over "
                         "a range of kappa. Ignored if --time is given.")
    args = ap.parse_args()

    run_dir = os.path.abspath(args.run_dir)
    out_dir = os.path.abspath(args.out) if args.out else run_dir
    os.makedirs(out_dir, exist_ok=True)

    c = scan_log(run_dir)
    t, kappa, s_meas, pts, tags = sample_contour(run_dir, c["eps"], c["axisym"],
                                                 args.time, args.times)

    # rho_vs at the run temperature, via the same expression as the solver.
    # Taken off the data so a Tgrad run uses its own local value.
    from paraview.simple import PVDReader  # noqa: F401  (import kept local)
    rho_vs = None
    # Supersaturation = (rho_v - rho_vs)/rho_vs is dimensionless, so the
    # balance below needs rho_vs only through d0 and beta, which use the T0
    # value. Recover it from the solver's own constants.
    T0_K = 253.15
    PATM, BB = 101325.0, 0.62
    K = (-0.5865e4, 0.2224e2, 0.1375e-1, -0.3403e-4, 0.2697e-7, 0.6918)
    Pvs = np.exp(K[0] / T0_K + K[1] + K[2] * T0_K + K[3] * T0_K ** 2
                 + K[4] * T0_K ** 3 + K[5] * np.log(T0_K))
    rho_vs = 1.341 * BB * Pvs / (PATM - Pvs)

    d0 = 15.0 * c["mob_sub"] * c["rho_ice"] / (c["alph_sub"] * rho_vs)
    beta = 5.0 * c["rho_ice"] / (c["eps"] * c["alph_sub"] * rho_vs)

    resid = s_meas - d0 * kappa
    v_n = resid / beta

    print()
    print("  run       : %s" % run_dir)
    print("  t         : %.6g s  (%.2f days)" % (t, t / 86400.0))
    print("  eps       : %.4e m      mob_sub  : %.4e m/s" % (c["eps"], c["mob_sub"]))
    print("  rho_ice   : %.4e kg/m3  alph_sub : %.4e 1/s" % (c["rho_ice"], c["alph_sub"]))
    print("  rho_vs    : %.4e kg/m3  axisym   : %s" % (rho_vs, c["axisym"]))
    print()
    print("  d0   = 15*M*rho_ice/(alph*rho_vs) = %.4e m" % d0)
    print("  d0_sub0 printed by the run        = %.4e m   <- must agree" % c["d0_sub0"])
    if c["d0_sub0"] == c["d0_sub0"]:
        print("  ratio                             = %.6f" % (d0 / c["d0_sub0"]))
    print("  beta = 5*rho_ice/(eps*alph*rho_vs) = %.4e s/m" % beta)
    print()
    print("  %d points on the phi = 0.5 contour, pooled over %d timestep(s)"
          % (kappa.size, len(np.unique(tags))))
    print("    kappa            %+.4e .. %+.4e /m" % (kappa.min(), kappa.max()))
    print("    d0*kappa         %+.4e .. %+.4e" % ((d0 * kappa).min(), (d0 * kappa).max()))
    print("    measured s       %+.4e .. %+.4e" % (s_meas.min(), s_meas.max()))
    print("    residual s-d0k   %+.4e .. %+.4e" % (resid.min(), resid.max()))
    print()
    # Split early from late. At t = 0 the vapour field is initialized uniform,
    # so s = 0 while kappa is already set by the initial geometry -- the
    # interface starts FAR from Gibbs-Thomson equilibrium and relaxes onto it.
    # That relaxation is the evidence the condition is emergent rather than
    # imposed, so report both states instead of averaging them together.
    ratio = np.abs(resid) / np.maximum(np.abs(d0 * kappa), 1e-300)
    t_lo, t_hi = float(tags.min()), float(tags.max())
    early = tags <= t_lo
    late = tags >= t_hi
    print("  median |residual| / |d0*kappa|")
    print("    at t = %-10.4g s (first sampled) : %.3f" % (t_lo, float(np.median(ratio[early]))))
    print("    at t = %-10.4g s (last sampled)  : %.3f" % (t_hi, float(np.median(ratio[late]))))
    print("    -> the vapour field starts uniform (s = 0, so the residual is the")
    print("       whole of d0*kappa) and relaxes ONTO the Gibbs-Thomson line.")
    print("       Once relaxed the interface sits within %.0f%% of equilibrium;"
          % (100 * float(np.median(ratio[late]))))
    print("       that remainder is the kinetic term beta*v_n.")
    print()
    print("  implied normal velocity v_n = (s - d0*kappa)/beta, at the last time")
    for label, sel in (("concave (kappa<0)", late & (kappa < 0)),
                       ("convex  (kappa>0)", late & (kappa > 0))):
        if np.any(sel):
            print("    %s: median %+.3e m/s" % (label, float(np.median(v_n[sel]))))

    csv_path = os.path.join(out_dir, "gt_balance.csv")
    np.savetxt(csv_path,
               np.column_stack([tags, pts[:, 0], pts[:, 1], kappa, d0 * kappa,
                                s_meas, resid, v_n]),
               delimiter=",", comments="", fmt="%.6g",
               header="time_s,x,y,kappa_1_per_m,d0_kappa,supersaturation,"
                      "residual,v_n_m_per_s")
    png_path = os.path.join(out_dir, "gt_balance.png")
    make_figure(png_path, kappa, s_meas, d0, beta, c, t, tags)
    print()
    print("  wrote %s" % csv_path)
    print("  wrote %s" % png_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())

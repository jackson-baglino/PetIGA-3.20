#!/usr/bin/env python3
"""Normal velocity of the two menisci of a channel-spanning ice bridge.

WHAT IS MEASURED

At the channel mid-plane the meniscus normal is exactly +/- x_hat, so the
interface's normal velocity there is just dx/dt of the phi = 0.5 crossing --
no projection, no curvature bookkeeping. That is the same trick
wedge_gt_velocity.py uses on the wedge centreline, and the reason both scripts
measure on a line rather than over the whole arc.

Sign convention: v_n > 0 means the interface advances INTO THE VAPOUR, i.e. the
ice grows. The left meniscus moving to smaller x and the right meniscus moving
to larger x therefore both count as positive.

Three things are reported per snapshot:
  - x_left, x_right  : meniscus positions at the mid-plane
  - v_left, v_right  : their normal velocities
  - area             : ice area, and its rate dA/dt

The area rate is the integral check on the two point velocities: if the bridge
kept a rectangular cross-section, dA/dt would be H*(v_left + v_right). It does
not exactly -- the menisci are curved -- so the two disagree by a shape factor,
and that disagreement is itself informative about how far the bridge is from its
equilibrium shape.

POSITION EXTRACTION

Interface positions are refined with pplib.refine_tanh rather than taken from
the raw phi = 0.5 crossing. v_n is a time derivative of this position, so a bias
that repeats with the sample grid appears as a periodic ripple in the velocity;
fitting atanh(2*phi - 1), which is exactly linear in x for this model's
half-normalised well, removes most of it.

Usage:
    python3 postprocess/meniscus_velocity.py --dir <run> [--save fig.png]
    python3 postprocess/meniscus_velocity.py --dir <legA> --dir <legB> ...
        (multiple --dir stitch a restart chain in time order)
"""
import argparse
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pplib import refine_tanh, read_opts, opt_float, step_times   # noqa: E402

CSV_NAME = "meniscus_velocity.csv"


def snapshots(run_dir, nu, nv):
    """Yield (step, x1d, y1d, phi[ix,iy]) over the patch's own knot range."""
    from igakit.io import PetIGA
    io = PetIGA()
    nrb = io.read(os.path.join(run_dir, "igasol.dat"))
    rng = lambda i: (float(nrb.knots[i][nrb.degree[i]]),
                     float(nrb.knots[i][-nrb.degree[i] - 1]))
    u = np.linspace(*rng(0), nu)
    v = np.linspace(*rng(1), nv)
    for f in sorted(x for x in os.listdir(run_dir)
                    if re.fullmatch(r"sol_\d+\.dat", x)):
        C, F = nrb(u, v, fields=io.read_vec(os.path.join(run_dir, f), nrb))
        yield (int(re.search(r"sol_(\d+)", f).group(1)),
               C[:, 0, 0].copy(), C[0, :, 1].copy(), F[:, :, 0].copy())


def measure(x1d, y1d, phi, eps, H):
    """Per-snapshot geometry of the bridge.

    Returns a dict, or None if there is no bridge.

    Three different velocities can be extracted from a meniscus and they are
    NOT interchangeable, which the first version of this script got wrong:

      wall   -- the CONTACT LINE, at y = 0 and y = H. This is where a wetting
                bridge grows: the contact line advances along the wall.
      mid    -- the mid-plane crossing. For a CONCAVE meniscus (theta < 90) this
                is the narrowest point of the bridge, so while the bridge is
                still relaxing toward its equilibrium shape the mid-plane
                RECEDES even as the ice grows. Measuring only here reports
                shrinkage for a growing bridge.
      mean   -- dA/dt divided by the total interface length. This is the
                honest single number: the mean normal velocity of the whole
                ice-vapour interface, sign-consistent with growth by
                construction.

    Once the shape has equilibrated all three converge, because the interface
    then translates without changing shape. Their spread is therefore a direct
    readout of how far the bridge still is from its equilibrium shape.
    """
    ny = y1d.size
    mid = ny // 2

    def crossings_in(col):
        c = np.flatnonzero((col[:-1] - 0.5) * (col[1:] - 0.5) < 0.0)
        if c.size < 2:
            return None
        def refine(i):
            f = (0.5 - col[i]) / (col[i + 1] - col[i])
            return refine_tanh(x1d, col, x1d[i] + f * (x1d[i + 1] - x1d[i]), eps)
        return refine(c[0]), refine(c[-1])

    m = crossings_in(phi[:, mid])
    w0 = crossings_in(phi[:, 0])
    w1 = crossings_in(phi[:, -1])
    if m is None or w0 is None or w1 is None:
        return None

    # ice area, and the length of the two menisci (the phi=0.5 contour,
    # excluding the wall segments -- those are ice/regolith, not ice/vapour)
    area = np.trapezoid(np.trapezoid(phi, y1d, axis=1), x1d)
    XL, XR = [], []
    for j in range(ny):
        c = crossings_in(phi[:, j])
        if c is not None:
            XL.append(c[0]); XR.append(c[1])
    XL, XR = np.array(XL), np.array(XR)
    yy = y1d[:len(XL)] if len(XL) == ny else np.linspace(y1d[0], y1d[-1], len(XL))
    arc = (np.sum(np.hypot(np.diff(XL), np.diff(yy)))
           + np.sum(np.hypot(np.diff(XR), np.diff(yy))))
    return dict(x_mid_l=m[0], x_mid_r=m[1],
                x_wall_l=0.5 * (w0[0] + w1[0]), x_wall_r=0.5 * (w0[1] + w1[1]),
                area=area, arc=arc)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", action="append", required=True,
                    help="run directory; repeat to stitch a restart chain")
    ap.add_argument("--save", default=None)
    ap.add_argument("--nu", type=int, default=1501)
    ap.add_argument("--nv", type=int, default=401)
    args = ap.parse_args()

    opts = read_opts(args.dir[0])
    eps = opt_float(opts, "-eps")
    H = opt_float(opts, "-Ly")
    if eps is None or H is None:
        raise SystemExit("could not read -eps and -Ly from the staged .opts")

    rows = []
    for run in args.dir:
        o = read_opts(run)
        times = step_times(run)
        for step, x1d, y1d, phi in snapshots(run, args.nu, args.nv):
            g = measure(x1d, y1d, phi, eps, H)
            if g is None:
                continue
            t = times.get(step, np.nan)
            if np.isfinite(t):
                g["t"] = t
                rows.append(g)
    if len(rows) < 3:
        raise SystemExit("need at least 3 measurable snapshots")

    rows.sort(key=lambda r: r["t"])
    t  = np.array([r["t"] for r in rows])
    A  = np.array([r["area"] for r in rows])
    L  = np.array([r["arc"] for r in rows])
    xwl= np.array([r["x_wall_l"] for r in rows]); xwr=np.array([r["x_wall_r"] for r in rows])
    xml= np.array([r["x_mid_l"] for r in rows]);  xmr=np.array([r["x_mid_r"] for r in rows])

    # v_n > 0 = advancing into the vapour. Left menisci move to smaller x.
    v_wall = 0.5 * (-np.gradient(xwl, t) + np.gradient(xwr, t))
    v_mid  = 0.5 * (-np.gradient(xml, t) + np.gradient(xmr, t))
    dA     = np.gradient(A, t)
    v_mean = dA / L                      # mean normal speed of the whole interface

    out = os.path.join(args.dir[-1], CSV_NAME)
    with open(out, "w") as fh:
        fh.write("# v_n > 0 = interface advances into the vapour (ice grows)\n")
        fh.write("# H = %.6e m, eps = %.6e m\n" % (H, eps))
        fh.write("# v_mean = (dA/dt)/L_interface is the honest growth velocity;\n")
        fh.write("# v_wall and v_mid are local and differ while the shape relaxes.\n")
        fh.write("time_s,time_d,area_m2,arc_len_m,x_wall_l,x_wall_r,x_mid_l,x_mid_r,"
                 "v_wall_m_s,v_mid_m_s,v_mean_m_s,dA_dt_m2_s\n")
        for i in range(len(t)):
            fh.write("%.6e,%.6f,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n"
                     % (t[i], t[i]/86400.0, A[i], L[i], xwl[i], xwr[i], xml[i], xmr[i],
                        v_wall[i], v_mid[i], v_mean[i], dA[i]))

    sel = t > 0.5 * t.max()
    print(f"\n  MENISCUS VELOCITY  ({os.path.basename(args.dir[-1].rstrip('/'))})")
    print(f"    span          {t.min()/86400:.1f} -> {t.max()/86400:.1f} d, {len(t)} snapshots")
    print(f"    ice area      {A[0]:.6e} -> {A[-1]:.6e} m^2  ({100*(A[-1]-A[0])/A[0]:+.3f} %)")
    print(f"    contact line  {xwr[0]-xwl[0]:.4e} -> {xwr[-1]-xwl[-1]:.4e} m wide")
    print(f"    mid-plane     {xmr[0]-xml[0]:.4e} -> {xmr[-1]-xml[-1]:.4e} m wide")
    print(f"\n    second half of the run:")
    print(f"      v_mean  (dA/dt / L)   {v_mean[sel].mean():+.4e} m/s   <-- the growth velocity")
    print(f"      v_wall  (contact line){v_wall[sel].mean():+.4e} m/s")
    print(f"      v_mid   (mid-plane)   {v_mid[sel].mean():+.4e} m/s")
    sp = abs(v_wall[sel].mean()-v_mid[sel].mean())/max(abs(v_mean[sel].mean()),1e-30)
    print(f"      spread / v_mean       {sp:.2f}   "
          f"({'shape still relaxing' if sp>0.5 else 'shape ~equilibrated'})")
    print(f"    -> {out}")

    if args.save:
        fig, ax = plt.subplots(3, 1, figsize=(7.6, 8.6), sharex=True,
                               constrained_layout=True)
        d = t / 86400.0
        ax[0].plot(d, A * 1e12, color="#2ca02c")
        ax[0].set_ylabel("ice area [µm²]"); ax[0].grid(alpha=.3)
        ax[0].set_title("ice area, and the three interface velocities")
        ax[1].plot(d, v_mean, lw=2, color="#1f77b4", label="$v_{mean}$ = (dA/dt)/L")
        ax[1].plot(d, v_wall, lw=1.2, ls="--", color="#ff7f0e", label="$v$ contact line")
        ax[1].plot(d, v_mid,  lw=1.2, ls=":",  color="#d62728", label="$v$ mid-plane")
        ax[1].axhline(0, color="k", lw=.8)
        ax[1].set_ylabel("$v_n$ [m/s]   (+ = ice grows)")
        ax[1].legend(fontsize=8); ax[1].grid(alpha=.3)
        ax[2].plot(d, (xwr - xwl) * 1e6, label="contact-line width")
        ax[2].plot(d, (xmr - xml) * 1e6, label="mid-plane width")
        ax[2].set_ylabel("width [µm]"); ax[2].set_xlabel("time [days]")
        ax[2].legend(fontsize=8); ax[2].grid(alpha=.3)
        fig.savefig(args.save, dpi=150)
        print(f"    -> {args.save}")


if __name__ == "__main__":
    main()

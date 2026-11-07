#!/usr/bin/env python3
"""Normal velocity of the two menisci of a channel-spanning ice bridge.

WHAT IS MEASURED

At the channel mid-plane the meniscus normal is exactly +/- x_hat, so the
interface's normal velocity there is just dx/dt of the phi = 0.5 crossing --
no projection, no curvature bookkeeping. That is the same trick
wedge_gt_velocity.py uses on the wedge centreline, and the reason both scripts
measure on a line rather than over the whole arc.

Sign convention: v_n > 0 means the interface advances INTO THE VAPOR, i.e. the
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
from pplib import node_weights, opt_float, read_opts, refine_tanh, step_times   # noqa: E402

CSV_NAME = "meniscus_velocity.csv"


def snapshots(run_dir, nu, nv):
    """Yield (step, x1d, Y, phi[ix,iy]) over the patch's own knot range.

    Y is the full (nu, nv) array of node heights: on a mapped (wedge) mesh the
    channel width changes with x, so no single column of y values describes
    the grid. Columns are assumed to sit at constant x (x1d), as on every
    channel and wedge mesh these runs use."""
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
               C[:, 0, 0].copy(), C[:, :, 1].copy(), F[:, :, 0].copy())


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
                ice-vapor interface, sign-consistent with growth by
                construction.

    Once the shape has equilibrated all three converge, because the interface
    then translates without changing shape. Their spread is therefore a direct
    readout of how far the bridge still is from its equilibrium shape.
    """
    Y = np.asarray(y1d, dtype=float)
    if Y.ndim == 1:                                    # rectangular mesh
        Y = np.broadcast_to(Y, (x1d.size, Y.size))
    ny = Y.shape[1]
    mid = ny // 2

    def all_crossings(col):
        c = np.flatnonzero((col[:-1] - 0.5) * (col[1:] - 0.5) < 0.0)
        out = []
        for i in c:
            f = (0.5 - col[i]) / (col[i + 1] - col[i])
            out.append(refine_tanh(x1d, col, x1d[i] + f * (x1d[i + 1] - x1d[i]), eps))
        return out

    m = all_crossings(phi[:, mid])
    w0 = all_crossings(phi[:, 0])
    w1 = all_crossings(phi[:, -1])
    if len(m) < 2 or len(w0) < 2 or len(w1) < 2:
        return None

    # Ice area from the 2-D node coordinates (exact for a mapped mesh too).
    X = np.broadcast_to(np.asarray(x1d, dtype=float)[:, None], Y.shape)
    area = float((node_weights(X, Y) * phi).sum())

    # Length of every ice-vapor face (the phi = 0.5 contour, excluding the
    # wall segments, which are ice/regolith). A bridge that has split into
    # two plugs has four faces; counting only the outer two would double
    # v_mean. Rows are grouped by how many crossings they have, and each
    # face k is the k-th crossing of the rows that have the usual count.
    rows = [(j, all_crossings(phi[:, j])) for j in range(ny)]
    counts = [len(c) for _, c in rows if c]
    n_faces = max(set(counts), key=counts.count) if counts else 0
    sel = [(j, c) for j, c in rows if len(c) == n_faces]
    arc = 0.0
    for k in range(n_faces):
        xs = np.array([c[k] for _, c in sel])
        ys = np.array([np.interp(c[k], x1d, Y[:, j]) for j, c in sel])
        arc += float(np.sum(np.hypot(np.diff(xs), np.diff(ys))))
    return dict(x_mid_l=m[0], x_mid_r=m[-1],
                x_wall_l=0.5 * (w0[0] + w1[0]), x_wall_r=0.5 * (w0[-1] + w1[-1]),
                area=area, arc=arc, n_faces=n_faces)


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
        times = step_times(run)
        for step, x1d, y1d, phi in snapshots(run, args.nu, args.nv):
            g = measure(x1d, y1d, phi, eps, H)
            if g is None:
                continue
            t = times.get(step, np.nan)
            if np.isfinite(t):
                g["t"] = t
                rows.append(g)
    # A restart leg rewrites the previous leg's last state as its step 0 at the
    # same time; keep one row per time (the later leg's), or np.gradient
    # divides by dt = 0.
    by_t = {}
    for r in rows:
        by_t[r["t"]] = r
    rows = [by_t[k] for k in sorted(by_t)]
    if len(rows) < 3:
        raise SystemExit("need at least 3 measurable snapshots")
    split = sorted({r["n_faces"] for r in rows if r["n_faces"] != 2})
    if split:
        print(f"  NOTE: some snapshots have {split} ice-vapor faces per row, not 2 (the bridge "
              "split or a face reached a wall); x_mid/x_wall use the outermost faces, "
              "v_mean uses every face.", file=sys.stderr)
    t  = np.array([r["t"] for r in rows])
    A  = np.array([r["area"] for r in rows])
    L  = np.array([r["arc"] for r in rows])
    xwl= np.array([r["x_wall_l"] for r in rows]); xwr=np.array([r["x_wall_r"] for r in rows])
    xml= np.array([r["x_mid_l"] for r in rows]);  xmr=np.array([r["x_mid_r"] for r in rows])

    # v_n > 0 = advancing into the vapor. Left menisci move to smaller x.
    v_wall = 0.5 * (-np.gradient(xwl, t) + np.gradient(xwr, t))
    v_mid  = 0.5 * (-np.gradient(xml, t) + np.gradient(xmr, t))
    dA     = np.gradient(A, t)
    v_mean = dA / L                      # mean normal speed of the whole interface

    out = os.path.join(args.dir[-1], CSV_NAME)
    with open(out, "w") as fh:
        fh.write("# v_n > 0 = interface advances into the vapor (ice grows)\n")
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
    print("\n    second half of the run:")
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

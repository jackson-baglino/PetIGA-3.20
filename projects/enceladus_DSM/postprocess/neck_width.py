#!/usr/bin/env python3
"""Extract neck width vs time from a two-grain sintering run.

Method: MINIMUM CROSS-SECTION, not curvature tracking.

For each snapshot, the ice body's vertical chord width w(x) is measured on
every mesh column (the span between the outermost phi = 0.5 crossings,
sub-cell interpolated). Restricted to the region strictly between the two
grain-center columns, w(x) has a single interior minimum — the neck plane —
and the neck width is min w(x). This definition:

  * tracks the neck plane automatically as it migrates toward the smaller
    grain (the minimum moves with it; nothing is anchored to the t = 0
    location);
  * matches the experimental convention (Molaro et al. report the neck
    WIDTH, i.e. the full waist of the bonded pair);
  * needs no derivatives of the contour. The proposed alternative — locate
    the two curvature maxima flanking the neck — identifies the same plane
    but requires second differences of a discretized contour, which is
    noise-amplifying at exactly the saddle points where the contour is
    sparsest. The waist minimum is the integral-robust version of the same
    idea. (If the rim points themselves are ever needed, e.g. for rim-radius
    measurements, add them as a separate diagnostic then.)

Caveat: if a run coarsens so far that the waist disappears (single convex
body), the interior minimum merges with a grain shoulder; the script reports
the minimum location so that regime is visible in the output.

Usage:
    python neck_width.py <run_dir> [--phi 0.5] [--out neck_width.csv]

Reads <run_dir>/vtkOut/solV_*.vts (snapshot index = timestep, -outp 1) and
<run_dir>/outp.txt (step -> time map from the monitor tables). Writes a CSV
(t, neck_width, x_neck) and a PNG plot next to it.
"""

import argparse
import glob
import math
import re
import sys
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np
import pplib
from pplib import step_times


def read_vts_phi(fn):
    """Return (phi[ny,nx], x[nx], y[ny]) from a solV_*.vts snapshot.

    Delegates to pplib.read_vts. This used to be a local copy of the parser
    that called a `decode` helper it never imported, so EVERY snapshot raised
    NameError -- and the caller's broad `except Exception` reported that as
    "empty/corrupt file, likely an incomplete transfer" and carried on to write
    an empty CSV. The script lived in postprocess/scratch/ at the time, so
    nothing exercised it after pplib absorbed the decoder.
    """
    fields, X, Y = pplib.read_vts(fn, want=("IcePhase",))
    if "IcePhase" not in fields:
        raise RuntimeError(f"IcePhase missing in {fn}")
    return fields["IcePhase"], X[0, :], Y[:, 0]


def chord_width(col, y, level):
    """Full span between outermost `level`-crossings of phi along a column."""
    above = col >= level
    if not above.any():
        return 0.0
    idx = np.flatnonzero(above)
    lo_i, hi_i = idx[0], idx[-1]
    # sub-cell interpolation at both ends
    y_lo = y[lo_i]
    if lo_i > 0:
        f = (level - col[lo_i - 1]) / (col[lo_i] - col[lo_i - 1])
        y_lo = y[lo_i - 1] + f * (y[lo_i] - y[lo_i - 1])
    y_hi = y[hi_i]
    if hi_i < len(y) - 1:
        f = (col[hi_i] - level) / (col[hi_i] - col[hi_i + 1])
        y_hi = y[hi_i] + f * (y[hi_i + 1] - y[hi_i])
    return y_hi - y_lo


def refine_min(w, x, jn):
    """Sub-grid neck: parabola vertex through the 3 columns bracketing w[jn].

    min(w) over discrete columns is a BIASED estimator -- a sampled minimum can
    only ever overshoot the true one, by roughly |w'|*dx/2 near a kink or
    w''*dx^2/8 near a smooth minimum. That matters here because w(x) at the neck
    is steep (slopes of 13-21 measured on the 2026-07-15 union series), so a
    fraction-of-a-cell misalignment becomes a ~1 um width error, and the
    misalignment differs per mesh -- manufacturing an apparent eps-dependence
    out of fields that actually agree.

    The vertex fit removes the smooth-minimum part of that bias exactly. It does
    NOT fix a genuine crease (a V has no parabola): on the union-IC t=0 crease
    it recovers only ~0.4 um of a 2.7 um spread, because ~2 um of that spread is
    a real O(eps) difference -- the sharp crease is rounded by the spline basis
    at the scale of eps. Once the fillet rounds (t >~ 1 min) the minimum IS
    smooth and this is accurate: on that series the refined and 40x-upsampled
    measurements agree to 0.01 um at t_final.
    """
    if jn <= 0 or jn >= len(w) - 1:
        return w[jn], x[jn]
    y0, y1, y2 = w[jn - 1], w[jn], w[jn + 1]
    den = y0 - 2.0 * y1 + y2
    if den <= 0.0:                      # not a convex minimum -> keep the sample
        return w[jn], x[jn]
    off = 0.5 * (y0 - y2) / den         # vertex offset in cells, |off| <= 0.5
    if abs(off) > 0.5:
        return w[jn], x[jn]
    dx = x[jn + 1] - x[jn - 1]
    return y1 - 0.25 * (y0 - y2) * off, x[jn] + off * 0.5 * dx


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir", type=Path)
    ap.add_argument("--phi", type=float, default=0.5, help="contour level")
    ap.add_argument("--out", type=Path, default=None, help="CSV path")
    ap.add_argument("--data", type=Path, default=None,
                    help="experimental neck-width CSV to overlay "
                         "(default: the repo's Molaro Fig. 11 table, if found)")
    ap.add_argument("--axisym", action="store_true",
                    help="axisymmetric r-z run: the grid's y is the radius and "
                         "the axis is y = 0, so the measured chord (axis to the "
                         "phi contour) is the neck RADIUS; report width = 2x.")
    args = ap.parse_args()

    files = sorted(glob.glob(str(args.run_dir / "vtkOut" / "solV_*.vts")))
    if not files:
        sys.exit(f"no solV_*.vts under {args.run_dir}/vtkOut")
    tmap = step_times(args.run_dir / "outp.txt")

    rows = []
    centers = None
    for fn in files:
        step = int(Path(fn).stem.split("_")[1])
        try:
            phi, x, y = read_vts_phi(fn)
        except (ET.ParseError, OSError, ValueError, RuntimeError) as e:
            # Deliberately NOT a bare `except Exception`. That is how a
            # NameError in the reader spent this script's whole life being
            # reported as a truncated download, once per snapshot, while the
            # run still exited 0 with an empty CSV. Anything outside this
            # tuple is a bug in the reader and should crash loudly.
            print(f"  WARNING: skipping {Path(fn).name} ({e}) — "
                  f"empty/corrupt file, likely an incomplete transfer",
                  file=sys.stderr)
            continue
        w = np.array([chord_width(phi[:, j], y, args.phi) for j in range(len(x))])
        if args.axisym:
            # chord runs from the axis (y = 0, inside the grain) up to the
            # contour: it IS the local radius; physical width is its double.
            w = 2.0 * w

        if centers is None:
            # Grain centers = the two prominent local maxima of w(x). A local
            # max is a column beating its +-5-column neighborhood and taller
            # than 0.5*max(w); the leftmost and rightmost such columns are
            # the two grain centers. (Do NOT use the global minimum to split:
            # the grain TIPS have the smallest nonzero w and would hijack it.)
            thr = 0.5 * w.max()
            peaks = [j for j in range(5, len(w) - 5)
                     if w[j] >= thr and w[j] == w[j - 5:j + 6].max()]
            if len(peaks) < 2:
                sys.exit("could not locate two grain-center peaks in w(x)")
            centers = (peaks[0], peaks[-1])

        lo, hi = centers
        interior = np.arange(lo + 1, hi)
        interior = interior[w[interior] > 0]
        if len(interior) == 0:
            neck, xneck = 0.0, np.nan
        else:
            jn = interior[np.argmin(w[interior])]
            neck, xneck = refine_min(w, x, jn)
        rows.append((tmap.get(step, np.nan), neck, xneck))

    if not rows:
        sys.exit(f"every snapshot under {args.run_dir}/vtkOut failed to read — "
                 f"see the warnings above; refusing to write an empty CSV")

    out = args.out or (args.run_dir / "neck_width.csv")
    with open(out, "w") as f:
        f.write("t_s,neck_width_m,x_neck_m\n")
        for t, nw, xn in rows:
            f.write(f"{t:.6e},{nw:.6e},{xn:.6e}\n")
    print(f"csv -> {out}")

    t = np.array([r[0] for r in rows]) / 60.0        # minutes
    nw = np.array([r[1] for r in rows]) * 1e6        # um
    print(f"neck width: {nw[0]:.2f} um (t=0) -> {nw[-1]:.2f} um "
          f"(t={t[-1]:.0f} min): {100*(nw[-1]/nw[0]-1):+.1f}%")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(t, nw, lw=2, color="#3d74d9", label="model")

    # --- power-law fit, for an at-a-glance check against r ~ t^(1/3) ---------
    #
    # Fitted on the WIDTH rather than the radius on purpose: a is invariant
    # under scaling, so the exponent is identical either way and there is no
    # factor-of-2 convention to get wrong here.
    #
    # Uses Demmenie's own form, r = C*(t+t0)^a with t0 free, because our clock
    # zero is not theirs: these runs start from a finite neck, so forcing t0=0
    # would report a spuriously small exponent. Restricted to the window where
    # the neck fillet (rho ~ r^2/2R) is actually resolved, r/R >= sqrt(12*eps/R)
    # -- the same floor fit_neck_growth.py uses. Below it the measurement is a
    # discretisation artifact and fitting through it is meaningless.
    try:
        from fit_neck_growth import fit_d_free
        opts = pplib.read_opts(str(args.run_dir))
        eps_m = pplib.opt_float(opts, "-eps")
        radii = opts.get("-ice_grain_R", "")
        R0_m = (float(np.mean([float(v) for v in radii.split(",")]))
                if radii else None)
        floor_um = (2.0 * math.sqrt(12.0 * eps_m / R0_m) * R0_m * 1e6
                    if (eps_m and R0_m) else 0.0)          # as a WIDTH, in um

        ts = np.array([r[0] for r in rows])                # seconds
        sel = (ts > 0) & (nw >= floor_um) & np.isfinite(nw)
        fit = fit_d_free(ts[sel], nw[sel]) if sel.sum() >= 4 else None
        if fit:
            tf = np.linspace(ts[sel][0], ts[sel][-1], 200)
            ax.plot(tf / 60.0,
                    fit["C"] * (tf + fit["t0"]) ** fit["a"],
                    lw=1.6, ls="--", color="#d1495b", zorder=5,
                    label=(f"fit $r\\propto(t+t_0)^{{a}}$: "
                           f"a = {fit['a']:.3f} ± {fit['a_ci']:.3f}"))
            # 1/3 reference through the same start point, for visual scale.
            a3, t0 = 1.0 / 3.0, fit["t0"]
            C3 = nw[sel][0] / (ts[sel][0] + t0) ** a3
            ax.plot(tf / 60.0, C3 * (tf + t0) ** a3, lw=1.2, ls=":",
                    color="#3f434a", zorder=4,
                    label="$t^{1/3}$ (sublimation–condensation)")
            if floor_um:
                ax.axhline(floor_um, lw=0.8, ls=(0, (1, 3)), color="#8a8a8a")
                ax.annotate(f"resolution floor  $r/R=\\sqrt{{12\\epsilon/R}}$",
                            (t[-1], floor_um), fontsize=7, color="#8a8a8a",
                            ha="right", va="bottom")
            print(f"power-law fit (t+t0)^a over {sel.sum()} pts above the "
                  f"{floor_um:.1f} um floor: a = {fit['a']:.4f} "
                  f"+-{fit['a_ci']:.4f}  (t0 = {fit['t0']:.0f} s, "
                  f"R2 = {fit['r2']:.4f});  1/3 = 0.3333")
        else:
            print("power-law fit: skipped — fewer than 4 samples above the "
                  f"{floor_um:.1f} um resolution floor")
    except Exception as e:                                   # noqa: BLE001
        print(f"power-law fit: skipped ({e})", file=sys.stderr)

    # Experimental overlay (Molaro et al. 2019 Fig. 11 table). Default file
    # is resolved against the repo (this script's parent's parent), so the
    # overlay also appears when analyzing runs outside the repo tree; a
    # postprocess/ snapshot copy inside a run dir won't find it — pass
    # --data explicitly there, or it is skipped with a note.
    data_csv = args.data
    if data_csv is None:
        cand = Path(__file__).resolve().parent.parent / "inputs" / "validation" \
               / "molaro2019_fig11_T-20.csv"
        data_csv = cand if cand.exists() else None
    if data_csv is not None and Path(data_csv).exists():
        dt_, dw_, ep_, em_ = [], [], [], []
        for line in open(data_csv):
            if line.startswith("#") or not line.strip():
                continue
            c = line.split(",")
            dt_.append(float(c[0])); dw_.append(float(c[1]))
            ep_.append(float(c[2])); em_.append(float(c[3]))
        ax.errorbar(dt_, dw_, yerr=[em_, ep_], fmt="o", ms=5, lw=1.2,
                    capsize=3, color="#3f434a", label="Molaro et al. (2019)")
        ax.legend(fontsize=9)
    else:
        print("note: experimental data CSV not found — plotting model only "
              "(pass --data <csv> to overlay)")

    ax.set_xlabel("time [min]")
    ax.set_ylabel("neck width [µm]")
    ax.set_title("Neck width vs time (minimum cross-section of φ ≥ 0.5 body)")
    ax.grid(alpha=0.25, lw=0.5)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    png = Path(str(out).replace(".csv", ".png"))
    fig.savefig(png, dpi=150)
    print(f"plot -> {png}")


if __name__ == "__main__":
    main()

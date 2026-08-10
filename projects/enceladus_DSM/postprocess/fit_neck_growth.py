#!/usr/bin/env python3
"""Fit sintering growth exponents to neck-radius-vs-time series.

WHY THIS EXISTS
---------------
Demmenie, Woutersen & Bonn (J. Phys. Chem. Lett. 2025, 16, 2104-2109,
doi:10.1021/acs.jpclett.5c00050) sintered pairs of ~1 mm ice spheres at -3 C in
a box held at the saturation vapour pressure of ice, and fitted

    r(t)/R0 = [C_m(T) * t]^(1/m)          reported as   r ~ t^alpha

getting alpha = 0.29 / 0.33 / 0.30 / 0.26 (+-0.01) over four runs, i.e. m = 3,
i.e. evaporation-condensation. Their Table 1 maps mechanism -> m as

    viscous flow 1   |   sublimation-condensation 3   |
    bulk diffusion 5 |   surface diffusion 7

(note they argue viscous flow is m = 1, not the textbook 2, because Frenkel's
derivation neither conserves mass nor localises dissipation at the neck).

THE POINT OF THIS SCRIPT: the fitted exponent is not a property of a curve, it
is a property of a curve PLUS A PROTOCOL. On the committed Molaro model runs,
the same neck curve yields a = 0.026 fitted over the solver's own dense output
and a = 0.065 fitted at the nine times Molaro et al. sampled -- a factor 2.5
from the sampling alone. So every number this script prints carries the fit
form and the fit window with it, and all four forms below are reported side by
side rather than one being blessed.

THE FOUR FITS
-------------
d_free   r = C*(t+t0)^a,  t0 free
         Demmenie's own protocol. Their t0 is the unknown delay between first
         contact and the first frame, fitted per run (their error bars come
         from shifting t0 by one 15 s frame). This is the ONLY form that
         compares like-for-like against their alpha = 0.26-0.33.

d_fixed  r = C*t^a,  t0 = 0
         Legitimate only when the clock zero is exact -- i.e. our own runs
         started from tangent contact. Meaningless for a run that opens with a
         mature neck, where it collapses toward a = 0 for the trivial reason
         that r barely changes over the sampled decades.

kucz     r^m - r0^m = K*t,  m and K free, r0 = first in-window radius
         The correct form when r0 > 0, and the only one worth quoting for the
         pre-necked Molaro geometries. m here IS the Kuczynski exponent, so it
         is directly comparable to Demmenie's m (not to their mechanism index
         n, which runs 1-4 -- an easy and consequential off-by-one).

slope    d ln r / d ln t, evaluated locally
         The honest diagnostic. A single power law is a claim that this is
         flat; the plot shows whether it ever is, and over which r/R0.

THE RESOLUTION FLOOR
--------------------
The fillet at the neck has radius rho ~ r^2/(2R). A diffuse interface only
resolves a feature once rho is a few times the 5%-95% band width 6*eps, so a
phase-field neck measurement is only trustworthy above

    r/R  >=  sqrt(12*eps/R)                    (--rmin default)

For the Molaro geometry (eps = 3.48e-7, R_ave = 85 um) that floor is r/R =
0.22, and the ENTIRE Molaro dataset (r/R = 0.19 -> 0.38) sits at or barely
above it -- which is why that geometry cannot resolve an exponent no matter
how the fit is done. For 1 mm spheres at comparable eps the floor drops to
0.087, and that is the actual reason the big-grain experiment is the one worth
simulating. The floor is drawn on both figures; points below it are plotted
faint and excluded from the fits unless --rmin is overridden.

Maeno & Ebinuma (1983) independently put the vapour/surface-diffusion
crossover at r/R0 ~ 0.08, so for 1 mm spheres the numerical floor and the
physical validity limit of a vapour-only model nearly coincide.

INPUT FORMATS (sniffed from the header)
---------------------------------------
model  t_s,neck_width_m,x_neck_m            <- postprocess/neck_width.py
data   time_min,neck_size_um,err_plus_um,err_minus_um,lg_diam_um,sm_diam_um
                                            <- inputs/validation/molaro2019_*.csv

Both are neck WIDTH; this script halves them to a RADIUS, because Demmenie
report a radius. R0 comes from the staged .opts (-ice_grain_R) for model
series and from the grain-diameter columns for data series; override with
--R0. eps likewise from -eps, override with --eps.

Usage
-----
    # the committed Molaro comparison, model arms vs the -20 C data
    python postprocess/fit_neck_growth.py \\
        ~/SimulationResults/enceladus_DSM/scratch/_neck_csv/*.csv \\
        inputs/validation/molaro2019_fig11_T-20.csv \\
        --results-base ~/SimulationResults/enceladus_DSM/scratch \\
        --out studies/sinter_exponent/verification

    # a single run against the Demmenie band
    python postprocess/fit_neck_growth.py <run_dir>/neck_width.csv --demmenie \\
        --out studies/sinter_exponent/verification

    # apples-to-apples: fit the model only where the experiment sampled it
    python postprocess/fit_neck_growth.py <model.csv> \\
        --resample-at inputs/validation/molaro2019_fig11_T-20.csv --out .
"""

from __future__ import annotations

import argparse
import csv
import glob
import math
import os
import re
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit

# Okabe-Ito, matching plot_neck_convergence.py.
COLORS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9"]
STYLES = ["-", "--", ":", "-."]
INK, MUTED, GRID = "#1a1a1a", "#5c5c5c", "#d8d8d8"

# Demmenie et al. (2025), four ice-sphere runs, Figure 3.
DEMMENIE_ALPHA = (0.29, 0.33, 0.30, 0.26)
DEMMENIE_ERR = 0.01

# Kuczynski exponents worth drawing as guide slopes, in Demmenie's Table 1
# convention (their m, NOT their mechanism index n).
GUIDE_M = {1: "viscous flow", 3: "sublim.-cond.", 5: "bulk diff.", 7: "surf. diff."}


# ---------------------------------------------------------------------------
# Series loading
# ---------------------------------------------------------------------------

class Series:
    """One neck-radius-vs-time curve plus the metadata needed to fit it.

    t   [s]  strictly increasing, > 0 (the t = 0 sample is dropped: ln t)
    r   [m]  neck RADIUS (half the measured width)
    R0  [m]  grain radius used for the r/R0 normalisation, or None
    eps [m]  interface parameter, or None (then no resolution floor)
    """

    def __init__(self, label, t, r, R0=None, eps=None, kind="model", path=None):
        keep = (t > 0) & np.isfinite(r) & (r > 0)
        self.t, self.r = t[keep], r[keep]
        self.label, self.R0, self.eps, self.kind, self.path = label, R0, eps, kind, path

    @property
    def rn(self):
        """r/R0, or r in metres when R0 is unknown."""
        return self.r / self.R0 if self.R0 else self.r

    @property
    def floor(self):
        """Resolution floor sqrt(12*eps/R) in r/R0 units, or None."""
        if not (self.eps and self.R0):
            return None
        return math.sqrt(12.0 * self.eps / self.R0)


def _opt_from_text(text, key):
    m = re.search(rf"^{re.escape(key)}\s+(\S+)", text, re.MULTILINE)
    return m.group(1) if m else None


def _opts_near(csv_path, results_base=None):
    """Merged text of the geometry .opts belonging to a neck CSV.

    Two layouts are supported. neck_width.py's default writes the CSV into the
    run folder, where the staged .opts sit alongside it. plot_neck_convergence.py
    instead collects arms into a shared _neck_csv/ directory and names each file
    "<geometry_stem>__<hash>.csv"; there the geometry stem is also the run
    folder's parent directory name under $RESULTS_BASE, so the opts are one
    lookup away.
    """
    d = Path(csv_path).parent
    cand = [f for f in d.glob("*.opts") if f.name != "solver.opts"]
    if not cand and results_base:
        stem = Path(csv_path).stem.split("__")[0]
        cand = [f for f in Path(results_base).glob(f"{stem}/*/*.opts")
                if f.name != "solver.opts"]
    return "\n".join(f.read_text(errors="replace") for f in cand)


def load_series(path, results_base=None, R0=None, eps=None, label=None):
    with open(path, errors="replace") as fh:
        lines = [ln for ln in fh if ln.strip()]
    header = next((ln for ln in lines if not ln.lstrip().startswith("#")), "")

    if header.startswith("t_s"):
        d = np.loadtxt(path, delimiter=",", skiprows=1, ndmin=2)
        t, r = d[:, 0], d[:, 1] / 2.0                      # width -> radius
        kind = "model"
        text = _opts_near(path, results_base)
        if R0 is None:
            radii = _opt_from_text(text, "-ice_grain_R")
            if radii:
                vals = [float(v) for v in radii.split(",")]
                R0 = float(np.mean(vals))
        if eps is None:
            e = _opt_from_text(text, "-eps")
            eps = float(e) if e else None
    else:
        d = np.loadtxt(path, delimiter=",", comments="#", ndmin=2)
        t, r = d[:, 0] * 60.0, d[:, 1] * 1e-6 / 2.0        # min->s, um width->m radius
        kind = "data"
        if R0 is None and d.shape[1] >= 6:
            R0 = float(np.mean(d[:, 4] + d[:, 5])) / 4.0 * 1e-6   # mean of two radii

    return Series(label or Path(path).stem, t, r, R0, eps, kind, path)


# ---------------------------------------------------------------------------
# Fits
# ---------------------------------------------------------------------------

def _ci(pcov, i):
    """95% half-width for parameter i, or nan if the covariance is degenerate."""
    try:
        v = pcov[i][i]
        return 1.96 * math.sqrt(v) if np.isfinite(v) and v >= 0 else float("nan")
    except (IndexError, TypeError):
        return float("nan")


def _r2(y, yhat):
    ss = np.sum((y - np.mean(y)) ** 2)
    return 1.0 - np.sum((y - yhat) ** 2) / ss if ss > 0 else float("nan")


def fit_d_fixed(t, r):
    """r = C*t^a. Ordinary least squares on (ln t, ln r) -- closed form."""
    if len(t) < 3:
        return None
    x, y = np.log(t), np.log(r)
    n = len(x)
    a, lnC = np.polyfit(x, y, 1)
    resid = y - (a * x + lnC)
    # OLS standard error of the slope.
    sxx = np.sum((x - x.mean()) ** 2)
    se = math.sqrt(np.sum(resid ** 2) / (n - 2) / sxx) if n > 2 and sxx > 0 else float("nan")
    return {"form": "d_fixed", "a": a, "a_ci": 1.96 * se, "C": math.exp(lnC),
            "t0": 0.0, "m": 1.0 / a if a else float("nan"),
            "r2": _r2(y, a * x + lnC)}


def fit_d_free(t, r):
    """r = C*(t+t0)^a with t0 free, fitted on ln r -- Demmenie's protocol.

    Fitted in log space because that is where they read the slope off Figure 3;
    a least-squares fit in the r domain would let the last decade dominate.
    t0 is bounded below by -min(t) (keeps t+t0 > 0) and initialised from the
    d_fixed residual curvature via a coarse scan, since a bad t0 seed sends
    curve_fit into the flat region where a -> 0.
    """
    if len(t) < 4:
        return None
    y = np.log(r)

    def model(tt, lnC, a, t0):
        return lnC + a * np.log(np.maximum(tt + t0, 1e-300))

    best, seed = None, None
    for t0 in np.concatenate(([0.0], np.geomspace(1e-3 * t.max(), 10 * t.max(), 40))):
        try:
            a, lnC = np.polyfit(np.log(t + t0), y, 1)
        except (ValueError, np.linalg.LinAlgError):
            continue
        sse = np.sum((y - (lnC + a * np.log(t + t0))) ** 2)
        if best is None or sse < best:
            best, seed = sse, (lnC, a, t0)
    if seed is None:
        return None

    try:
        p, pcov = curve_fit(model, t, y, p0=seed, maxfev=200000,
                            bounds=([-np.inf, -np.inf, -0.999 * t.min()],
                                    [np.inf, np.inf, np.inf]))
    except (RuntimeError, ValueError):
        return None
    lnC, a, t0 = p
    return {"form": "d_free", "a": a, "a_ci": _ci(pcov, 1), "C": math.exp(lnC),
            "t0": t0, "m": 1.0 / a if a else float("nan"),
            "r2": _r2(y, model(t, *p))}


def fit_kuczynski(t, r):
    """r^m - r0^m = K*t, with r0 pinned to the first in-window radius.

    Parameterised as (m, ln K) and fitted on u = r/r0 rather than on r itself.
    Both rescalings matter. u >= 1 keeps u**m finite for the large m a nearly
    flat curve drives the fit toward (raw r**m in metres underflows to zero
    around m ~ 13 and the fit silently stops being a fit), and a unit-scale
    response keeps the residuals above curve_fit's absolute ftol -- fitted in
    metres the cost is O(1e-14) at the seed, curve_fit declares convergence on
    the first iteration, and every series comes back with m exactly equal to
    whichever seed won. That failure is silent and looks like a real answer.
    """
    if len(t) < 4:
        return None
    r0 = r[0]
    u = r / r0                                  # O(1) response, u >= 1

    def model(tt, m, lnK):
        return np.power(np.maximum(1.0 + math.exp(lnK) * tt, 1e-300), 1.0 / m)

    best = None
    for m0 in (2.0, 3.0, 5.0, 7.0, 11.0):
        # K such that the model passes through the last point at this m.
        K0 = (u[-1] ** m0 - 1.0) / t[-1] if t[-1] > 0 else 1.0
        if K0 <= 0:
            continue
        try:
            p, pcov = curve_fit(model, t, u, p0=[m0, math.log(K0)], maxfev=200000,
                                bounds=([0.2, -np.inf], [60.0, np.inf]))
        except (RuntimeError, ValueError):
            continue
        sse = np.sum((u - model(t, *p)) ** 2)
        if best is None or sse < best[0]:
            best = (sse, p, pcov)
    if best is None:
        return None
    _, p, pcov = best
    m, lnK = p
    return {"form": "kucz", "a": 1.0 / m, "a_ci": _ci(pcov, 0) / m ** 2,
            "C": math.exp(lnK), "t0": 0.0, "m": m,
            "r2": _r2(u, model(t, *p)), "m_ci": _ci(pcov, 0)}


def local_slope(t, r, window=9):
    """d ln r / d ln t from a sliding least-squares line in log-log space.

    A bare np.gradient on log-log is dominated by the measurement jitter in r
    (the neck is a sub-grid-refined minimum, so consecutive samples differ by
    much less than their own noise). A short regression window is the cheapest
    estimator that is not noise; window is in samples, forced odd.
    """
    x, y = np.log(t), np.log(r)
    n = len(x)
    w = max(3, window | 1)
    out = np.full(n, np.nan)
    for i in range(n):
        lo, hi = max(0, i - w // 2), min(n, i + w // 2 + 1)
        if hi - lo >= 3:
            out[i] = np.polyfit(x[lo:hi], y[lo:hi], 1)[0]
    return out


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def window_mask(s, rmin, rmax, tmin, tmax):
    m = np.ones(len(s.t), bool)
    rn = s.rn
    lo = rmin if rmin is not None else (s.floor or 0.0)
    if lo:
        m &= rn >= lo
    if rmax is not None:
        m &= rn <= rmax
    if tmin is not None:
        m &= s.t >= tmin
    if tmax is not None:
        m &= s.t <= tmax
    return m, lo


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv", nargs="+", type=Path,
                    help="neck CSVs (model and/or experimental; format sniffed)")
    ap.add_argument("--out", type=Path, required=True, help="output directory")
    ap.add_argument("--results-base", type=Path, default=None,
                    help="$RESULTS_BASE, so a shared _neck_csv/ file can find "
                         "its geometry .opts (for R0 and eps)")
    ap.add_argument("--labels", default=None,
                    help="comma-separated series labels, in argument order")
    ap.add_argument("--R0", type=float, default=None,
                    help="grain radius [m] override, applied to every series")
    ap.add_argument("--eps", type=float, default=None,
                    help="interface parameter [m] override, applied to every series")
    ap.add_argument("--rmin", type=float, default=None,
                    help="fit-window floor in r/R0 (default: sqrt(12*eps/R) per series)")
    ap.add_argument("--rmax", type=float, default=None, help="fit-window ceiling in r/R0")
    ap.add_argument("--tmin", type=float, default=None, help="fit-window start [s]")
    ap.add_argument("--tmax", type=float, default=None, help="fit-window end [s]")
    ap.add_argument("--resample-at", type=Path, default=None,
                    help="resample every MODEL series onto this series' sample "
                         "times before fitting (log-interpolated). Exposes how "
                         "much of a fitted exponent is the sampling.")
    ap.add_argument("--slope-window", type=int, default=9,
                    help="samples per local-slope regression window")
    ap.add_argument("--demmenie", action="store_true",
                    help="draw the Demmenie et al. (2025) alpha = 0.26-0.33 band")
    ap.add_argument("--prefix", default="neck", help="output filename prefix")
    args = ap.parse_args()

    paths = []
    for p in args.csv:
        hits = sorted(glob.glob(str(p)))
        paths.extend(hits if hits else [str(p)])
    labels = args.labels.split(",") if args.labels else [None] * len(paths)
    if len(labels) != len(paths):
        sys.exit(f"--labels has {len(labels)} entries for {len(paths)} series")

    series = [load_series(p, args.results_base, args.R0, args.eps, lb)
              for p, lb in zip(paths, labels)]
    series = [s for s in series if len(s.t) >= 3]
    if not series:
        sys.exit("no usable series (need >= 3 samples with t > 0)")

    if args.resample_at:
        ref = load_series(args.resample_at, args.results_base)
        for s in series:
            if s.kind != "model":
                continue
            tt = ref.t[(ref.t >= s.t[0]) & (ref.t <= s.t[-1])]
            if len(tt) >= 3:
                s.t, s.r = tt, np.exp(np.interp(np.log(tt), np.log(s.t), np.log(s.r)))
                s.label += " @data-t"

    args.out.mkdir(parents=True, exist_ok=True)
    rows, fits_by_series = [], {}

    for s in series:
        mask, lo = window_mask(s, args.rmin, args.rmax, args.tmin, args.tmax)
        t, r = s.t[mask], s.r[mask]
        fits = [f for f in (fit_d_free(t, r), fit_d_fixed(t, r), fit_kuczynski(t, r))
                if f is not None]
        fits_by_series[s.label] = (mask, fits)
        if not fits:
            # An empty window is a RESULT, not a skip: it means the whole
            # series lies below sqrt(12*eps/R), i.e. this geometry cannot
            # resolve a neck exponent at this eps. Say so instead of quietly
            # dropping the series from the table.
            rn = s.rn
            print(f"NOTE: '{s.label}' has no fittable window — "
                  f"{mask.sum()}/{len(s.t)} samples pass "
                  f"r/R0 >= {lo:.3f} (series spans {rn.min():.3f}-{rn.max():.3f}). "
                  f"Plotted faint, excluded from the fit table.", file=sys.stderr)
        for f in fits:
            rows.append({
                "series": s.label, "kind": s.kind, "form": f["form"],
                "a": f["a"], "a_ci95": f["a_ci"], "m": f["m"],
                "C": f["C"], "t0_s": f["t0"], "r2": f["r2"],
                "n_pts": len(t), "t_lo_s": t[0] if len(t) else float("nan"),
                "t_hi_s": t[-1] if len(t) else float("nan"),
                "rn_lo": (r[0] / s.R0) if (len(r) and s.R0) else float("nan"),
                "rn_hi": (r[-1] / s.R0) if (len(r) and s.R0) else float("nan"),
                "R0_m": s.R0 or float("nan"), "eps_m": s.eps or float("nan"),
                "rn_floor": s.floor if s.floor else float("nan"),
                "rmin_used": lo,
            })

    csv_out = args.out / f"{args.prefix}_fits.csv"
    with open(csv_out, "w", newline="") as fh:
        # lineterminator: csv defaults to CRLF, which git flags on every commit.
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    # ---- console table ----------------------------------------------------
    print(f"\n{'series':<38} {'form':<8} {'a':>8} {'+-95%':>8} {'m':>7} "
          f"{'t0 [s]':>10} {'R2':>7} {'window r/R0':>16} {'n':>4}")
    print("-" * 118)
    for r_ in rows:
        win = (f"{r_['rn_lo']:.3f}-{r_['rn_hi']:.3f}"
               if np.isfinite(r_["rn_lo"]) else "n/a")
        print(f"{r_['series'][:38]:<38} {r_['form']:<8} {r_['a']:>8.4f} "
              f"{r_['a_ci95']:>8.4f} {r_['m']:>7.2f} {r_['t0_s']:>10.1f} "
              f"{r_['r2']:>7.4f} {win:>16} {r_['n_pts']:>4}")
    print(f"\nDemmenie et al. (2025) ice spheres: alpha = "
          f"{', '.join(f'{a:.2f}' for a in DEMMENIE_ALPHA)} (+-{DEMMENIE_ERR})"
          f"  -> m = 3, sublimation-condensation")
    print(f"csv -> {csv_out}")

    # ---- figures ----------------------------------------------------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7.6, 5.4))
    for i, s in enumerate(series):
        c, st = COLORS[i % len(COLORS)], STYLES[i % len(STYLES)]
        mask, fits = fits_by_series[s.label]
        # Out-of-window samples stay visible but faint: the reader should see
        # what was excluded, not just what survived.
        ax.plot(s.t[~mask], s.rn[~mask], st, color=c, lw=1.0, alpha=0.25, zorder=2)
        if s.kind == "data":
            ax.plot(s.t[mask], s.rn[mask], "o", color=c, ms=5, zorder=5, label=s.label)
        else:
            ax.plot(s.t[mask], s.rn[mask], st, color=c, lw=2.0, zorder=4, label=s.label)
        if s.floor:
            ax.axhline(s.floor, color=c, lw=0.8, ls=(0, (1, 3)), alpha=0.6, zorder=1)

    # Framing and guide slopes are both driven by the IN-WINDOW cloud, not by
    # everything plotted. A phase-field run's first accepted steps sit at
    # t ~ 1e-4 s, seven decades before anything physical happens; letting them
    # set the axes squashes the entire fitted region into one tick and fans
    # the guide slopes across the full height of the figure.
    inw = [(s.t[m_], s.rn[m_]) for s, (m_, _) in
           ((s, fits_by_series[s.label]) for s in series) if m_.any()]
    if not inw:
        inw = [(s.t, s.rn) for s in series]
    tw = np.concatenate([a for a, _ in inw])
    rw = np.concatenate([b for _, b in inw])

    ta, ra = math.sqrt(tw.min() * tw.max()), math.sqrt(rw.min() * rw.max())
    xlo, xhi = tw.min() / 30.0, tw.max() * 3.0
    tg = np.geomspace(xlo, xhi, 60)
    for m, name in GUIDE_M.items():
        ax.plot(tg, ra * (tg / ta) ** (1.0 / m), color=MUTED, lw=0.7,
                ls=(0, (4, 3)), alpha=0.55, zorder=0)
        ax.annotate(f"1/{m}  {name}", (tg[-1], ra * (tg[-1] / ta) ** (1.0 / m)),
                    fontsize=7, color=MUTED, ha="right", va="bottom")

    if args.demmenie:
        lo, hi = min(DEMMENIE_ALPHA) - DEMMENIE_ERR, max(DEMMENIE_ALPHA) + DEMMENIE_ERR
        ax.fill_between(tg, ra * (tg / ta) ** lo, ra * (tg / ta) ** hi,
                        color="#0072B2", alpha=0.10, zorder=0,
                        label=f"Demmenie 2025  $\\alpha$={lo:.2f}-{hi:.2f}")

    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(xlo, xhi)
    # Vertical extent from whatever is actually visible in that x range,
    # so the excluded-but-plotted samples do not blow the axis open either.
    vis = np.concatenate([s.rn[(s.t >= xlo) & (s.t <= xhi)] for s in series]
                         + [rw])
    ax.set_ylim(vis.min() / 1.8, vis.max() * 1.8)
    ax.set_xlabel("time [s]")
    ax.set_ylabel("$r/R_0$" if all(s.R0 for s in series) else "neck radius [m]")
    ax.set_title("Neck growth: $r = C\\,t^{a}$  (dotted = $\\sqrt{12\\epsilon/R}$ "
                 "resolution floor)", fontsize=10)
    ax.grid(alpha=0.25, lw=0.5, which="both", color=GRID)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    p1 = args.out / f"{args.prefix}_growth_loglog.png"
    fig.savefig(p1, dpi=150)
    print(f"plot -> {p1}")

    fig, ax = plt.subplots(figsize=(7.6, 4.6))
    slopes = []
    for i, s in enumerate(series):
        c, st = COLORS[i % len(COLORS)], STYLES[i % len(STYLES)]
        sl = local_slope(s.t, s.r, args.slope_window)
        slopes.append(sl)
        ax.plot(s.rn, sl, st, color=c, lw=1.8, label=s.label)
        if s.floor:
            ax.axvline(s.floor, color=c, lw=0.8, ls=(0, (1, 3)), alpha=0.6)
    if args.demmenie:
        ax.axhspan(min(DEMMENIE_ALPHA) - DEMMENIE_ERR, max(DEMMENIE_ALPHA) + DEMMENIE_ERR,
                   color="#0072B2", alpha=0.10, label="Demmenie 2025")
    ax.set_xscale("log")
    ax.set_xlabel("$r/R_0$" if all(s.R0 for s in series) else "neck radius [m]")
    ax.set_ylabel("$d\\ln r\\,/\\,d\\ln t$")
    hi = np.nanmax(np.concatenate(slopes)) if slopes else 0.5
    top = min(1.2, max(0.4, (hi if np.isfinite(hi) else 0.5) * 1.2))
    ax.set_ylim(0, top)
    # Guide lines are drawn after the limits are set, and only the ones that
    # land inside them: an annotation parked off-axis defeats tight_layout.
    for m, name in GUIDE_M.items():
        if 1.0 / m > top:
            continue
        ax.axhline(1.0 / m, color=MUTED, lw=0.7, ls=(0, (4, 3)), alpha=0.55)
        ax.annotate(f"1/{m}  {name}", (0.995, 1.0 / m), xycoords=("axes fraction", "data"),
                    fontsize=7, color=MUTED, ha="right", va="bottom")
    ax.set_title("Local growth exponent (dotted = resolution floor)", fontsize=10)
    ax.grid(alpha=0.25, lw=0.5, which="both", color=GRID)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    p2 = args.out / f"{args.prefix}_local_slope.png"
    fig.savefig(p2, dpi=150)
    print(f"plot -> {p2}\n")


if __name__ == "__main__":
    main()

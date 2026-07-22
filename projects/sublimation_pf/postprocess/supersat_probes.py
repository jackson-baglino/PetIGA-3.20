#!/usr/bin/env python3
"""Supersaturation probe time series for two-grain sintering runs.

Tracks sigma = (rho_v - rho_vs(T))/rho_vs(T) at six moving probes:

  neck        just outside the phi contour at the waist (the deposition sink)
  top-sm/lg   above each grain center at r = R(t) + delta ("top of the grains,
              in line with the centers" — halfway position along the surface
              between neck and far pole)
  pole-sm/lg  on the axis just beyond each grain's far tip (furthest points
              from the neck)
  farfield    fixed point far from all ice (0.9*Ly off-axis, mid-z)

All ice-adjacent probes float with the geometry (recomputed per snapshot from
the phi = 0.5 contour, offset delta = 3*eps into the vapor), so the series
stay physical as the neck migrates and the grains shrink.

REGIME DISCRIMINANT (printed and annotated on the plot)
-------------------------------------------------------
Two-sphere sintering by the vapor route has two limits (Kingery & Berg):

  kinetic-limited   x^3 - x0^3 = A3 * (d0*R/beta_sub) * t          (n = 3)
  diffusion-limited x^5 - x0^5 = A5 * D_v*(rho_vs/rho_i)*d0*R^2 * t (n = 5)

with A3, A5 O(1) geometry constants. The crossover is set by the kinetic
screening length L_k = D_v * beta_HK (the distance over which vapor diffusion
resistance equals attachment resistance):

  Lambda(t) = L_k / x(t)     Lambda >> 1  -> attachment-limited (n -> 3)
                             Lambda << 1  -> vapor-diffusion-limited (n -> 5)

The probes give the same answer non-asymptotically: the transport share of
the total Gibbs-Thomson driving is

  f_T(t) = (sigma_pole - sigma_neck) / (d0 * (kappa_grain - kappa_neck))

f_T -> 0 kinetic-limited (vapor field flat, all driving dropped at the
interface); f_T -> 1 diffusion-limited (interfaces at local equilibrium, all
driving dropped across the vapor field).

Usage:
    python supersat_probes.py <run_dir> [--axisym] [--eps 4.64e-7]
                              [--beta-hk B] [--d0 D0]
"""

import argparse
import base64
import csv
import glob
import math
import re
import struct
import sys
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np

# Thermodynamics inlined (mirrors material_properties.c / comp_eps.py) so
# this script also runs from the postprocess/ snapshot copy inside a run
# directory, where ../preprocess does not exist.
_KJ = [-0.5865e4, 0.2224e2, 0.1375e-1, -0.3403e-4, 0.2697e-7, 0.6918]
_RHO_AIR, _BB, _PATM = 1.341, 0.62, 101325.0
_DV0 = 2.178e-5


def rho_vs_sat(T_C):
    """Saturation vapor density over ice [kg/m^3] (ASHRAE, as in RhoVS_I)."""
    T = T_C + 273.15
    Pvs = math.exp(_KJ[0] / T + _KJ[1] + _KJ[2] * T + _KJ[3] * T ** 2
                   + _KJ[4] * T ** 3 + _KJ[5] * math.log(T))
    return _RHO_AIR * _BB * Pvs / (_PATM - Pvs)


def Dv_T(T_C):
    """Vapor diffusivity D_v(T) [m^2/s]."""
    return _DV0 * ((T_C + 273.15) / 273.15) ** 1.81


def read_vts(fn):
    root = ET.parse(fn).getroot()
    grid = root.find(".//StructuredGrid")
    ext = [int(v) for v in grid.get("WholeExtent").split()]
    nx, ny = ext[1] - ext[0] + 1, ext[3] - ext[2] + 1

    def dec(da):
        raw = base64.b64decode("".join(da.text.split()))
        n = struct.unpack("<Q", raw[:8])[0]
        return np.frombuffer(raw[8:8 + n], dtype=np.float64)

    pts = dec(root.find(".//Points/DataArray")).reshape(ny, nx, 3)
    out = {}
    for da in root.findall(".//PointData/DataArray"):
        out[da.get("Name")] = dec(da).reshape(ny, nx)
    return out, pts[0, :, 0], pts[:, 0, 1]


def step_times(outp):
    tmap = {}
    pat = re.compile(r"^\s+(\d+)\s+\|\s+([0-9.eE+-]+)\s+\|")
    for line in open(outp, errors="replace"):
        m = pat.match(line)
        if m and line.count("|") == 8:
            tmap[int(m.group(1))] = float(m.group(2))
    return tmap


def contour_r(col, y, level=0.5):
    """Outermost phi=level crossing height along a column (0 if no ice)."""
    a = np.flatnonzero(col >= level)
    if a.size == 0:
        return 0.0
    hi = a[-1]
    if hi < len(y) - 1:
        f = (col[hi] - level) / (col[hi] - col[hi + 1])
        return y[hi] + f * (y[hi + 1] - y[hi])
    return y[hi]


def sample(field_dict, x, y, xp, yp):
    """Bilinear sample of sigma at physical point (xp, yp)."""
    T = field_dict["Temperature"]; rv = field_dict["VaporDensity"]
    j = np.clip(np.searchsorted(x, xp) - 1, 0, len(x) - 2)
    i = np.clip(np.searchsorted(y, yp) - 1, 0, len(y) - 2)
    fx = (xp - x[j]) / (x[j + 1] - x[j]); fy = (yp - y[i]) / (y[i + 1] - y[i])
    def bl(F):
        return ((1 - fx) * (1 - fy) * F[i, j] + fx * (1 - fy) * F[i, j + 1]
                + (1 - fx) * fy * F[i + 1, j] + fx * fy * F[i + 1, j + 1])
    Tv, rvv = bl(T), bl(rv)
    rvs = rho_vs_sat(Tv)
    return (rvv - rvs) / rvs


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir", type=Path)
    ap.add_argument("--axisym", action="store_true",
                    help="(accepted for symmetry with neck_width.py; probe "
                         "logic is identical — the axis is y=0 either way "
                         "for on-axis grains)")
    ap.add_argument("--eps", type=float, default=4.64e-7, help="interface eps [m]")
    ap.add_argument("--beta-hk", type=float, default=None,
                    help="beta_HK [s/m] for Lambda; parsed from outp.txt if omitted")
    ap.add_argument("--d0", type=float, default=1.0166e-9, help="capillary length [m]")
    args = ap.parse_args()

    files = sorted(glob.glob(str(args.run_dir / "vtkOut" / "solV_*.vts")))
    if not files:
        sys.exit(f"no solV_*.vts under {args.run_dir}/vtkOut")
    tmap = step_times(args.run_dir / "outp.txt")

    beta_hk = args.beta_hk
    if beta_hk is None:
        for line in open(args.run_dir / "outp.txt", errors="replace"):
            m = re.search(r"beta_sub \(SCALED.*?=\s*([0-9.eE+-]+)", line)
            if m:
                beta_hk = float(m.group(1)); break
        # SCALED beta printed as "beta_sub (SCALED = ...) = X s/m"; fall back:
        if beta_hk is None:
            beta_hk = 0.733
    Dv = Dv_T(-20.0)
    delta = 3.0 * args.eps

    names = ["neck", "top-sm", "top-lg", "pole-sm", "pole-lg", "farfield"]
    series = {k: [] for k in names}; times = []; lam = []

    for fn in files:
        step = int(Path(fn).stem.split("_")[1])
        try:
            d, x, y = read_vts(fn)
        except Exception:
            continue
        phi = d["IcePhase"]
        rprof = np.array([contour_r(phi[:, j], y) for j in range(len(x))])

        occ = np.flatnonzero(rprof > 0)
        thr = 0.5 * rprof.max()
        peaks = [j for j in range(5, len(rprof) - 5)
                 if rprof[j] >= thr and rprof[j] == rprof[j - 5:j + 6].max()]
        if len(peaks) < 2:
            continue
        c_sm, c_lg = peaks[0], peaks[-1]
        if rprof[c_sm] > rprof[c_lg]:
            c_sm, c_lg = c_lg, c_sm
        interior = np.arange(min(c_sm, c_lg) + 1, max(c_sm, c_lg))
        interior = interior[rprof[interior] > 0]
        jn = interior[np.argmin(rprof[interior])]

        x_neck = rprof[jn]                       # neck radius
        probes = {
            "neck":     (x[jn], rprof[jn] + delta),
            "top-sm":   (x[c_sm], rprof[c_sm] + delta),
            "top-lg":   (x[c_lg], rprof[c_lg] + delta),
            "pole-sm":  (x[min(occ[0], occ[-1])] - delta if occ[0] < jn else x[occ[0]] - delta, y[1]),
            "pole-lg":  (x[occ[-1]] + delta, y[1]),
            "farfield": (0.5 * x[-1], 0.9 * y[-1]),
        }
        # pole-sm: leftmost ice extent minus delta (small grain is left by
        # construction in our geometries; peaks ordering already handles size)
        probes["pole-sm"] = (x[occ[0]] - delta, y[1])

        for k, (px, py) in probes.items():
            series[k].append(sample(d, x, y, px, py))
        times.append(tmap.get(step, np.nan))
        lam.append(Dv * beta_hk / x_neck)

    t = np.array(times) / 60.0
    out = args.run_dir / "supersat_probes.csv"
    with open(out, "w") as f:
        wcsv = csv.writer(f)
        wcsv.writerow(["t_min"] + names + ["Lambda"])
        for i in range(len(t)):
            wcsv.writerow([f"{t[i]:.4f}"] + [f"{series[k][i]:.6e}" for k in names]
                          + [f"{lam[i]:.4f}"])
    print(f"csv -> {out}")

    # transport share f_T from the probes (avg of both grains' poles)
    sig_pole = 0.5 * (np.array(series["pole-sm"]) + np.array(series["pole-lg"]))
    sig_neck = np.array(series["neck"])
    print(f"Lambda = D_v*beta_HK/x_neck: {lam[0]:.2f} (t0) -> {lam[-1]:.2f} (end)"
          f"   [>>1 kinetic-limited, <<1 diffusion-limited]")
    print(f"transport driving split (sigma_pole - sigma_neck): "
          f"{(sig_pole-sig_neck)[1]:.2e} (early) -> {(sig_pole-sig_neck)[-1]:.2e} (end)")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    colors = {"neck": "#3d74d9", "top-sm": "#e8883a", "top-lg": "#b05cc7",
              "pole-sm": "#4aa66c", "pole-lg": "#2b9aa8", "farfield": "#8a8f98"}
    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(8, 7), sharex=True,
                                  height_ratios=[2.2, 1])
    for k in names:
        ax.plot(t, np.array(series[k]) * 1e4, lw=2, color=colors[k], label=k)
    ax.axhline(0, color="k", lw=0.6)
    ax.set_ylabel(r"supersaturation $\sigma \times 10^4$")
    ax.legend(fontsize=8, ncol=3)
    ax.grid(alpha=0.25, lw=0.5)
    ax.set_title("Supersaturation probes (moving with the interface)")
    ax2.plot(t, lam, lw=2, color="#3d74d9")
    ax2.axhline(1.0, color="k", lw=0.8, ls="--")
    ax2.set_yscale("log")
    ax2.set_ylabel(r"$\Lambda = D_v \beta_{HK}/x_{neck}$")
    ax2.set_xlabel("time [min]")
    ax2.grid(alpha=0.25, lw=0.5)
    ax2.text(0.02, 0.8, "kinetic-limited above dashes; diffusion-limited below",
             transform=ax2.transAxes, fontsize=8)
    for a in (ax, ax2):
        a.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    png = args.run_dir / "supersat_probes.png"
    fig.savefig(png, dpi=150)
    print(f"plot -> {png}")


if __name__ == "__main__":
    main()

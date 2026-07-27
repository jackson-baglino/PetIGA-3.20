#!/usr/bin/env python3
"""mass_from_vts.py — mass-conservation plot integrated straight from the
vtkOut/solV_*.vts snapshots, versus STEP number.

Why this exists: plot_mass.py takes its TIME axis from SSA_evo.dat / outp.txt.
Those are single append-files that a mid-run download often leaves stale (the
sol_/solV_ snapshots have fresh filenames each output and sync current, while
SSA_evo.dat keeps the same name and lags). When that happens the mass plot
stops at an early step even though the run is far ahead. This script ignores
the stale time-series files and integrates the conserved quantities directly
from the snapshots that ARE current, so the curve reaches the latest step.

It matches monitoring.c's integrands (assembly.c Stats S[]):
    tot_ice  = ∫ φ dA
    tot_air  = ∫ (1-φ) dA
    tot_rhov = ∫ ρv·(1-φ) dA
    tot_mass = ρ_ice·tot_ice + tot_rhov          (the conserved quantity)
integrated on the true (deformed, bumpy) cells via the shoelace area of each
quad. Planar 2D only (no r-weight).

CAVEAT: x-axis is STEP, not physical time — the step→time map lives only in
the (stale) SSA_evo.dat/outp.txt. Re-sync those and use plot_mass.py for a
true time axis; use this for a quick current look mid-run.

Usage:
    python mass_from_vts.py <run_dir> [--out FILE.png] [--rho-ice 919]
"""

import argparse
import base64
import glob
import re
import struct
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_vts(fn):
    root = ET.parse(fn).getroot()
    grid = root.find(".//StructuredGrid")
    ext = [int(v) for v in grid.get("WholeExtent").split()]
    nx, ny = ext[1] - ext[0] + 1, ext[3] - ext[2] + 1

    def decode(da):
        raw = base64.b64decode("".join(da.text.split()))
        n = struct.unpack("<Q", raw[:8])[0]
        return np.frombuffer(raw[8:8 + n], dtype=np.float64)

    pts = None
    for da in root.findall(".//Points/DataArray"):
        pts = decode(da).reshape(ny, nx, 3)
    f = {}
    for da in root.findall(".//PointData/DataArray"):
        f[da.get("Name")] = decode(da).reshape(ny, nx)
    return f, pts[:, :, 0], pts[:, :, 1]


def step_of(fn):
    return int(re.search(r"solV_(\d+)\.vts", fn).group(1))


def integrate(fn, rho_ice):
    """Return (tot_ice, tot_air, tot_rhov, tot_mass) over the deformed grid."""
    f, X, Y = read_vts(fn)
    phi, rhov = f["IcePhase"], f["VaporDensity"]
    # quad corner arrays (ny-1, nx-1)
    x00, x10 = X[:-1, :-1], X[:-1, 1:]
    x11, x01 = X[1:, 1:], X[1:, :-1]
    y00, y10 = Y[:-1, :-1], Y[:-1, 1:]
    y11, y01 = Y[1:, 1:], Y[1:, :-1]
    # shoelace area of each quad (00->10->11->01)
    area = 0.5 * np.abs(
        x00 * y10 - x10 * y00 + x10 * y11 - x11 * y10
        + x11 * y01 - x01 * y11 + x01 * y00 - x00 * y01)
    pc = 0.25 * (phi[:-1, :-1] + phi[:-1, 1:] + phi[1:, 1:] + phi[1:, :-1])
    rc = 0.25 * (rhov[:-1, :-1] + rhov[:-1, 1:] + rhov[1:, 1:] + rhov[1:, :-1])
    tot_ice = float(np.sum(area * pc))
    tot_air = float(np.sum(area * (1.0 - pc)))
    tot_rhov = float(np.sum(area * rc * (1.0 - pc)))
    return tot_ice, tot_air, tot_rhov, rho_ice * tot_ice + tot_rhov


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir", type=Path)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--rho-ice", type=float, default=919.0)
    args = ap.parse_args()

    files = sorted(glob.glob(str(args.run_dir / "vtkOut" / "solV_*.vts")), key=step_of)
    if not files:
        sys.exit(f"no solV_*.vts under {args.run_dir}/vtkOut")

    steps, ice, air, rhov, mass = [], [], [], [], []
    for k, fn in enumerate(files):
        try:
            ti, ta, tr, tm = integrate(fn, args.rho_ice)
        except Exception as e:
            print(f"  skip {Path(fn).name}: {e}", file=sys.stderr); continue
        steps.append(step_of(fn)); ice.append(ti); air.append(ta)
        rhov.append(tr); mass.append(tm)
        if k % 50 == 0:
            print(f"  {k}/{len(files)}  step {step_of(fn)}", flush=True)
    steps = np.array(steps)
    ice, air, rhov, mass = map(np.array, (ice, air, rhov, mass))

    def pct(a): return (a - a[0]) / a[0] * 100.0 if a[0] else np.zeros_like(a)

    fig, ax = plt.subplots(2, 2, figsize=(13, 8))
    fig.suptitle(f"Mass conservation vs STEP (integrated from vtkOut snapshots)\n"
                 f"{args.run_dir.name}   —   latest step {steps[-1]}", fontsize=11)
    ax[0, 0].plot(steps, pct(mass), color="#3a8f8a")
    ax[0, 0].set_title("total mass drift  ρ_ice·∫φ + ∫ρv(1-φ)")
    ax[0, 0].set_ylabel("% change from step 0"); ax[0, 0].axhline(0, color="k", lw=0.5)
    ax[0, 1].plot(steps, ice * 1e6, color="#3d74d9")
    ax[0, 1].set_title(r"total ice  $\int\phi\,dA$  [×10$^{-6}$ m²]")
    ax[1, 0].plot(steps, rhov, color="#e8883a")
    ax[1, 0].set_title(r"total vapor  $\int\rho_v(1-\phi)\,dA$")
    ax[1, 1].plot(steps, pct(ice), label="ice", color="#3d74d9")
    ax[1, 1].plot(steps, pct(air), label="air", color="#8a8f3a")
    ax[1, 1].plot(steps, pct(rhov), label="vapor", color="#e8883a")
    ax[1, 1].set_title("% change from step 0"); ax[1, 1].legend(); ax[1, 1].axhline(0, color="k", lw=0.5)
    for a in ax.flat:
        a.set_xlabel("step"); a.grid(alpha=0.3)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = args.out or (args.run_dir / "mass_from_vts.png")
    fig.savefig(out, dpi=140)
    print(f"\nmass drift over run: {pct(mass)[-1]:+.4f}%   (step 0 -> {steps[-1]})")
    print(f"plot -> {out}")


if __name__ == "__main__":
    main()

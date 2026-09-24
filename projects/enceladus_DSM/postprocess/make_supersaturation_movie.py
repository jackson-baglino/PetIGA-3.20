#!/usr/bin/env python3
"""Movie of the vapour supersaturation in the pore space.

    venv_enceladus/bin/python postprocess/make_supersaturation_movie.py --dir <run>

WHAT IS PLOTTED

    sigma(x,t) = ( rho_v(x,t) - rho_vs(T(x,t)) ) / rho_vs(T(x,t))

the local supersaturation, dimensionless: positive where vapour is supersaturated
(net deposition), negative where it is undersaturated (net sublimation). Shown
in the AIR PHASE only -- ice is masked to a neutral grey -- because rho_v inside
a grain is a phase-field continuation, not a physical pore vapour density.

rho_vs mirrors `rho_vs_sat` in preprocess/comp_eps.py and `RhoVS_I` in
src/material_properties.c; the constants are imported from comp_eps rather than
restated, so the three cannot drift apart.

THE COLOUR SCALE

cmocean `balance`, symmetric about zero. A diverging map whose midpoint is not
the field's zero is a lie about the sign, so the scale is always +/- one bound.

Which bound is the question. The field is strongly asymmetric -- on the
2026-09-16 pilot, over all 55 frames, sigma spans -2.76e-4 to +4.72e-5, i.e.
the undersaturated extreme is 5.9x the oversaturated one. So:

  --bound-from smaller  (default)  +/- min(|min|,|max|). The less extreme side
                                   uses the full half-range and the dominant
                                   side clips. This is "take the smallest
                                   extreme, not the largest".
  --bound-from under               +/- |min sigma|
  --bound-from over                +/- |max sigma|
  --bound-from larger              +/- max(|min|,|max|)
  --bound VALUE                    set it explicitly

`--percentile P` (default 100 = the true extreme) computes the bound from the
P-th percentile instead, which is worth trying if a handful of cells at a neck
are setting the scale for the whole field.

Bounds are computed ONCE over every frame, so the colour of a given sigma does
not change as the movie plays. A per-frame rescale makes a movie in which
nothing appears to happen.
"""
from __future__ import annotations

import argparse
import glob
import os
import re
import shutil
import subprocess
import sys
import tempfile

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cmocean
import vtk
from vtk.util.numpy_support import vtk_to_numpy

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "preprocess"))
from comp_eps import _KJ, _RHO_AIR, _BB, _PATM     # noqa: E402

DAY = 86400.0


def rho_vs(T_C):
    """Saturation vapour density over ice [kg/m^3], vectorised."""
    T = np.asarray(T_C, dtype=float) + 273.15
    Pvs = np.exp(_KJ[0] / T + _KJ[1] + _KJ[2] * T + _KJ[3] * T**2
                 + _KJ[4] * T**3 + _KJ[5] * np.log(T))
    return _RHO_AIR * _BB * Pvs / (_PATM - Pvs)


def read_vts(path):
    """(sigma, ice, nx, ny, extent_m) from one .vts."""
    r = vtk.vtkXMLStructuredGridReader()
    r.SetFileName(path)
    r.Update()
    out = r.GetOutput()
    dims = [0, 0, 0]
    out.GetDimensions(dims)
    nx, ny = dims[0], dims[1]
    pd = out.GetPointData()

    def arr(name):
        a = pd.GetArray(name)
        if a is None:
            raise KeyError(f"{path}: no array '{name}'")
        return vtk_to_numpy(a)

    phi = arr("IcePhase")
    rvs = rho_vs(arr("Temperature"))
    sig = (arr("VaporDensity") - rvs) / rvs
    b = out.GetBounds()
    # VTK structured points vary x fastest.
    return (sig.reshape(ny, nx), phi.reshape(ny, nx), nx, ny,
            (b[0], b[1], b[2], b[3]))


def time_map(run_dir):
    """{step: time_s} from SSA_evo.dat -- see plot_fields.py::_load_time_map_ssa."""
    p = os.path.join(run_dir, "SSA_evo.dat")
    m = {}
    if os.path.isfile(p):
        for line in open(p):
            f = line.split()
            if len(f) >= 4:
                try:
                    m[int(float(f[3]))] = float(f[2])
                except ValueError:
                    pass
    return m


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", required=True, help="run directory (holds vtkOut/)")
    ap.add_argument("--out", default=None, help="output .mp4 (default <dir>/supersaturation.mp4)")
    ap.add_argument("--fps", type=int, default=8)
    ap.add_argument("--dpi", type=int, default=150)
    ap.add_argument("--bound-from", default="smaller",
                    choices=("smaller", "larger", "under", "over"))
    ap.add_argument("--bound", type=float, default=None, help="explicit bound, overrides --bound-from")
    ap.add_argument("--percentile", type=float, default=100.0,
                    help="use this percentile instead of the true extreme")
    ap.add_argument("--ice-threshold", type=float, default=0.5)
    ap.add_argument("--keep-frames", action="store_true")
    args = ap.parse_args()

    run = os.path.abspath(args.dir)
    files = sorted(glob.glob(os.path.join(run, "vtkOut", "solV_*.vts")))
    if not files:
        print(f"no .vts files in {run}/vtkOut", file=sys.stderr)
        return 1
    tmap = time_map(run)
    steps = [int(re.search(r"_(\d+)\.vts", f).group(1)) for f in files]
    times = [tmap.get(s, float(s)) for s in steps]
    if not all(b > a for a, b in zip(times, times[1:])):
        print("⚠ frame times are not increasing — check SSA_evo.dat", file=sys.stderr)

    # -- pass 1: the colour bounds, over every frame ----------------------
    lo = hi = 0.0
    for f in files:
        sig, ice, *_ = read_vts(f)
        s = sig[ice < args.ice_threshold]
        if s.size == 0:
            continue
        if args.percentile >= 100.0:
            lo, hi = min(lo, s.min()), max(hi, s.max())
        else:
            p = args.percentile
            lo = min(lo, np.percentile(s, 100.0 - p))
            hi = max(hi, np.percentile(s, p))
    print(f"sigma over {len(files)} frames, air phase only (ice < {args.ice_threshold}):")
    print(f"  undersaturated extreme {lo:+.4e}   oversaturated extreme {hi:+.4e}")
    print(f"  |under|/|over| = {abs(lo)/abs(hi):.1f}x" if hi else "")

    if args.bound is not None:
        bound = abs(args.bound)
    else:
        bound = {"smaller": min(abs(lo), abs(hi)), "larger": max(abs(lo), abs(hi)),
                 "under": abs(lo), "over": abs(hi)}[args.bound_from]
    frac_clip = None
    print(f"  colour bounds: +/- {bound:.4e}   (--bound-from {args.bound_from})")

    # -- pass 2: render ----------------------------------------------------
    frames = tempfile.mkdtemp(prefix="sigma_frames_")
    cmap = cmocean.cm.balance.copy()
    cmap.set_bad("0.72")                      # ice
    nclip = ntot = 0
    for i, (f, t) in enumerate(zip(files, times)):
        sig, ice, nx, ny, ext = read_vts(f)
        m = np.ma.masked_where(ice >= args.ice_threshold, sig)
        nclip += int(np.sum(np.abs(m.compressed()) > bound))
        ntot += m.count()

        fig, ax = plt.subplots(figsize=(6.6, 6.0))
        im = ax.imshow(m, origin="lower", cmap=cmap, vmin=-bound, vmax=bound,
                       extent=[ext[0]*1e3, ext[1]*1e3, ext[2]*1e3, ext[3]*1e3],
                       interpolation="nearest")
        ax.set_xlabel("x [mm]"); ax.set_ylabel("y [mm]")
        ax.set_title(f"supersaturation in the pore space\n"
                     f"t = {t/DAY:6.2f} d     step {steps[i]}", fontsize=11)
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03, extend="both")
        cb.set_label(r"$\sigma=(\rho_v-\rho_{vs})/\rho_{vs}$")
        cb.formatter.set_powerlimits((0, 0))
        fig.tight_layout()
        fig.savefig(os.path.join(frames, f"f{i:05d}.png"), dpi=args.dpi)
        plt.close(fig)
    print(f"  clipped cells: {100.0*nclip/max(ntot,1):.2f}% of the pore space")

    out = args.out or os.path.join(run, "supersaturation.mp4")
    cmd = ["ffmpeg", "-y", "-framerate", str(args.fps),
           "-i", os.path.join(frames, "f%05d.png"),
           "-c:v", "libx264", "-pix_fmt", "yuv420p",
           "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2", out]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stderr[-2000:], file=sys.stderr)
        print(f"ffmpeg failed; frames left in {frames}", file=sys.stderr)
        return 1
    if args.keep_frames:
        print(f"  frames: {frames}")
    else:
        shutil.rmtree(frames, ignore_errors=True)
    print(f"\nwrote {out}  ({len(files)} frames @ {args.fps} fps)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

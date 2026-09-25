#!/usr/bin/env python3
"""make_packing_movie.py — animate a grain-PACKING sintering run.

The packing counterpart of make_neck_movie.py, and it deliberately shares that
script's visual language so the two read as one set:

  background   supersaturation sigma = rho_v/rho_vs(T) - 1, cmocean `balance`,
               re-sampled by centered_cmap so the pale middle lands on
               sigma = 0 -- blue undersaturated, red supersaturated, pale
               saturation over flat ice.
  ice          painted ON TOP as one opaque layer, cmocean `ice`, transparent
               below phi = 0.5 (ice_alpha_cmap). Every pixel is covered by the
               vapour base, so no pixel is left unowned at the phi = 0.5
               boundary.
  scale        AsinhNorm: log-like over the decades, linear through zero.

All four helpers -- SIGMA_SCALE, ice_alpha_cmap, centered_cmap, sigma_ticks --
are IMPORTED from make_neck_movie rather than reimplemented, so the two movies
cannot drift apart in convention. sigma itself comes from pplib.supersaturation,
which mirrors the solver's RhoVS_I, and uses the local Temperature array rather
than the -temp option because a run with a gradient has no single rho_vs.

WHY ASINH AND NOT SYMMETRIC LIMITS. A packing's sigma field is strongly
one-sided: on the 2026-09-16 pilot it spans about -2.8e-4 to +4.7e-5, i.e. the
undersaturated extreme is ~6x the oversaturated one. Symmetric limits at the
larger extreme wash the structure out; at the smaller one they clip a fifth of
the pore space. asinh keeps both -- the decades where the pore lives and the
sign change at zero -- and centered_cmap keeps the diverging map's promise that
pale means neutral. `--symmetric` forces the older behaviour, which puts
sigma = 0 at the midpoint of the BAR at the cost of range.

A second panel tracks the run's own k_eff(t) with a dot at the current frame,
read from whichever k_eff CSVs the run directory holds -- so a frame can be
pointed at while saying "this is where the conductivity is".

Usage:
    python make_packing_movie.py <run_dir> [--out FILE.mp4] [--fps 10]
        [--stride N] [--dpi 150] [--frame-png STEP] [--sat-clip P]
        [--symmetric] [--no-keff] [--cmap NAME]
"""
from __future__ import annotations

import argparse
import csv
import glob
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import AsinhNorm
import cmocean

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
from pplib import read_vts, step_of, step_times          # noqa: E402
from make_neck_movie import (SIGMA_SCALE, ice_alpha_cmap,  # noqa: E402
                             centered_cmap, sigma_ticks)
import pplib                                              # noqa: E402

DAY = 86400.0
WANT = ("IcePhase", "VaporDensity", "Temperature")


def load_keff(run: Path):
    """{law: (t_days, k_iso)} from every k_eff*.csv in the run directory."""
    out = {}
    for p in sorted(run.glob("k_eff*.csv")):
        m = re.fullmatch(r"k_eff(?:_(\w+))?\.csv", p.name)
        if not m:
            continue
        law = m.group(1) or "arith"
        try:
            rows = list(csv.DictReader(p.open()))
            t = np.array([float(r["time"]) for r in rows]) / DAY
            k = np.array([float(r["k_iso"]) for r in rows])
        except Exception:                                  # noqa: BLE001
            continue
        if t.size:
            out[law] = (t, k)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir", type=Path)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--fps", type=int, default=10)
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--dpi", type=int, default=150)
    ap.add_argument("--frame-png", type=int, default=None,
                    help="render only this step, to PNG, and stop")
    ap.add_argument("--sat-clip", type=float, default=0.1,
                    help="percentile clipped off each end when setting the "
                         "colour range (default 0.1)")
    ap.add_argument("--symmetric", action="store_true",
                    help="symmetric limits about zero instead of asinh over "
                         "the clipped range")
    ap.add_argument("--cmap", default="balance")
    ap.add_argument("--no-keff", action="store_true")
    ap.add_argument("--keep-frames", action="store_true")
    args = ap.parse_args()

    run = args.run_dir.resolve()
    files = sorted(glob.glob(str(run / "vtkOut" / "solV_*.vts")), key=step_of)
    if not files:
        print(f"no solV_*.vts under {run}/vtkOut", file=sys.stderr)
        return 1
    files = files[::max(1, args.stride)]
    tmap = step_times(str(run))
    steps = [step_of(f) for f in files]
    times = [tmap.get(s, float(s)) for s in steps]

    # -- pass 1: the colour range, over every frame -----------------------
    lo = hi = None
    raw_lo = raw_hi = None
    for f in files:
        fl, _, _ = read_vts(f, want=WANT)
        if "VaporDensity" not in fl or "Temperature" not in fl:
            print("  snapshots carry no vapour field; nothing to animate",
                  file=sys.stderr)
            return 1
        s = SIGMA_SCALE * pplib.supersaturation(fl["VaporDensity"], fl["Temperature"])
        s = s[fl["IcePhase"] < 0.5]                # pore only sets the scale
        if s.size == 0:
            continue
        a, b = np.percentile(s, [args.sat_clip, 100.0 - args.sat_clip])
        lo = a if lo is None else min(lo, a)
        hi = b if hi is None else max(hi, b)
        raw_lo = s.min() if raw_lo is None else min(raw_lo, s.min())
        raw_hi = s.max() if raw_hi is None else max(raw_hi, s.max())

    if args.symmetric:
        v = max(abs(lo), abs(hi))
        lo, hi = -v, v
    norm = AsinhNorm(linear_width=max(max(abs(lo), abs(hi)) / 300.0, 1e-12),
                     vmin=lo, vmax=hi)
    base = getattr(cmocean.cm, args.cmap, None) or plt.get_cmap(args.cmap)
    vapcm = centered_cmap(base, norm)
    icecm = ice_alpha_cmap()
    print(f"  sigma x{SIGMA_SCALE:g}: raw {raw_lo:+.4g} .. {raw_hi:+.4g}, "
          f"p{args.sat_clip:g} {lo:+.4g} .. {hi:+.4g}")

    keff = {} if args.no_keff else load_keff(run)
    if keff:
        print("  k_eff panel: " + ", ".join(f"{k} ({len(v[0])} pts)"
                                            for k, v in keff.items()))

    # -- pass 2: render ----------------------------------------------------
    frames = tempfile.mkdtemp(prefix="packing_frames_")
    ncol = 1 if not keff else 2
    n_out = 0
    for i, (fn, t) in enumerate(zip(files, times)):
        if args.frame_png is not None and steps[i] != args.frame_png:
            continue
        fl, X, Y = read_vts(fn, want=WANT)
        sig = SIGMA_SCALE * pplib.supersaturation(fl["VaporDensity"], fl["Temperature"])
        phi = fl["IcePhase"]

        fig = plt.figure(figsize=(7.0 if ncol == 1 else 11.4, 6.2))
        gs = fig.add_gridspec(1, ncol, width_ratios=[1] if ncol == 1 else [1, 0.82])
        ax = fig.add_subplot(gs[0, 0])
        XX, YY = X * 1e3, Y * 1e3
        vap = ax.pcolormesh(XX, YY, sig, cmap=vapcm, norm=norm, shading="gouraud")
        ax.pcolormesh(XX, YY, phi, cmap=icecm, vmin=0.0, vmax=1.0, shading="gouraud")
        ax.set_aspect("equal")
        ax.set_xlabel("x [mm]"); ax.set_ylabel("y [mm]")
        ax.set_title(f"t = {t/DAY:6.2f} d      step {steps[i]}", fontsize=10)
        cb = fig.colorbar(vap, ax=ax, fraction=0.046, pad=0.03, extend="both",
                          ticks=sigma_ticks(norm))
        # Same label and 3-sig-fig formatter as make_neck_movie, so a packing
        # frame and a neck frame can sit on one slide without the reader
        # having to re-learn the bar.
        cb.set_label(r"supersaturation  $\sigma = \rho_v/\rho_{vs}-1$   "
                     r"[$\times 10^{-4}$]", fontsize=9)
        cb.ax.yaxis.set_major_formatter(plt.FuncFormatter(
            lambda v, _p: f"{v:.3g}"))
        cb.ax.tick_params(labelsize=8)

        if keff:
            axk = fig.add_subplot(gs[0, 1])
            for law, (tk, kk) in sorted(keff.items()):
                c = "#c0392b" if law == "arith" else "#1f6fb4"
                axk.plot(tk, kk, color=c, lw=1.4,
                         label="arith (legacy)" if law == "arith" else law)
                j = int(np.argmin(np.abs(tk - t / DAY)))
                axk.plot([tk[j]], [kk[j]], "o", color=c, ms=6, zorder=5)
            axk.axvline(t / DAY, color="0.6", lw=0.8, ls="--")
            axk.set_xlabel("time [days]")
            axk.set_ylabel(r"$k_{\mathrm{eff}}$  [W m$^{-1}$ K$^{-1}$]")
            axk.set_title("effective conductivity", fontsize=10)
            axk.grid(alpha=0.25, lw=0.6)
            axk.legend(frameon=False, fontsize=9, loc="lower right")

        fig.tight_layout()
        if args.frame_png is not None:
            out = args.out or run / f"frame_{steps[i]:05d}.png"
            fig.savefig(out, dpi=args.dpi, bbox_inches="tight")
            plt.close(fig)
            print(f"wrote {out}")
            shutil.rmtree(frames, ignore_errors=True)
            return 0
        fig.savefig(os.path.join(frames, f"f{n_out:05d}.png"), dpi=args.dpi)
        plt.close(fig)
        n_out += 1

    out = args.out or run / "packing_sintering.mp4"
    cmd = ["ffmpeg", "-y", "-framerate", str(args.fps),
           "-i", os.path.join(frames, "f%05d.png"),
           "-c:v", "libx264", "-pix_fmt", "yuv420p",
           "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2", str(out)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stderr[-2000:], file=sys.stderr)
        print(f"ffmpeg failed; frames kept in {frames}", file=sys.stderr)
        return 1
    if args.keep_frames:
        print(f"  frames: {frames}")
    else:
        shutil.rmtree(frames, ignore_errors=True)
    print(f"\nwrote {out}  ({n_out} frames @ {args.fps} fps)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

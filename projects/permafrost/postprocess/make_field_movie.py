#!/usr/bin/env python3
"""make_field_movie.py — render an ice/vapor animation from solV_*.vts, in pure
Python (matplotlib + cmocean), no ParaView.

WHAT IT DRAWS (per the user's ParaView recipe, reproduced faithfully)
  * ICE region  (IcePhase >= 0.5): coloured by IcePhase through cmocean 'ice'.
  * AIR region  (IcePhase <  0.5): coloured by VaporDensity through cmocean 'amp'.
The phi = 0.5 contour is the partition, exactly as the two ParaView iso-volumes
[0.5, 1.01] and [-0.01, 0.5] intended.

WHY NOT TWO MASKED LAYERS (and why ParaView dropped pixels): thresholding into
two disjoint volumes leaves the phi=0.5 boundary cells owned by neither at the
export resolution -> the "missing pixels inside the grains" on animation export.
Here the VAPOR field is drawn as a FULL-DOMAIN base layer and the ICE layer is
drawn OPAQUE on top, masked to phi>=0.5. Every pixel is covered by the base, so
nothing can be dropped; the ice simply paints over the air where it exists.

CURVILINEAR MESH: the bumpy floor makes the grid deformed (the bottom row of
nodes follows the bumps). We plot with the ACTUAL node coordinates via
pcolormesh, not an imshow on a regular grid -- otherwise the floor bumps are
misplaced. The mesh is static (Eulerian), so the QuadMeshes are built once and
only their colour arrays are updated each frame.

VAPOR IS NEARLY FLAT: in a saturated cell VaporDensity varies only in its ~5th
significant digit (e.g. 8.4868e-4..8.4873e-4). The colour range is therefore
auto-scaled to robust percentiles of the vapour in the AIR region, sampled
GLOBALLY across the run so the scale is fixed and comparable frame-to-frame.
Override with --vmin-vapor/--vmax-vapor. The colourbar prints the true values.

Usage:
    python make_field_movie.py <run_dir> [--out FILE.mp4] [--stride N]
        [--fps 24] [--dpi 150] [--vmin-vapor V --vmax-vapor V]
        [--frame-png STEP]   # render one frame to PNG and exit (a preview)
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
from matplotlib.animation import FFMpegWriter
from matplotlib.cm import ScalarMappable
from matplotlib.colors import ListedColormap, Normalize
import cmocean


def ice_alpha_cmap(n=256):
    """cmocean 'ice' for phi in [0.5,1], fully TRANSPARENT below 0.5.

    Lets the ice layer be drawn with gouraud shading (smooth, anti-aliased
    grain edges -- the fix for the 'pixelated grains') WITHOUT a masked array,
    which gouraud pcolormesh cannot handle. The transparent rows carry the
    same RGB as the phi=0.5 ice colour (not black), so gouraud interpolation
    across the interface ramps only alpha -- no dark fringe. Render with
    vmin=0, vmax=1; show a clean 0.5..1 ice colourbar via a separate mappable.
    """
    ice = cmocean.cm.ice
    edge = ice(0.0)
    lut = np.zeros((n, 4))
    for i in range(n):
        phi = i / (n - 1)
        if phi < 0.5:
            lut[i] = (edge[0], edge[1], edge[2], 0.0)
        else:
            c = ice((phi - 0.5) / 0.5)
            lut[i] = (c[0], c[1], c[2], 1.0)
    return ListedColormap(lut)


def read_vts(fn, want=("IcePhase", "VaporDensity")):
    """Return (fields dict, X[ny,nx], Y[ny,nx]) from a solV_*.vts snapshot."""
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
    fields = {}
    for da in root.findall(".//PointData/DataArray"):
        if da.get("Name") in want:
            fields[da.get("Name")] = decode(da).reshape(ny, nx)
    return fields, pts[:, :, 0], pts[:, :, 1]


def step_of(fn):
    return int(re.search(r"solV_(\d+)\.vts", fn).group(1))


def step_times(run_dir):
    """step -> time[s] from the monitor tables in outp.txt (8-pipe rows)."""
    tmap = {}
    outp = Path(run_dir) / "outp.txt"
    if not outp.exists():
        return tmap
    pat = re.compile(r"^\s+(\d+)\s+\|\s+([0-9.eE+-]+)\s+\|")
    for line in open(outp, errors="replace"):
        if line.count("|") == 8:
            m = pat.match(line)
            if m:
                tmap[int(m.group(1))] = float(m.group(2))
    return tmap


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir", type=Path)
    ap.add_argument("--out", type=Path, default=None,
                    help="output mp4 (default: <run_dir>/ice_vapor_movie.mp4)")
    ap.add_argument("--stride", type=int, default=1,
                    help="use every Nth snapshot (default 1 = all)")
    ap.add_argument("--fps", type=int, default=24)
    ap.add_argument("--dpi", type=int, default=200)
    ap.add_argument("--frames-dir", type=Path, default=None,
                    help="where to save per-frame PNGs (default: "
                         "<run_dir>/movie_frames/)")
    ap.add_argument("--no-frames", action="store_true",
                    help="render only the mp4, skip the per-frame PNG dump")
    ap.add_argument("--vmin-vapor", type=float, default=None)
    ap.add_argument("--vmax-vapor", type=float, default=None)
    ap.add_argument("--frame-png", type=int, default=None,
                    help="render only the snapshot with this step index to a PNG "
                         "(preview) and exit")
    ap.add_argument("--title", type=str, default=None,
                    help="static title line (default: auto from -temp/-humidity "
                         "in the run's opts). Pass '' for no static line.")
    args = ap.parse_args()

    files = sorted(glob.glob(str(args.run_dir / "vtkOut" / "solV_*.vts")),
                   key=step_of)
    if not files:
        sys.exit(f"no solV_*.vts under {args.run_dir}/vtkOut")
    tmap = step_times(args.run_dir)

    # Static title: describe the PHYSICS, not the folder. The folder name is a
    # filing label, not a plot caption. Auto-build from the run's opts unless
    # the user supplies --title.
    if args.title is None:
        opt = {}
        for f in glob.glob(str(args.run_dir / "*.opts")):
            for line in open(f, errors="replace"):
                line = line.split("#")[0].split()
                if len(line) >= 2 and line[0] in ("-temp", "-humidity"):
                    opt[line[0]] = line[1]
        bits = []
        if "-temp" in opt:
            bits.append(f"T = {float(opt['-temp']):g} °C")
        if "-humidity" in opt:
            bits.append(f"RH = {float(opt['-humidity'])*100:g}%")
        static_title = "   ".join(bits)
    else:
        static_title = args.title

    # Static mesh: read coordinates once from the first snapshot.
    _, X, Y = read_vts(files[0])

    # Global vapour colour range: robust percentiles over the AIR region,
    # sampled across the run so the scale is fixed frame-to-frame.
    if args.vmin_vapor is None or args.vmax_vapor is None:
        sample = files[:: max(1, len(files) // 40)]
        los, his = [], []
        for fn in sample:
            f, _, _ = read_vts(fn)
            air = f["IcePhase"] < 0.5
            if air.any():
                lo, hi = np.percentile(f["VaporDensity"][air], [0.5, 99.5])
                los.append(lo); his.append(hi)
        vmin = args.vmin_vapor if args.vmin_vapor is not None else min(los)
        vmax = args.vmax_vapor if args.vmax_vapor is not None else max(his)
    else:
        vmin, vmax = args.vmin_vapor, args.vmax_vapor
    if vmax <= vmin:
        vmax = vmin + 1e-30

    # ---- Figure: wide domain (aspect ~5:1) with colourbars BELOW ---------
    # Vertical colourbars collide on a short-wide plot; horizontal bars under
    # the axes have room. Axes rect leaves space top (title) and bottom (bars).
    Lx, Ly = X.max(), Y.max()
    fig_w = 16.0
    ax_frac_h = 0.56
    fig_h = fig_w * 0.90 * (Ly / Lx) / ax_frac_h + 1.6
    fig = plt.figure(figsize=(fig_w, fig_h))
    ax = fig.add_axes([0.05, 0.30, 0.90, ax_frac_h])
    ax.set_aspect("equal")
    ax.set_xlim(X.min(), X.max()); ax.set_ylim(Y.min(), Y.max())
    ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]")

    f0, _, _ = read_vts(files[0])
    phi0 = f0["IcePhase"]; vap0 = f0["VaporDensity"]

    # Base: vapour over the WHOLE domain, gouraud (smooth). Guarantees coverage.
    base = ax.pcolormesh(X, Y, vap0, cmap=cmocean.cm.amp,
                         vmin=vmin, vmax=vmax, shading="gouraud", rasterized=True)
    # Ice on top: gouraud with the alpha-cutoff cmap (transparent below 0.5).
    # gouraud interpolates between nodes -> anti-aliased grain edges (the fix
    # for the pixelation) without resampling the .vts. No masked array, so no
    # gouraud+mask crash.
    ice_cmap = ice_alpha_cmap()
    ice = ax.pcolormesh(X, Y, phi0, cmap=ice_cmap,
                        vmin=0.0, vmax=1.0, shading="gouraud", rasterized=True)

    # Two horizontal colourbars, side by side under the axes -- no overlap.
    cax_i = fig.add_axes([0.10, 0.15, 0.35, 0.03])
    cax_v = fig.add_axes([0.57, 0.15, 0.35, 0.03])
    sm_i = ScalarMappable(norm=Normalize(0.5, 1.0), cmap=cmocean.cm.ice)
    cb_i = fig.colorbar(sm_i, cax=cax_i, orientation="horizontal")
    cb_i.set_label(r"IcePhase $\phi_i$ (ice region)")
    cb_v = fig.colorbar(base, cax=cax_v, orientation="horizontal")
    cb_v.set_label(r"VaporDensity [kg/m$^3$] (air region)")
    cb_v.formatter.set_powerlimits((0, 0)); cb_v.update_ticks()

    title = ax.set_title("")

    def draw(fn):
        f, _, _ = read_vts(fn)
        phi = f["IcePhase"]; vap = f["VaporDensity"]
        base.set_array(vap.ravel())
        ice.set_array(phi.ravel())
        st = step_of(fn)
        t = tmap.get(st)
        tstr = f"t = {t/86400:.2f} days" if t is not None else f"step {st}"
        ice_frac = 100.0 * (phi >= 0.5).mean()
        dyn = f"{tstr}     ice area {ice_frac:.1f}%"
        title.set_text(f"{static_title}\n{dyn}" if static_title else dyn)

    # ---- Preview single frame -------------------------------------------
    if args.frame_png is not None:
        fn = min(files, key=lambda p: abs(step_of(p) - args.frame_png))
        draw(fn)
        out = args.run_dir / f"frame_{step_of(fn):05d}.png"
        fig.savefig(out, dpi=args.dpi)
        print(f"preview -> {out}  (vapor range {vmin:.6e}..{vmax:.6e})")
        return

    frames = files[:: args.stride]
    out = args.out or (args.run_dir / "ice_vapor_movie.mp4")
    # Screenshots: save each rendered frame as a PNG too, unless --no-frames.
    frames_dir = None
    if not args.no_frames:
        frames_dir = args.frames_dir or (args.run_dir / "movie_frames")
        frames_dir.mkdir(parents=True, exist_ok=True)

    writer = FFMpegWriter(fps=args.fps, bitrate=-1,
                          metadata={"title": args.run_dir.name})
    print(f"rendering {len(frames)} frames -> {out}")
    print(f"  vapor colour range: {vmin:.6e} .. {vmax:.6e} kg/m^3")
    if frames_dir:
        print(f"  screenshots -> {frames_dir}/frame_*.png")
    with writer.saving(fig, str(out), dpi=args.dpi):
        for i, fn in enumerate(frames):
            draw(fn)
            writer.grab_frame()
            if frames_dir:
                fig.savefig(frames_dir / f"frame_{step_of(fn):05d}.png",
                            dpi=args.dpi)
            if i % 50 == 0:
                print(f"  frame {i}/{len(frames)}  ({step_of(fn)})", flush=True)
    print(f"movie -> {out}")
    if frames_dir:
        print(f"screenshots -> {frames_dir}")


if __name__ == "__main__":
    main()

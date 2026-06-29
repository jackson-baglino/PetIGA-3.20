#!/usr/bin/env python3
"""
make_movie.py — render a DSM ice/vapor movie from a PVD time series.

Typical ParaView GUI workflow:
  - "ice volume": IsoVolume of IcePhase in [0.5, 1.1], colored by IcePhase
    with the cmocean "ice" colormap.
  - "air volume": IsoVolume of IcePhase in [-0.1, 0.5],
    colored by VaporDensity.

Must be run with ParaView's own Python, not a normal venv:

    /Applications/ParaView-5.12.0.app/Contents/bin/pvpython postprocess/make_movie.py --dir <rundir>

Before first use, generate the colormap preset once with the regular
project Python (it has cmocean; pvpython does not):

    python3 postprocess/make_cmocean_preset.py ice

Each run also needs a dense, true-NURBS-interpolated dsm_highres.pvd
generated first (same regular project Python, not pvpython -- it needs igakit):

    python3 postprocess/plotDSM.py --dir <rundir>

------------------------------------------------------------------------
Linear-time playback
------------------------------------------------------------------------
PetIGA's adaptive dt means snapshots are dense early and sparse late.
This script samples --n-frames frames at EVENLY SPACED simulated times
across [t-start, t-end] (default: the full run), independent of the
underlying snapshot positions. Use --interpolate (on by default) to
linearly blend between sparse late snapshots instead of holding them.

------------------------------------------------------------------------
Usage examples
------------------------------------------------------------------------
  # Full run, 600 frames linear in simulated time, 30 fps
  pvpython postprocess/make_movie.py --dir /path/to/run

  # Fixed vapor colorbar bounds
  pvpython postprocess/make_movie.py --dir /path/to/run --vapor-range 1e-4 5e-4

  # First 2 days only
  pvpython postprocess/make_movie.py --dir /path/to/run --t-end 172800

  # 2x supersampled frames
  pvpython postprocess/make_movie.py --dir /path/to/run --supersample 2
"""

import argparse
import json
import os
import shutil
import subprocess

from paraview.simple import (
    OpenDataFile, IsoVolume, TemporalInterpolator, Show, ColorBy,
    GetColorTransferFunction, GetActiveViewOrCreate, GetAnimationScene,
    Render, SaveScreenshot,
)
from paraview import servermanager
from vtk.numpy_interface import dataset_adapter as dsa
import numpy as np


def save_vector_colorbar(lut, label, out_path):
    """Export a standalone vertical colorbar as SVG."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, Normalize
    from matplotlib.colorbar import ColorbarBase

    pts = list(lut.RGBPoints)
    n = len(pts) // 4
    values = [pts[4 * i] for i in range(n)]
    colors = [(pts[4 * i + 1], pts[4 * i + 2], pts[4 * i + 3]) for i in range(n)]
    lo, hi = values[0], values[-1]
    positions = [(v - lo) / (hi - lo) if hi > lo else 0.0 for v in values]
    cmap = LinearSegmentedColormap.from_list("custom", list(zip(positions, colors)), N=256)

    fig, ax = plt.subplots(figsize=(0.5, 6))
    cb = ColorbarBase(ax, cmap=cmap, norm=Normalize(vmin=lo, vmax=hi), orientation="vertical")
    cb.set_label(label, fontsize=14, fontweight="bold", color="white")
    cb.ax.tick_params(labelsize=12, colors="white")
    for tick_label in cb.ax.get_yticklabels():
        tick_label.set_fontweight("bold")
    fig.savefig(out_path, format="svg", bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"  wrote {out_path}  (range [{lo:.4g}, {hi:.4g}])")


def find_pvd(run_dir: str) -> str:
    """Require dsm_highres.pvd; fall back to dsm.pvd with a warning."""
    highres_path = os.path.join(run_dir, "dsm_highres.pvd")
    if os.path.isfile(highres_path):
        return highres_path
    coarse_path = os.path.join(run_dir, "dsm.pvd")
    if os.path.isfile(coarse_path):
        raise FileNotFoundError(
            f"No dsm_highres.pvd in {run_dir} (found only the coarse dsm.pvd). "
            f"Generate the dense, true-NURBS interpolated version first:\n"
            f"    python3 postprocess/plotDSM.py --dir {run_dir}\n"
            f"then re-run make_movie.py. Pass --pvd {coarse_path} explicitly "
            f"if you really want the coarse mesh instead.")
    raise FileNotFoundError(
        f"No dsm_highres.pvd or dsm.pvd found in {run_dir}")


def auto_vapor_range(air_volume, timestep_values, n_samples, lo_pct, hi_pct):
    """Sample n_samples timesteps, pool VaporDensity, return percentile range."""
    n = len(timestep_values)
    idx = np.unique(np.linspace(0, n - 1, min(n_samples, n)).astype(int))
    pooled = []
    for i in idx:
        air_volume.UpdatePipeline(timestep_values[i])
        data = dsa.WrapDataObject(servermanager.Fetch(air_volume))
        arr = np.asarray(data.PointData["VaporDensity"]).ravel()
        if arr.size:
            pooled.append(arr)
    pooled = np.concatenate(pooled)
    lo, hi = np.percentile(pooled, [lo_pct, hi_pct])
    return float(lo), float(hi)


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dir", default=".", help="run directory (default: .)")
    p.add_argument("--pvd", default=None,
                   help="override path to the .pvd file (default: auto-detect dsm_highres.pvd)")
    p.add_argument("--out", default=None,
                   help="output movie path (default: <dir>/movie.mp4)")
    p.add_argument("--n-frames", type=int, default=600,
                   help="number of frames, evenly spaced in simulated time (default: 600)")
    p.add_argument("--fps", type=int, default=30, help="playback frame rate (default: 30)")
    p.add_argument("--t-start", type=float, default=None)
    p.add_argument("--t-end",   type=float, default=None)
    p.add_argument("--no-interpolate", action="store_true",
                   help="hold nearest snapshot instead of blending between sparse ones")
    p.add_argument("--ice-colormap", default=None,
                   help="path to ParaView JSON preset for ice volume "
                        "(default: postprocess/colormaps/cmocean_ice.json)")
    p.add_argument("--vapor-colormap", default=None,
                   help="ParaView preset name or path for VaporDensity colormap "
                        "(default: postprocess/colormaps/cmocean_balance.json)")
    p.add_argument("--vapor-range", type=float, nargs=2, default=None,
                   help="fixed [min, max] for the vapor-density colorbar; "
                        "if omitted, auto-computed from sampled timesteps")
    p.add_argument("--vapor-percentile", type=float, nargs=2, default=[2.0, 98.0])
    p.add_argument("--vapor-range-samples", type=int, default=40)
    p.add_argument("--resolution", default=None, help="WxH override (default: native)")
    p.add_argument("--supersample", type=float, default=1.0)
    p.add_argument("--min-height", type=int, default=1080)
    p.add_argument("--max-width",  type=int, default=4096)
    p.add_argument("--no-colorbars", action="store_true")
    p.add_argument("--keep-frames",  action="store_true")
    args = p.parse_args()

    pvd_path = args.pvd or find_pvd(args.dir)
    out_path = args.out or os.path.join(args.dir, "movie.mp4")
    script_dir = os.path.dirname(os.path.abspath(__file__))
    ice_preset_path = args.ice_colormap or os.path.join(script_dir, "colormaps", "cmocean_ice.json")
    vapor_colormap  = args.vapor_colormap or os.path.join(script_dir, "colormaps", "cmocean_balance.json")

    print(f"Reading {pvd_path}")
    reader = OpenDataFile(pvd_path)
    reader.UpdatePipeline(0.0)
    timestep_values = list(reader.TimestepValues)
    if not timestep_values:
        raise RuntimeError("No timesteps found in the PVD file")

    t_start = args.t_start if args.t_start is not None else timestep_values[0]
    t_end   = args.t_end   if args.t_end   is not None else timestep_values[-1]
    print(f"Simulated time window: [{t_start:.6g}, {t_end:.6g}] "
          f"({len(timestep_values)} source snapshots, {args.n_frames} output frames)")

    source = reader
    if not args.no_interpolate:
        source = TemporalInterpolator(Input=reader)

    ice_volume = IsoVolume(Input=source, InputScalars=["POINTS", "IcePhase"],
                            ThresholdRange=[0.5, 1.1])
    air_volume = IsoVolume(Input=source, InputScalars=["POINTS", "IcePhase"],
                            ThresholdRange=[-0.1, 0.5])

    if args.vapor_range is not None:
        vmin, vmax = args.vapor_range
    else:
        vmin, vmax = auto_vapor_range(
            air_volume, timestep_values, args.vapor_range_samples, *args.vapor_percentile)
        print(f"Auto vapor-density range: [{vmin:.4g}, {vmax:.4g}]")

    view = GetActiveViewOrCreate("RenderView")
    view.InteractionMode = "2D"
    view.OrientationAxesVisibility = 0

    ice_display = Show(ice_volume, view)
    ColorBy(ice_display, ("POINTS", "IcePhase"))
    presets = servermanager.vtkSMTransferFunctionPresets.GetInstance()
    if not presets.HasPreset("cmocean_ice"):
        if not presets.ImportPresets(ice_preset_path):
            raise RuntimeError(f"Failed to import ice colormap preset: {ice_preset_path}")
    ice_lut = GetColorTransferFunction("IcePhase")
    ice_lut.ApplyPreset("cmocean_ice", True)
    ice_lut.RescaleTransferFunction(0.0, 1.0)

    air_display = Show(air_volume, view)
    ColorBy(air_display, ("POINTS", "VaporDensity"))
    vapor_lut = GetColorTransferFunction("VaporDensity")
    if os.path.isfile(vapor_colormap):
        preset_name = json.load(open(vapor_colormap))[0]["Name"]
        if not presets.HasPreset(preset_name):
            presets.ImportPresets(vapor_colormap)
        vapor_lut.ApplyPreset(preset_name, True)
    else:
        vapor_lut.ApplyPreset(vapor_colormap, True)
    vapor_lut.RescaleTransferFunction(vmin, vmax)

    ice_display.SetScalarBarVisibility(view, False)
    air_display.SetScalarBarVisibility(view, False)

    if args.resolution:
        w, h = (int(v) for v in args.resolution.split("x"))
    else:
        ext = reader.GetDataInformation().GetExtent()
        nx, ny = ext[1] - ext[0] + 1, ext[3] - ext[2] + 1
        w = int(round(nx * args.supersample))
        h = int(round(ny * args.supersample))
    if args.min_height and h < args.min_height:
        scale = args.min_height / h
        w, h = int(round(w * scale)), args.min_height
        if w > args.max_width:
            args.max_width = w
    if w > args.max_width:
        scale = args.max_width / w
        w, h = args.max_width, max(1, int(round(h * scale)))
    w, h = w + (w % 2), h + (h % 2)
    view.ViewSize = [w, h]
    print(f"Render resolution: {w}x{h}")

    Render(view)
    b = reader.GetDataInformation().GetBounds()
    xmid, ymid = 0.5 * (b[0] + b[1]), 0.5 * (b[2] + b[3])
    cam = view.GetActiveCamera()
    cam.SetParallelProjection(1)
    cam.SetFocalPoint(xmid, ymid, 0.0)
    cam.SetPosition(xmid, ymid, 1.0)
    cam.SetViewUp(0.0, 1.0, 0.0)
    cam.SetParallelScale(0.5 * (b[3] - b[2]))

    if not args.no_colorbars:
        base, _ = os.path.splitext(out_path)
        save_vector_colorbar(ice_lut, "Ice phase", base + "_ice_colorbar.svg")
        save_vector_colorbar(vapor_lut, "Vapor density [kg/m³]", base + "_vapor_colorbar.svg")

    frame_dir = out_path + "_frames"
    if os.path.isdir(frame_dir):
        shutil.rmtree(frame_dir)
    os.makedirs(frame_dir)

    scene = GetAnimationScene()
    scene.UpdateAnimationUsingDataTimeSteps()

    n = args.n_frames
    for i in range(n):
        t = t_start + (t_end - t_start) * (i / (n - 1) if n > 1 else 0.0)
        scene.AnimationTime = t
        Render(view)
        frame_path = os.path.join(frame_dir, f"frame_{i:05d}.png")
        SaveScreenshot(frame_path, view, ImageResolution=[w, h])
        if i % 50 == 0 or i == n - 1:
            print(f"  frame {i+1}/{n}  t={t:.6g}")

    print(f"Encoding {out_path} at {args.fps} fps")
    subprocess.run([
        "ffmpeg", "-y", "-framerate", str(args.fps),
        "-i", os.path.join(frame_dir, "frame_%05d.png"),
        "-c:v", "libx264", "-pix_fmt", "yuv420p", out_path,
    ], check=True)

    if not args.keep_frames:
        shutil.rmtree(frame_dir)

    print(f"Done: {out_path}")


if __name__ == "__main__":
    main()

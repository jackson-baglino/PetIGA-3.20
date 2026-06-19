#!/usr/bin/env python3
"""
make_movie.py — render the ice/air two-phase movie from a PetIGA run's PVD
time series, the way it's normally built by hand in the ParaView GUI:

  - "ice volume": IsoVolume of IcePhase in [0.5, 1.1], colored by IcePhase
    with the cmocean "ice" colormap.
  - "air volume": IsoVolume of IcePhase in [-0.1, 0.5] (i.e. air = 1-ice),
    colored by VaporDensity.

Must be run with ParaView's own Python, not this project's venv:

    /Applications/ParaView-5.12.0.app/Contents/bin/pvpython postprocess/make_movie.py --dir <rundir>

(adjust the pvpython path for your ParaView install/platform).

Before first use, generate the ice colormap preset once with the regular
project Python (it has cmocean; pvpython does not):

    python3 postprocess/make_cmocean_preset.py ice

------------------------------------------------------------------------
Linear-time playback
------------------------------------------------------------------------
PetIGA's adaptive dt means snapshots are dense early (small dt) and sparse
late (large dt). A naive "one video frame per snapshot" movie is therefore
*not* linear in simulated time -- it looks like it crawls early and
speeds up late, because each frame early covers a tiny dt and each frame
late covers a huge one, but all frames get equal screen time.

This script instead samples --n-frames frames at EVENLY SPACED simulated
times across [t-start, t-end] (default: the full run), independent of
where the underlying snapshots fall. Two consequences of that, both
inherent to "linear time" and not fully avoidable without changing what
"linear" means:

  - Early on, many original snapshots can fall inside a single video
    frame's time slot -- those in-between snapshots are simply not
    individually visible (the fast dynamics get compressed, not skipped:
    the frame shown is whichever moment in that slot you sample, by
    default the nearest snapshot).
  - Late on, one snapshot may need to cover many video frames' worth of
    time -- by default those frames just repeat the same (held) data.

--interpolate (on by default) softens the second issue: it inserts a
TemporalInterpolator filter so frames that fall *between* two sparse late
snapshots are linearly blended instead of held static. It does nothing
for the first issue (compression of dense early data is unavoidable
under strictly linear time) and isn't real new physics -- just a smoother
visual transition between the two real frames that bracket each video
frame.

If the early dynamics are the point of the movie, consider --t-end well
short of the full run (a separate "early" movie at linear time over just
that window), or ask about a variable-speed two-segment mode instead of
truly linear time -- that's a bigger script change, not implemented here.

------------------------------------------------------------------------
Usage examples
------------------------------------------------------------------------
  # Full run, 600 frames linear in simulated time, 30 fps, auto vapor range
  pvpython postprocess/make_movie.py --dir /path/to/run

  # Fixed vapor colorbar bounds instead of auto-detecting
  pvpython postprocess/make_movie.py --dir /path/to/run --vapor-range 1e-4 5e-4

  # Only the first 2 days, denser sampling, no temporal interpolation
  pvpython postprocess/make_movie.py --dir /path/to/run --t-end 172800 \\
      --n-frames 900 --no-interpolate
"""

import argparse
import json
import os
import shutil
import subprocess

from paraview.simple import (
    OpenDataFile, IsoVolume, TemporalInterpolator, Show, ColorBy,
    GetColorTransferFunction, GetActiveViewOrCreate, GetAnimationScene,
    Render, SaveScreenshot, ResetCamera,
)
from paraview import servermanager
from vtk.numpy_interface import dataset_adapter as dsa
import numpy as np


def find_pvd(run_dir: str) -> str:
    for name in ("permafrost_highres.pvd", "permafrost.pvd"):
        path = os.path.join(run_dir, name)
        if os.path.isfile(path):
            return path
    raise FileNotFoundError(
        f"No permafrost_highres.pvd or permafrost.pvd found in {run_dir}")


def auto_vapor_range(reader, timestep_values, n_samples, lo_pct, hi_pct):
    """Sample n_samples evenly-spaced timesteps, pool VaporDensity across
    them, and return a global [lo_pct, hi_pct] percentile range -- a fixed
    colorbar range that stays meaningful across the whole movie instead of
    auto-rescaling (and thus changing meaning) every frame."""
    n = len(timestep_values)
    idx = np.unique(np.linspace(0, n - 1, min(n_samples, n)).astype(int))
    pooled = []
    for i in idx:
        reader.UpdatePipeline(timestep_values[i])
        data = dsa.WrapDataObject(servermanager.Fetch(reader))
        arr = np.asarray(data.PointData["VaporDensity"]).ravel()
        pooled.append(arr)
    pooled = np.concatenate(pooled)
    lo, hi = np.percentile(pooled, [lo_pct, hi_pct])
    return float(lo), float(hi)


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dir", default=".", help="run directory (default: .)")
    p.add_argument("--pvd", default=None,
                   help="override path to the .pvd file (default: auto-detect "
                        "permafrost_highres.pvd, falling back to permafrost.pvd)")
    p.add_argument("--out", default=None,
                   help="output movie path (default: <dir>/movie.mp4)")
    p.add_argument("--n-frames", type=int, default=600,
                   help="number of frames, evenly spaced in simulated time "
                        "(default: 600)")
    p.add_argument("--fps", type=int, default=30, help="playback frame rate (default: 30)")
    p.add_argument("--t-start", type=float, default=None,
                   help="start of the simulated-time window (default: first timestep)")
    p.add_argument("--t-end", type=float, default=None,
                   help="end of the simulated-time window (default: last timestep)")
    p.add_argument("--no-interpolate", action="store_true",
                   help="disable TemporalInterpolator (frames hold the nearest "
                        "snapshot instead of blending between sparse late ones)")
    p.add_argument("--ice-colormap", default=None,
                   help="path to a ParaView JSON preset for the ice volume "
                        "(default: postprocess/colormaps/cmocean_ice.json, "
                        "generated by make_cmocean_preset.py)")
    p.add_argument("--vapor-colormap", default="Viridis (matplotlib)",
                   help="ParaView preset name (or path to a JSON preset) for "
                        "the air/vapor-density volume (default: 'Viridis (matplotlib)')")
    p.add_argument("--vapor-range", type=float, nargs=2, default=None,
                   help="fixed [min, max] for the vapor-density colorbar; if "
                        "omitted, auto-computed from --vapor-percentile over "
                        "--vapor-range-samples sampled timesteps")
    p.add_argument("--vapor-percentile", type=float, nargs=2, default=[2.0, 98.0],
                   help="percentile clipping used for auto vapor range (default: 2 98)")
    p.add_argument("--vapor-range-samples", type=int, default=40,
                   help="number of timesteps sampled for auto vapor range (default: 40)")
    p.add_argument("--resolution", default=None,
                   help="WxH output resolution (default: auto from data aspect ratio, "
                        "capped at 1600 wide)")
    p.add_argument("--keep-frames", action="store_true",
                   help="keep the rendered PNG frame sequence after muxing")
    args = p.parse_args()

    pvd_path = args.pvd or find_pvd(args.dir)
    out_path = args.out or os.path.join(args.dir, "movie.mp4")
    ice_preset_path = args.ice_colormap or os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "colormaps", "cmocean_ice.json")

    print(f"Reading {pvd_path}")
    reader = OpenDataFile(pvd_path)
    timestep_values = list(reader.TimestepValues)
    if not timestep_values:
        raise RuntimeError("No timesteps found in the PVD file")

    t_start = args.t_start if args.t_start is not None else timestep_values[0]
    t_end = args.t_end if args.t_end is not None else timestep_values[-1]
    print(f"Simulated time window: [{t_start:.6g}, {t_end:.6g}] "
          f"({len(timestep_values)} source snapshots, {args.n_frames} output frames)")

    # ---- vapor colorbar range (fixed across the whole movie) --------------
    if args.vapor_range is not None:
        vmin, vmax = args.vapor_range
    else:
        vmin, vmax = auto_vapor_range(
            reader, timestep_values, args.vapor_range_samples, *args.vapor_percentile)
        print(f"Auto vapor-density range ({args.vapor_percentile[0]:.0f}-"
              f"{args.vapor_percentile[1]:.0f} percentile over "
              f"{args.vapor_range_samples} samples): [{vmin:.4g}, {vmax:.4g}]")

    # ---- pipeline -----------------------------------------------------------
    source = reader
    if not args.no_interpolate:
        source = TemporalInterpolator(Input=reader)

    ice_volume = IsoVolume(Input=source, InputScalars=["POINTS", "IcePhase"],
                            ThresholdRange=[0.5, 1.1])
    air_volume = IsoVolume(Input=source, InputScalars=["POINTS", "IcePhase"],
                            ThresholdRange=[-0.1, 0.5])

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
    ice_display.SetScalarBarVisibility(view, True)

    air_display = Show(air_volume, view)
    ColorBy(air_display, ("POINTS", "VaporDensity"))
    vapor_lut = GetColorTransferFunction("VaporDensity")
    if os.path.isfile(args.vapor_colormap):
        preset_name = json.load(open(args.vapor_colormap))[0]["Name"]
        if not presets.HasPreset(preset_name):
            presets.ImportPresets(args.vapor_colormap)
        vapor_lut.ApplyPreset(preset_name, True)
    else:
        vapor_lut.ApplyPreset(args.vapor_colormap, True)
    vapor_lut.RescaleTransferFunction(vmin, vmax)
    air_display.SetScalarBarVisibility(view, True)

    ResetCamera(view)

    if args.resolution:
        w, h = (int(v) for v in args.resolution.split("x"))
    else:
        b = reader.GetDataInformation().GetBounds()
        lx, ly = b[1] - b[0], b[3] - b[2]
        w = 1600
        h = max(200, int(round(w * ly / lx)))
    view.ViewSize = [w, h]
    print(f"Render resolution: {w}x{h}")

    # ---- frame loop: evenly spaced in simulated time, not in snapshot index --
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

    # ---- mux with ffmpeg -----------------------------------------------------
    print(f"Encoding {out_path} at {args.fps} fps")
    cmd = [
        "ffmpeg", "-y", "-framerate", str(args.fps),
        "-i", os.path.join(frame_dir, "frame_%05d.png"),
        "-c:v", "libx264", "-pix_fmt", "yuv420p", out_path,
    ]
    subprocess.run(cmd, check=True)

    if not args.keep_frames:
        shutil.rmtree(frame_dir)

    print(f"Done: {out_path}")


if __name__ == "__main__":
    main()

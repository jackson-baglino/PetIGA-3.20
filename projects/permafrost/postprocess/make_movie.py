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

STREAMING (default): each frame's dense, true-NURBS .vts is generated on the
fly (by shelling out to the igakit venv running plot_permafrost_highres.py),
rendered, and then DELETED -- so at most one dense file is on disk at a time.
This avoids materializing the entire vtkOut_highres/ series, which at
--n-per-elem 4 can be hundreds of GB for a long run. The generator interpreter
is autodetected (venv_pf311 / venv_DSM); override with --gen-python.

    /path/to/pvpython postprocess/make_movie.py --dir <rundir>

Streaming frames are the actual sol_*.dat snapshots (optionally down-sampled to
--n-frames), not points evenly spaced in simulated time, and --interpolate is
unavailable (temporal blending needs the whole series loaded at once).

--no-stream uses a PRE-BUILT permafrost_highres.pvd instead (all dense .vts
generated up front) -- needed only for strict linear-time playback. Generate it
first with the venv Python (not pvpython -- it needs igakit):

    python3 postprocess/plot_permafrost_highres.py --dir <rundir>

find_pvd() (the --no-stream path) requires that file and refuses to silently
fall back to the coarse control-point permafrost.pvd (the raw B-spline
control-point grid, not the actual field shape -- faceted/blocky).

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
Resolution, cropping, and colorbars
------------------------------------------------------------------------
Frames render at the underlying data's own point resolution by default
(read straight from the reader -- whatever density plot_permafrost_highres.py
sampled at, e.g. ~2432x488 for a 608x122-element mesh at n_per_elem=4), not
some fixed/capped video size. Use --supersample for extra antialiasing or
--resolution WxH to override outright.

Most of these domains are wide and short (Lx >> Ly), so native height is
often well under 1080 (e.g. 488px in the example above) even when native
width is already huge -- "render at native resolution" alone does not
imply HD. --min-height (default 1080) scales the whole frame UP
proportionally whenever native/supersampled height would fall short, so
movies are never sub-HD; --max-width is still a downscale safety cap on
top of that, but yields to --min-height if the two conflict (a printed
note explains when this happens). Pass --min-height 0 to disable and get
the old "purely native/supersampled, no HD floor" behavior back.

The camera is set explicitly to the data's bounding box with zero margin
(not ResetCamera's default padding), so frames are cropped tight to the
geometry -- no background border. In-frame scalar bars are hidden (they
were eating into that crop and don't antialias/scale well baked into video
anyway); use --no-colorbars to skip exporting the standalone ones below.

Each run also exports two standalone vector colorbars next to the movie
(<out>_ice_colorbar.svg, <out>_vapor_colorbar.svg) reconstructed directly
from the actual ParaView transfer functions used for rendering (so they're
guaranteed to match), as SVG for easy resizing/relabeling in Inkscape.

------------------------------------------------------------------------
Usage examples
------------------------------------------------------------------------
  # Full run, 600 frames linear in simulated time, 30 fps, auto vapor range,
  # native resolution, separate vector colorbars
  pvpython postprocess/make_movie.py --dir /path/to/run

  # Fixed vapor colorbar bounds instead of auto-detecting
  pvpython postprocess/make_movie.py --dir /path/to/run --vapor-range 1e-4 5e-4

  # Only the first 2 days, denser sampling, no temporal interpolation
  pvpython postprocess/make_movie.py --dir /path/to/run --t-end 172800 \\
      --n-frames 900 --no-interpolate

  # 2x supersampled frames for a crisper video
  pvpython postprocess/make_movie.py --dir /path/to/run --supersample 2
"""

import argparse
import glob
import json
import os
import re
import shutil
import subprocess
import sys

from paraview.simple import (
    OpenDataFile, XMLStructuredGridReader, IsoVolume, TemporalInterpolator,
    Show, ColorBy, GetColorTransferFunction, GetActiveViewOrCreate,
    GetAnimationScene, Render, SaveScreenshot,
)
from paraview import servermanager
from vtk.numpy_interface import dataset_adapter as dsa
import numpy as np


def save_vector_colorbar(lut, label, out_path):
    """Export a standalone vertical colorbar as SVG, reconstructed from the
    actual ParaView transfer function (lut.RGBPoints) so it's guaranteed to
    match what's rendered in the frames -- regardless of whether the LUT
    came from a named built-in preset or a custom JSON import."""
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
    """Require the dense, true-NURBS-interpolated permafrost_highres.pvd
    (written by plot_permafrost_highres.py) -- never silently fall back to
    the coarse control-point permafrost.pvd. The coarse mesh is the raw
    B-spline control-point grid, not the actual field shape; rendering it
    directly looks faceted/blocky and is not the same data plot_
    permafrost_highres.py's dense .vts files show. If a movie ever needs to
    render that coarse mesh on purpose, pass --pvd explicitly."""
    highres_path = os.path.join(run_dir, "permafrost_highres.pvd")
    if os.path.isfile(highres_path):
        return highres_path
    coarse_path = os.path.join(run_dir, "permafrost.pvd")
    if os.path.isfile(coarse_path):
        raise FileNotFoundError(
            f"No permafrost_highres.pvd in {run_dir} (found only the coarse "
            f"control-point permafrost.pvd). Generate the dense, true-NURBS "
            f"interpolated version first, with this project's regular "
            f"Python (not pvpython):\n"
            f"    python3 postprocess/plot_permafrost_highres.py --dir {run_dir}\n"
            f"then re-run make_movie.py. Pass --pvd {coarse_path} explicitly "
            f"if you really want the coarse mesh instead.")
    raise FileNotFoundError(
        f"No permafrost_highres.pvd or permafrost.pvd found in {run_dir}")


def auto_vapor_range(air_volume, timestep_values, n_samples, lo_pct, hi_pct):
    """Sample n_samples evenly-spaced timesteps, pool VaporDensity across
    them, and return a global [lo_pct, hi_pct] percentile range -- a fixed
    colorbar range that stays meaningful across the whole movie instead of
    auto-rescaling (and thus changing meaning) every frame.

    Samples from air_volume (the IsoVolume already clipped to the air
    region, IcePhase in [-0.1, 0.5]) rather than the raw reader. VaporDensity
    inside the ICE region sits near the (much higher, more uniform)
    saturation density rho_vs -- pooling that in would dominate the
    percentile range and wash out the actual variation in the air region,
    which is the only place this colormap is ever rendered."""
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


def find_gen_python(explicit):
    """Locate the venv Python that has igakit (for the high-res generator).
    pvpython itself does NOT have igakit, so generation must run in the venv."""
    if explicit:
        return explicit
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    for cand in (os.path.join(root, "venv_pf311", "bin", "python3"),
                 os.path.expanduser("~/venvs/venv_DSM/bin/python3"),
                 os.path.join(root, "venv_permafrost", "bin", "python3")):
        if os.path.isfile(cand):
            return cand
    return "python3"


def sol_steps(run_dir):
    """Step indices from sol_*.dat, dropping the stray sol_<jobid>.dat outlier."""
    steps = []
    for f in glob.glob(os.path.join(run_dir, "sol_*.dat")):
        m = re.search(r"sol_(\d+)\.dat", os.path.basename(f))
        if m:
            s = int(m.group(1))
            if s < 10_000_000:            # exclude sol_<jobid>.dat (~65e6)
                steps.append(s)
    return sorted(set(steps))


def gen_highres(run_dir, step, gen_python, n_per_elem):
    """Generate ONE dense true-NURBS .vts for `step` via plot_permafrost_highres.py
    (run in the igakit venv). Returns the .vts path. Caller deletes it."""
    script = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "plot_permafrost_highres.py")
    out = os.path.join(run_dir, "vtkOut_highres", f"solV_{step:05d}.vts")
    r = subprocess.run(
        [gen_python, script, "--dir", run_dir, "--steps", str(step),
         "--n-per-elem", str(n_per_elem), "--force"],
        capture_output=True, text=True)
    if r.returncode != 0 or not os.path.isfile(out) or os.path.getsize(out) == 0:
        raise RuntimeError(
            f"high-res generation failed for step {step}:\n{r.stdout[-500:]}\n{r.stderr[-1500:]}")
    return out


def run_stream(args):
    """Streaming movie: for each frame, generate that step's dense .vts, render
    it, then DELETE the .vts -- so at most one dense file is on disk at a time
    (avoids the ~500 GB of a full pre-generated vtkOut_highres/).

    Differs from the pre-built-PVD path (--no-stream): frames are the actual
    sol_*.dat snapshots (optionally down-sampled to --n-frames), NOT points
    evenly spaced in simulated time, and --interpolate is unavailable (temporal
    interpolation needs the whole series loaded at once). Use --no-stream with a
    pre-built permafrost_highres.pvd if you need strict linear-time playback.
    """
    run_dir = args.dir
    gen_python = find_gen_python(args.gen_python)
    out_path = args.out or os.path.join(run_dir, "movie.mp4")
    ice_preset_path = args.ice_colormap or os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "colormaps", "cmocean_ice.json")
    vapor_colormap = args.vapor_colormap or os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "colormaps", "cmocean_balance.json")

    steps = sol_steps(run_dir)
    if not steps:
        raise RuntimeError(f"no sol_*.dat snapshots in {run_dir}")
    # frame selection: all snapshots, or evenly sampled down to --n-frames
    if args.n_frames and args.n_frames < len(steps):
        idx = [round(i * (len(steps) - 1) / (args.n_frames - 1))
               for i in range(args.n_frames)]
        frame_steps = [steps[j] for j in sorted(set(idx))]
    else:
        frame_steps = steps
    print(f"Streaming {len(frame_steps)} frames from {len(steps)} snapshots "
          f"(gen: {gen_python}, n-per-elem {args.n_per_elem})")

    outdir = os.path.join(run_dir, "vtkOut_highres")
    os.makedirs(outdir, exist_ok=True)

    # ---- build the pipeline ONCE on a re-pointable single-file reader --------
    first_vts = gen_highres(run_dir, frame_steps[0], gen_python, args.n_per_elem)
    reader = XMLStructuredGridReader(FileName=[first_vts])
    reader.UpdatePipeline()

    ice_volume = IsoVolume(Input=reader, InputScalars=["POINTS", "IcePhase"],
                            ThresholdRange=[0.5, 1.1])
    air_volume = IsoVolume(Input=reader, InputScalars=["POINTS", "IcePhase"],
                            ThresholdRange=[-0.1, 0.5])

    # ---- vapor range: sample a subset of steps (gen -> read -> delete) -------
    if args.vapor_range is not None:
        vmin, vmax = args.vapor_range
    else:
        ns = min(args.vapor_range_samples, len(frame_steps))
        samp = [frame_steps[round(i * (len(frame_steps) - 1) / max(1, ns - 1))]
                for i in range(ns)]
        import numpy as np
        vals = []
        for s in sorted(set(samp)):
            vts = first_vts if s == frame_steps[0] else gen_highres(
                run_dir, s, gen_python, args.n_per_elem)
            reader.FileName = [vts]; reader.UpdatePipeline()
            air_volume.UpdatePipeline()
            d = dsa.WrapDataObject(servermanager.Fetch(air_volume))
            arr = d.PointData.GetArray("VaporDensity")
            if arr is not None and len(arr):
                vals.append(np.asarray(arr))
            if s != frame_steps[0] and os.path.isfile(vts):
                os.remove(vts)
        allv = np.concatenate(vals) if vals else np.array([0.0, 1.0])
        vmin, vmax = (float(np.percentile(allv, args.vapor_percentile[0])),
                      float(np.percentile(allv, args.vapor_percentile[1])))
        # restore the reader to the first frame's file for pipeline setup
        reader.FileName = [first_vts]; reader.UpdatePipeline()
        print(f"Auto vapor range ({args.vapor_percentile[0]:.0f}-"
              f"{args.vapor_percentile[1]:.0f}%, air region): [{vmin:.4g}, {vmax:.4g}]")

    view, w, h, regolith_bg = _setup_scene(
        reader, ice_volume, air_volume, args, vmin, vmax,
        ice_preset_path, vapor_colormap, out_path)

    # ---- frame loop: gen -> render -> delete --------------------------------
    frame_dir = out_path + "_frames"
    if os.path.isdir(frame_dir):
        shutil.rmtree(frame_dir)
    os.makedirs(frame_dir)

    n = len(frame_steps)
    for i, step in enumerate(frame_steps):
        vts = first_vts if step == frame_steps[0] else gen_highres(
            run_dir, step, gen_python, args.n_per_elem)
        reader.FileName = [vts]; reader.UpdatePipeline()
        Render(view)
        frame_path = os.path.join(frame_dir, f"frame_{i:05d}.png")
        if regolith_bg is not None:
            from PIL import Image
            SaveScreenshot(frame_path, view, ImageResolution=[w, h], TransparentBackground=1)
            fg = Image.open(frame_path).convert("RGBA")
            Image.alpha_composite(regolith_bg, fg).convert("RGB").save(frame_path)
        else:
            SaveScreenshot(frame_path, view, ImageResolution=[w, h])
        if os.path.isfile(vts):          # delete this step's dense file
            os.remove(vts)
        if i % 25 == 0 or i == n - 1:
            print(f"  frame {i+1}/{n}  step {step}", flush=True)

    # tidy the (now-empty) high-res dir and the stale pvd the generator wrote
    for junk in (os.path.join(run_dir, "permafrost_highres.pvd"),):
        if os.path.isfile(junk):
            os.remove(junk)
    try:
        os.rmdir(outdir)
    except OSError:
        pass

    print(f"Encoding {out_path} at {args.fps} fps")
    subprocess.run(["ffmpeg", "-y", "-framerate", str(args.fps),
                    "-i", os.path.join(frame_dir, "frame_%05d.png"),
                    "-c:v", "libx264", "-pix_fmt", "yuv420p", out_path], check=True)
    if not args.keep_frames:
        shutil.rmtree(frame_dir)
    print(f"Done: {out_path}")


def _setup_scene(reader, ice_volume, air_volume, args, vmin, vmax,
                 ice_preset_path, vapor_colormap, out_path):
    """Shared display/colormap/camera/resolution setup. Returns (view,w,h,bg).
    Mirrors the inline setup in main() so both paths render identically."""
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
        w = int(round(nx * args.supersample)); h = int(round(ny * args.supersample))
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
        save_vector_colorbar(vapor_lut, "Vapor density", base + "_vapor_colorbar.svg")

    regolith_bg = None
    if not args.no_sediment_texture:
        sediment_texture = args.sediment_texture or os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "textures", "lunar_regolith.png")
        if os.path.isfile(sediment_texture):
            from PIL import Image
            regolith_bg = Image.open(sediment_texture).convert("RGBA").resize((w, h))
    return view, w, h, regolith_bg


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
    p.add_argument("--vapor-colormap", default=None,
                   help="ParaView preset name (or path to a JSON preset) for "
                        "the air/vapor-density volume (default: "
                        "postprocess/colormaps/cmocean_balance.json, generated by "
                        "make_cmocean_preset.py)")
    p.add_argument("--vapor-range", type=float, nargs=2, default=None,
                   help="fixed [min, max] for the vapor-density colorbar; if "
                        "omitted, auto-computed from --vapor-percentile over "
                        "--vapor-range-samples sampled timesteps")
    p.add_argument("--vapor-percentile", type=float, nargs=2, default=[2.0, 98.0],
                   help="percentile clipping used for auto vapor range (default: 2 98)")
    p.add_argument("--vapor-range-samples", type=int, default=40,
                   help="number of timesteps sampled for auto vapor range (default: 40)")
    p.add_argument("--resolution", default=None,
                   help="WxH output resolution (default: native -- the reader's own "
                        "point grid resolution, e.g. matching plot_permafrost_highres.py's "
                        "dense sampling 1:1)")
    p.add_argument("--supersample", type=float, default=1.0,
                   help="multiply the native resolution by this factor (default: 1.0)")
    p.add_argument("--min-height", type=int, default=1080,
                   help="floor on output height in pixels (the '1080' in '1080p'); "
                        "native/supersampled resolution is scaled UP proportionally "
                        "if it would fall short, so the movie is never sub-HD even for "
                        "short, low-resolution domains. Set to 0 to disable. (default: 1080)")
    p.add_argument("--max-width", type=int, default=4096,
                   help="safety cap on output width in pixels; resolution is downscaled "
                        "proportionally if it would exceed this -- unless that would also "
                        "violate --min-height, in which case --min-height wins and "
                        "--max-width is raised to match, with a printed note (default: 4096)")
    p.add_argument("--no-colorbars", action="store_true",
                   help="skip exporting the standalone SVG colorbars")
    p.add_argument("--keep-frames", action="store_true",
                   help="keep the rendered PNG frame sequence after muxing")
    p.add_argument("--sediment-texture", default=None,
                   help="path to a background texture image composited behind "
                        "the ice/air rendering, showing through wherever the "
                        "domain excludes a sediment bump (default: "
                        "postprocess/textures/lunar_regolith.png, generated by "
                        "make_regolith_texture.py)")
    p.add_argument("--no-sediment-texture", action="store_true",
                   help="disable the sediment-region background texture "
                        "entirely (plain background instead)")
    p.add_argument("--no-stream", action="store_true",
                   help="use a PRE-BUILT permafrost_highres.pvd (all dense .vts "
                        "generated up front by plot_permafrost_highres.py -- can be "
                        "hundreds of GB). Default is streaming: generate each "
                        "frame's dense .vts on the fly, render it, delete it, so "
                        "at most one dense file exists at a time.")
    p.add_argument("--gen-python", default=None,
                   help="Python interpreter (with igakit) used to generate each "
                        "frame's dense .vts in streaming mode. Default: autodetect "
                        "venv_pf311 / venv_DSM. Ignored with --no-stream.")
    p.add_argument("--n-per-elem", type=int, default=4,
                   help="dense sample points per element per direction for the "
                        "on-the-fly high-res .vts (streaming mode; default 4)")
    args = p.parse_args()

    # Streaming is the default: generate/render/delete one dense .vts per frame.
    if not args.no_stream:
        run_stream(args)
        return

    pvd_path = args.pvd or find_pvd(args.dir)
    out_path = args.out or os.path.join(args.dir, "movie.mp4")
    ice_preset_path = args.ice_colormap or os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "colormaps", "cmocean_ice.json")
    vapor_colormap = args.vapor_colormap or os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "colormaps", "cmocean_balance.json")
    sediment_texture = None
    if not args.no_sediment_texture:
        sediment_texture = args.sediment_texture or os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "textures", "lunar_regolith.png")
        if not os.path.isfile(sediment_texture):
            print(f"  (sediment texture not found at {sediment_texture} -- "
                  f"skipping; run postprocess/make_regolith_texture.py to generate it)")
            sediment_texture = None

    print(f"Reading {pvd_path}")
    reader = OpenDataFile(pvd_path)
    reader.UpdatePipeline(0.0)
    timestep_values = list(reader.TimestepValues)
    if not timestep_values:
        raise RuntimeError("No timesteps found in the PVD file")

    t_start = args.t_start if args.t_start is not None else timestep_values[0]
    t_end = args.t_end if args.t_end is not None else timestep_values[-1]
    print(f"Simulated time window: [{t_start:.6g}, {t_end:.6g}] "
          f"({len(timestep_values)} source snapshots, {args.n_frames} output frames)")

    # ---- pipeline -----------------------------------------------------------
    source = reader
    if not args.no_interpolate:
        source = TemporalInterpolator(Input=reader)

    ice_volume = IsoVolume(Input=source, InputScalars=["POINTS", "IcePhase"],
                            ThresholdRange=[0.5, 1.1])
    air_volume = IsoVolume(Input=source, InputScalars=["POINTS", "IcePhase"],
                            ThresholdRange=[-0.1, 0.5])

    # ---- vapor colorbar range (fixed across the whole movie) --------------
    if args.vapor_range is not None:
        vmin, vmax = args.vapor_range
    else:
        vmin, vmax = auto_vapor_range(
            air_volume, timestep_values, args.vapor_range_samples, *args.vapor_percentile)
        print(f"Auto vapor-density range ({args.vapor_percentile[0]:.0f}-"
              f"{args.vapor_percentile[1]:.0f} percentile over "
              f"{args.vapor_range_samples} samples, air region only): "
              f"[{vmin:.4g}, {vmax:.4g}]")

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
    if os.path.isfile(vapor_colormap):
        preset_name = json.load(open(vapor_colormap))[0]["Name"]
        if not presets.HasPreset(preset_name):
            presets.ImportPresets(vapor_colormap)
        vapor_lut.ApplyPreset(preset_name, True)
    else:
        vapor_lut.ApplyPreset(vapor_colormap, True)
    vapor_lut.RescaleTransferFunction(vmin, vmax)

    # In-frame legends off: they eat into the tight crop below and don't
    # scale/antialias well baked into video. Standalone vector colorbars
    # (matching these exact LUTs) are exported separately further down.
    ice_display.SetScalarBarVisibility(view, False)
    air_display.SetScalarBarVisibility(view, False)

    # ---- resolution: native data point grid by default, not a fixed cap ----
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
        print(f"  (upscaled to meet --min-height {args.min_height}: now {w}x{h})")
        if w > args.max_width:
            print(f"  (--min-height requires width {w} > --max-width {args.max_width} "
                  f"-- raising --max-width to match; min-height wins)")
            args.max_width = w
    if w > args.max_width:
        scale = args.max_width / w
        w, h = args.max_width, max(1, int(round(h * scale)))
        print(f"  (downscaled to stay under --max-width {args.max_width})")
    w, h = w + (w % 2), h + (h % 2)  # libx264/yuv420p needs even dimensions
    view.ViewSize = [w, h]
    print(f"Render resolution: {w}x{h}")

    # ---- tight crop: explicit camera on the data bounds, zero margin -------
    # ParaView auto-resets the camera (fit-to-data, with padding) on the
    # FIRST render after new representations are shown -- render once now
    # to absorb that reset, then override the camera explicitly. Setting
    # the camera before this first Render() is silently undone.
    Render(view)
    b = reader.GetDataInformation().GetBounds()
    xmid, ymid = 0.5 * (b[0] + b[1]), 0.5 * (b[2] + b[3])
    cam = view.GetActiveCamera()
    cam.SetParallelProjection(1)
    cam.SetFocalPoint(xmid, ymid, 0.0)
    cam.SetPosition(xmid, ymid, 1.0)
    cam.SetViewUp(0.0, 1.0, 0.0)
    cam.SetParallelScale(0.5 * (b[3] - b[2]))

    # ---- standalone vector colorbars (match the LUTs exactly) --------------
    if not args.no_colorbars:
        base, _ = os.path.splitext(out_path)
        save_vector_colorbar(ice_lut, "Ice phase", base + "_ice_colorbar.svg")
        save_vector_colorbar(vapor_lut, "Vapor density", base + "_vapor_colorbar.svg")

    # ---- frame loop: evenly spaced in simulated time, not in snapshot index --
    frame_dir = out_path + "_frames"
    if os.path.isdir(frame_dir):
        shutil.rmtree(frame_dir)
    os.makedirs(frame_dir)

    scene = GetAnimationScene()
    scene.UpdateAnimationUsingDataTimeSteps()

    # ---- sediment background: render frames transparent where the domain
    # excludes a bump, then alpha-composite onto a static regolith-colored
    # texture -- ParaView's own background-texture/environment properties
    # don't render in this 2D parallel-projection view (tried BackgroundTexture
    # and EnvironmentalBGTexture; both no-op here), so this is done as a
    # per-frame post-composite instead.
    regolith_bg = None
    if sediment_texture:
        from PIL import Image
        regolith_bg = Image.open(sediment_texture).convert("RGBA").resize((w, h))
        print(f"  Sediment background texture: {sediment_texture}")

    n = args.n_frames
    for i in range(n):
        t = t_start + (t_end - t_start) * (i / (n - 1) if n > 1 else 0.0)
        scene.AnimationTime = t
        Render(view)
        frame_path = os.path.join(frame_dir, f"frame_{i:05d}.png")
        if regolith_bg is not None:
            SaveScreenshot(frame_path, view, ImageResolution=[w, h], TransparentBackground=1)
            fg = Image.open(frame_path).convert("RGBA")
            Image.alpha_composite(regolith_bg, fg).convert("RGB").save(frame_path)
        else:
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

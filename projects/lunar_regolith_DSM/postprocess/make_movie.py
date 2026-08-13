#!/usr/bin/env python3
"""
make_movie.py — render the ice/air two-phase movie from a PetIGA run's PVD
time series, the way it's normally built by hand in the ParaView GUI:

  - "ice volume": IsoVolume of IcePhase in [0.5, 1.1], colored by IcePhase
    with the cmocean "ice" colormap.
  - "air volume": IsoVolume of IcePhase in [-0.1, 0.5] (i.e. air = 1-ice),
    colored by supersaturation sigma = (rhov - rho_vs(T)) / rho_vs(T) with
    the diverging cmocean "balance" colormap, on a range forced symmetric
    about zero so white always marks the saturated state (--air-field
    VaporDensity restores the old raw-density coloring).

    sigma, not rhov, is the field that shows the physics: rhov is dominated
    by the strong, smooth temperature dependence of rho_vs itself, so under
    a gradient it mostly redraws the background saturation profile. sigma
    divides that out, leaving the departure from equilibrium -- which is the
    actual driving force for deposition (sigma > 0) and sublimation
    (sigma < 0). It is written into each .vts by plot_fields_highres.py.

Must be run with ParaView's own Python, not this project's venv:

    /Applications/ParaView-5.12.0.app/Contents/bin/pvpython postprocess/make_movie.py --dir <rundir>

(adjust the pvpython path for your ParaView install/platform).

Before first use, generate the ice colormap preset once with the regular
project Python (it has cmocean; pvpython does not):

    python3 postprocess/make_cmocean_preset.py ice

STREAMING (default): each frame's dense, true-NURBS .vts is generated on the
fly (by shelling out to the igakit venv running plot_fields_highres.py),
rendered, and then DELETED -- so at most one dense file is on disk at a time.
This avoids materializing the entire vtkOut_highres/ series, which at
--n-per-elem 4 can be hundreds of GB for a long run. The generator interpreter
is autodetected (venv_pf311 / venv_DSM); override with --gen-python.

    /path/to/pvpython postprocess/make_movie.py --dir <rundir>

Streaming frames are the actual sol_*.dat snapshots (optionally down-sampled to
--n-frames), not points evenly spaced in simulated time, and --interpolate is
unavailable (temporal blending needs the whole series loaded at once).

--no-stream uses a PRE-BUILT pf_highres.pvd instead (all dense .vts
generated up front) -- needed only for strict linear-time playback. Generate it
first with the venv Python (not pvpython -- it needs igakit):

    python3 postprocess/plot_fields_highres.py --dir <rundir>

find_pvd() (the --no-stream path) requires that file and refuses to silently
fall back to the coarse control-point pf.pvd (the raw B-spline
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
(read straight from the reader -- whatever density plot_fields_highres.py
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

Each run also exports standalone vector colorbars next to the movie,
reconstructed directly from the actual ParaView transfer functions used for
rendering (so they're guaranteed to match), as SVG for easy
resizing/relabeling in Inkscape. For each field, BOTH ink colours are
written -- identical bars differing only in the colour of the label, ticks
and frame:

    <out>_ice_colorbar_black.svg    <out>_ice_colorbar_white.svg
    <out>_sigma_colorbar_black.svg  <out>_sigma_colorbar_white.svg

The SVG background is transparent, so the ink is what carries the contrast:
black for a white page or light slide, white for a dark slide or when laid
over the movie itself. Writing both costs nothing and saves recolouring an
SVG by hand later. Re-running on a run whose frames are already cached
re-exports the colorbars without re-rendering anything, so a change to how
they're drawn reaches an existing movie cheaply.

------------------------------------------------------------------------
Frame cache and movie length
------------------------------------------------------------------------
Rendered PNGs are KEPT in <out>.mp4_frames/ and, in streaming mode, named
by their SOURCE STEP (frame_01234.png) rather than by position in the
video. Rendering is the expensive part -- each frame costs a dense .vts
generation -- so a rerun reuses every cached PNG it can and renders only
what is missing. Alongside them sits render_meta.json recording what the
cache depicts (field, colour range, resolution, colormaps, n_per_elem);
if the current invocation would draw a different picture, the cache is
discarded rather than silently mixed into one movie. --force-frames
re-renders unconditionally; --delete-frames cleans up afterwards.

--duration sets the movie length in seconds, as frame count = duration *
--fps (so playback stays at a normal frame rate and a shorter movie is
proportionally cheaper). Because the frame set is a subsample of the same
snapshot list, changing --duration on a rerun reuses whatever frames it
lands on and renders only the rest; re-encoding an unchanged frame set at
a new length is a pure re-mux with no ParaView work at all.

The first snapshot (step 0, the raw initial condition) is SKIPPED by
default -- it is not a solution, and its uniform initial vapor field is a
flat sigma that appears nowhere else in the run and drags the pooled
percentile range used to set the colorbar. The solver always writes step 1
as well (OutputMonitor in src/monitoring.c), so the movie still opens
essentially at t=0, just on real solved data. --include-first restores it.

------------------------------------------------------------------------
Usage examples
------------------------------------------------------------------------
  # Full run, 600 frames linear in simulated time, 30 fps, auto sigma range,
  # native resolution, separate vector colorbars
  pvpython postprocess/make_movie.py --dir /path/to/run

  # A 5-second movie (150 frames at 30 fps)
  pvpython postprocess/make_movie.py --dir /path/to/run --duration 5

  # Re-cut the frames already on disk to 10 seconds -- no re-rendering
  pvpython postprocess/make_movie.py --dir /path/to/run --duration 10

  # Fixed colorbar bounds instead of auto-detecting
  pvpython postprocess/make_movie.py --dir /path/to/run --vapor-range -0.02 0.02

  # Colour by raw vapor density, as before
  pvpython postprocess/make_movie.py --dir /path/to/run \\
      --air-field VaporDensity --no-symmetric-range

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
import time

from paraview.simple import (
    OpenDataFile, XMLStructuredGridReader, IsoVolume, TemporalInterpolator,
    Show, ColorBy, GetColorTransferFunction, GetActiveViewOrCreate,
    GetAnimationScene, Render, SaveScreenshot,
)
from paraview import servermanager
from vtk.numpy_interface import dataset_adapter as dsa
import numpy as np


COLORBAR_INKS = ("black", "white")


def save_vector_colorbar(lut, label, out_base):
    """Export standalone vertical colorbars as SVG, reconstructed from the
    actual ParaView transfer function (lut.RGBPoints) so they're guaranteed
    to match what's rendered in the frames -- regardless of whether the LUT
    came from a named built-in preset or a custom JSON import.

    TWO files are written per field, identical but for the ink colour:
    <out_base>_black.svg and <out_base>_white.svg. The SVG background is
    transparent, so the label, ticks and frame have to carry the contrast
    themselves: black reads on a white page or light slide, white on a dark
    slide or when laid over the movie itself. Which one a given figure needs
    is not knowable from here, and recolouring an SVG by hand afterwards is
    exactly the fiddling this script exists to avoid -- so write both, every
    time, and let the figure pick."""
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

    for ink in COLORBAR_INKS:
        fig, ax = plt.subplots(figsize=(0.5, 6))
        cb = ColorbarBase(ax, cmap=cmap, norm=Normalize(vmin=lo, vmax=hi),
                          orientation="vertical")
        cb.set_label(label, fontsize=14, fontweight="bold", color=ink)
        cb.ax.tick_params(labelsize=12, colors=ink)
        # The frame around the ramp is drawn in the default black; left alone
        # it would vanish against a dark slide in the white variant.
        cb.outline.set_edgecolor(ink)
        for tick_label in cb.ax.get_yticklabels():
            tick_label.set_fontweight("bold")
        fig.savefig(f"{out_base}_{ink}.svg", format="svg",
                    bbox_inches="tight", transparent=True)
        plt.close(fig)
    print(f"  wrote {os.path.basename(out_base)}_{{{','.join(COLORBAR_INKS)}}}.svg"
          f"  (range [{lo:.4g}, {hi:.4g}])")


def build_lut(array_name, colormap, vmin, vmax):
    """Fetch and configure the transfer function for `array_name`.

    Depends only on the preset and the value range -- no reader, no data --
    which is what lets the cached re-mux path re-export colorbars without
    rebuilding the render pipeline."""
    presets = servermanager.vtkSMTransferFunctionPresets.GetInstance()
    lut = GetColorTransferFunction(array_name)
    if os.path.isfile(colormap):
        preset_name = json.load(open(colormap))[0]["Name"]
        if not presets.HasPreset(preset_name):
            if not presets.ImportPresets(colormap):
                raise RuntimeError(f"Failed to import colormap preset: {colormap}")
        lut.ApplyPreset(preset_name, True)
    else:
        lut.ApplyPreset(colormap, True)
    lut.RescaleTransferFunction(vmin, vmax)
    return lut


def export_colorbars(args, out_path, ice_lut, air_lut):
    base, _ = os.path.splitext(out_path)
    save_vector_colorbar(ice_lut, "Ice phase", base + "_ice_colorbar")
    save_vector_colorbar(air_lut, air_field_label(args.air_field),
                         base + "_" + air_field_slug(args.air_field) + "_colorbar")


def find_pvd(run_dir: str) -> str:
    """Require the dense, true-NURBS-interpolated pf_highres.pvd
    (written by plot_fields_highres.py) -- never silently fall back to
    the coarse control-point pf.pvd. The coarse mesh is the raw
    B-spline control-point grid, not the actual field shape; rendering it
    directly looks faceted/blocky and is not the same data plot_
    pf_highres.py's dense .vts files show. If a movie ever needs to
    render that coarse mesh on purpose, pass --pvd explicitly."""
    highres_path = os.path.join(run_dir, "pf_highres.pvd")
    if os.path.isfile(highres_path):
        return highres_path
    coarse_path = os.path.join(run_dir, "pf.pvd")
    if os.path.isfile(coarse_path):
        raise FileNotFoundError(
            f"No pf_highres.pvd in {run_dir} (found only the coarse "
            f"control-point pf.pvd). Generate the dense, true-NURBS "
            f"interpolated version first, with this project's regular "
            f"Python (not pvpython):\n"
            f"    python3 postprocess/plot_fields_highres.py --dir {run_dir}\n"
            f"then re-run make_movie.py. Pass --pvd {coarse_path} explicitly "
            f"if you really want the coarse mesh instead.")
    raise FileNotFoundError(
        f"No pf_highres.pvd or pf.pvd found in {run_dir}")


def symmetrize(lo, hi):
    """Widen [lo, hi] to the symmetric range about zero that contains it.

    The air field's colormap (cmocean balance) is DIVERGING: its white
    midpoint only means anything if it sits at zero. For supersaturation
    that midpoint is exactly the saturated state (sigma = 0, the boundary
    between deposition and sublimation), so an asymmetric range would put
    white at some arbitrary non-equilibrium value and silently mislabel
    which regions are growing versus shrinking."""
    m = max(abs(lo), abs(hi))
    return -m, m


def auto_air_range(air_volume, timestep_values, n_samples, lo_pct, hi_pct,
                   field):
    """Sample n_samples evenly-spaced timesteps, pool `field` across them,
    and return a global [lo_pct, hi_pct] percentile range -- a fixed
    colorbar range that stays meaningful across the whole movie instead of
    auto-rescaling (and thus changing meaning) every frame.

    Samples from air_volume (the IsoVolume already clipped to the air
    region, IcePhase in [-0.1, 0.5]) rather than the raw reader. Inside the
    ICE region the vapor field sits pinned near saturation -- pooling that
    in would dominate the percentile range and wash out the actual
    variation in the air region, which is the only place this colormap is
    ever rendered."""
    n = len(timestep_values)
    idx = np.unique(np.linspace(0, n - 1, min(n_samples, n)).astype(int))
    pooled = []
    for i in idx:
        air_volume.UpdatePipeline(timestep_values[i])
        data = dsa.WrapDataObject(servermanager.Fetch(air_volume))
        arr = np.asarray(data.PointData[field]).ravel()
        if arr.size:
            pooled.append(arr)
    pooled = np.concatenate(pooled)
    lo, hi = np.percentile(pooled, [lo_pct, hi_pct])
    return float(lo), float(hi)


AIR_FIELD_LABELS = {
    "Supersaturation": "Supersaturation  $\\sigma$",
    "VaporDensity":    "Vapor density",
    "Temperature":     "Temperature",
}


def air_field_label(field):
    return AIR_FIELD_LABELS.get(field, field)


def air_field_slug(field):
    """Filename stem for this field's standalone colorbar."""
    return {"Supersaturation": "sigma", "VaporDensity": "vapor"}.get(
        field, field.lower())


def resolve_n_frames(args, n_available=None):
    """Reconcile --duration / --n-frames / --fps into a frame count.

    --duration sets the frame COUNT at the given --fps (rather than raising
    fps and keeping every frame), so playback stays at a normal frame rate
    and a shorter movie is proportionally cheaper to render. To instead keep
    every frame and just play them faster, set --fps directly.

    n_available caps the count at the number of distinct source snapshots;
    pass None where frames are synthesised at arbitrary times (the
    --no-stream, temporally-interpolated path) and no such cap applies."""
    if args.duration:
        n = max(2, int(round(args.duration * args.fps)))
        if n_available is not None and n > n_available:
            # Not enough snapshots to fill the requested length at --fps;
            # drop fps so the movie still lasts --duration seconds.
            args.fps = max(1, int(round(n_available / args.duration)))
            print(f"  (only {n_available} snapshots available; lowering fps to "
                  f"{args.fps} to still fill {args.duration:g}s)")
            n = n_available
        print(f"Target duration {args.duration:g}s at {args.fps} fps "
              f"-> {n} frames")
        return n
    return args.n_frames


def render_key(args, w, h, vmin, vmax):
    """Everything that changes what a rendered PNG looks like.

    Stored alongside the cached frames; a mismatch means the cache is of a
    DIFFERENT picture (other field, other colormap, other resolution) and
    must be discarded rather than silently re-muxed into a movie that mixes
    two rendering conventions. Frame COUNT and fps are deliberately absent:
    those change which cached frames are used and how fast they play, not
    what any one of them looks like."""
    return {
        "air_field": args.air_field,
        "range": [round(vmin, 12), round(vmax, 12)],
        "resolution": [w, h],
        "ice_colormap": os.path.basename(args.ice_colormap or "cmocean_ice.json"),
        "air_colormap": os.path.basename(args.vapor_colormap or "cmocean_balance.json"),
        "n_per_elem": args.n_per_elem,
        "sediment_texture": (not args.no_sediment_texture),
    }


def load_frame_meta(frame_dir):
    path = os.path.join(frame_dir, "render_meta.json")
    if not os.path.isfile(path):
        return None
    try:
        with open(path) as fh:
            return json.load(fh)
    except (ValueError, OSError):
        return None


def save_frame_meta(frame_dir, key):
    with open(os.path.join(frame_dir, "render_meta.json"), "w") as fh:
        json.dump(key, fh, indent=2, sort_keys=True)


PNG_IEND = b"IEND\xaeB`\x82"


def _await_png(path, timeout=60.0):
    """Block until `path` is a complete, decodable PNG.

    SaveScreenshot can return BEFORE ParaView has finished flushing the file
    (observed on ParaView 6.1.1 with ~10 MB frames): reading it immediately
    fails with UnidentifiedImageError on a half-written header. Wait for the
    terminating IEND chunk AND for a decode to actually succeed."""
    from PIL import Image
    deadline = time.time() + timeout
    while True:
        try:
            with open(path, "rb") as fh:
                fh.seek(-len(PNG_IEND), os.SEEK_END)
                if fh.read(len(PNG_IEND)) == PNG_IEND:
                    im = Image.open(path)
                    im.load()
                    return im
        except (OSError, ValueError):
            pass                       # truncated / not yet written
        if time.time() > deadline:
            raise RuntimeError(
                f"ParaView never finished writing {path} (waited {timeout:.0f}s)")
        time.sleep(0.05)


def save_frame(view, w, h, path, background):
    """Render one frame to `path`, compositing onto `background` if given.

    Written via a temporary file and renamed into place, so `path` either
    does not exist or is a finished frame. The frame cache decides what to
    reuse by filename alone, so a run interrupted mid-write must not leave
    behind something that *looks* cached -- a half-written or (worse) an
    uncomposited transparent frame would be silently muxed into the next
    movie."""
    tmp = path + ".part.png"
    if background is not None:
        from PIL import Image
        SaveScreenshot(tmp, view, ImageResolution=[w, h], TransparentBackground=1)
        fg = _await_png(tmp)
        Image.alpha_composite(background, fg.convert("RGBA")
                              ).convert("RGB").save(tmp)
    else:
        SaveScreenshot(tmp, view, ImageResolution=[w, h])
        _await_png(tmp).close()
    os.replace(tmp, path)


def encode(frame_paths, out_path, fps):
    """Mux an arbitrary ordered list of PNGs into out_path at `fps`.

    The frames are cached under their SOURCE STEP number (so a rerun at a
    different length can reuse them), which is not the gapless 0,1,2,...
    sequence ffmpeg's numbered-input pattern needs. Rather than renaming
    (and so destroying) the cache, stage a directory of symlinks in the
    sequential order this particular movie wants."""
    seq_dir = os.path.join(os.path.dirname(frame_paths[0]), "_seq")
    if os.path.isdir(seq_dir):
        shutil.rmtree(seq_dir)
    os.makedirs(seq_dir)
    for i, src in enumerate(frame_paths):
        os.symlink(os.path.abspath(src), os.path.join(seq_dir, f"f_{i:05d}.png"))
    print(f"Encoding {out_path}: {len(frame_paths)} frames at {fps} fps "
          f"({len(frame_paths) / fps:.2f}s)")
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error",
                    "-framerate", str(fps),
                    "-i", os.path.join(seq_dir, "f_%05d.png"),
                    "-c:v", "libx264", "-pix_fmt", "yuv420p", out_path],
                   check=True)
    shutil.rmtree(seq_dir)


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
    """Generate ONE dense true-NURBS .vts for `step` via plot_fields_highres.py
    (run in the igakit venv). Returns the .vts path. Caller deletes it."""
    script = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "plot_fields_highres.py")
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
    pre-built pf_highres.pvd if you need strict linear-time playback.
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

    # ---- drop the initial condition ----------------------------------------
    # Step 0 is the raw IC, not a solution: its vapor field is a uniform
    # hum0*rho_vs, so its supersaturation is a flat constant that appears
    # nowhere else in the run. Kept as frame 0 it both opens the movie on a
    # non-physical still and, worse, drags the pooled percentile range used
    # for the colorbar. The solver now always writes step 1 as well (see
    # OutputMonitor), so dropping step 0 still leaves a genuine first frame
    # essentially at t=0.
    if not args.include_first and len(steps) > 1:
        print(f"Skipping the initial condition (step {steps[0]}); "
              f"first frame is step {steps[1]}")
        steps = steps[1:]

    n_frames = resolve_n_frames(args, len(steps))
    # frame selection: all snapshots, or evenly sampled down to n_frames
    if n_frames and n_frames < len(steps):
        idx = [round(i * (len(steps) - 1) / (n_frames - 1))
               for i in range(n_frames)]
        frame_steps = [steps[j] for j in sorted(set(idx))]
    else:
        frame_steps = steps

    frame_dir = out_path + "_frames"
    os.makedirs(frame_dir, exist_ok=True)

    def frame_path(step):
        return os.path.join(frame_dir, f"frame_{step:05d}.png")

    # ---- fast path: every requested frame is already rendered ---------------
    # Reuse is keyed on render_key(); a cached meta that matches means those
    # PNGs are the same picture this invocation would draw, so a rerun that
    # only changes --duration/--fps/--n-frames is a pure re-mux -- no
    # ParaView, no .vts generation, no percentile sampling.
    cached = None if args.force_frames else load_frame_meta(frame_dir)
    if cached is not None and all(os.path.isfile(frame_path(s))
                                  for s in frame_steps):
        # The cache supplies resolution and colour range (both are otherwise
        # only known after the expensive setup), so rebuild the key around
        # the cached values -- but honour them only where this invocation did
        # not ASK for something different, or "--vapor-range 0 1" would be
        # silently ignored in favour of whatever the cache happened to hold.
        c_res, c_rng = cached.get("resolution", [0, 0]), cached.get("range", [0.0, 0.0])
        want = None
        if args.resolution:
            want = [int(v) for v in args.resolution.split("x")]
        req_rng = args.vapor_range
        if req_rng and args.symmetric_range:
            req_rng = list(symmetrize(*req_rng))
        if ((want is None or want == c_res)
                and (req_rng is None
                     or [round(v, 12) for v in req_rng] == c_rng)
                and render_key(args, *c_res, *c_rng) == cached):
            print(f"Reusing {len(frame_steps)} cached frames in {frame_dir} "
                  f"(air field {cached['air_field']}, range "
                  f"[{c_rng[0]:.4g}, {c_rng[1]:.4g}])")
            # Colorbars are re-exported even here. They depend only on the
            # presets and the cached range -- no data, no pipeline -- so this
            # costs nothing, and it means a change to how they are DRAWN
            # (ink colour, labels) reaches an existing run without paying for
            # a full re-render of frames that would come out identical.
            if not args.no_colorbars:
                export_colorbars(
                    args, out_path,
                    build_lut("IcePhase", ice_preset_path, 0.0, 1.0),
                    build_lut(args.air_field, vapor_colormap, *c_rng))
            encode([frame_path(s) for s in frame_steps], out_path, args.fps)
            print(f"Done: {out_path}")
            return
        print("Cached frames were rendered with different settings "
              "-- re-rendering.")

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
        vals = []
        for s in sorted(set(samp)):
            vts = first_vts if s == frame_steps[0] else gen_highres(
                run_dir, s, gen_python, args.n_per_elem)
            reader.FileName = [vts]; reader.UpdatePipeline()
            air_volume.UpdatePipeline()
            d = dsa.WrapDataObject(servermanager.Fetch(air_volume))
            arr = d.PointData.GetArray(args.air_field)
            if arr is not None and len(arr):
                vals.append(np.asarray(arr))
            if s != frame_steps[0] and os.path.isfile(vts):
                os.remove(vts)
        if not vals:
            raise RuntimeError(
                f"field '{args.air_field}' not present in the generated .vts "
                f"(available via plot_fields_highres.py: IcePhase, Temperature, "
                f"VaporDensity, Supersaturation). Pass --air-field to pick "
                f"another, or regenerate with an updated plot_fields_highres.py.")
        allv = np.concatenate(vals)
        vmin, vmax = (float(np.percentile(allv, args.vapor_percentile[0])),
                      float(np.percentile(allv, args.vapor_percentile[1])))
        if args.symmetric_range:
            vmin, vmax = symmetrize(vmin, vmax)
        # restore the reader to the first frame's file for pipeline setup
        reader.FileName = [first_vts]; reader.UpdatePipeline()
        print(f"Auto {args.air_field} range ({args.vapor_percentile[0]:.0f}-"
              f"{args.vapor_percentile[1]:.0f}%, air region): "
              f"[{vmin:.4g}, {vmax:.4g}]")

    view, w, h, regolith_bg = _setup_scene(
        reader, ice_volume, air_volume, args, vmin, vmax,
        ice_preset_path, vapor_colormap, out_path)

    # ---- frame loop: gen -> render -> delete --------------------------------
    # Frames are cached under their SOURCE STEP, and any left over from a
    # previous run with identical settings are kept: a rerun that only widens
    # or re-samples the frame set pays only for the steps it does not have.
    key = render_key(args, w, h, vmin, vmax)
    stale = load_frame_meta(frame_dir) != key
    if stale or args.force_frames:
        for old in glob.glob(os.path.join(frame_dir, "frame_*.png")):
            os.remove(old)
    save_frame_meta(frame_dir, key)

    todo = [s for s in frame_steps if not os.path.isfile(frame_path(s))]
    if len(todo) < len(frame_steps):
        print(f"  reusing {len(frame_steps) - len(todo)} cached frames; "
              f"rendering {len(todo)}")

    n = len(todo)
    for i, step in enumerate(todo):
        vts = first_vts if step == frame_steps[0] else gen_highres(
            run_dir, step, gen_python, args.n_per_elem)
        reader.FileName = [vts]; reader.UpdatePipeline()
        Render(view)
        save_frame(view, w, h, frame_path(step), regolith_bg)
        if os.path.isfile(vts):          # delete this step's dense file
            os.remove(vts)
        if i % 25 == 0 or i == n - 1:
            print(f"  frame {i+1}/{n}  step {step}", flush=True)

    # tidy the (now-empty) high-res dir and the stale pvd the generator wrote.
    # first_vts is only deleted by the loop if its step was actually rendered;
    # when frame_steps[0] came from the cache it is still sitting there.
    for junk in (os.path.join(run_dir, "pf_highres.pvd"), first_vts):
        if os.path.isfile(junk):
            os.remove(junk)
    try:
        os.rmdir(outdir)
    except OSError:
        pass

    encode([frame_path(s) for s in frame_steps], out_path, args.fps)
    if args.delete_frames:
        shutil.rmtree(frame_dir)
    print(f"Done: {out_path}")


def _setup_scene(reader, ice_volume, air_volume, args, vmin, vmax,
                 ice_preset_path, vapor_colormap, out_path):
    """Shared display/colormap/camera/resolution setup. Returns (view,w,h,bg).
    The single implementation for BOTH the streaming and --no-stream paths."""
    view = GetActiveViewOrCreate("RenderView")
    view.InteractionMode = "2D"
    view.OrientationAxesVisibility = 0

    ice_display = Show(ice_volume, view)
    ColorBy(ice_display, ("POINTS", "IcePhase"))
    ice_lut = build_lut("IcePhase", ice_preset_path, 0.0, 1.0)

    air_display = Show(air_volume, view)
    ColorBy(air_display, ("POINTS", args.air_field))
    vapor_lut = build_lut(args.air_field, vapor_colormap, vmin, vmax)

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
        print(f"  (upscaled to meet --min-height {args.min_height}: now {w}x{h})")
        if w > args.max_width:
            print(f"  (--min-height requires width {w} > --max-width "
                  f"{args.max_width} -- raising --max-width to match)")
            args.max_width = w
    if w > args.max_width:
        scale = args.max_width / w
        w, h = args.max_width, max(1, int(round(h * scale)))
        print(f"  (downscaled to stay under --max-width {args.max_width})")
    w, h = w + (w % 2), h + (h % 2)  # libx264/yuv420p needs even dimensions
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
        export_colorbars(args, out_path, ice_lut, vapor_lut)

    # ---- sediment background: render frames transparent where the domain
    # excludes a bump, then alpha-composite onto a static regolith-colored
    # texture -- ParaView's own background-texture/environment properties
    # don't render in this 2D parallel-projection view (tried BackgroundTexture
    # and EnvironmentalBGTexture; both no-op here), so this is done as a
    # per-frame post-composite instead.
    regolith_bg = None
    if not args.no_sediment_texture:
        sediment_texture = args.sediment_texture or os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "textures", "lunar_regolith.png")
        if os.path.isfile(sediment_texture):
            from PIL import Image
            regolith_bg = Image.open(sediment_texture).convert("RGBA").resize((w, h))
            print(f"  Sediment background texture: {sediment_texture}")
        else:
            print(f"  (sediment texture not found at {sediment_texture} -- "
                  f"skipping; run postprocess/make_regolith_texture.py to "
                  f"generate it)")
    return view, w, h, regolith_bg


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dir", default=".", help="run directory (default: .)")
    p.add_argument("--pvd", default=None,
                   help="override path to the .pvd file (default: auto-detect "
                        "pf_highres.pvd, falling back to pf.pvd)")
    p.add_argument("--out", default=None,
                   help="output movie path (default: <dir>/movie.mp4)")
    p.add_argument("--n-frames", type=int, default=600,
                   help="number of frames, evenly spaced in simulated time "
                        "(default: 600). Ignored when --duration is given.")
    p.add_argument("--duration", type=float, default=None,
                   help="target movie length in SECONDS. Overrides --n-frames: "
                        "the frame count becomes duration*fps, so playback "
                        "stays at --fps and a shorter movie is proportionally "
                        "cheaper to render. (To instead keep every frame and "
                        "play them faster, set --fps and leave this off.) If "
                        "there are fewer snapshots than duration*fps, fps is "
                        "lowered instead so the movie still runs this long.")
    p.add_argument("--fps", type=int, default=30, help="playback frame rate (default: 30)")
    p.add_argument("--air-field", default="Supersaturation",
                   help="point field used to color the air region (default: "
                        "Supersaturation, sigma = (rhov - rho_vs)/rho_vs -- the "
                        "actual driving force for deposition/sublimation, and "
                        "signed about zero so the diverging colormap reads "
                        "correctly). Other options written by "
                        "plot_fields_highres.py: VaporDensity, Temperature.")
    p.add_argument("--no-symmetric-range", dest="symmetric_range",
                   action="store_false",
                   help="do not force the air-field color range to be symmetric "
                        "about zero. Symmetric is the default because the "
                        "colormap is diverging and its midpoint is only "
                        "meaningful at sigma = 0 (the saturated state).")
    p.set_defaults(symmetric_range=True)
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
                   help="fixed [min, max] for the air-field colorbar; if "
                        "omitted, auto-computed from --vapor-percentile over "
                        "--vapor-range-samples sampled timesteps")
    p.add_argument("--vapor-percentile", type=float, nargs=2, default=[2.0, 98.0],
                   help="percentile clipping used for auto vapor range (default: 2 98)")
    p.add_argument("--vapor-range-samples", type=int, default=40,
                   help="number of timesteps sampled for auto vapor range (default: 40)")
    p.add_argument("--resolution", default=None,
                   help="WxH output resolution (default: native -- the reader's own "
                        "point grid resolution, e.g. matching plot_fields_highres.py's "
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
    p.add_argument("--delete-frames", action="store_true",
                   help="delete the rendered PNG frame sequence after muxing. "
                        "Frames are KEPT by default (in <out>.mp4_frames/, "
                        "named by source step) so a rerun that only changes "
                        "--duration/--fps/--n-frames re-muxes them instead of "
                        "re-rendering -- which is the expensive part.")
    p.add_argument("--force-frames", action="store_true",
                   help="re-render every frame even if a matching cached PNG "
                        "already exists")
    p.add_argument("--include-first", action="store_true",
                   help="include the first snapshot (step 0, the raw initial "
                        "condition). Skipped by default: it is not a solution, "
                        "and its uniform vapor field skews the pooled "
                        "percentile range used for the colorbar.")
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
                   help="use a PRE-BUILT pf_highres.pvd (all dense .vts "
                        "generated up front by plot_fields_highres.py -- can be "
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

    print(f"Reading {pvd_path}")
    reader = OpenDataFile(pvd_path)
    reader.UpdatePipeline(0.0)
    timestep_values = list(reader.TimestepValues)
    if not timestep_values:
        raise RuntimeError("No timesteps found in the PVD file")

    if not args.include_first and len(timestep_values) > 1:
        # Start the window at the first SOLVED snapshot, not the raw IC --
        # see the same skip in run_stream() for why.
        print(f"Skipping the initial condition (t = {timestep_values[0]:.6g})")
        timestep_values = timestep_values[1:]

    t_start = args.t_start if args.t_start is not None else timestep_values[0]
    t_end = args.t_end if args.t_end is not None else timestep_values[-1]
    n_frames = resolve_n_frames(args)
    print(f"Simulated time window: [{t_start:.6g}, {t_end:.6g}] "
          f"({len(timestep_values)} source snapshots, {n_frames} output frames)")

    # ---- pipeline -----------------------------------------------------------
    source = reader
    if not args.no_interpolate:
        source = TemporalInterpolator(Input=reader)

    ice_volume = IsoVolume(Input=source, InputScalars=["POINTS", "IcePhase"],
                            ThresholdRange=[0.5, 1.1])
    air_volume = IsoVolume(Input=source, InputScalars=["POINTS", "IcePhase"],
                            ThresholdRange=[-0.1, 0.5])

    # ---- air-field colorbar range (fixed across the whole movie) ----------
    if args.vapor_range is not None:
        vmin, vmax = args.vapor_range
        if args.symmetric_range:
            vmin, vmax = symmetrize(vmin, vmax)
    else:
        vmin, vmax = auto_air_range(
            air_volume, timestep_values, args.vapor_range_samples,
            *args.vapor_percentile, field=args.air_field)
        if args.symmetric_range:
            vmin, vmax = symmetrize(vmin, vmax)
        print(f"Auto {args.air_field} range ({args.vapor_percentile[0]:.0f}-"
              f"{args.vapor_percentile[1]:.0f} percentile over "
              f"{args.vapor_range_samples} samples, air region only): "
              f"[{vmin:.4g}, {vmax:.4g}]")

    # Display/colormap/camera/resolution setup is shared with the streaming
    # path -- this used to be a second, hand-kept copy of it, which is exactly
    # how the two paths end up rendering subtly different movies.
    view, w, h, regolith_bg = _setup_scene(
        reader, ice_volume, air_volume, args, vmin, vmax,
        ice_preset_path, vapor_colormap, out_path)

    # ---- frame loop: evenly spaced in simulated time, not in snapshot index --
    # Frames here are keyed by output INDEX, not by source step: they are
    # synthesised at arbitrary interpolated times, so "which snapshot is this"
    # has no answer and a cache would only be reusable at an identical frame
    # count. Reuse therefore requires the same count as well as the same key.
    frame_dir = out_path + "_frames"
    os.makedirs(frame_dir, exist_ok=True)
    key = dict(render_key(args, w, h, vmin, vmax),
               mode="pvd", n_frames=n_frames,
               window=[round(t_start, 9), round(t_end, 9)],
               interpolate=(not args.no_interpolate))
    frames = [os.path.join(frame_dir, f"frame_{i:05d}.png") for i in range(n_frames)]

    if (not args.force_frames and load_frame_meta(frame_dir) == key
            and all(os.path.isfile(f) for f in frames)):
        print(f"Reusing {n_frames} cached frames in {frame_dir}")
        encode(frames, out_path, args.fps)
        if args.delete_frames:
            shutil.rmtree(frame_dir)
        print(f"Done: {out_path}")
        return

    for old in glob.glob(os.path.join(frame_dir, "frame_*.png")):
        os.remove(old)
    save_frame_meta(frame_dir, key)

    scene = GetAnimationScene()
    scene.UpdateAnimationUsingDataTimeSteps()

    n = n_frames
    for i in range(n):
        t = t_start + (t_end - t_start) * (i / (n - 1) if n > 1 else 0.0)
        scene.AnimationTime = t
        Render(view)
        save_frame(view, w, h, frames[i], regolith_bg)
        if i % 50 == 0 or i == n - 1:
            print(f"  frame {i+1}/{n}  t={t:.6g}")

    # ---- mux with ffmpeg -----------------------------------------------------
    encode(frames, out_path, args.fps)

    if args.delete_frames:
        shutil.rmtree(frame_dir)

    print(f"Done: {out_path}")


if __name__ == "__main__":
    main()

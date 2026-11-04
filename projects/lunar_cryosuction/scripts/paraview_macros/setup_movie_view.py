"""
setup_movie_view.py — ParaView macro: build the make_movie.py scene in the GUI,
framed tight to the data so File -> Save Animation writes frames with no
background border.

Why this exists
---------------
postprocess/make_movie.py renders every frame itself, and in streaming mode it
regenerates each frame's dense .vts by shelling out to igakit first -- which is
what makes it slow. ParaView's own Save Animation is far faster because the data
is already loaded and it renders in-process. What you lose going through the GUI
is make_movie.py's explicit camera: ResetCamera pads around the data, and if the
saved image's aspect ratio differs from the domain's you get background bars.

This macro reproduces the script's scene setup on whatever source is selected:

  - ice IsoVolume, IcePhase in [0.5, 1.1], coloured by IcePhase (cmocean ice,
    range 0..1)
  - air IsoVolume, IcePhase in [-0.1, 0.5], coloured by VaporDensity
    (cmocean balance, fixed range -- see VAPOR_RANGE below)
  - scalar bars and orientation axes off, 2D interaction
  - parallel-projection camera locked to the data bounding box with ZERO margin
  - the layout resized to the domain's exact aspect ratio

The two IsoVolumes tile the whole domain between them, so once the camera and
the aspect ratio both match the data, every pixel of the frame is data and the
background never shows.

Install
-------
ParaView -> Macros -> Add new macro... -> pick this file.

Use
---
1. Open the run's pf_highres.pvd (or a solV_*.vts series) and select it in the
   Pipeline Browser.
2. Click the macro. It prints the frame resolution it locked in.
3. File -> Save Animation..., keep the pre-filled resolution (or any resolution
   with the SAME aspect ratio), set the frame rate, save.

   Save Animation stretches the rendered view into whatever resolution it is
   given rather than re-rendering at it, so the crop stays tight either way --
   but a different aspect ratio comes out geometrically distorted.

Colormaps
---------
The cmocean presets live in <project>/postprocess/colormaps/*.json. ParaView
COPIES a macro into its own Macros directory, so this file usually cannot find
the project by relative path -- set PV_DSM_ROOT to the project root before
launching ParaView, or edit PROJECT_ROOT below. Once imported, the presets are
remembered in ParaView's settings and the JSON is no longer needed. Without
them the macro falls back to built-in presets and says so.

Tunables are the constants directly below.
"""

import os

import numpy as np

from paraview import servermanager
from paraview.simple import (
    ColorBy,
    CreateRenderView,
    GetActiveSource,
    GetActiveView,
    GetAnimationScene,
    GetColorTransferFunction,
    GetLayout,
    GetViews,
    Hide,
    IsoVolume,
    RenameSource,
    Render,
    SetActiveView,
    Show,
    TemporalInterpolator,
)
from vtk.numpy_interface import dataset_adapter as dsa

# --- scene ------------------------------------------------------------------
ICE_RANGE = [0.5, 1.1]       # IsoVolume window for the ice region
AIR_RANGE = [-0.1, 0.5]      # IsoVolume window for the air region

# Fixed [min, max] for the vapor-density colormap. None => auto-detect from
# VAPOR_PERCENTILE over VAPOR_RANGE_SAMPLES timesteps, which costs one data
# fetch per sample; set it explicitly (e.g. (1e-4, 5e-4)) to skip that entirely.
VAPOR_RANGE = None
VAPOR_PERCENTILE = (2.0, 98.0)
VAPOR_RANGE_SAMPLES = 5      # make_movie.py uses 40; fewer keeps the GUI snappy

# --- animation --------------------------------------------------------------
# Sequence mode samples N_FRAMES times evenly in SIMULATED time, so playback is
# linear even though PetIGA's adaptive dt makes snapshots dense early and sparse
# late. INTERPOLATE blends between sparse late snapshots instead of holding them.
N_FRAMES = 600
INTERPOLATE = True

# --- frame size -------------------------------------------------------------
# Native = the data's own point grid. The aspect ratio is what makes the crop
# tight; MIN_HEIGHT/MAX_WIDTH only scale it.
SUPERSAMPLE = 1.0
MIN_HEIGHT = 1080            # 0 disables the HD floor
MAX_WIDTH = 4096

# Shows only where the domain excludes geometry (e.g. a sediment bump); the two
# IsoVolumes cover everything else. make_movie.py composites a regolith texture
# there, which the GUI's Save Animation cannot do per frame.
BACKGROUND_RGB = (0.0, 0.0, 0.0)

# Project root holding postprocess/colormaps/. Env var wins; edit as a fallback.
PROJECT_ROOT = os.environ.get("PV_DSM_ROOT", "")


def _colormap_dir():
    """Locate <project>/postprocess/colormaps, or "" if it can't be found."""
    candidates = []
    if PROJECT_ROOT:
        candidates.append(os.path.join(PROJECT_ROOT, "postprocess", "colormaps"))
    # Works only when the macro is run from its place in the repo, not from the
    # copy ParaView keeps in its Macros directory.
    here = os.path.dirname(os.path.abspath(__file__))
    candidates.append(os.path.join(os.path.dirname(os.path.dirname(here)),
                                   "postprocess", "colormaps"))
    candidates.append(os.path.join(os.getcwd(), "postprocess", "colormaps"))
    for cand in candidates:
        if os.path.isdir(cand):
            return cand
    return ""


def _apply_preset(lut, preset_name, json_name, fallbacks):
    """Apply preset_name, importing its JSON if ParaView doesn't know it yet.

    Presets imported once persist in ParaView's settings, so the JSON is only
    needed the first time. Falls back to a built-in when neither is available.
    """
    presets = servermanager.vtkSMTransferFunctionPresets.GetInstance()
    if not presets.HasPreset(preset_name):
        cmap_dir = _colormap_dir()
        path = os.path.join(cmap_dir, json_name) if cmap_dir else ""
        if path and os.path.isfile(path):
            presets.ImportPresets(path)
    if presets.HasPreset(preset_name):
        lut.ApplyPreset(preset_name, True)
        return preset_name
    for name in fallbacks:
        try:
            lut.ApplyPreset(name, True)
            print("  (%s preset not found -- set PV_DSM_ROOT to the project "
                  "root to get it; using built-in '%s')" % (preset_name, name))
            return name
        except RuntimeError:
            continue
    return None


def _render_view():
    """A render view to build the scene in (never a chart/spreadsheet view)."""
    view = GetActiveView()
    if view is not None and view.GetXMLName() == "RenderView":
        return view
    for candidate in GetViews():
        if candidate.GetXMLName() == "RenderView":
            SetActiveView(candidate)
            return candidate
    view = CreateRenderView()
    SetActiveView(view)
    return view


def _auto_vapor_range(air_volume, timesteps):
    """Percentile range of VaporDensity pooled over sampled timesteps.

    Sampled from the AIR IsoVolume, not the raw source: inside the ice the
    vapor density sits near the much higher saturation value, and pooling that
    in would swamp the percentiles and wash out the air-region variation that
    is the only place this colormap is ever drawn.
    """
    if timesteps:
        n = min(VAPOR_RANGE_SAMPLES, len(timesteps))
        sample_times = [timesteps[round(i * (len(timesteps) - 1) / max(1, n - 1))]
                        for i in range(n)]
    else:
        sample_times = [None]

    pooled = []
    for t in sorted(set(sample_times)):
        if t is None:
            air_volume.UpdatePipeline()
        else:
            air_volume.UpdatePipeline(t)
        data = dsa.WrapDataObject(servermanager.Fetch(air_volume))
        arr = data.PointData.GetArray("VaporDensity")
        if arr is not None and len(arr):
            pooled.append(np.asarray(arr).ravel())
    if not pooled:
        raise RuntimeError(
            "setup_movie_view: no VaporDensity values in the air region -- "
            "set VAPOR_RANGE explicitly at the top of the macro.")
    pooled = np.concatenate(pooled)
    lo, hi = np.percentile(pooled, list(VAPOR_PERCENTILE))
    return float(lo), float(hi)


def _frame_size(src, bounds):
    """Frame resolution: native point grid, scaled to honour MIN_HEIGHT/MAX_WIDTH.

    Width is always derived from the DATA's aspect ratio, never from the point
    counts directly: a uniform nx-by-ny grid spans (nx-1) by (ny-1) cells, so
    nx/ny is not the domain's aspect (120/40 = 3.000 where the data is 3.051).
    Sizing off the point counts leaves a sliver of background down one edge.
    Point counts still set the SCALE, so frames stay at native resolution.
    """
    aspect = (bounds[1] - bounds[0]) / max(bounds[3] - bounds[2], 1e-30)

    ext = src.GetDataInformation().GetExtent()
    ny = ext[3] - ext[2] + 1
    if ny > 1:
        h = int(round(ny * SUPERSAMPLE))
    else:
        # Unstructured input has no point grid to inherit.
        h = int(MIN_HEIGHT or 1080)

    if MIN_HEIGHT and h < MIN_HEIGHT:
        h = MIN_HEIGHT
    if int(round(h * aspect)) > MAX_WIDTH:
        h = max(2, int(MAX_WIDTH / aspect))
    h += h % 2
    w = int(round(h * aspect))
    return w + (w % 2), h


def main():
    src = GetActiveSource()
    if src is None:
        raise RuntimeError(
            "setup_movie_view: no active source. Open the run's pf_highres.pvd "
            "(or a solV_*.vts series) and select it in the Pipeline Browser.")

    src.UpdatePipeline()
    timesteps = list(getattr(src, "TimestepValues", []) or [])

    view = _render_view()
    view.InteractionMode = "2D"
    view.OrientationAxesVisibility = 0
    view.UseColorPaletteForBackground = 0
    view.BackgroundColorMode = "Single Color"
    view.Background = list(BACKGROUND_RGB)

    # TemporalInterpolator blends frames that fall between two snapshots; with
    # Sequence play mode below, that is what makes late, sparse snapshots ease
    # into each other instead of visibly holding.
    pipeline_src = src
    if INTERPOLATE and len(timesteps) > 1:
        pipeline_src = TemporalInterpolator(Input=src)

    ice_volume = IsoVolume(Input=pipeline_src, InputScalars=["POINTS", "IcePhase"],
                           ThresholdRange=ICE_RANGE)
    air_volume = IsoVolume(Input=pipeline_src, InputScalars=["POINTS", "IcePhase"],
                           ThresholdRange=AIR_RANGE)
    RenameSource("ice volume", ice_volume)
    RenameSource("air volume", air_volume)

    if VAPOR_RANGE is not None:
        vmin, vmax = VAPOR_RANGE
    else:
        vmin, vmax = _auto_vapor_range(air_volume, timesteps)
        print("Auto vapor range (%.0f-%.0f%%, air region, %d samples): "
              "[%.4g, %.4g]" % (VAPOR_PERCENTILE[0], VAPOR_PERCENTILE[1],
                                VAPOR_RANGE_SAMPLES, vmin, vmax))

    Hide(src, view)
    ice_display = Show(ice_volume, view)
    ColorBy(ice_display, ("POINTS", "IcePhase"))
    ice_lut = GetColorTransferFunction("IcePhase")
    _apply_preset(ice_lut, "cmocean_ice", "cmocean_ice.json", ("Blues", "Cool to Warm"))
    ice_lut.RescaleTransferFunction(0.0, 1.0)

    air_display = Show(air_volume, view)
    ColorBy(air_display, ("POINTS", "VaporDensity"))
    vapor_lut = GetColorTransferFunction("VaporDensity")
    _apply_preset(vapor_lut, "cmocean_balance", "cmocean_balance.json",
                  ("Cool to Warm (Extended)", "Cool to Warm"))
    vapor_lut.RescaleTransferFunction(vmin, vmax)

    # In-frame legends would eat into the tight crop and don't scale well baked
    # into video; make_movie.py exports standalone SVG colorbars instead.
    ice_display.SetScalarBarVisibility(view, False)
    air_display.SetScalarBarVisibility(view, False)

    bounds = src.GetDataInformation().GetBounds()
    w, h = _frame_size(src, bounds)
    layout = GetLayout(view)
    if layout is not None:
        layout.SetSize(w, h)     # locks the GUI view to the data's aspect ratio
    else:
        view.ViewSize = [w, h]

    # ParaView auto-resets the camera on the FIRST render after new
    # representations appear, which re-adds its padding -- so render once to
    # absorb that reset, THEN override the camera. Doing it in the other order
    # is silently undone.
    Render(view)

    # Read the size BACK rather than trusting the request: a layout can refuse
    # to grow past the ParaView window. What the camera must match is the size
    # the view actually has, or the framing is wrong by exactly that error.
    actual = list(view.ViewSize)
    w = max(2, int(actual[0]))
    h = max(2, int(actual[1]))

    xmid = 0.5 * (bounds[0] + bounds[1])
    ymid = 0.5 * (bounds[2] + bounds[3])
    cam = view.GetActiveCamera()
    cam.SetParallelProjection(1)
    cam.SetFocalPoint(xmid, ymid, 0.0)
    cam.SetPosition(xmid, ymid, 1.0)
    cam.SetViewUp(0.0, 1.0, 0.0)
    # ParallelScale is the visible HALF-HEIGHT; the visible half-width follows
    # as scale*(w/h). Take whichever of the two constraints is tighter so the
    # data always over-fills the view: rounding to integer pixels leaves a
    # sub-pixel aspect error, and this makes that error crop a sliver of data
    # rather than expose a line of background.
    cam.SetParallelScale(min(0.5 * (bounds[3] - bounds[2]),
                             0.5 * (bounds[1] - bounds[0]) * h / w))
    Render(view)

    scene = GetAnimationScene()
    scene.UpdateAnimationUsingDataTimeSteps()
    if len(timesteps) > 1:
        scene.PlayMode = "Sequence"
        scene.StartTime = timesteps[0]
        scene.EndTime = timesteps[-1]
        scene.NumberOfFrames = N_FRAMES

    print("Frame size locked to %dx%d (data aspect %.4f)"
          % (w, h, (bounds[1] - bounds[0]) / max(bounds[3] - bounds[2], 1e-30)))
    if len(timesteps) > 1:
        print("Animation: %d frames, linear in simulated time over [%.6g, %.6g]%s"
              % (N_FRAMES, timesteps[0], timesteps[-1],
                 ", temporally interpolated" if INTERPOLATE else ""))
    # Save Animation STRETCHES the rendered view to whatever resolution it is
    # given (verified in ParaView 6.1: a sphere saved at 2:1 comes out an
    # ellipse) -- it does not re-render at that aspect. So the framing is
    # already tight whatever is typed there; what a mismatched aspect ratio
    # costs is distorted geometry, not a border.
    print("Now: File -> Save Animation... -- keep the pre-filled %dx%d, or any "
          "resolution with the SAME aspect ratio (%.4f). ParaView stretches the "
          "view into whatever size it is given, so a different aspect ratio "
          "comes out squashed or stretched." % (w, h, w / h))

    return view


main()

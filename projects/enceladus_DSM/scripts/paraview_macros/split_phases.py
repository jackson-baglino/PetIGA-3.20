"""
split_phases.py -- ParaView macro: split the active source into an ice
IsoVolume and an air IsoVolume on the phase field.

  ice  : ICE_RANGE = [0.5, 1.01]
  air  : AIR_RANGE = [-0.01, 0.5]

The two windows share the phi = 0.5 boundary and are padded past 0 and 1 so
that the overshoot the phase field picks up near a sharp interface still lands
inside one of them -- together they tile the whole domain with nothing dropped.

This is the isovolume part of setup_movie_view.py on its own, without the
camera, colormaps or animation setup, for when you just want the two phases as
separate pipeline objects to look at or run filters on.

Install
-------
ParaView -> Macros -> Add new macro... -> pick this file.
(ParaView copies the file into its own Macros directory, so re-add it after
editing this one.)

Use
---
Select a source with a phase array in the Pipeline Browser and click the macro.
Works on the raw .vts/.pvd reader or on anything downstream of it.

In 2D the two isovolumes tile the plane side by side, so both are shown. In 3D
the air region encloses the ice and would hide it, so air is created but left
invisible -- toggle its eye icon, or set SHOW_AIR_IN_3D below.
"""

from paraview.simple import (
    ColorBy,
    GetActiveSource,
    GetActiveViewOrCreate,
    GetDisplayProperties,
    Hide,
    IsoVolume,
    RenameSource,
    Render,
    SetActiveSource,
    Show,
)

# --- tunables ---------------------------------------------------------------
ICE_RANGE = [0.5, 1.01]
AIR_RANGE = [-0.01, 0.5]

# Phase array to threshold on. The first name present on the input wins; add
# yours to the front if you run a differently-named field.
PHASE_CANDIDATES = ["IcePhase", "phaseice", "phi", "Phase"]

ICE_COLOR = [0.75, 0.87, 0.95]   # pale ice blue
AIR_COLOR = [0.25, 0.25, 0.28]   # dark neutral, reads as "not ice"
AIR_OPACITY = 1.0

SHOW_AIR_IN_3D = False           # air encloses the ice and hides it


def _phase_array(source):
    """Name of the phase array on `source`, or raise with what was available."""
    point_data = source.GetPointDataInformation()
    present = [point_data.GetArray(i).GetName()
               for i in range(point_data.GetNumberOfArrays())]
    for name in PHASE_CANDIDATES:
        if name in present:
            return name
    raise RuntimeError(
        "split_phases: no phase array on the selected source. Looked for %s, "
        "found %s. Add the right name to PHASE_CANDIDATES at the top of the "
        "macro." % (PHASE_CANDIDATES, present or "no point arrays"))


def _is_3d(source):
    """True if the input has a non-degenerate extent in all three directions."""
    bounds = source.GetDataInformation().GetBounds()
    return all(bounds[2 * i + 1] > bounds[2 * i] for i in range(3))


def _isovolume(source, array, value_range, name, color, opacity, visible, view):
    volume = IsoVolume(Input=source, registrationName=name)
    volume.InputScalars = ["POINTS", array]
    volume.ThresholdRange = value_range
    RenameSource(name, volume)

    display = Show(volume, view)
    ColorBy(display, None)                    # solid colour, no scalar mapping
    display.AmbientColor = color
    display.DiffuseColor = color
    display.Opacity = opacity
    if not visible:
        Hide(volume, view)
    return volume


def main():
    source = GetActiveSource()
    if source is None:
        raise RuntimeError(
            "split_phases: nothing selected -- pick a source in the Pipeline "
            "Browser first.")

    array = _phase_array(source)
    view = GetActiveViewOrCreate("RenderView")
    three_d = _is_3d(source)

    ice = _isovolume(source, array, ICE_RANGE, "IcePhaseVolume",
                     ICE_COLOR, 1.0, True, view)
    air = _isovolume(source, array, AIR_RANGE, "AirPhaseVolume",
                     AIR_COLOR, AIR_OPACITY,
                     (not three_d) or SHOW_AIR_IN_3D, view)

    # The input itself would sit on top of both; the isovolumes replace it.
    Hide(source, view)
    GetDisplayProperties(ice, view)           # realise the reps before render
    GetDisplayProperties(air, view)
    SetActiveSource(ice)
    Render()

    print("split_phases: thresholded '%s' on %s" % (array, source.GetXMLLabel()))
    print("  IcePhaseVolume  %s in [%g, %g]" % (array, ICE_RANGE[0], ICE_RANGE[1]))
    print("  AirPhaseVolume  %s in [%g, %g]%s"
          % (array, AIR_RANGE[0], AIR_RANGE[1],
             "  (hidden -- 3D input, it would enclose the ice)"
             if three_d and not SHOW_AIR_IN_3D else ""))


main()

#!/usr/bin/env pvpython
"""
verify_curvature.py — analytic gate on the Gibbs–Thomson curvature that
plot_rhovsI.py computes.

Nothing in a solver output tells you whether a curvature field is right, so
this builds phase fields whose curvature is known in closed form, pushes them
through plot_rhovsI.py's OWN filter script (not a reimplementation), and checks
kappa at the phi = 0.5 contour against the exact value.

Cases
-----
  planar convex   ice disc, radius R                  kappa = +1/R
  planar concave  air bubble of radius R inside ice   kappa = -1/R
  axisym convex   ice sphere on the axis, radius R    kappa = +2/R
  axisym concave  air bubble on the axis, radius R    kappa = -2/R

The concave cases are the sign check: kappa = -div(grad phi/|grad phi|) with
grad phi pointing into the ice must come out POSITIVE on a convex grain, so
that small grains get the higher rho_vs^eff and vapor flows small -> large.
Flip that sign and the model runs Ostwald ripening backwards.

The axisym cases are the factor-of-2 check: the planar formula on an r–z grid
gives a sphere 1/R instead of 2/R.

Run
---
    /Applications/ParaView-6.1.1.app/Contents/bin/pvpython \
        scripts/paraview_macros/verify_curvature.py

Exits non-zero if any case misses its analytic value by more than TOL_PCT.
"""

import builtins
import os
import sys
import tempfile

import numpy as np
from vtkmodules.util import numpy_support as _numpy_support
from vtkmodules.vtkCommonCore import vtkPoints as _vtkPoints
from vtkmodules.vtkCommonDataModel import vtkStructuredGrid as _vtkStructuredGrid
from vtkmodules.vtkIOXML import (
    vtkXMLStructuredGridWriter as _vtkStructuredGridWriter,
)

# ParaView's Programmable Filter execs its script in __main__'s namespace, so
# when THIS file is __main__ the filter's preamble writes into these globals:
# it star-imports numpy_interface.algorithms (rebinding max/min/sum/round to
# array reductions taking `axis` second) and rebinds `vtk` to paraview's
# partial module. Either one breaks the SECOND case, after the first filter has
# run -- a confusing failure to chase. A leading underscore is the fix: star
# imports skip those names, so anything bound as _foo survives. Hence the
# _vtk* aliases above and the builtins below.
_max, _round, _abs = builtins.max, builtins.round, builtins.abs

HERE = os.path.dirname(os.path.abspath(__file__))
MACRO = os.path.join(HERE, "plot_rhovsI.py")

# Geometry of the synthetic fields. eps is a 20th of the radius and the mesh
# runs 3 points per eps, so the discrete gradients are well inside their
# asymptotic regime and the residual error is the scheme's, not the mesh's.
R_M = 50.0e-6            # feature radius [m]
EPS_M = R_M / 20.0       # phase-field decay length [m]
DX_M = EPS_M / 3.0       # grid spacing [m]
TEMP_C = -20.0

# Points within this much of phi = 0.5 are where kappa is compared.
CONTOUR_HALF_WIDTH = 0.05
TOL_PCT = 5.0

# Axisym only: how many cells off y = 0 get their own reported column.
AXIS_REPORT_CELLS = 2.0


def load_macro():
    """Exec plot_rhovsI.py without firing its main(), and hand back its namespace.

    The macro ends in a bare main() call so ParaView runs it on click; strip
    exactly that trailing call so the module-level definitions can be reused
    here. Importing the real thing is the point — a copy of the formula in this
    file would verify nothing.
    """
    with open(MACRO) as handle:
        source = handle.read()
    marker = "\nmain()\n"
    if not source.endswith(marker):
        raise SystemExit("verify_curvature: %s no longer ends in a bare main() "
                         "call; update load_macro()." % MACRO)
    namespace = {"__name__": "plot_rhovsI_undertest", "__file__": MACRO}
    exec(compile(source[: -len(marker)], MACRO, "exec"), namespace)
    return namespace


def build_grid(axisym, concave):
    """A vtkStructuredGrid carrying a logistic phase field around a circle.

    Planar: the disc sits in the middle of a square domain. Axisym: y is the
    radius and the domain starts at y = 0, so the circle centred on y = 0 is a
    sphere of revolution.
    """
    if axisym:
        x0, x1 = -2.0 * R_M, 2.0 * R_M
        y0, y1 = 0.0, 2.0 * R_M
        xc, yc = 0.0, 0.0
    else:
        x0, x1 = -2.0 * R_M, 2.0 * R_M
        y0, y1 = -2.0 * R_M, 2.0 * R_M
        xc, yc = 0.0, 0.0

    nx = int(_round((x1 - x0) / DX_M)) + 1
    ny = int(_round((y1 - y0) / DX_M)) + 1
    xs = np.linspace(x0, x1, nx)
    ys = np.linspace(y0, y1, ny)
    # VTK structured grids vary the first index fastest.
    gx, gy = np.meshgrid(xs, ys, indexing="ij")
    gx, gy = gx.ravel(order="F"), gy.ravel(order="F")

    rad = np.sqrt((gx - xc) ** 2 + (gy - yc) ** 2)
    # signed distance, positive inside the ice
    sdf = (rad - R_M) if concave else (R_M - rad)
    phi = 1.0 / (1.0 + np.exp(-sdf / EPS_M))

    points = _vtkPoints()
    coords = np.column_stack([gx, gy, np.zeros_like(gx)])
    points.SetData(_numpy_support.numpy_to_vtk(np.ascontiguousarray(coords),
                                               deep=1))

    grid = _vtkStructuredGrid()
    grid.SetDimensions(nx, ny, 1)
    grid.SetPoints(points)

    def add(name, values):
        arr = _numpy_support.numpy_to_vtk(np.ascontiguousarray(values), deep=1)
        arr.SetName(name)
        grid.GetPointData().AddArray(arr)

    add("IcePhase", phi)
    add("Temperature", np.full_like(phi, TEMP_C))
    add("VaporDensity", np.full_like(phi, 8.4868e-4))
    return grid, rad


def run_case(namespace, axisym, concave):
    """Apply the macro's filter to a synthetic case; return (measured, exact)."""
    from paraview.simple import (
        Delete, OpenDataFile, ProgrammableFilter, servermanager,
    )
    from vtk.numpy_interface import dataset_adapter as dsa

    grid, rad = build_grid(axisym, concave)

    handle, path = tempfile.mkstemp(suffix=".vts")
    os.close(handle)
    writer = _vtkStructuredGridWriter()
    writer.SetFileName(path)
    writer.SetInputData(grid)
    writer.Write()

    try:
        reader = OpenDataFile(path)
        reader.UpdatePipeline()

        pf = ProgrammableFilter(Input=reader)
        pf.OutputDataSetType = "Same as Input"
        # EPS is passed explicitly: a scratch .vts has no staged .opts beside
        # it, which is exactly the fallback path, but the gate should test the
        # formula rather than the auto-detection.
        pf.Script = namespace["_filter_script"](EPS_M, axisym)
        pf.RequestInformationScript = ""
        pf.RequestUpdateExtentScript = ""
        pf.UpdatePipeline()

        data = dsa.WrapDataObject(servermanager.Fetch(pf))
        if "Curvature" not in data.PointData.keys():
            # The filter raises inside VTK, which only prints the traceback and
            # hands back the input untouched. Turn that into a real failure.
            raise SystemExit(
                "verify_curvature: the filter produced no Curvature array for "
                "axisym=%s concave=%s -- see the traceback above."
                % (axisym, concave))
        phi = np.asarray(data.PointData["IcePhase"])
        kappa = np.asarray(data.PointData["Curvature"])
        ys = np.asarray(data.Points)[:, 1]
        Delete(pf)
        Delete(reader)
    finally:
        os.unlink(path)

    # Compare only on the phi = 0.5 contour, where the interface actually is.
    # Off-contour the logistic profile's own level sets have their own radii,
    # so kappa legitimately differs from 1/R there.
    on_contour = np.abs(phi - 0.5) < CONTOUR_HALF_WIDTH

    dim_factor = 2.0 if axisym else 1.0
    exact = (-1.0 if concave else 1.0) * dim_factor / R_M

    # Nothing is excluded from the score, the axis rows included: they are
    # where the azimuthal term is a 0/0 limit AND where the mesh edge corrupts
    # the second derivative, so they are the part most worth testing -- a
    # sintering neck sits exactly there. The near-axis column reports them
    # separately as well, to catch a regression there specifically.
    near_axis = np.zeros_like(on_contour)
    if axisym:
        near_axis = on_contour & (ys < AXIS_REPORT_CELLS * DX_M)

    sample = kappa[on_contour]
    if sample.size == 0:
        raise SystemExit("verify_curvature: no points on the phi=0.5 contour")

    axis_err = float("nan")
    if near_axis.any():
        axis_err = 100.0 * abs(float(np.median(kappa[near_axis])) - exact) / abs(exact)

    return float(np.median(sample)), exact, sample, axis_err


def main():
    namespace = load_macro()
    rows = []
    worst = 0.0

    for axisym in (False, True):
        for concave in (False, True):
            measured, exact, sample, axis_err = run_case(namespace, axisym, concave)
            err = 100.0 * _abs(measured - exact) / _abs(exact)
            worst = _max(worst, err)
            rows.append((
                "%s %s" % ("axisym " if axisym else "planar ",
                           "concave" if concave else "convex "),
                exact, measured, err,
                100.0 * float(np.percentile(np.abs(sample - exact), 90)) / abs(exact),
                sample.size, axis_err,
            ))

    print()
    print("  case              exact [1/m]   median [1/m]   err     p90 err   pts"
          "    near-axis err")
    print("  " + "-" * 88)
    for name, exact, measured, err, p90, npts, axis_err in rows:
        axis_txt = "     --" if axis_err != axis_err else "%6.2f%%" % axis_err
        print("  %-16s %+12.4g  %+13.4g  %5.2f%%  %6.2f%%  %6d    %s"
              % (name, exact, measured, err, p90, npts, axis_txt))
    print()
    print("  near-axis = the first %g cells off y = 0, also counted in the score."
          % AXIS_REPORT_CELLS)
    print("  Those rows carry both the 0/0 azimuthal limit and the mesh-edge")
    print("  second-derivative repair, and a sintering neck sits on the axis.")
    print()
    print("  radius R = %.4g m, eps = %.4g m, dx = %.4g m, tol = %.1f%%"
          % (R_M, EPS_M, DX_M, TOL_PCT))

    if worst > TOL_PCT:
        print("\n  FAIL: worst case is %.2f%% off (tolerance %.1f%%)"
              % (worst, TOL_PCT))
        return 1
    print("\n  PASS: all four cases within %.1f%% (worst %.2f%%)"
          % (TOL_PCT, worst))
    return 0


if __name__ == "__main__":
    sys.exit(main())

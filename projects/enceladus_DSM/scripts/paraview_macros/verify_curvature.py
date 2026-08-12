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

    PV_DENOM=gradient ... verify_curvature.py     # the other denominator

Exits non-zero if any case misses its analytic value by more than TOL_PCT.

On these synthetic fields the "gradient" denominator scores slightly BETTER
than the "equipartition" default (max 3.05% vs 4.73%), which is expected and
not an argument for switching: these profiles are exactly at equilibrium and
noise-free, so the measured |grad phi| is the exact denominator while
equipartition is only as good as the discretization. Equipartition earns its
place on real runs, where |grad phi| is a differentiated quantity that decays
to zero at the band edges and carries noise into the divide. The macro prints
the measured equipartition ratio per run so the assumption stays checkable.
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

# Band mask, matching the macro's own localizer: 16*(phi*(1-phi))^2 >= this.
LOC_FLOOR = 0.01
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


def build_grid(axisym, kind):
    """A vtkStructuredGrid carrying a logistic phase field.

    kind "flat" gives one flat interface through x = 0 (ice at x > 0); "convex"
    an ice disc of radius R; "concave" an air bubble of radius R inside ice.
    Planar puts the feature in the middle of a square domain. Axisym makes y
    the radius with the domain starting at y = 0, so a circle centred on y = 0
    is a sphere of revolution.
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
    if kind == "flat":
        sdf = gx - xc
    elif kind == "concave":
        sdf = rad - R_M
    else:
        sdf = R_M - rad
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


def run_case(namespace, axisym, kind):
    """Apply the macro's filter to a synthetic case.

    Returns (err_pct over the band, npts, near-axis err, profile) where the
    profile is a list of (phi, median kappa, median exact) rows for printing.
    """
    from paraview.simple import (
        Delete, OpenDataFile, ProgrammableFilter, servermanager,
    )
    from vtk.numpy_interface import dataset_adapter as dsa

    grid, rad = build_grid(axisym, kind)

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
                "axisym=%s kind=%s -- see the traceback above."
                % (axisym, kind))
        phi = np.asarray(data.PointData["IcePhase"])
        kappa = np.asarray(data.PointData["Curvature"])
        ys = np.asarray(data.Points)[:, 1]
        Delete(pf)
        Delete(reader)
    finally:
        os.unlink(path)

    # Score the WHOLE band, using the macro's own localizer. Each level set of
    # the logistic profile is its own circle, so the exact kappa follows the
    # point's actual radius rather than R -- which is why `rad` comes back from
    # build_grid. For the flat case the exact answer is 0 at every phi.
    loc = 16.0 * (phi * (1.0 - phi)) ** 2
    band = loc >= LOC_FLOOR

    dim_factor = 2.0 if axisym else 1.0
    if kind == "flat":
        exact = np.zeros_like(kappa)
        # kappa = 0 has no relative scale, so measure the artifact against a
        # real curvature: that of an R_M-radius feature.
        scale = np.full_like(kappa, 1.0 / R_M)
    else:
        sign = -1.0 if kind == "concave" else 1.0
        exact = sign * dim_factor / np.maximum(rad, 1e-30)
        scale = np.abs(exact)

    if not band.any():
        raise SystemExit("verify_curvature: no points in the band")

    err = 100.0 * np.abs(kappa[band] - exact[band]) / scale[band]

    axis_err = float("nan")
    if axisym:
        near = band & (ys < AXIS_REPORT_CELLS * DX_M)
        if near.any():
            axis_err = float(np.median(
                100.0 * np.abs(kappa[near] - exact[near]) / scale[near]))

    # A phi-binned profile: this is what makes a regularization residue
    # legible -- it is antisymmetric about phi = 0.5 and worst at the edges.
    profile = []
    for lo_phi, hi_phi in ((0.02, 0.10), (0.10, 0.30), (0.30, 0.70),
                           (0.70, 0.90), (0.90, 0.98)):
        sel = band & (phi >= lo_phi) & (phi < hi_phi)
        if sel.any():
            profile.append((0.5 * (lo_phi + hi_phi),
                            float(np.median(kappa[sel])),
                            float(np.median(exact[sel])),
                            float(np.median(100.0 * np.abs(kappa[sel] - exact[sel])
                                            / scale[sel]))))

    return err, int(band.sum()), axis_err, profile


def main():
    namespace = load_macro()
    rows = []
    worst = 0.0
    profiles = {}

    cases = [(False, "flat"), (False, "convex"), (False, "concave"),
             (True, "convex"), (True, "concave")]
    for axisym, kind in cases:
        err, npts, axis_err, profile = run_case(namespace, axisym, kind)
        name = "%s %s" % ("axisym" if axisym else "planar", kind)
        worst = _max(worst, float(np.max(err)))
        rows.append((name, float(np.median(err)), float(np.percentile(err, 90)),
                     float(np.max(err)), npts, axis_err))
        profiles[name] = profile

    print()
    print("  Error vs the analytic kappa, over the whole diffuse band")
    print()
    print("  case              median     p90       MAX      pts    near-axis")
    print("  " + "-" * 68)
    for name, med, p90, mx, npts, axis_err in rows:
        axis_txt = "      --" if axis_err != axis_err else "%7.2f%%" % axis_err
        print("  %-16s %7.2f%%  %7.2f%%  %8.2f%%  %6d  %s"
              % (name, med, p90, mx, npts, axis_txt))

    print()
    print("  Profile through the band (median per phi bin) -- a regularization")
    print("  residue shows up here as a sign flip about phi = 0.5:")
    print()
    print("    case              phi    kappa [1/m]   exact [1/m]     err")
    print("    " + "-" * 62)
    for name in ("planar flat", "planar convex"):
        for phi_mid, k_med, k_exact, err_med in profiles[name]:
            print("    %-16s %4.2f  %+11.4g  %+12.4g  %7.2f%%"
                  % (name, phi_mid, k_med, k_exact, err_med))
        print()

    print("  near-axis = the first %g cells off y = 0, also counted in the score."
          % AXIS_REPORT_CELLS)
    print("  radius R = %.4g m, eps = %.4g m, dx = %.4g m, band = 16(phi(1-phi))^2 "
          ">= %g" % (R_M, EPS_M, DX_M, LOC_FLOOR))
    print("  tolerance %.1f%% on the MAX error over the band." % TOL_PCT)

    if worst > TOL_PCT:
        print("\n  FAIL: worst point is %.2f%% off (tolerance %.1f%%)"
              % (worst, TOL_PCT))
        return 1
    print("\n  PASS: all %d cases within %.1f%% across the band (worst %.2f%%)"
          % (len(rows), TOL_PCT, worst))
    return 0


if __name__ == "__main__":
    sys.exit(main())

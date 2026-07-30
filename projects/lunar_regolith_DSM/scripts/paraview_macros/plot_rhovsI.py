"""
plot_rhovsI.py — ParaView macro: saturation vapor density over ice, rho_vs^I(T).

Adds a Programmable Filter on top of the active source that evaluates the exact
same expression as RhoVS_I() in src/material_properties.c:

    temK = T[C] + 273.15
    Pvs  = exp(K0/temK + K1 + K2*temK + K3*temK^2 + K4*temK^3 + K5*ln(temK))
    rho_vs^I = rho_air * bb * Pvs / (Patm - Pvs)          [kg/m^3]

and appends two point arrays:

    RhoVS_I          saturation vapor density [kg/m^3]
    Supersaturation  (rho_v - rho_vs^I) / rho_vs^I   [-]   (dimensionless;
                     positive => deposition, negative => sublimation)

Supersaturation is only added when the input carries a VaporDensity array.

Install
-------
ParaView -> Macros -> Add new macro... -> pick this file.
(Or copy it into ~/.config/ParaView/Macros/ ; on macOS
 ~/Library/Applications/ParaView/Macros/ depending on version.)

Use
---
Select a solV_*.vts reader (or any filter carrying a "Temperature" array) in
the Pipeline Browser, then click the macro button. The new filter is shown in a
render view and coloured by RhoVS_I; switch the colouring to Supersaturation in
the toolbar to see the driving force.

Notes
-----
* Temperature is expected in DEGREES CELSIUS, which is what the solver writes
  to DOF 1 and what scripts/plot_fields.py labels "Temperature".
* rho_air is the solver's -rho_air (AppCtx default 1.341 kg/m^3). If a run used
  a different value, change RHO_AIR below or set the environment variable
  PV_RHO_AIR before launching ParaView.
"""

import os

from paraview.simple import (
    ColorBy,
    CreateRenderView,
    GetActiveSource,
    GetActiveView,
    GetColorTransferFunction,
    GetViews,
    ProgrammableFilter,
    RenameSource,
    Render,
    SetActiveView,
    Show,
)

# --- solver constants (material_properties.c / enceladus_main.c) -------------
RHO_AIR = float(os.environ.get("PV_RHO_AIR", 1.341))  # kg/m^3, AppCtx default

# Candidate names for the input arrays, in order of preference.
TEMP_ARRAY_CANDIDATES = ("Temperature", "temp", "T")
RHOV_ARRAY_CANDIDATES = ("VaporDensity", "rho_v", "rhov")


# The body below runs inside ParaView's Programmable Filter, so it must be
# self-contained: it cannot see anything defined in this module except what is
# textually substituted in (RHO_AIR, the candidate name tuples).
FILTER_SCRIPT = '''
import numpy as np

RHO_AIR         = {rho_air!r}
TEMP_CANDIDATES = {temp_candidates!r}
RHOV_CANDIDATES = {rhov_candidates!r}

# Constants copied verbatim from RhoVS_I() in src/material_properties.c
PATM = 101325.0
BB   = 0.62
K0, K1, K2 = -0.5865e4, 0.2224e2, 0.1375e-1
K3, K4, K5 = -0.3403e-4, 0.2697e-7, 0.6918

inp = inputs[0]
output.ShallowCopy(inp.VTKObject)

available = list(inp.PointData.keys())


def pick(candidates):
    for cand in candidates:
        if cand in available:
            return cand
    return None


temp_name = pick(TEMP_CANDIDATES)
if temp_name is None:
    raise RuntimeError(
        "plot_rhovsI: no temperature point array found (looked for %s); "
        "available arrays: %s" % (list(TEMP_CANDIDATES), available)
    )

tem  = np.asarray(inp.PointData[temp_name], dtype=np.float64)   # deg Celsius
temK = tem + 273.15

Pvs = np.exp(K0 / temK + K1 + K2 * temK
             + K3 * temK ** 2 + K4 * temK ** 3
             + K5 * np.log(temK))

# Pvs -> Patm only near the boiling point; guard so the macro never divides by
# zero on a stray out-of-range temperature.
denom = np.maximum(PATM - Pvs, 1.0e-12)

rho_vs = RHO_AIR * BB * Pvs / denom
output.PointData.append(rho_vs, "RhoVS_I")

# Supersaturation, the driving force for deposition/sublimation.
rhov_name = pick(RHOV_CANDIDATES)
if rhov_name is not None:
    rho_v = np.asarray(inp.PointData[rhov_name], dtype=np.float64)
    output.PointData.append((rho_v - rho_vs) / rho_vs, "Supersaturation")
else:
    print("plot_rhovsI: no vapor-density array (looked for %s), "
          "skipping Supersaturation" % (list(RHOV_CANDIDATES),))
'''.format(rho_air=RHO_AIR,
           temp_candidates=TEMP_ARRAY_CANDIDATES,
           rhov_candidates=RHOV_ARRAY_CANDIDATES)


def _render_view():
    """Return a render view to display in.

    The active view may well be a chart view (Line Chart, Spreadsheet, ...),
    which has no colour map — colouring a representation there raises
    AttributeError on UseSeparateColorMap. So prefer the active view only when
    it is a render view, then any existing render view, then make one.
    """
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


def _apply_preset(lut):
    """Apply the first colour preset that this ParaView build knows about.

    ParaView 6 renamed "Viridis (matplotlib)" to plain "Viridis", and
    ApplyPreset raises on an unknown name — so try in order and give up
    quietly rather than fail the whole macro over cosmetics.
    """
    for name in ("Viridis", "Viridis (matplotlib)", "Cool to Warm"):
        try:
            lut.ApplyPreset(name, True)
            return
        except RuntimeError:
            continue


def main():
    src = GetActiveSource()
    if src is None:
        raise RuntimeError(
            "plot_rhovsI: no active source. Select a solV_*.vts reader (or any "
            "filter with a Temperature array) in the Pipeline Browser first."
        )

    pf = ProgrammableFilter(Input=src)
    pf.OutputDataSetType = "Same as Input"
    pf.Script = FILTER_SCRIPT
    pf.RequestInformationScript = ""
    pf.RequestUpdateExtentScript = ""
    RenameSource("RhoVS_I", pf)
    pf.UpdatePipeline()

    view = _render_view()
    display = Show(pf, view)
    ColorBy(display, ("POINTS", "RhoVS_I"))
    display.RescaleTransferFunctionToDataRange(True, False)
    display.SetScalarBarVisibility(view, True)
    _apply_preset(GetColorTransferFunction("RhoVS_I"))
    Render(view)

    return pf


main()

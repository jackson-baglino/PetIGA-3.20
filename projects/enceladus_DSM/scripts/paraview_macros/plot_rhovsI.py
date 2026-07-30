"""
plot_rhovsI.py — ParaView macro: saturation vapor density over ice, rho_vs^I(T).

Adds a Programmable Filter on top of the active source that evaluates the exact
same expression as RhoVS_I() in src/material_properties.c:

    temK = T[C] + 273.15
    Pvs  = exp(K0/temK + K1 + K2*temK + K3*temK^2 + K4*temK^3 + K5*ln(temK))
    rho_vs^I = rho_air * bb * Pvs / (Patm - Pvs)          [kg/m^3]

and appends it as a point array "RhoVS_I".

Install
-------
ParaView -> Macros -> Add new macro... -> pick this file.
(Or copy it into ~/.config/ParaView/Macros/ ; on macOS
 ~/Library/Applications/ParaView/Macros/ depending on version.)

Use
---
Select a solV_*.vts reader (or any filter carrying a "Temperature" array) in
the Pipeline Browser, then click the macro button. The new filter is shown and
coloured by RhoVS_I.

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
    GetActiveSource,
    GetActiveView,
    GetColorTransferFunction,
    ProgrammableFilter,
    RenameSource,
    Render,
    Show,
)

# --- solver constants (material_properties.c / enceladus_main.c) -------------
RHO_AIR = float(os.environ.get("PV_RHO_AIR", 1.341))  # kg/m^3, AppCtx default

# Candidate names for the temperature array, in order of preference.
TEMP_ARRAY_CANDIDATES = ("Temperature", "temp", "T")


# The body below runs inside ParaView's Programmable Filter, so it must be
# self-contained: it cannot see anything defined in this module except what is
# textually substituted in (RHO_AIR, TEMP_ARRAY_CANDIDATES).
FILTER_SCRIPT = '''
import numpy as np

RHO_AIR = {rho_air!r}
CANDIDATES = {candidates!r}

# Constants copied verbatim from RhoVS_I() in src/material_properties.c
PATM = 101325.0
BB   = 0.62
K0, K1, K2 = -0.5865e4, 0.2224e2, 0.1375e-1
K3, K4, K5 = -0.3403e-4, 0.2697e-7, 0.6918

inp = inputs[0]
output.ShallowCopy(inp.VTKObject)

name = None
for cand in CANDIDATES:
    if cand in inp.PointData.keys():
        name = cand
        break
if name is None:
    raise RuntimeError(
        "plot_rhovsI: no temperature point array found (looked for %s); "
        "available arrays: %s" % (list(CANDIDATES), list(inp.PointData.keys()))
    )

tem  = np.asarray(inp.PointData[name], dtype=np.float64)   # degrees Celsius
temK = tem + 273.15

Pvs = np.exp(K0 / temK + K1 + K2 * temK
             + K3 * temK ** 2 + K4 * temK ** 3
             + K5 * np.log(temK))

# Pvs -> Patm only near the boiling point; guard so the macro never divides by
# zero on a stray out-of-range temperature.
denom = np.maximum(PATM - Pvs, 1.0e-12)

rho_vs = RHO_AIR * BB * Pvs / denom
output.PointData.append(rho_vs, "RhoVS_I")
'''.format(rho_air=RHO_AIR, candidates=TEMP_ARRAY_CANDIDATES)


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
    pf.CopyArrays = 1
    RenameSource("RhoVS_I", pf)
    pf.UpdatePipeline()

    view = GetActiveView()
    display = Show(pf, view)
    ColorBy(display, ("POINTS", "RhoVS_I"))
    display.RescaleTransferFunctionToDataRange(True, False)
    display.SetScalarBarVisibility(view, True)
    GetColorTransferFunction("RhoVS_I").ApplyPreset("Viridis (matplotlib)", True)
    Render()

    return pf


main()

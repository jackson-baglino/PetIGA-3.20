r"""
plot_growth_rate.py — ParaView macro: interface normal velocity v_n, split into
its curvature and phase-change parts.

Runs ON TOP of the filter that plot_rhovsI.py creates. Select that filter (it is
named "RhoVS_I" in the pipeline) and click this macro. It needs the `Curvature`
and `Supersaturation` arrays that one produces; it deliberately does not
recompute kappa, so there is exactly one curvature implementation in the repo.

What it computes
----------------
docs/curvature_driven_growth.md derives, by projecting the ice residual in
src/assembly.c onto the translation mode f' = dphi/ds, that the interface moves at

    v_n  =  -3*M*eps*kappa  +  (eps*alph_sub/(5*rho_ice)) * (rho_v - rho_vs)
            \____________/     \_______________________________________/
             curvature                     phase change

with M = mob_sub. **v_n > 0 means the ice is GROWING** — the interface advances
into the air. The 5 is the ratio of the two profile integrals,
(int f'^2 ds)/(int loc*f' ds) = (1/6eps)/(1/30); it is the same number as the
hard-coded `a1 = 5.0` in enceladus_main.c.

Appended point arrays, all masked to the diffuse band (where kappa is defined):

    V_normal            v_n  [m/s]          <- the headline field
    V_curvature         the Allen-Cahn part [m/s]
    V_phasechange       the sublimation part [m/s]
    V_normal_um_per_day v_n in um/day, because 5e-13 m/s is unreadable

Why the split is the interesting part
-------------------------------------
On the wedge run the two terms are individually ~1e-11 m/s and cancel to within
7-12%, leaving v_n ~ 5e-13. Sublimation really is happening at the growing
interface — it is simply outrun by the curvature term. Plotting only v_n hides
that; plotting V_curvature and V_phasechange beside it is what shows a
near-cancellation rather than a single dominant process.

This is also the same quantity as `Supersaturation_GT / beta`: the identity
v_n = (sigma - d0*kappa)/beta is just this equation rearranged, with
d0 = 15*M*rho_ice/(alph_sub*rho_vs) and beta = 5*rho_ice/(eps*alph_sub*rho_vs).
The two-term form is used here because it takes the run's own mob_sub and
alph_sub directly and gives the decomposition for free.

Where the constants come from
-----------------------------
mob_sub, alph_sub, rho_ice and eps are read from the run's own `outp.txt`
startup banner, found by walking up from the data file. That banner is the
authority: mob_sub and alph_sub are DERIVED at startup from d0_sub0/beta_sub0
and usually appear in no .opts file at all. Override any of them with the
constants below or the PV_* environment variables; the macro always prints what
it used and where each value came from.

Install
-------
ParaView -> Macros -> Add new macro... -> pick this file. (ParaView copies it,
so re-add after editing.) Open View -> Python Shell first — the provenance and
the term-balance summary are printed there.

Use
---
1. Open the run's pf.pvd.
2. Click `plot_rhovsI` to build the Curvature / Supersaturation fields.
3. With that filter selected, click `plot_growth_rate`.
4. To read v_n on the interface itself rather than smeared across the band,
   apply Contour on the result: Contour By = IcePhase, value 0.5, colour by
   V_normal. That is where the sharp-interface v_n actually lives.

Caveats
-------
* v_n is an INTERFACE quantity. The band values are the same formula evaluated
  off-contour and are only approximately meaningful; read the contour.
* The derivation assumes eps*kappa << 1 and an equilibrium profile. Check the
  equipartition ratio plot_rhovsI.py prints (1.0 = at equilibrium) before
  trusting the numbers at a tight neck.
* Leading order in the matched asymptotics: it omits the a2 thin-interface
  correction, so absolute magnitudes carry ~20% systematic (see
  docs/curvature_driven_growth.md section 7). Signs and ratios are solid.
"""

import builtins
import os
import re

# Clicking plot_rhovsI first runs a Programmable Filter, whose preamble
# star-imports numpy_interface.algorithms into __main__ and rebinds max/min/sum
# to array reductions taking `axis` second. Macros share that namespace, so a
# bare max(a, b) below would fail with "cannot be interpreted as an integer"
# — but only when this macro is used the way it is meant to be, after
# plot_rhovsI. Bind the real ones up front.
_max, _abs = builtins.max, builtins.abs

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

# --- constants; None => read from the run's outp.txt -------------------------
MOB_SUB = None      # M   [m/s]      PV_MOB_SUB
ALPH_SUB = None     # alpha_sub [1/s] PV_ALPH_SUB
RHO_ICE = None      # [kg/m^3]        PV_RHO_ICE
EPS = None          # [m]             PV_EPS

# (int f'^2 ds)/(int loc*f' ds) in units of 1/eps -- see the module docstring.
# Same number as `a1` in enceladus_main.c. Only change this if the double well
# or the phase-change localizer changes.
A1 = 5.0

CURVATURE_ARRAY = "Curvature"
SUPERSAT_ARRAY = "Supersaturation"
RHOVS_ARRAY = "RhoVS_I"


FILTER_BODY = r'''
import numpy as np
import builtins   # ParaView shadows max/min/sum here; see verify_curvature.py

inp = inputs[0]
output.ShallowCopy(inp.VTKObject)

available = list(inp.PointData.keys())
for needed in (CURV_NAME, SUPERSAT_NAME, RHOVS_NAME):
    if needed not in available:
        raise RuntimeError(
            "plot_growth_rate: no '%s' array on the selected source. Run the "
            "plot_rhovsI macro first and select the filter it creates (named "
            "RhoVS_I). Available: %s" % (needed, available))

kappa  = np.asarray(inp.PointData[CURV_NAME],     dtype=np.float64)
sigma  = np.asarray(inp.PointData[SUPERSAT_NAME], dtype=np.float64)
rho_vs = np.asarray(inp.PointData[RHOVS_NAME],    dtype=np.float64)

# (rho_v - rho_vs) rebuilt from the dimensionless supersaturation, so this
# macro stays agnostic about the vapour array's name.
drho = sigma * rho_vs

v_curv = -3.0 * MOB_SUB * EPS * kappa
v_pc = (EPS * ALPH_SUB / (A1 * RHO_ICE)) * drho

# kappa is already zero outside the diffuse band (plot_rhovsI masks it there),
# which would leave v_pc alone painting a velocity across the bulk where there
# is no interface at all. Mask both parts to the same band.
band = kappa != 0.0
v_curv = np.where(band, v_curv, 0.0)
v_pc = np.where(band, v_pc, 0.0)
v_n = v_curv + v_pc

output.PointData.append(v_curv, "V_curvature")
output.PointData.append(v_pc, "V_phasechange")
output.PointData.append(v_n, "V_normal")
output.PointData.append(v_n * 86400.0 * 1.0e6, "V_normal_um_per_day")

n_band = int(band.sum())
print("plot_growth_rate: v_n on %d band points" % n_band)
if n_band:
    grow = int((v_n[band] > 0).sum())
    kb, cb, pb = kappa[band], v_curv[band], v_pc[band]
    print("  V_curvature    median %+.4e m/s   |max| %.4e" % (np.median(cb), np.abs(cb).max()))
    print("  V_phasechange  median %+.4e m/s   |max| %.4e" % (np.median(pb), np.abs(pb).max()))
    print("  V_normal       median %+.4e m/s   |max| %.4e" % (np.median(v_n[band]), np.abs(v_n[band]).max()))
    print("  %d of %d band points growing (v_n > 0), %d sublimating"
          % (grow, n_band, n_band - grow))
    # How completely do the two terms cancel? Near 100% means the net motion is
    # a small residual of two much larger opposing rates -- the situation on the
    # wedge run, and the reason v_n alone is a misleading thing to look at.
    denom = np.maximum(np.abs(cb), 1.0e-300)
    cancel = 100.0 * (1.0 - np.abs(cb + pb) / denom)
    print("  the two terms cancel by a median of %.1f%% "
          "(100%% = perfect standoff)" % np.median(cancel))
    for name, sel in (("concave (kappa<0)", kb < 0.0), ("convex  (kappa>0)", kb > 0.0)):
        if np.any(sel):
            print("    %s: median v_n = %+.4e m/s (%+.4g um/day)"
                  % (name, np.median((cb + pb)[sel]),
                     np.median((cb + pb)[sel]) * 86400.0 * 1.0e6))
'''


def _filter_script(mob_sub, alph_sub, rho_ice, eps):
    header = (
        "MOB_SUB  = %r\n"
        "ALPH_SUB = %r\n"
        "RHO_ICE  = %r\n"
        "EPS      = %r\n"
        "A1       = %r\n"
        "CURV_NAME     = %r\n"
        "SUPERSAT_NAME = %r\n"
        "RHOVS_NAME    = %r\n"
    ) % (mob_sub, alph_sub, rho_ice, eps, A1,
         CURVATURE_ARRAY, SUPERSAT_ARRAY, RHOVS_ARRAY)
    return header + FILTER_BODY


def _walk_to_reader(src):
    """Follow Input links upstream until something exposes a FileName."""
    seen = set()
    stack = [src]
    while stack:
        node = stack.pop(0)
        key = id(node)
        if node is None or key in seen:
            continue
        seen.add(key)
        for prop in ("FileName", "FileNames"):
            try:
                value = getattr(node, prop)
            except AttributeError:
                continue
            if not value:
                continue
            if isinstance(value, str):
                return value
            try:
                return str(value[0])
            except (TypeError, IndexError):
                continue
        try:
            upstream = node.Input
        except AttributeError:
            continue
        if upstream is None:
            continue
        try:
            stack.extend(list(upstream))
        except TypeError:
            stack.append(upstream)
    return None


def _scan_log(src):
    """mob_sub, alph_sub, rho_ice, eps from the run's outp.txt startup banner.

    Returns a dict of whatever was found. The banner is the authority because
    mob_sub and alph_sub are derived at startup from d0_sub0/beta_sub0 and
    typically appear in no .opts file.
    """
    path = _walk_to_reader(src)
    if not path:
        return {}
    directory = os.path.dirname(os.path.abspath(path))
    patterns = {
        "eps": r"eps\s*=\s*([0-9.eE+-]+)\s*m",
        "mob_sub": r"mob_sub\s*=\s*([0-9.eE+-]+)",
        "alph_sub": r"alph_sub\s*=\s*([0-9.eE+-]+)",
        "rho_ice": r"rho_ice\s*=\s*([0-9.eE+-]+)",
    }
    for _ in range(4):                       # <run>/vtkOut/x.vts -> <run> and up
        candidate = os.path.join(directory, "outp.txt")
        if os.path.isfile(candidate):
            try:
                with open(candidate, errors="replace") as handle:
                    text = handle.read()
            except OSError:
                return {}
            found = {}
            for key, pattern in patterns.items():
                match = re.search(pattern, text)
                if match:
                    try:
                        found[key] = float(match.group(1))
                    except ValueError:
                        pass
            return found
        parent = os.path.dirname(directory)
        if parent == directory:
            break
        directory = parent
    return {}


def _resolve(name, env_key, constant, from_log):
    """PV_* env var, then the module constant, then the run's outp.txt."""
    if env_key in os.environ:
        return float(os.environ[env_key]), "the %s env var" % env_key
    if constant is not None:
        return float(constant), "the %s constant in the macro" % name
    if from_log is not None:
        return float(from_log), "outp.txt"
    return None, "NOT FOUND"


def _render_view():
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


def _apply_preset(lut, names):
    for name in names:
        try:
            lut.ApplyPreset(name, True)
            return
        except RuntimeError:
            continue


def main():
    src = GetActiveSource()
    if src is None:
        raise RuntimeError(
            "plot_growth_rate: no active source. Run the plot_rhovsI macro "
            "first, then select the filter it creates (named RhoVS_I).")

    present = list(src.PointData.keys())
    if CURVATURE_ARRAY not in present:
        raise RuntimeError(
            "plot_growth_rate: the selected source has no '%s' array, so there "
            "is nothing to build v_n from.\n  Run the plot_rhovsI macro first "
            "and select the filter it creates (named RhoVS_I).\n  Arrays here: "
            "%s" % (CURVATURE_ARRAY, present))

    log = _scan_log(src)
    values, sources, missing = {}, {}, []
    for name, env_key, constant in (("MOB_SUB", "PV_MOB_SUB", MOB_SUB),
                                    ("ALPH_SUB", "PV_ALPH_SUB", ALPH_SUB),
                                    ("RHO_ICE", "PV_RHO_ICE", RHO_ICE),
                                    ("EPS", "PV_EPS", EPS)):
        key = name.lower()
        value, origin = _resolve(name, env_key, constant, log.get(key))
        values[key], sources[key] = value, origin
        if value is None:
            missing.append((name, env_key))

    if missing:
        raise RuntimeError(
            "plot_growth_rate: could not determine %s.\n"
            "  These set the velocity scale, so the macro will not guess. The "
            "run's outp.txt was not found next to the data (looked up to 4 "
            "levels above it).\n  Either point ParaView at the data inside its "
            "run directory, or set %s before launching ParaView, or edit the "
            "constants at the top of the macro."
            % (", ".join(n for n, _ in missing),
               "/".join(e for _, e in missing)))

    print("plot_growth_rate: constants")
    for key, label in (("mob_sub", "mob_sub  [m/s]  "),
                       ("alph_sub", "alph_sub [1/s]  "),
                       ("rho_ice", "rho_ice  [kg/m3]"),
                       ("eps", "eps      [m]    ")):
        print("  %s = %-12.4e  (%s)" % (label, values[key], sources[key]))

    pf = ProgrammableFilter(Input=src)
    pf.OutputDataSetType = "Same as Input"
    pf.Script = _filter_script(values["mob_sub"], values["alph_sub"],
                               values["rho_ice"], values["eps"])
    pf.RequestInformationScript = ""
    pf.RequestUpdateExtentScript = ""
    RenameSource("GrowthRate", pf)
    pf.UpdatePipeline()

    view = _render_view()
    display = Show(pf, view)
    ColorBy(display, ("POINTS", "V_normal"))
    display.RescaleTransferFunctionToDataRange(True, False)
    display.SetScalarBarVisibility(view, True)

    # v_n is signed and the sign is the whole message (growing vs sublimating),
    # so it gets a diverging map centred on zero. Auto-rescaling to an
    # asymmetric range would put white off zero and misreport which parts of
    # the interface are advancing.
    lut = GetColorTransferFunction("V_normal")
    _apply_preset(lut, ("Cool to Warm (Extended)", "Cool to Warm"))
    lo, hi = pf.PointData["V_normal"].GetRange()
    half = _max(_abs(lo), _abs(hi))
    if half > 0.0:
        lut.RescaleTransferFunction(-half, half)
    Render(view)

    print("plot_growth_rate: added V_normal, V_curvature, V_phasechange, "
          "V_normal_um_per_day")
    print("  Contour the result at IcePhase = 0.5 and colour by V_normal to "
          "read v_n on the interface itself.")
    return pf


main()

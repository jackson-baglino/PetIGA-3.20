"""
plot_rhovsI.py — ParaView macro: saturation vapor density over ice, rho_vs^I(T),
flat-interface and with the Gibbs–Thomson curvature correction.

Adds a Programmable Filter on top of the active source that evaluates the exact
same expression as RhoVS_I() in src/material_properties.c:

    temK = T[C] + 273.15
    Pvs  = exp(K0/temK + K1 + K2*temK + K3*temK^2 + K4*temK^3 + K5*ln(temK))
    rho_vs^I = rho_air * bb * Pvs / (Patm - Pvs)          [kg/m^3]

and then re-applies the Gibbs–Thomson correction that used to live in the
solver (removed 2026-07-21 along with d0_GT/Curvature/Hessian):

    rho_vs^eff = rho_vs^I * (1 + d0 * kappa),   kappa = -div(grad phi / |grad phi|)
    d0(T)      = gamma_iv * V_m / (R_g * temK)            [m]

Appended point arrays:

    Curvature           kappa [1/m], masked to the diffuse interface (see below)
    RhoVS_I             flat-interface saturation vapor density [kg/m^3]
    RhoVS_I_GT          curvature-corrected saturation vapor density [kg/m^3]
    Supersaturation     (rho_v - rho_vs^I)   / rho_vs^I      [-]
    Supersaturation_GT  (rho_v - rho_vs^eff) / rho_vs^eff    [-]

Positive supersaturation => deposition, negative => sublimation. The two
Supersaturation arrays are only added when the input carries a VaporDensity
array; the GT pair is only added when it carries a phase array.

Sign convention: grad phi points INTO the ice, so kappa = +(d-1)/R on a convex
ice grain. Small grains therefore get a HIGHER rho_vs^eff than large ones —
the Lifshitz–Slyozov–Wagner driving force that makes vapor flow small -> large.

How kappa is computed
---------------------
kappa = -div(n) with n = grad phi/|grad phi|, expanded as

    kappa = -(laplacian(phi) - n.H.n) / G

where H is the Hessian of phi. grad phi and H come from VTK's mesh-aware
gradient operator applied once and twice. Two deliberate choices:

* **The bracket, not the two-term form.** The solver's deleted Curvature()
  used the algebraically identical -L/G + (g.H.g)/G^3 with a Tikhonov
  G^2 = |grad phi|^2 + eps_reg^2. That regularization does NOT cancel between
  the two terms: on a flat interface, where kappa must be zero, it leaves
  -eps_reg^2 * phi'' / G^3 -- zero at phi = 0.5, sign-flipped either side of
  it, and divergent toward the band edges (~1e5 /m at phi = 0.03 for
  eps = 0.86 um, far larger than any real curvature in the run). Factoring G
  out kills it identically: a flat interface has L = n.H.n exactly.

* **Direction and magnitude from different places.** n = grad phi/|grad phi|
  is a unit vector, so gradient noise can tilt it but cannot rescale kappa;
  and it points toward increasing phi, i.e. into the ice, automatically -- a
  convex ice grain gives kappa > 0 with no sign convention to impose. The
  DENOMINATOR instead comes from equipartition, |grad phi| = phi(1-phi)/eps,
  an analytic function of phi alone. That is the smoothing: the measured
  |grad phi| is a differentiated quantity that decays to zero at the band
  edges, and dividing by it there amplifies noise into exactly the region
  where the interface is least resolved. Set DENOM = "gradient" (or
  PV_DENOM=gradient) to divide by the measured |grad phi| instead -- exact
  where the profile is at equilibrium, noisier at the edges.

  The macro reports how well equipartition actually holds, as the median and
  spread of |grad phi| * eps / (phi(1-phi)) over the band. That ratio is 1 at
  equilibrium; far from 1 means the profile is being driven or under-resolved,
  and the two DENOM choices will disagree by that factor.

**Away from the interface kappa is meaningless** -- grad phi -> 0 there and
the profile relation does not apply. So kappa is zeroed wherever the localizer
16*phi^2*(1-phi)^2 (the solver's ice^2*air^2, scaled to peak at 1) drops below
LOC_FLOOR. At the default 0.01 that keeps phi in [0.026, 0.974], about the
visibly diffuse band. Outside it kappa = 0 makes RhoVS_I_GT fall back exactly
to RhoVS_I, which is the right bulk limit.

Axisymmetric runs
-----------------
In an r–z run the grid is 2D but the surface is a surface of revolution, and
the planar formula misses the azimuthal mode. With AXISYM on, the macro adds it:

    kappa_3D = kappa_planar - n_y / y,        n = grad phi / |grad phi|

using the project's convention that **y is the radius and the axis is y = 0**
(same as `postprocess/neck_width.py --axisym`). Without it a sphere reads 1/R
instead of 2/R — every curvature, and so every GT correction, comes out half.

AXISYM is auto-detected from `-axisym` in the staged .opts and can be forced
with PV_AXISYM=1/0. The macro always prints which mode it used — check that
line, since nothing in the .vts itself says whether a run was axisymmetric.

Install
-------
ParaView -> Macros -> Add new macro... -> pick this file.
(Or copy it into ~/.config/ParaView/Macros/ ; on macOS
 ~/Library/Applications/ParaView/Macros/ depending on version.)

Use
---
Select a solV_*.vts reader (or any filter carrying a "Temperature" array) in
the Pipeline Browser, then click the macro button. The new filter is shown in a
render view coloured by Curvature; switch the colouring in the toolbar to
Supersaturation_GT to see the driving force, or to Supersaturation to see what
the run actually solved (the solver has no GT term — these fields are a
diagnostic laid over the result, not a re-solve).

Notes
-----
* Temperature is expected in DEGREES CELSIUS, which is what the solver writes
  to DOF 1 and what scripts/plot_fields.py labels "Temperature".
* rho_air is the solver's -rho_air (AppCtx default 1.341 kg/m^3). If a run used
  a different value, change RHO_AIR below or set the environment variable
  PV_RHO_AIR before launching ParaView.
* GAMMA_VM_OVER_R is preprocess/comp_eps.py's gamma*V_m/R = 2.574e-7 m*K,
  giving d0 = 1.02e-9 m at -20 C. The legacy monitoring.c hardcoded 2.548e-7
  instead (~1% lower); comp_eps.py flags the discrepancy as unreconciled.
"""

import builtins
import os
import re

builtins_max = builtins.max

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

# gamma_iv * V_m / R_gas [m*K]; d0(T) = this / temK. See comp_eps.py.
GAMMA_VM_OVER_R = float(os.environ.get("PV_GAMMA_VM_OVER_R", 2.574e-7))

# Interface localizer floor: kappa is kept where 16*phi^2*(1-phi)^2 >= this,
# zeroed elsewhere. 0.01 => phi in [0.026, 0.974].
LOC_FLOOR = float(os.environ.get("PV_LOC_FLOOR", 0.01))

# eps [m] -- the phase-field decay length. NOT optional book-keeping: it sets
# the scale of the curvature denominator, so kappa scales as 1/eps if it is
# wrong, and it changes from run to run. Resolved in this order:
#
#   1. the PV_EPS environment variable
#   2. the EPS constant just below
#   3. `-eps` in a .opts the run script staged next to the data
#   4. estimated from the data as 1/(4 max|grad phi|)   <- last resort
#
# The macro prints which one it used every time. Set EPS here (or export
# PV_EPS) whenever the data is not sitting next to its .opts.
EPS = None

# Denominator for kappa: "equipartition" uses phi(1-phi)/eps, "gradient" uses
# the measured |grad phi|. See the block comment in the filter body.
DENOM = os.environ.get("PV_DENOM", "equipartition")

# Axisymmetric r-z mode. None => auto-detect from the staged .opts.
AXISYM = None

# Candidate names for the input arrays, in order of preference.
TEMP_ARRAY_CANDIDATES = ("Temperature", "temp", "T")
RHOV_ARRAY_CANDIDATES = ("VaporDensity", "rho_v", "rhov")
PHASE_ARRAY_CANDIDATES = ("IcePhase", "phaseice", "phi", "Phase")


# The body below runs inside ParaView's Programmable Filter, so it must be
# self-contained: it cannot see anything defined in this module. Constants are
# prepended as a separate header rather than substituted into the body with
# .format(), so the body can use dicts, sets and format strings freely without
# every brace in it having to be escaped.
FILTER_BODY = r'''
import numpy as np
from vtk.numpy_interface import algorithms as algs

# ParaView star-imports numpy_interface.algorithms into this namespace, which
# SHADOWS the builtin max/min/sum/round with array reductions whose second
# positional argument is `axis`. So a plain max(n, 1) here silently becomes
# algs.max(n, axis=1) and dies with "axis 1 is out of bounds". Reach for the
# real ones through `builtins`.
import builtins

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

# --- Gibbs-Thomson curvature correction -------------------------------------
# kappa = -div(grad phi/|grad phi|) = -(L - n.H.n)/G, with the unit normal from
# the numerical gradient and G from equipartition. This is NOT the deleted
# Curvature()'s two-term regularized form -- see the block comment below and
# the module docstring for why that form cannot be used pointwise.
phase_name = pick(PHASE_CANDIDATES)
rho_vs_eff = None

if phase_name is None:
    print("plot_rhovsI: no phase array (looked for %s), skipping the "
          "Gibbs-Thomson fields" % (list(PHASE_CANDIDATES),))
else:
    phi = inp.PointData[phase_name]

    # Mesh-aware gradients: g is (N,3), H is (N,3,3) with H[:,k,l] = d g_k/d x_l.
    g = algs.gradient(phi)
    H = algs.gradient(g)
    g = np.asarray(g, dtype=np.float64)
    H = np.asarray(H, dtype=np.float64).reshape(-1, 3, 3)

    phi_np = np.asarray(phi, dtype=np.float64)
    gmag = np.sqrt(np.sum(g * g, axis=1))
    gmax = float(gmag.max())

    # eps sets the denominator scale, so it is a real input now, not a
    # regularization knob: kappa scales as 1/eps if you get it wrong. Order of
    # preference is PV_EPS, the EPS constant, the run's staged .opts, and only
    # then a guess off the data (the equilibrium profile peaks at 1/(4 eps)).
    if EPS is not None and EPS > 0.0:
        EPS_USED = float(EPS)
        eps_src = EPS_ORIGIN
    elif gmax > 0.0:
        EPS_USED = 1.0 / (4.0 * gmax)
        eps_src = ("ESTIMATED as 1/(4 max|grad phi|) -- no eps given and none "
                   "found in a .opts; pass one via PV_EPS if this is wrong")
    else:
        raise RuntimeError("plot_rhovsI: phase field is constant; cannot "
                           "determine eps. Set PV_EPS or the EPS constant.")
    print("plot_rhovsI: eps = %.4g m (%s)" % (EPS_USED, eps_src))

    # --- direction and magnitude are taken from DIFFERENT places -------------
    #
    # kappa = -div(n) with n = grad phi/|grad phi| expands two ways that are
    # algebraically identical but numerically are not:
    #
    #   (a)  -L/G + (g.H.g)/G^3           <- one G per term
    #   (b)  -(L - n.H.n)/G               <- G factored out
    #
    # Form (a) is what the solver's Curvature() used, with a Tikhonov G^2 =
    # |g|^2 + eps_reg^2 to survive the bulk. That regularization does not
    # cancel between the two terms: on a FLAT interface, where kappa must be
    # 0, (a) leaves -eps_reg^2 * phi'' / G^3. That residue is zero at
    # phi = 0.5, flips sign across it, and diverges toward the band edges --
    # ~1e5 /m at phi = 0.03 for eps = 0.86 um, swamping any real curvature.
    #
    # Form (b) cannot do that: for a flat interface L = n.H.n exactly, so the
    # bracket vanishes identically whatever the denominator is. Use (b).
    #
    # The two inputs then come from different places, on purpose:
    #   * DIRECTION  n = g/|g| -- from the numerical gradient. It is a unit
    #     vector, so gradient noise tilts it slightly but cannot scale kappa.
    #     n automatically points toward increasing phi, i.e. INTO the ice.
    #   * MAGNITUDE  G -- from equipartition, |grad phi| = phi(1-phi)/eps, an
    #     analytic function of phi alone. This is the smoothing: the true
    #     |grad phi| is a differentiated quantity that goes to zero at the band
    #     edges, and dividing by it there amplifies noise. phi(1-phi)/eps has
    #     the same profile without the noise. Set DENOM = "gradient" to divide
    #     by the measured |grad phi| instead (exact, but noisier at the edges).
    L = H[:, 0, 0] + H[:, 1, 1] + H[:, 2, 2]

    # Floor only guards the bulk, which is masked out below anyway.
    dir_floor = 1.0e-6 * gmax if gmax > 0.0 else 1.0
    n_hat = g / np.maximum(gmag, dir_floor)[:, None]
    nHn = np.einsum("nk,nkl,nl->n", n_hat, H, n_hat)

    if DENOM == "gradient":
        G = gmag
        denom_src = "measured |grad phi|"
    else:
        G = phi_np * (1.0 - phi_np) / EPS_USED
        denom_src = "equipartition phi(1-phi)/eps, eps = %.4g m" % EPS_USED
    # phi(1-phi) bottoms out at 0.0253/eps at the band edge vs 0.25/eps at the
    # centre; the floor is well below both and only bites outside the band.
    G = np.maximum(G, 1.0e-4 * (0.25 / EPS_USED))

    kappa = -(L - nHn) / G

    # Axisymmetric r-z: add the azimuthal mode. y is the radius, axis at y = 0.
    # kappa_3D = kappa_planar - n_y/y.
    #
    # On the axis that quotient is 0/0. Flooring the denominator is what you
    # reach for first, but it is wrong by ~13% over the first cells, and a
    # sintering neck sits exactly there. Use the limit instead: n_y = g_y/G, so
    #
    #   lim_{y->0} n_y/y = d n_y/dy = H_yy/G - g_y (dG/dy)/G^2 = H_yy/G
    #
    # since g_y vanishes on the axis by symmetry. H_yy and G are already in
    # hand, so the exact limit costs nothing. Switch to it within one cell.

    # Grid geometry, needed below for the axis switch, the boundary repair and
    # the resolution check. GetExtent covers the structured types;
    # vtkStructuredGrid's GetDimensions() wants an out parameter and raises
    # TypeError when called bare, so don't reach for it.
    bnds = inp.VTKObject.GetBounds()
    try:
        ext = inp.VTKObject.GetExtent()
        nx_p = int(ext[1]) - int(ext[0]) + 1
        ny_p = int(ext[3]) - int(ext[2]) + 1
        nz_p = int(ext[5]) - int(ext[4]) + 1
    except (AttributeError, TypeError):
        nx_p = ny_p = nz_p = 0

    if AXISYM:
        pts = np.asarray(inp.Points, dtype=np.float64)
        y = pts[:, 1]
        ny = n_hat[:, 1]
        ny_rows = ny_p if ny_p > 1 else builtins.max(int(len(y) ** 0.5), 2)
        dy = (bnds[3] - bnds[2]) / (ny_rows - 1)
        dy = builtins.max(dy, 1.0e-30)
        ay = np.abs(y)
        azimuthal = np.where(ay > dy,
                             ny / np.maximum(ay, dy),
                             H[:, 1, 1] / G)
        kappa = kappa - azimuthal

    # Repair the outer TWO layers of points. kappa takes a second derivative,
    # so boundary error reaches two layers in: the outermost layer gets a
    # one-sided first pass, and the next one in is central but reads that bad
    # outer gradient in the second pass. On an analytic sphere both come out
    # ~39% and ~26% low, in planar runs as much as axisymmetric ones -- this is
    # the mesh edge, not the axis. Replicate the first clean layer outward
    # instead of publishing those. Needs the structured layout (point index =
    # i + nx*(j + ny*k)); unstructured input keeps the raw boundary values.
    edge = 2
    if nx_p * ny_p * nz_p == kappa.size and kappa.size > 0:
        block = kappa.reshape(nz_p, ny_p, nx_p).copy()
        for axis, n_axis in ((2, nx_p), (1, ny_p), (0, nz_p)):
            if n_axis <= 2 * edge:      # too thin to have a clean interior
                continue
            lo = [slice(None)] * 3
            hi = [slice(None)] * 3
            for layer in range(edge):
                lo[axis] = layer
                hi[axis] = n_axis - 1 - layer
                src_lo, src_hi = [slice(None)] * 3, [slice(None)] * 3
                src_lo[axis] = edge
                src_hi[axis] = n_axis - 1 - edge
                block[tuple(lo)] = block[tuple(src_lo)]
                block[tuple(hi)] = block[tuple(src_hi)]
        kappa = block.reshape(-1)
    else:
        print("plot_rhovsI: unstructured input, leaving the boundary ring of "
              "kappa uncorrected (it reads low; see the module docstring)")

    # Zero kappa away from the diffuse band, where it is a regularization
    # artifact rather than geometry. Localizer = solver's ice^2*air^2, scaled
    # to peak at 1.
    loc = 16.0 * (phi_np * (1.0 - phi_np)) ** 2
    band = loc >= LOC_FLOOR
    kappa = np.where(band, kappa, 0.0)
    output.PointData.append(kappa, "Curvature")

    d0 = GAMMA_VM_OVER_R / temK                      # capillary length [m]
    rho_vs_eff = rho_vs * (1.0 + d0 * kappa)
    output.PointData.append(rho_vs_eff, "RhoVS_I_GT")

    n_band = int(band.sum())
    print("plot_rhovsI: curvature on %d of %d points (%.1f%% in the band), "
          "denominator = %s%s"
          % (n_band, len(kappa), 100.0 * n_band / builtins.max(len(kappa), 1),
             denom_src,
             ", AXISYMMETRIC (y = radius)" if AXISYM else ", planar"))
    if n_band:
        kb = np.abs(kappa[band])
        rel = np.abs(d0[band] * kappa[band])
        print("  |kappa| in band: median %.4g, p99 %.4g /m  "
              "(radius of curvature %.4g m at the median)"
              % (np.median(kb), np.percentile(kb, 99),
                 1.0 / builtins.max(np.median(kb), 1e-300)))
        print("  |d0*kappa|: median %.3e, max %.3e  "
              "<- relative size of the GT correction"
              % (np.median(rel), rel.max()))

        # Does equipartition actually hold here? The ratio is 1 for a profile
        # at equilibrium. A ratio far from 1 means the two DENOM choices
        # disagree by that factor, and kappa carries the same bias.
        eq = gmag[band] * EPS_USED / np.maximum(
            phi_np[band] * (1.0 - phi_np[band]), 1.0e-12)
        print("  equipartition |grad phi|*eps/(phi(1-phi)): median %.3f, "
              "p10 %.3f, p90 %.3f  (1.0 = at equilibrium)"
              % (np.median(eq), np.percentile(eq, 10), np.percentile(eq, 90)))

        # A curvature radius under one cell is not a measurement, it is the
        # mesh. Grain-contact cusps in a t = 0 packing produce a few of these
        # legitimately; a large fraction means kappa is being read off features
        # the grid cannot resolve.
        if nx_p > 1 and ny_p > 1:
            spacing = builtins.min((bnds[1] - bnds[0]) / (nx_p - 1),
                                   (bnds[3] - bnds[2]) / (ny_p - 1))
            unresolved = int((kb * spacing > 1.0).sum())
            if unresolved:
                print("  %d of %d band points (%.2f%%) have a curvature radius "
                      "below the %.3g m cell -- unresolved features "
                      "(grain-contact cusps at t = 0, typically)"
                      % (unresolved, n_band, 100.0 * unresolved / n_band,
                         spacing))

# --- Supersaturation, the driving force for deposition/sublimation -----------
rhov_name = pick(RHOV_CANDIDATES)
if rhov_name is None:
    print("plot_rhovsI: no vapor-density array (looked for %s), "
          "skipping Supersaturation" % (list(RHOV_CANDIDATES),))
else:
    rho_v = np.asarray(inp.PointData[rhov_name], dtype=np.float64)
    output.PointData.append((rho_v - rho_vs) / rho_vs, "Supersaturation")
    if rho_vs_eff is not None:
        output.PointData.append((rho_v - rho_vs_eff) / rho_vs_eff,
                                "Supersaturation_GT")
'''


def _filter_script(eps, axisym, eps_origin="given", denom=None):
    """FILTER_BODY with a constants header prepended."""
    header = (
        "RHO_AIR         = %r\n"
        "GAMMA_VM_OVER_R = %r\n"
        "LOC_FLOOR       = %r\n"
        "EPS             = %r\n"
        "EPS_ORIGIN      = %r\n"
        "DENOM           = %r\n"
        "AXISYM          = %r\n"
        "TEMP_CANDIDATES = %r\n"
        "RHOV_CANDIDATES = %r\n"
        "PHASE_CANDIDATES = %r\n"
    ) % (RHO_AIR, GAMMA_VM_OVER_R, LOC_FLOOR, eps, eps_origin,
         denom if denom is not None else DENOM, bool(axisym),
         TEMP_ARRAY_CANDIDATES, RHOV_ARRAY_CANDIDATES, PHASE_ARRAY_CANDIDATES)
    return header + FILTER_BODY


def _source_files(src):
    """Filenames behind a reader, or [] for a filter that has none."""
    for prop in ("FileName", "FileNames"):
        try:
            value = getattr(src, prop)
        except AttributeError:
            continue
        if value is None:
            continue
        if isinstance(value, str):
            return [value]
        try:
            return [str(v) for v in value]
        except TypeError:
            return [str(value)]
    return []


def _scan_opts(src):
    """(eps, axisym) parsed from the .opts the run script staged by the data.

    Returns (None, None) when the source is a filter, or when no run directory
    with .opts files turns up within a few levels of the data file. Comment
    lines are skipped so a '# -axisym 1' in prose never counts.
    """
    files = _source_files(src)
    if not files:
        return None, None

    directory = os.path.dirname(os.path.abspath(files[0]))
    eps, axisym = None, None
    eps_re = re.compile(r"^\s*-eps\s+([0-9eE.+-]+)")
    axisym_re = re.compile(r"^\s*-axisym(?:\s+(\S+))?\s*(?:#.*)?$")

    for _ in range(4):                      # <run>/vtkOut/x.vts -> <run> and up
        try:
            names = sorted(os.listdir(directory))
        except OSError:
            break
        for name in names:
            if not name.endswith(".opts"):
                continue
            try:
                with open(os.path.join(directory, name)) as handle:
                    lines = handle.readlines()
            except OSError:
                continue
            for line in lines:
                if eps is None:
                    match = eps_re.match(line)
                    if match:
                        try:
                            eps = float(match.group(1))
                        except ValueError:
                            pass
                if axisym is None:
                    match = axisym_re.match(line)
                    if match:
                        arg = (match.group(1) or "1").lower()
                        axisym = arg in ("1", "true", "yes", "on")
        if eps is not None or axisym is not None:
            break
        parent = os.path.dirname(directory)
        if parent == directory:
            break
        directory = parent

    return eps, axisym


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


def _apply_preset(lut, names):
    """Apply the first colour preset that this ParaView build knows about.

    ParaView 6 renamed "Viridis (matplotlib)" to plain "Viridis", and
    ApplyPreset raises on an unknown name — so try in order and give up
    quietly rather than fail the whole macro over cosmetics.
    """
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
            "plot_rhovsI: no active source. Select a solV_*.vts reader (or any "
            "filter with a Temperature array) in the Pipeline Browser first."
        )

    opts_eps, opts_axisym = _scan_opts(src)
    if "PV_EPS" in os.environ:
        eps, eps_origin = float(os.environ["PV_EPS"]), "the PV_EPS env var"
    elif EPS is not None:
        eps, eps_origin = float(EPS), "the EPS constant in the macro"
    elif opts_eps is not None:
        eps, eps_origin = opts_eps, "-eps in the run's staged .opts"
    else:
        eps, eps_origin = None, "unavailable"
    if "PV_AXISYM" in os.environ:
        axisym = os.environ["PV_AXISYM"].strip().lower() in ("1", "true", "yes", "on")
        axisym_src = "PV_AXISYM"
    elif AXISYM is not None:
        axisym, axisym_src = bool(AXISYM), "the AXISYM constant"
    elif opts_axisym is not None:
        axisym, axisym_src = opts_axisym, "the staged .opts"
    else:
        axisym, axisym_src = False, "nothing found — DEFAULTED to planar"

    print("plot_rhovsI: axisym = %s (from %s)" % (axisym, axisym_src))
    if eps is None:
        print("plot_rhovsI: WARNING -- no eps from PV_EPS, the EPS constant or "
              "any .opts.\n  Falling back to an estimate off the data. eps sets "
              "the curvature scale,\n  so set PV_EPS=<value> or edit EPS if this "
              "run's eps is known.")

    pf = ProgrammableFilter(Input=src)
    pf.OutputDataSetType = "Same as Input"
    pf.Script = _filter_script(eps, axisym, eps_origin)
    pf.RequestInformationScript = ""
    pf.RequestUpdateExtentScript = ""
    RenameSource("RhoVS_I", pf)
    pf.UpdatePipeline()

    view = _render_view()
    display = Show(pf, view)

    # Curvature is the array whose correctness is worth eyeballing first, and
    # it is signed, so it gets a diverging map. Switch in the toolbar for the
    # others. (Change COLOR_BY below if you always want a different default.)
    color_by = "Curvature"
    if color_by not in pf.PointData.keys():
        color_by = "RhoVS_I"
    ColorBy(display, ("POINTS", color_by))
    display.RescaleTransferFunctionToDataRange(True, False)
    display.SetScalarBarVisibility(view, True)
    lut = GetColorTransferFunction(color_by)
    if color_by == "Curvature":
        _apply_preset(lut, ("Cool to Warm (Extended)", "Cool to Warm"))
        # Centre the diverging map on zero. Auto-rescaling to an asymmetric
        # range puts white at (min+max)/2, so a field that is mostly one sign
        # reads as if it straddled zero -- the opposite of what a diverging
        # map is for, and easy to misread as a sign error in the data.
        info = pf.PointData["Curvature"]
        lo, hi = info.GetRange()
        half = builtins_max(abs(lo), abs(hi))
        if half > 0.0:
            lut.RescaleTransferFunction(-half, half)
    else:
        _apply_preset(lut, ("Viridis", "Viridis (matplotlib)", "Cool to Warm"))
    Render(view)

    print("plot_rhovsI: added %s"
          % ", ".join(n for n in ("Curvature", "RhoVS_I", "RhoVS_I_GT",
                                  "Supersaturation", "Supersaturation_GT")
                      if n in pf.PointData.keys()))

    return pf


main()

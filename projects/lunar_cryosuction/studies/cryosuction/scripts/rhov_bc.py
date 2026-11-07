"""Vapor reservoir values for -rhovfix_lo / -rhovfix_hi (simulation plan tool T4).

The solver pins rho_v on a reservoir face to  frac * rho_vs(temp0)
(src/cryosuction_main.c, the IGASetBoundaryValue call near line 1354): frac is
a MULTIPLE of the saturation density at the MEAN temperature temp0, not a
humidity at the face's own temperature. In a temperature gradient the IC is

    T(x) = temp0 + G (x - Lx/2)          (src/initial_conditions.c:252)

so a face at x = 0 sits at temp0 - G Lx/2. This script converts what we mean
physically into the number the solver wants:

  --local-sat        each face holds the saturation density of its own
                     temperature: frac = rho_vs(T_face) / rho_vs(temp0)
  --h-face H         each face holds a humidity H at its own temperature:
                     frac = H rho_vs(T_face) / rho_vs(temp0)
  --h-open H         the collector of the plan: frac = H on the open face
                     (defined at the mid-domain temperature)
  --chi CHI          multiply by the Gibbs-Thomson factor (1 + d0 chi) so the
                     reservoir is in equilibrium with an ice face of curvature
                     chi [1/m] (positive for convex ice; a 2D circle of radius
                     R has chi = 1/R). d0 is the solver's: computed once from
                     temp0 (K&P Eq. 13), as in the run.

The saturation law is the solver's (ASHRAE polynomial over ice, mirrored in
preprocess/mesh_params/physics.py). Its fit is published for -100 C to 0 C;
outside that range the script prints a warning, because a frac computed there
is extrapolated (the lunar targets near 100 K = -173 C are outside).

Every number is printed with 10 significant digits: a grain's own
Gibbs-Thomson excess is ~1e-4 for a 10 um radius, so rounding frac to 1e-6
is safe but a humidity-style 0.99 swamps it (see the solver's -rhovfix_lo help).

    venv_cryosuction/bin/python studies/cryosuction/scripts/rhov_bc.py \
        --T0 -20 --G 100 --Lx 2e-3 --local-sat [--chi 1e5]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from preprocess.mesh_params import physics  # noqa: E402

T_FIT = (-100.0, 0.0)   # ASHRAE ice saturation fit range [C]


def face_temps(T0: float, G: float, Lx: float) -> tuple[float, float]:
    """T at x = 0 and x = Lx for the solver's IC T = T0 + G (x - Lx/2)."""
    return T0 - 0.5 * G * Lx, T0 + 0.5 * G * Lx


def fracs(T0, G, Lx, mode, h=1.0, chi=0.0):
    """(frac_lo, frac_hi, T_lo, T_hi, d0) for one mode."""
    Tlo, Thi = face_temps(T0, G, Lx)
    r0 = physics.rho_vs_sat_kg_per_m3(T0)
    d0 = physics.capillary_length_m(T0)
    gt = 1.0 + d0 * chi
    if mode == "local-sat":
        lo, hi = physics.rho_vs_sat_kg_per_m3(Tlo) / r0, physics.rho_vs_sat_kg_per_m3(Thi) / r0
    elif mode == "h-face":
        lo, hi = h * physics.rho_vs_sat_kg_per_m3(Tlo) / r0, h * physics.rho_vs_sat_kg_per_m3(Thi) / r0
    elif mode == "h-open":
        lo = hi = h
    else:
        raise ValueError(mode)
    return lo * gt, hi * gt, Tlo, Thi, d0


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--T0", type=float, required=True, help="temp0, the mean temperature [C]")
    p.add_argument("--G", type=float, default=0.0, help="grad_temp0 along x [K/m]")
    p.add_argument("--Lx", type=float, default=0.0, help="domain length along x [m]")
    m = p.add_mutually_exclusive_group(required=True)
    m.add_argument("--local-sat", action="store_true")
    m.add_argument("--h-face", type=float)
    m.add_argument("--h-open", type=float)
    p.add_argument("--chi", type=float, default=0.0, help="ice face curvature [1/m], convex > 0")
    a = p.parse_args(argv)

    if a.T0 <= -273.15:
        p.error("--T0 must be above absolute zero")
    if a.Lx < 0:
        p.error("--Lx must be >= 0")
    if a.G != 0.0 and a.Lx == 0.0:
        p.error("--G needs --Lx (the face temperatures depend on the domain length)")
    for name, h in (("--h-face", a.h_face), ("--h-open", a.h_open)):
        if h is not None and h <= 0:
            p.error(f"{name} must be positive")
    mode, h = ("local-sat", 1.0) if a.local_sat else (("h-face", a.h_face) if a.h_face is not None
                                                      else ("h-open", a.h_open))
    lo, hi, Tlo, Thi, d0 = fracs(a.T0, a.G, a.Lx, mode, h, a.chi)
    if Tlo <= -273.15:
        p.error(f"the x = 0 face would be at {Tlo:.2f} C, below absolute zero")
    print(f"# temp0 {a.T0:g} C, G {a.G:g} K/m, Lx {a.Lx:g} m -> T(x=0) {Tlo:.6f} C, T(x=Lx) {Thi:.6f} C")
    print(f"# rho_vs(temp0) {physics.rho_vs_sat_kg_per_m3(a.T0):.10e} kg/m3, d0(temp0) {d0:.6e} m, "
          f"chi {a.chi:g} 1/m -> Gibbs-Thomson factor 1 + d0 chi = {1 + d0 * a.chi:.10f}")
    print(f"# mode {mode}" + (f", h {h:g}" if mode != "local-sat" else ""))
    for T in (a.T0, Tlo, Thi):
        if not (T_FIT[0] <= T <= T_FIT[1]):
            print(f"# WARNING: {T:.2f} C is outside the saturation fit range {T_FIT} C; extrapolated")
            break
    if a.chi != 0.0 and abs(d0 * a.chi) > 0.05:
        print(f"# WARNING: d0 chi = {d0 * a.chi:.3g}; the linear Gibbs-Thomson law assumes d0 chi << 1")
    print(f"-rhovfix_lo {lo:.10f}")
    print(f"-rhovfix_hi {hi:.10f}")
    if mode == "h-open":
        print("# (--h-open sets one value for both; pin only the open face with -rhovfix_faces)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

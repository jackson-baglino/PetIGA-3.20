#!/usr/bin/env python3
"""
generate_calibration_opts.py — inputs for the Phase-4 safety-factor / xi_v
calibration.

Goal: find the LOOSEST (cheapest) discretisation that still reproduces the
observables of the finest run, so the 200-run production matrix is not paid for
at unnecessary resolution.

Two knobs, and they push in opposite directions:

  safety   eps = safety * min(K&P bounds). With R_feat = 2 um the kinetic bound
           binds, so eps = 2*safety um exactly. Cost ~ 1/eps^2.

  xi_v     The model solves (1/xi_v)*d(phi_a*rhov)/dt = div(D*grad rhov)
           + rho_ice*phi_a_t (assembly.c:21), i.e. 1/xi_v AMPLIFIES storage.
           SMALLER xi_v therefore slows vapour relaxation and permits a larger
           dt -- it is the cheap direction, not the expensive one. xi_v = 1 is
           the unscaled physical case.

           Validity is quasi-steadiness: tau_vap = L^2/(xi_v*D_v) must stay far
           below the ice-evolution time. At -20 C over 28 days:
                        1 mm            2 mm
             xi_v=1e-4  527 s (4.6e3x)  2108 s (1.1e3x)
             xi_v=1e-3   53 s (4.6e4x)   211 s (1.1e4x)
           tau_vap grows as L^2, so the 2 mm production domain is 4x less
           quasi-steady than a 1 mm test at the same xi_v. That is the reason
           to check it rather than assume.

           Note K&P Eq. 48 (xi_v <= rho_vs/rho_ice ~ 1e-6) is a DIFFERENT
           condition -- it is what you need to resolve the true vapour
           transient. Deliberately exceeding it is the whole point of the
           temporal scaling; comp_eps.py's warning wording is misleading here.

Three independent observables, so a setting has to survive all of them:

  A. two-grain sintering   neck width vs time      (existing Molaro axisym
                           geometries, neck held fixed across eps by
                           calibrate_neck_geometry.py -- reused, not regenerated)
  B. Ostwald ripening      large-grain growth / small-grain loss, ice mass
                           conservation. Grain radii are set to the STUDY's
                           45/90 um rather than the legacy 9/19 um so eps/R
                           matches production.
  C. packing evolution     SSA(t) and k_eff(t) on a 1 mm packing -- same grain
                           size, gap and porosity as production, quarter the
                           domain so the sweep is affordable.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))
from comp_eps import compute_eps, rho_vs_sat, beta_HK, capillary_length, _RHO_ICE  # noqa: E402

T0 = -20.0                 # the validated anchor temperature
ALPHA_C = 1.341e-2
R_FEAT = 2.0e-6
SAFETYS = {"s025": 0.25, "s050": 0.50, "s075": 0.75, "s100": 1.00}


def eps_for(safety, Lx, Ly, Rave):
    rho_rat = rho_vs_sat(T0) / _RHO_ICE
    beta_uns = beta_HK(T0, ALPHA_C) / rho_rat
    vn = capillary_length(T0) / (beta_uns * R_FEAT)
    return compute_eps(Lx=Lx, Ly=Ly, T0_C=T0, alpha_c=ALPHA_C,
                       Rave=Rave, safety=safety, v_n=vn)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", type=Path, default=PROJECT_ROOT / "inputs")
    ap.add_argument("--packing", type=Path, required=True,
                    help="1 mm calibration packing directory")
    args = ap.parse_args(argv)

    geo = args.out_dir / "geometry"
    geo.mkdir(parents=True, exist_ok=True)

    pk = json.loads((args.packing / "metadata.json").read_text())
    grains = (args.packing / "grains.dat").resolve().relative_to(PROJECT_ROOT)

    print(f"{'tag':>22} {'safety':>7} {'eps[um]':>9} {'Nx':>6} {'Ny':>6} "
          f"{'DOF':>11} {'cores':>6}")

    # ---- Test B: Ostwald ripening, study-scale grains ---------------------
    # Semicircles on the x=0 and x=Lx walls; Ly must clear the larger grain.
    Bx, By, R0, R1 = 5.0e-4, 2.5e-4, 45e-6, 90e-6
    for tag, s in SAFETYS.items():
        p = eps_for(s, Bx, By, R0)
        f = geo / f"calib_ripen_{tag}.opts"
        f.write_text(f"""\
# Ostwald ripening calibration, safety = {s:g}  (generated)
#
# Unequal semicircles centred on the x=0 and x=Lx walls. The large grain should
# grow at the small one's expense; the transfer rate and the ice-mass drift are
# the observables. Radii are the STUDY's grain scale (45 / 90 um), not the
# legacy 9 / 19 um, so eps/R matches production.
#
#   eps/R_small = {p['eps']/R0:.3f}

-ic_type two_ice_grains_boundary
-dim 2
-Lx {Bx:.6e}
-Ly {By:.6e}
-Nx {p['Nx']}
-Ny {p['Ny']}
-RCice0 {R0:.6e}                # small grain, x = 0
-RCice1 {R1:.6e}                # large grain, x = Lx
-eps {p['eps']:.6e}             # safety={s:g}, binding: {p['binding']}
-eps_valid_temp {T0:g}
-periodic 0
-delt_t 1.0e-4
""")
        dof = 3 * p["Nx"] * p["Ny"]
        print(f"{'calib_ripen_'+tag:>22} {s:>7.2f} {p['eps']*1e6:>9.3f} "
              f"{p['Nx']:>6} {p['Ny']:>6} {dof:>11,} {math.ceil(dof/50000):>6}")

    # ---- Test C: packing evolution ---------------------------------------
    for tag, s in SAFETYS.items():
        p = eps_for(s, pk["Lx"], pk["Ly"], 45e-6)
        f = geo / f"calib_pack_{tag}.opts"
        f.write_text(f"""\
# Packing calibration, safety = {s:g}  (generated)
#
# Same grain size, contact gap and porosity as production, on a 1 mm domain so
# the sweep is affordable. Observables: SSA(t) from SSA_evo.dat and k_eff(t)
# from effective_thermal_cond.
#
#   porosity = {pk['porosity_achieved']:.4f}   grains = {pk['n_grains']}
#   contact gap = {pk['contact_gap_m']*1e6:.2f} um   Z = {pk['coordination_number']:.2f}
#   diffuse band 6*eps = {6*p['eps']*1e6:.2f} um
#   NOTE tau_vap ~ L^2, so this 1 mm test is 4x MORE quasi-steady than the
#   2 mm production domain at the same xi_v -- read the xi_v result as a lower
#   bound on what production needs.

-ic_type multi_grains_file
-grains_file {grains}
-ic_grain_union 1
-dim 2
-Lx {pk['Lx']:.6e}
-Ly {pk['Ly']:.6e}
-Nx {p['Nx']}
-Ny {p['Ny']}
-eps {p['eps']:.6e}             # safety={s:g}, binding: {p['binding']}
-eps_valid_temp {T0:g}
-delt_t 1.0e-4
""")
        dof = 3 * p["Nx"] * p["Ny"]
        print(f"{'calib_pack_'+tag:>22} {s:>7.2f} {p['eps']*1e6:>9.3f} "
              f"{p['Nx']:>6} {p['Ny']:>6} {dof:>11,} {math.ceil(dof/50000):>6}")

    # ---- Shared experiment: shorter than production ----------------------
    exp = args.out_dir / "experiment"
    exp.mkdir(parents=True, exist_ok=True)
    p = eps_for(0.5, 1e-3, 1e-3, 45e-6)
    (exp / "calib_T-20_7d.opts").write_text(f"""\
# Calibration experiment: -20 C, 7 days  (generated)
#
# Shorter than the 28-day production runs: the convergence question is whether
# the observables AGREE between discretisations, which is visible well before
# the microstructure fully coarsens.
#
# xi_v is NOT set here. Sweep it on the command line so one experiment file
# serves every arm:
#     ./scripts/HPC/submit_enceladus.sh <geom> calib_T-20_h1.00_7d <tag> -- -xi_v 1e-2
# Unset, the solver default (1e-3) applies.

-temp {T0:g}
-humidity 1.00
-t_final 6.048e5                  # 7 days

-beta_sub0 {p['beta_uns']:.6e}
-d0_sub0 {p['d0']:.6e}

-dtCFL 1
-dtCFL_dphimax 0.2
-dtmax 2.0e2
""")
    print(f"\n  experiment -> {exp/'calib_T-20_7d.opts'}")
    print(f"  geometry   -> {geo}/calib_*.opts")
    return 0


if __name__ == "__main__":
    sys.exit(main())

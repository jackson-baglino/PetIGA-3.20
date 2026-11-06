"""Order-of-magnitude numbers for the paper-2 simulation plan (2026-11-06).

Every number quoted in the simulation proposal's "scales" boxes comes from
this script. Constants come from preprocess/mesh_params (which mirrors the C
solver), so rho_vs, D_v(T), d0 and beta0 are the solver's values at T0.
These are estimates for choosing parameters and run lengths, not results.

Sections
  1. B1 series-law slope: predicted with D_v at T0 and with D_v at 0 C, against
     the measured 5.73e12 s/m^2 (studies/cryosuction/RUNS.md, row 0c-B1).
  2. Ripening pair: Gibbs-Thomson drive, the gradient that cancels it, and a
     rough vanishing time (2D quasi-steady diffusion in series with kinetics).
  3. Gradient-driven migration speed of ice (and K&P Eq. 45 resolution number
     Pi = beta0 v eps / d0 at --eps; the design limit is 0.5), v ~ D_v * d(rho_vs)/dT * G / rho_ice
     (no conductivity-contrast factor; ice conducts ~100x better than air, so
     the real local gradient in the air differs by an O(1) factor).
  4. Harvest time of a filled channel, front moving from l0 to l1 under
     suction sigma: t = [beta0 (l1 - l0) + K (l1^2 - l0^2) / 2] * f / sigma,
     K = rho_ice / (D_v rho_vs), f = ice fraction of the pore cross-section.

    venv_cryosuction/bin/python studies/cryosuction/scripts/plan_estimates.py
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from preprocess.mesh_params import physics  # noqa: E402
from preprocess.mesh_params.constants import CONSTANTS as C  # noqa: E402

DAY = 86400.0


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--T0", type=float, default=-20.0, help="temperature [C]")
    p.add_argument("--alpha", type=float, default=1e-2, help="alpha_c")
    p.add_argument("--b1-slope", type=float, default=5.73e12, help="measured B1 slope of 1/v vs l [s/m^2]")
    p.add_argument("--b1-sigma", type=float, default=0.01, help="B1 suction 1 - h")
    p.add_argument("--pair", type=float, nargs=3, default=[20e-6, 50e-6, 200e-6],
                   help="ripening pair: R_small R_large center distance [m]")
    p.add_argument("--G", type=float, nargs="+", default=[1, 10, 100, 1000], help="gradients [K/m]")
    p.add_argument("--harvest", type=float, nargs="+", default=[0.99, 0.9, 0.5], help="open-end h values")
    p.add_argument("--l0", type=float, default=1e-4, help="harvest: start distance [m]")
    p.add_argument("--l1", type=float, nargs="+", default=[0.85e-3, 2e-3, 5e-3], help="harvest: end distances [m]")
    p.add_argument("--fill", type=float, default=1.0, help="harvest: ice fraction of the pore cross-section")
    p.add_argument("--eps", type=float, default=8.584e-7, help="eps for the resolution number Pi [m]")
    a = p.parse_args()

    T = a.T0
    rho_vs = physics.rho_vs_sat_kg_per_m3(T)
    D = physics.dv_m2_per_s(T)
    D0 = C.dv0_m2_per_s
    d0 = physics.capillary_length_m(T)
    beta0 = physics.beta_hk_s_per_m(T, a.alpha) * C.rho_ice_kg_per_m3 / rho_vs  # unscaled K&P beta0, as front_speed.py
    rho_i = C.rho_ice_kg_per_m3
    dT = 0.01
    s = (math.log(physics.rho_vs_sat_kg_per_m3(T + dT)) - math.log(physics.rho_vs_sat_kg_per_m3(T - dT))) / (2 * dT)
    K = rho_i / (D * rho_vs)
    print("# inputs")
    print(f"T0 = {T} C, alpha_c = {a.alpha}")
    print(f"rho_vs = {rho_vs:.5e} kg/m3, D_v(T0) = {D:.5e} m2/s, D_v(0 C) = {D0:.5e} m2/s, rho_ice = {rho_i}")
    print(f"d0 = {d0:.5e} m, beta0 = {beta0:.5e} s/m, s = dln(rho_vs)/dT = {s:.5e} 1/K, K = {K:.5e} s/m2")

    print("\n# 1. B1 series-law slope, 1/v = (beta0 + K l)/sigma")
    sl_T = K / a.b1_sigma
    sl_0 = rho_i / (D0 * rho_vs) / a.b1_sigma
    print(f"slope with D_v(T0) = {sl_T:.4e} s/m2 -> measured/predicted = {a.b1_slope / sl_T:.4f}")
    print(f"slope with D_v(0C) = {sl_0:.4e} s/m2 -> measured/predicted = {a.b1_slope / sl_0:.4f}")
    print(f"intercept beta0/sigma = {beta0 / a.b1_sigma:.4e} s/m;"
          f" a 4.2e7 intercept equals a shift in l of {(beta0 / a.b1_sigma - 4.2e7) / sl_T * 1e6:.2f} um")

    print("\n# 2. Ripening pair (2D discs)")
    Rs, Rl, L = a.pair
    dsig = d0 * (1 / Rs - 1 / Rl)
    gap = L - Rs - Rl
    print(f"R_small = {Rs:.2e}, R_large = {Rl:.2e}, center distance {L:.2e}, gap {gap:.2e} m")
    print(f"drive d0 (1/Rs - 1/Rl) = {dsig:.4e}; mean-curvature h for the pair = 1 + {d0 * 0.5 * (1 / Rs + 1 / Rl):.4e}")
    print(f"gradient that cancels it over the center distance: G* = {dsig / (s * L):.3f} K/m"
          f" (temperature difference {dsig / s * 1e3:.4f} mK)")
    v_dif = D * rho_vs * dsig / (rho_i * Rs * math.log(L / Rs))
    v_kin = dsig / beta0
    v = 1 / (1 / v_dif + 1 / v_kin)
    print(f"small-grain shrink speed: diffusion {v_dif:.3e}, kinetics {v_kin:.3e}, series {v:.3e} m/s"
          f" = {v * DAY * 1e6:.3f} um/day; Rs/v = {Rs / v / DAY:.1f} days (rough vanishing time)")

    print("\n# 3. Gradient-driven migration speed v ~ D rho_vs s G / rho_ice")
    for G in a.G:
        vm = D * rho_vs * s * G / rho_i
        print(f"G = {G:7.1f} K/m: v = {vm:.3e} m/s = {vm * DAY * 1e6:8.3f} um/day;"
              f" saturation change over 100 um = {s * G * 1e-4:.3e};"
              f" Pi = beta0 v eps/d0 = {beta0 * vm * a.eps / d0:.3f} at eps {a.eps:.3e} m")

    print(f"\n# 4. Harvest time, front from l0 = {a.l0:.2e} m, ice fill fraction {a.fill}")
    for h in a.harvest:
        sig = 1 - h
        for l1 in a.l1:
            t = (beta0 * (l1 - a.l0) + K * (l1 ** 2 - a.l0 ** 2) / 2) * a.fill / sig
            v_start = sig / (beta0 + K * a.l0)
            print(f"h = {h:5.3f}  l1 = {l1:.2e} m: t = {t / DAY:9.1f} days;"
                  f" starting front speed {v_start:.3e} m/s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

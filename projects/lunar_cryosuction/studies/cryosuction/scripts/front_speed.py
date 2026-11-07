"""Expected front speed of a sublimating ice front for comp_eps.py --vn.

Series law (kinetic attachment in series with vapor diffusion along a path of
length l from the front to the open end), as used in the velocity study
(studies/contact_angle/velocity_plan_2026-10-01/README.md, where it matched
the solver to 1-3 %):

    v = |sigma| / (beta0 + K*l),   K = rho_ice / (D_v * rho_vs)

beta0 is the UNSCALED Hertz-Knudsen coefficient (K&P beta0, the solver's
-beta_sub0): beta0 = (1/alpha_c) sqrt(2 pi m/(k_B T)) * rho_ice/rho_vs.
sigma = h - 1 at the open end. The largest v (shortest l) sets the K&P Eq. 45
bound, so pass the smallest front-to-open-end distance of the geometry.
Constants come from preprocess/mesh_params (mirrors the C solver). Prints its
inputs and the comp_eps.py command to run.

    venv_cryosuction/bin/python studies/cryosuction/scripts/front_speed.py \
        --T0 -20 --alpha 1e-2 --h 0.99 --l 9.4e-5 [--Lx .. --Ly .. --Rave ..]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from preprocess.mesh_params import physics  # noqa: E402
from preprocess.mesh_params.constants import CONSTANTS as C  # noqa: E402


def speed(T0: float, alpha: float, h: float, l_m: float) -> tuple[float, float, float]:
    rvs = physics.rho_vs_sat_kg_per_m3(T0)
    beta0 = physics.beta_hk_s_per_m(T0, alpha) * C.rho_ice_kg_per_m3 / rvs
    K = C.rho_ice_kg_per_m3 / (physics.dv_m2_per_s(T0) * rvs)
    return abs(h - 1.0) / (beta0 + K * l_m), beta0, K


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--T0", type=float, default=-20.0)
    p.add_argument("--alpha", type=float, required=True)
    p.add_argument("--h", type=float, required=True, help="open-end humidity (multiple of rho_vs)")
    p.add_argument("--l", type=float, required=True, help="smallest front-to-open-end distance [m]")
    p.add_argument("--Lx", type=float)
    p.add_argument("--Ly", type=float)
    p.add_argument("--Rave", type=float)
    a = p.parse_args()
    # Same input range as comp_eps (preprocess/mesh_params/tool.py), so the
    # two never disagree about what is a valid case.
    if not (-150.0 < a.T0 < -0.001):
        p.error(f"--T0 must be in (-150, -0.001) C, where the ice formulas apply (got {a.T0})")
    if not (0.0 < a.alpha <= 1.0):
        p.error(f"--alpha must be in (0, 1] (got {a.alpha})")
    if not (a.h > 0.0):
        p.error(f"--h must be positive (got {a.h})")
    if not (a.l >= 0.0):
        p.error(f"--l must be >= 0 (got {a.l})")
    if a.T0 < -100.0:
        print(f"# WARNING: {a.T0} C is below the -100 C limit of the saturation fit; "
              "rho_vs is extrapolated, and D_v is the air value", file=sys.stderr)
    v, beta0, K = speed(a.T0, a.alpha, a.h, a.l)
    print(f"# T0 {a.T0} C, alpha_c {a.alpha:g}, h {a.h}, l {a.l:g} m")
    print(f"beta0 = {beta0:.4e} s/m, K = {K:.4e} s/m^2, K*l = {K * a.l:.4e} s/m")
    print(f"v_n = {v:.4e} m/s   (kinetic-limit upper bound |sigma|/beta0 = {abs(a.h - 1) / beta0:.4e})")
    if a.Lx:
        dims = f"--Lx {a.Lx:g}" + (f" --Ly {a.Ly:g}" if a.Ly else "")
        rave = f" --Rave {a.Rave:g}" if a.Rave else ""
        print(f"comp_eps: venv_cryosuction/bin/python preprocess/comp_eps.py {dims}{rave} "
              f"--T0 {a.T0:g} --alpha {a.alpha:g} --vn {v:.4e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""Compare vapor deviations from saturation: imposed humidity vs curvature.

Prints every input it used. Formulas and constants come from
preprocess/mesh_params (which mirrors the C solver), so the numbers match the
solver's d0 (Kaempfer and Plapp 2009, Eq. 13) and rho_vs (ASHRAE polynomial).

    sigma_h   = h - 1                     imposed reservoir undersaturation
    sigma_GT  = d0 * kappa                Gibbs-Thomson (Kelvin, linearized)
        grain:    kappa = 1/R (2D disc)  or  2/R (3D sphere)
        meniscus: kappa = 2 cos(theta) / H (2D channel of width H)
    K&P Eq. 48 check: rho_vs / (xi_v * rho_ice) must be << 1.
    Quasi-steady check: vapor relaxation time L^2 / (xi_v * D_v).

Run from the project root:
    venv_cryosuction/bin/python studies/cryosuction/scripts/supersat_compare.py \
        --temp -20 --radius 50e-6 --width 125e-6 --theta 30 60 90 120 150 \
        --humidity 0.998 0.995 0.99 0.9 0.5 --xi-v 1e-2 1e-3 1e-4 1e-5 --length 375e-6
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from preprocess.mesh_params import physics  # noqa: E402
from preprocess.mesh_params.constants import CONSTANTS as C  # noqa: E402


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--temp", type=float, default=-20.0, help="temperature [C]")
    p.add_argument("--radius", type=float, default=50e-6, help="grain radius of curvature [m]")
    p.add_argument("--width", type=float, nargs="+", default=[125e-6], help="channel widths H [m]")
    p.add_argument("--theta", type=float, nargs="+", default=[30, 60, 90, 120, 150], help="contact angles [deg]")
    p.add_argument("--humidity", type=float, nargs="+", default=[0.998, 0.995, 0.99, 0.9, 0.5])
    p.add_argument("--xi-v", type=float, nargs="+", default=[1e-2, 1e-3, 1e-4, 1e-5])
    p.add_argument("--length", type=float, default=375e-6, help="domain length for the relaxation time [m]")
    a = p.parse_args()

    d0 = physics.capillary_length_m(a.temp)
    rvs = physics.rho_vs_sat_kg_per_m3(a.temp)
    dv = physics.dv_m2_per_s(a.temp)

    print("# inputs")
    print(f"T = {a.temp} C, R = {a.radius:g} m, H = {a.width}, theta = {a.theta} deg")
    print(f"h = {a.humidity}, xi_v = {a.xi_v}, L = {a.length:g} m")
    print(f"d0 = {d0:.5e} m, rho_vs = {rvs:.5e} kg/m3, D_v = {dv:.5e} m2/s, rho_ice = {C.rho_ice_kg_per_m3}")

    print("\n# curvature (Gibbs-Thomson), sigma = d0*kappa")
    s2 = d0 / a.radius
    s3 = 2.0 * d0 / a.radius
    print(f"grain R = {a.radius:g} m: 2D disc {s2:.3e}, 3D sphere {s3:.3e}")
    for H in a.width:
        for th in a.theta:
            sm = d0 * 2.0 * math.cos(math.radians(th)) / H
            print(f"meniscus H = {H:g} m, theta = {th:5.1f}: {sm:+.3e}")

    print("\n# imposed undersaturation sigma_h = h - 1, and ratio to the 3D grain term")
    for h in a.humidity:
        sh = h - 1.0
        print(f"h = {h:6.3f}: sigma_h = {sh:+.3e}, |sigma_h|/sigma_GT(3D) = {abs(sh) / s3:9.1f}")

    print("\n# K&P scaling checks")
    for xi in a.xi_v:
        ratio = rvs / (xi * C.rho_ice_kg_per_m3)
        t_rel = a.length**2 / (xi * dv)
        print(f"xi_v = {xi:.0e}: rho_vs/(xi rho_ice) = {ratio:.3e}, L^2/(xi D_v) = {t_rel:.3e} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

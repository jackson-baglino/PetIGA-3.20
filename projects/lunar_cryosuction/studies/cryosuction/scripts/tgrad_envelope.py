"""Back-of-envelope numbers: how an applied temperature gradient competes with
open-end suction when ice is harvested as vapor.

Geometry (1D): the open end (reservoir or vacuum, vapor density
rho_open = h_open * rho_vs(T_open)) is at the top; a dry layer of thickness l
separates it from the ice front at T_f; ice continues below the front. Heat
applied at the top gives T rising toward the top with gradient G [K/m] (warm
top, cold depth).

1. Warm open end against suction. Vapor leaves the front only if
   rho_open < rho_vs(T_f), i.e. h_open < h* = rho_vs(T_f) / rho_vs(T_open).
2. Cold-trapping below the front. The saturated ice below the front carries a
   vapor gradient d(rho_vs)/dz = rho_vs * (dln rho_vs/dT) * G, so vapor also
   diffuses DOWN into the colder ice and redeposits there. With the same
   transport coefficient D in both directions (it cancels):
       F_down / F_up = l / l_T,   l_T = 1 / (G * dln rho_vs/dT)   (h_open = 0)
   Once the dry layer is thicker than l_T, more vapor is pumped down than is
   harvested.
3. Equivalent saturation drive of a temperature difference: dT changes the
   saturation level by about dT * dln rho_vs/dT; compared with the suction
   1 - h_open in our channels.
4. Rate gain from warming the ice itself: rho_vs(T + 10 K) / rho_vs(T).

Vapor pressure:
  * -20 C analog: the solver's RhoVS_I (src/material_properties.c), through
    postprocess/pplib.rho_vs, with rho_air 1.341 kg/m^3.
  * Lunar temperatures: Murphy & Koop (2005), QJRMS 131, 1539, Eq. 7 (ice,
    stated valid T > 110 K; values below 110 K are extrapolated), with
    rho = p M / (R T) (ideal gas, no air).
Transport law: Fick diffusion of rho_v (as in the model). In lunar vacuum the
pores are in the Knudsen regime; there the driving quantity is p/sqrt(T), which
changes dln/dT by 1/(2T) against L/(R_v T^2): printed below as a correction.

    venv_cryosuction/bin/python studies/cryosuction/scripts/tgrad_envelope.py
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "postprocess"))
from pplib import rho_vs as rho_vs_solver  # noqa: E402

R = 8.314462618   # J/(mol K)
M_W = 0.01801528  # kg/mol
RHO_AIR = 1.341   # kg/m^3, inputs/solver.opts default used by vapor_movie


def p_mk(T: float) -> float:
    """Murphy & Koop 2005 Eq. 7, vapor pressure over hexagonal ice [Pa]."""
    return math.exp(9.550426 - 5723.265 / T + 3.53068 * math.log(T) - 0.00728332 * T)


def rho_mk(T: float) -> float:
    return p_mk(T) * M_W / (R * T)


def rho_solver(T: float) -> float:
    return float(rho_vs_solver(T - 273.15, RHO_AIR))


def dln(f, T: float, h: float = 1e-3) -> float:
    return (math.log(f(T + h)) - math.log(f(T - h))) / (2 * h)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--T-analog", type=float, default=253.15, help="analog ice temperature [K] (-20 C)")
    ap.add_argument("--T-lunar", type=float, nargs="+", default=[100.0, 120.0, 150.0, 200.0])
    ap.add_argument("--G", type=float, nargs="+", default=[10.0, 100.0, 1e3, 1e4], help="gradients [K/m]")
    ap.add_argument("--dT-open", type=float, nargs="+", default=[0.1, 1.0, 5.0, 10.0, 20.0],
                    help="open end warmer than the ice by [K]")
    ap.add_argument("--Lx", type=float, default=8.5e-4, help="simulation channel length [m]")
    ap.add_argument("--h-open", type=float, default=0.99, help="open-end saturation in the runs")
    a = ap.parse_args()
    if not (a.T_analog > 0.0) or any(not (t > 0.0) for t in a.T_lunar):
        ap.error("temperatures are in kelvin and must be > 0")
    if any(not (g > 0.0) for g in a.G):
        ap.error("--G values must be > 0 (the lengths divide by G)")
    if not (0.0 < a.h_open < 1.0):
        ap.error("--h-open must be in (0, 1): at h >= 1 there is no suction")
    if not (a.Lx > 0.0):
        ap.error("--Lx must be positive")

    print("# inputs:", vars(a))
    print(f"# check Murphy-Koop at 273.16 K: {p_mk(273.16):.2f} Pa (expected ~611.7)")
    print(f"# check solver rho_vs at -20 C: {rho_solver(253.15)*1e3:.4f} g/m^3; "
          f"Murphy-Koop ideal gas: {rho_mk(253.15)*1e3:.4f} g/m^3")

    Ta = a.T_analog
    print(f"\n== 1. Warm open end against suction (ice at {Ta - 273.15:.1f} C, solver rho_vs) ==")
    print("  open end warmer by dT [K] -> h* (vapor leaves the ice only if h_open < h*)")
    for dT in a.dT_open:
        print(f"    dT {dT:5.1f}   h* = {rho_solver(Ta) / rho_solver(Ta + dT):.4f}")

    print("\n== 2. Cold-trapping length l_T = 1/(G dln rho_vs/dT); F_down/F_up = l/l_T ==")
    rows = [("analog (solver)", Ta, rho_solver)] + [(f"lunar (M&K{'*' if T < 110 else ''})", T, rho_mk)
                                                    for T in a.T_lunar]
    print("  " + f"{'case':18s} {'T [K]':>6s} {'dln/dT [1/K]':>13s} {'1/(..) [K]':>10s} {'Knudsen corr.':>13s}  "
          + "  ".join(f"l_T@G={g:g}" for g in a.G))
    for name, T, f in rows:
        k = dln(f, T)
        corr = 1 / (2 * T) / k
        lts = "  ".join(f"{1 / (g * k):>11.3e}" for g in a.G)
        print(f"  {name:18s} {T:6.1f} {k:13.4e} {1 / k:10.3f} {corr:12.2%}   {lts}   [m]")
    print("  (* below the 110 K validity limit of M&K Eq. 7: extrapolated)")

    print(f"\n== 3. Equivalent drive across the simulation channel (Lx = {a.Lx*1e6:.0f} um, analog) ==")
    k = dln(rho_solver, Ta)
    print(f"  suction drive 1 - h_open = {1 - a.h_open:.3g}")
    for g in a.G:
        dT = g * a.Lx
        print(f"    G {g:8.0f} K/m: dT across channel {dT:.3e} K -> saturation change {dT * k:.3e} "
              f"({dT * k / (1 - a.h_open):.2f} x the suction); l/l_T = {a.Lx * g * k:.3e}")

    print("\n== 4. Warming the ice by 10 K multiplies rho_vs (and the harvest rate at fixed geometry) by ==")
    print(f"    analog {Ta:.0f} K: {rho_solver(Ta + 10) / rho_solver(Ta):.2f}")
    for T in a.T_lunar:
        print(f"    lunar  {T:.0f} K: {rho_mk(T + 10) / rho_mk(T):.3g}   (rho_vs({T:.0f} K) = {rho_mk(T):.3e} kg/m^3)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

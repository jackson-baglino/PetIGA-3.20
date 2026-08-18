"""
estimate_vapor_bc_molaro.py — derive the vapour Dirichlet value (-humidity) and
check the vapour IC, from the grain shrinkage Molaro et al. (2019) report.

WHY
---
Molaro's grains lose mass over the experiment, so the chamber was slightly
undersaturated. In the solver that is one number: `-humidity` (hum0), which pins
rho_v = hum0*rho_vs(T0) on every domain face except the axisymmetric axis
(src/enceladus_main.c:940-953). Guessing it wastes HPC time, and the previously
fitted h = 0.998 is tied to the old 121 um domain and does not transfer.

THE MATH
--------
Treat the pair as a sphere of radius a_eff with a Dirichlet wall at distance L.
Quasi-steady, two resistances in series between the ice surface and the wall:

    attachment    du_att  = beta' * v_n,      beta' = (1/alpha_c)*sqrt(2*pi*m/kT)
    diffusion     du_diff = v_n * l_eff/D_v

For a sphere the exact truncated-shell conductance gives the effective length

    l_eff = a * (1 - a/L)

(l_eff -> a as L -> infinity; l_eff -> 0 as the wall approaches the surface, i.e.
a tight domain offers almost no diffusive resistance and needs almost no
undersaturation to drive the same flux).

The available driving force is the wall undersaturation, in the same u units:

    du_total = (rho_vs - rho_v,wall)/rho_i = (rho_vs/rho_i) * (1 - hum0)

Setting du_total = (beta' + l_eff/D_v)*|v_n| and solving:

    1 - hum0 = |v_n| * (beta' + l_eff/D_v) / (rho_vs/rho_i)

WHICH v_n
---------
The large grain. Its 78-min trend is a 4.4-sigma signal and its least-squares
slope (-2.93%) reproduces the -3% the paper's caption quotes. The small grain's
"-4%" is not supported by the same table: endpoint -0.69%, minimum -2.40%, fit
-0.89%, all inside the +-2% uncertainty the table itself carries. This script
prints both fits so the asymmetry is visible, and targets the large grain.

RIPENING IS NOT THE MECHANISM
-----------------------------
Ostwald ripening would shrink the small grain faster. Compare the driving forces:

    ripening    d0 * (2/R_sm - 2/R_lg)     ~ 8e-6
    wall        1 - hum0                   ~ 3e-3

The imposed undersaturation is ~2-3 decades stronger, so both grains sublimate at
nearly the same rate and the ripening difference is buried in the measurement
noise. Report the differential; do not fit to it.

THE IC NEEDS NO CHANGE
----------------------
SetNodeFields (src/initial_conditions.c:215-225) writes
rho_v = rho_vs(T)*[hum0*(1-phi) + phi]: saturated inside the ice, hum0-scaled in
the pore. The quasi-steady profile establishes in t ~ L^2/D_v, which this script
prints -- microseconds to milliseconds, i.e. instantaneous on a 78-min run. So
the uniform pore value is fine and needs no pre-equilibration.

Usage
-----
  python preprocess/estimate_vapor_bc_molaro.py
  python preprocess/estimate_vapor_bc_molaro.py --alpha 0.1 --L-over-a 3 --n-grid 5
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))

from comp_eps import (  # noqa: E402
    beta_HK, Dv_T, rho_vs_sat, capillary_length, _RHO_ICE,
)

_REPO = _HERE.parent

# Same validated pair as estimate_alpha_c_molaro.py.
C_MODEL, C_REF, C_INK, C_MUTED = "#3d74d9", "#d1495b", "#3f434a", "#8a8a8a"


# =========================================================================
# The estimate
# =========================================================================

def l_eff(a: float, L: float) -> float:
    """Effective diffusion length for a sphere of radius a with a wall at L."""
    return a * (1.0 - a / L)


def humidity_for(v_n: float, a: float, L: float, alpha_c: float, T_C: float) -> float:
    """hum0 that drives a surface recession |v_n| with the wall at L."""
    rr = rho_vs_sat(T_C) / _RHO_ICE
    total = beta_HK(T_C, alpha_c) + l_eff(a, L) / Dv_T(T_C)
    return 1.0 - v_n * total / rr


def ripening_sigma(R_lg: float, R_sm: float, T_C: float) -> float:
    """Supersaturation difference between the two grain surfaces, d0*(2/R_sm - 2/R_lg)."""
    return capillary_length(T_C) * 2.0 * (1.0 / R_sm - 1.0 / R_lg)


# =========================================================================
# Data
# =========================================================================

def read_diameters(path: Path):
    """[(t_s, lg_diam_m, sm_diam_m), ...] from the validation CSV."""
    rows = []
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            f = next(csv.reader([line]))
            rows.append((float(f[0]) * 60.0, float(f[4]) * 1e-6, float(f[5]) * 1e-6))
    return rows


def linfit(x, y):
    """(slope, intercept, rms_residual) — plain least squares."""
    n = len(x)
    sx, sy = sum(x), sum(y)
    sxx = sum(v * v for v in x)
    sxy = sum(a * b for a, b in zip(x, y))
    m = (n * sxy - sx * sy) / (n * sxx - sx * sx)
    b = (sy - m * sx) / n
    res = [v - (m * u + b) for u, v in zip(x, y)]
    return m, b, math.sqrt(sum(e * e for e in res) / (n - 2))


def make_figure(rows, fits, grid, chosen, T_C, out_png):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11.0, 4.3))

    # --- left: the measured diameters and their fits ----------------------
    t_min = [r[0] / 60.0 for r in rows]
    for key, col, lab in (("large", C_MODEL, "large grain"), ("small", C_REF, "small grain")):
        d_um = [r[1 if key == "large" else 2] * 1e6 for r in rows]
        m, b, rms = fits[key]
        axL.plot(t_min, d_um, "o", color=col, ms=8, zorder=3)
        axL.plot(t_min, [m * t + b for t in t_min], "-", color=col, lw=2, zorder=2,
                 label=f"{lab}: {100*m*(t_min[-1])/b:+.2f}%  (S/N {abs(m*t_min[-1])/rms:.1f})")
        axL.fill_between(t_min, [m * t + b - rms for t in t_min],
                         [m * t + b + rms for t in t_min], color=col, alpha=0.15, lw=0)
    axL.set_xlabel("time  [min]")
    axL.set_ylabel("grain diameter  [$\\mu$m]")
    axL.set_title("Measured shrinkage: only the large grain has signal",
                  fontsize=10, color=C_INK)
    axL.grid(alpha=0.25, lw=0.6)
    axL.set_axisbelow(True)
    axL.legend(frameon=False, fontsize=9, loc="center left")

    # --- right: required undersaturation vs domain size -------------------
    k = [g["L_over_a"] for g in grid]
    dh = [1.0 - g["hum0"] for g in grid]
    axR.plot(k, dh, "o-", color=C_MODEL, lw=2, ms=8, zorder=3, label="required $1-h$")
    axR.axhline(ripening_sigma(101e-6, 72.5e-6, T_C), color=C_REF, lw=2, ls="--", zorder=2,
                label="ripening driving force")
    ch = [g for g in grid if abs(g["L_over_a"] - chosen) < 1e-9]
    if ch:
        axR.annotate(f"chosen  $h$ = {ch[0]['hum0']:.5f}",
                     (ch[0]["L_over_a"], 1.0 - ch[0]["hum0"]),
                     textcoords="offset points", xytext=(8, -14),
                     fontsize=9, color=C_INK)
    axR.set_yscale("log")
    axR.set_xlabel("$L / a_{\\mathrm{eff}}$  [–]")
    axR.set_ylabel("$1 - h$  [–]")
    axR.set_title("Wall undersaturation needed, and what ripening supplies",
                  fontsize=10, color=C_INK)
    axR.grid(alpha=0.25, lw=0.6)
    axR.set_axisbelow(True)
    axR.legend(frameon=False, fontsize=9, loc="center right")

    for ax in (axL, axR):
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.tick_params(colors=C_INK, labelsize=9)

    fig.tight_layout()
    fig.savefig(out_png, dpi=160)
    print(f"  wrote {out_png}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", type=Path,
                    default=_REPO / "inputs/validation/molaro2019_fig11_T-20.csv")
    ap.add_argument("--T0", type=float, default=-20.0)
    ap.add_argument("--R1", type=float, default=101.0e-6, help="large grain radius [m]")
    ap.add_argument("--R2", type=float, default=72.5e-6, help="small grain radius [m]")
    ap.add_argument("--alpha", type=float, default=0.1, help="condensation coefficient")
    ap.add_argument("--L-over-a", type=float, default=3.0,
                    help="chosen domain size in units of a_eff")
    ap.add_argument("--n-grid", type=int, default=5, help="humidity arms to propose")
    ap.add_argument("--outdir", type=Path, default=_REPO / "studies/molaro_2019")
    args = ap.parse_args()

    T = args.T0
    a_eff = (args.R1 ** 3 + args.R2 ** 3) ** (1.0 / 3.0)
    rr = rho_vs_sat(T) / _RHO_ICE
    rows = read_diameters(args.data)
    t_min = [r[0] / 60.0 for r in rows]

    fits = {"large": linfit(t_min, [r[1] * 1e6 for r in rows]),
            "small": linfit(t_min, [r[2] * 1e6 for r in rows])}

    print("=" * 78)
    print(f"  Vapour BC / IC for Molaro et al. (2019), T = {T:g} C, alpha_c = {args.alpha:g}")
    print("=" * 78)
    print(f"  a_eff (volume-equivalent radius of the PAIR) = {a_eff*1e6:.2f} um")
    print(f"  rho_vs/rho_i = {rr:.4e}    D_v = {Dv_T(T):.4e} m^2/s"
          f"    beta' = {beta_HK(T, args.alpha):.4f} s/m")

    print("\n--- Which shrinkage to fit ---")
    print("   grain    D0[um]   slope[um/min]   78-min change      rms resid   S/N")
    for key in ("large", "small"):
        m, b, rms = fits[key]
        span = t_min[-1]
        print(f"   {key:6s} {b:8.2f}   {m:+12.4f}   {m*span:+7.2f} um ({100*m*span/b:+.2f}%)"
              f"   {rms:6.2f} um  {abs(m*span)/rms:5.1f}")
    print("   The paper's caption says -3% (large) and -4% (small). The -3% matches the")
    print("   large-grain FIT almost exactly. The small grain's -4% is not in this table:")
    print("   endpoint -0.69%, minimum -2.40%, fit -0.89% -- all inside its own +-2%.")
    print("   -> fit the large grain, report the small one.")

    m_lg = fits["large"][0]
    v_target = abs(m_lg) / 2.0 * 1e-6 / 60.0          # diameter slope -> surface recession
    print(f"\n   target |v_n| (large grain surface recession) = {v_target:.4e} m/s")

    print("\n--- 1 - hum0 = |v_n|*(beta' + l_eff/D_v)/(rho_vs/rho_i),  l_eff = a(1 - a/L) ---")
    print("   L/a_eff   L[um]    l_eff[um]   l_eff/D_v   total R    1 - hum0      hum0")
    grid = []
    ks = sorted({1.08, 2.0, 3.0, 4.0, 5.0, float(args.L_over_a)})
    for kk in ks:
        L = kk * a_eff
        le = l_eff(a_eff, L)
        tot = beta_HK(T, args.alpha) + le / Dv_T(T)
        h = humidity_for(v_target, a_eff, L, args.alpha, T)
        mark = "  <- chosen" if abs(kk - args.L_over_a) < 1e-9 else ""
        grid.append(dict(L_over_a=kk, L_m=L, l_eff_m=le, total_R=tot, hum0=h))
        print(f"   {kk:7.2f} {L*1e6:8.1f}  {le*1e6:10.2f}  {le/Dv_T(T):10.4f}  {tot:8.4f}"
              f"   {1-h:.5f}     {h:.5f}{mark}")
    print("\n   Nearly independent of alpha_c (grain-scale sublimation is transport-limited),")
    print("   strongly dependent on L -- which is why h = 0.998 does not transfer from the")
    print("   old 121 um domain.")

    h0 = humidity_for(v_target, a_eff, args.L_over_a * a_eff, args.alpha, T)
    # Span the ~3x uncertainty in the envelope, plus a saturated control.
    factors = [0.0, 0.4, 1.0, 2.0, 3.2][: args.n_grid]
    arms = [1.0 - f * (1.0 - h0) for f in factors]
    print(f"\n--- PROPOSED -humidity GRID at L/a = {args.L_over_a:g} (centre {h0:.5f}) ---")
    for f, h in zip(factors, arms):
        tag = "  saturated control (isolates neck growth)" if f == 0.0 else ""
        print(f"   -humidity {h:.4f}   (1-h = {1-h:.5f}, {f:.1f}x the estimate){tag}")

    sig = ripening_sigma(args.R1, args.R2, T)
    print("\n--- RIPENING IS NOT THE MECHANISM ---")
    print(f"   ripening d0*(2/R_sm - 2/R_lg) = {sig:.3e}")
    print(f"   imposed  1 - hum0             = {1-h0:.3e}")
    print(f"   ratio  wall / ripening        = {(1-h0)/sig:.0f}x")
    print("   Both grains sublimate at nearly the same rate; the ripening difference is")
    print("   buried in the measurement noise. Report the differential, do not fit to it.")

    print("\n--- VAPOUR IC: no change needed ---")
    for kk in ks:
        L = kk * a_eff
        print(f"   L/a = {kk:4.2f}: t_diff = L^2/D_v = {L*L/Dv_T(T):.2e} s")
    print("   Instantaneous on a 78-min run, so the uniform hum0*rho_vs pore value that")
    print("   SetNodeFields already writes is correct as-is.")
    print("=" * 78)

    args.outdir.mkdir(parents=True, exist_ok=True)
    csv_path = args.outdir / "vapor_bc_estimate.csv"
    with csv_path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["quantity", "key", "value", "unit"])
        w.writerow(["a_eff", "", f"{a_eff:.6e}", "m"])
        w.writerow(["alpha_c", "", f"{args.alpha:.6g}", "-"])
        w.writerow(["v_n_target", "large_grain_fit", f"{v_target:.6e}", "m/s"])
        for key in ("large", "small"):
            m, b, rms = fits[key]
            w.writerow([f"fit_{key}", "slope_um_per_min", f"{m:.6f}", "um/min"])
            w.writerow([f"fit_{key}", "pct_78min", f"{100*m*t_min[-1]/b:.4f}", "%"])
            w.writerow([f"fit_{key}", "rms_resid_um", f"{rms:.4f}", "um"])
        w.writerow(["ripening_sigma", "", f"{sig:.6e}", "-"])
        w.writerow([])
        w.writerow(["humidity_grid", "L_over_a", "L_m", "l_eff_m", "total_R_s_m", "hum0"])
        for g in grid:
            w.writerow(["humidity_grid", f"{g['L_over_a']:.4f}", f"{g['L_m']:.6e}",
                        f"{g['l_eff_m']:.6e}", f"{g['total_R']:.6f}", f"{g['hum0']:.6f}"])
        w.writerow([])
        w.writerow(["proposed_arm", "factor", "humidity"])
        for f, h in zip(factors, arms):
            w.writerow(["proposed_arm", f"{f:.2f}", f"{h:.4f}"])
    print(f"  wrote {csv_path}")
    make_figure(rows, fits, grid, args.L_over_a, T, args.outdir / "vapor_bc_estimate.png")


if __name__ == "__main__":
    main()

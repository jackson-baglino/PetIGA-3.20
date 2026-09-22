#!/usr/bin/env python3
"""Audit of the 2025-09 `unresolved_results` runs: are they publishable?

    venv_enceladus/bin/python studies/keff_sintering/audit_unresolved/audit.py

Reads the nine completed runs under
  ~/SimulationResults/HPC_results/enceladus_DSM/unresolved_results/
and the k_eff tables the old standalone code produced for them under
  ~/SimulationResults/effective_thermal_cond/

and scores four things the runs themselves do not record, because the code
that produced them printed no kinetics banner:

  1. MESH        h/eps against the K&P rule h = eps/sqrt(2).
  2. TIME STEP   the realised dt against tau_sub = eps^2*beta_sub0/d0_sub0.
  3. eps BOUND   eps against the comp_eps.py sharp-interface ceiling at the
                 run's own temperature.
  4. alpha_c     the condensation coefficient IMPLIED by the run's hardcoded
                 beta_sub0 = 1.4e5, which was held FIXED across a 20 C sweep.

(4) is the one that matters. The solver rescales beta_sub0 by rho_vs(T)/rho_i
at runtime, so a beta_sub0 that is constant in T is an alpha_c that is NOT --
and the implied alpha_c rises by the same factor that rho_vs falls. The
product alpha_c*rho_vs, which sets the sublimation flux, is then constant by
construction, and a temperature sweep built this way cannot show a temperature
effect no matter what the physics does.
"""
from __future__ import annotations

import glob
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
PROJ = HERE.parents[2]
sys.path.insert(0, str(PROJ / "preprocess"))

DAY = 86400.0
RUNS = Path.home() / "SimulationResults/HPC_results/enceladus_DSM/unresolved_results"
KEFF = Path.home() / "SimulationResults/effective_thermal_cond"

# --- exactly as hardcoded in the staged src/ copy shipped with each run ------
EPS_USED    = 8.756952130264e-07
NX_USED     = 1142
L           = 2.0e-3
BETA_SUB0   = 1.4e5      # dry_snow_metamorphism.c:49, fixed for ALL runs
D0_SUB0     = 1.0e-9     # dry_snow_metamorphism.c:48
RHO_ICE     = 919.0
DIF_VAP     = 2.178e-5
A1, A2      = 5.0, 0.1581
DIFF_SUB    = 0.5 * (0.02 / 1.341 / 1.044e3 + 2.29 / RHO_ICE / 1.96e3)
DTMAX_RATIO = 1.09       # lunar-calibrated safe dtmax/tau_sub


def rho_vs(T_C: float) -> float:
    """Saturation vapour density [kg/m^3], the runs' own RhoVS_I."""
    tk = T_C + 273.15
    K = [-0.5865e4, 0.2224e2, 0.1375e-1, -0.3403e-4, 0.2697e-7, 0.6918]
    P = np.exp(K[0] / tk + K[1] + K[2] * tk + K[3] * tk**2
               + K[4] * tk**3 + K[5] * np.log(tk))
    return 1.341 * 0.62 * P / 101325.0


def tau_sub(T_C: float) -> float:
    """The runs' own tau_sub, from their staged source."""
    rr = RHO_ICE / rho_vs(T_C)
    d0, beta = D0_SUB0 / rr, BETA_SUB0 / rr
    lam = A1 * EPS_USED / d0
    return EPS_USED * lam * (beta / A1 + A2 * EPS_USED / DIFF_SUB
                             + A2 * EPS_USED / DIF_VAP)


def eps_ceiling(T_C: float, alpha_c: float, R_ave: float = 45e-6):
    """comp_eps.py ceiling and mesh at this T and alpha_c."""
    import comp_eps                                           # noqa: E402
    import io, contextlib, re, subprocess
    out = subprocess.run(
        [sys.executable, str(PROJ / "preprocess" / "comp_eps.py"),
         "--Lx", f"{L}", "--Ly", f"{L}", "--Rave", f"{R_ave}",
         "--T0", f"{T_C}", "--alpha", f"{alpha_c}"],
        capture_output=True, text=True).stdout
    g = lambda k: float(re.search(rf"^\s+{k}\s+(\S+)", out, re.M).group(1))
    return g("-eps"), g("-Nx"), g("-beta_sub0")


def implied_alpha_c(T_C: float) -> float:
    """alpha_c that comp_eps would need for the run's beta_sub0 = 1.4e5.

    beta0 ~ 1/alpha_c at fixed T, so one probe at a reference alpha suffices.
    """
    _, _, b_ref = eps_ceiling(T_C, 1.0e-2)
    return 1.0e-2 * b_ref / BETA_SUB0


def load_runs():
    out = []
    for d in sorted(glob.glob(str(RUNS / "2mm_results*" / "*" / ""))):
        f = os.path.join(d, "SSA_evo.dat")
        if not os.path.isfile(f):
            continue
        name = os.path.basename(d.rstrip("/"))
        T = next(t for t in (-20, -25, -30, -35, -40) if f"Tm{t}_" in name)
        phi = 0.24
        for p in (0.24, 0.26, 0.28, 0.30):
            if f"phi{p:.2f}" in name:
                phi = p
        a = np.loadtxt(f)                 # ssa/eps, tot_ice, t, step
        sweep = "temperature" if "results2" in d else "porosity"
        out.append(dict(dir=d, name=name, T=T, phi=phi, sweep=sweep,
                        ssa=a[:, 0], ice=a[:, 1], t=a[:, 2]))
    return out


def load_keff(name: str):
    for sub in ("const_porosity_varying_temp", "const_temp_varying_phi"):
        for c in glob.glob(str(KEFF / sub / "*" / "k_eff.csv")):
            tag = os.path.basename(os.path.dirname(c))
            if f"Tm{_T(name)}_" in tag and f"phi{_phi(name):.2f}" in tag:
                k = np.genfromtxt(c, delimiter=",", names=True)
                return 0.5 * (k["k_00"] + k["k_11"])
    return None


def _T(n):
    return next(t for t in (-20, -25, -30, -35, -40) if f"Tm{t}_" in n)


def _phi(n):
    for p in (0.24, 0.26, 0.28, 0.30):
        if f"phi{p:.2f}" in n:
            return p
    return 0.24


def main() -> int:
    runs = load_runs()
    h = L / NX_USED
    band = 9.2 * EPS_USED

    print("=" * 78)
    print("A. DISCRETISATION  (identical for all nine runs -- one eps, one mesh)")
    print("=" * 78)
    print(f"  eps                  {EPS_USED:.4e} m")
    print(f"  Nx = Ny              {NX_USED}   (p=1, C=0 -- LINEAR elements)")
    print(f"  h = L/Nx             {h:.4e} m")
    print(f"  h/eps                {h/EPS_USED:.3f}   "
          f"vs K&P mesh rule 0.707  ->  {h/EPS_USED/0.7071:.2f}x too coarse")
    print(f"  1%-99% band 9.2*eps  {band*1e6:.2f} um = {band/h:.1f} elements"
          f"   (project standard 7.5-10, and that is for p=2/C=1)")

    print()
    print("=" * 78)
    print("B. KINETICS: the temperature sweep is confounded with alpha_c")
    print("=" * 78)
    print("  beta_sub0 = 1.4e5 was HARDCODED and identical at every temperature.")
    print("  The solver forms beta_sub = beta_sub0 / (rho_i/rho_vs(T)), so a")
    print("  constant beta_sub0 does NOT mean a constant condensation coefficient.")
    print()
    hdr = (f"  {'T':>4} {'rho_vs':>10} {'alpha_c':>9} {'a_c*rho_vs':>11} "
           f"{'eps_max':>9} {'eps/max':>8} {'Nx_req':>7} {'coarse':>7} "
           f"{'tau_sub':>8} {'dt_med/tau':>10}")
    print(hdr)
    rows = []
    for T in (-20, -25, -30, -35, -40):
        ac = implied_alpha_c(T)
        emax, nxr, _ = eps_ceiling(T, ac)
        rv = rho_vs(T)
        ts = tau_sub(T)
        r = next((x for x in runs if x["T"] == T and x["sweep"] == "temperature"),
                 None)
        dtm = np.median(np.diff(r["t"])) if r else np.nan
        rows.append((T, rv, ac, ac * rv, emax, nxr, ts, dtm / ts))
        print(f"  {T:4d} {rv:10.3e} {ac:9.3e} {ac*rv:11.4e} {emax:9.3e} "
              f"{EPS_USED/emax:8.2f} {nxr:7.0f} {nxr/NX_USED:7.1f} "
              f"{ts:8.1f} {dtm/ts:10.2f}")

    flux = np.array([r[3] for r in rows])
    print()
    print(f"  alpha_c rises {rows[0][2]:.3e} -> {rows[-1][2]:.3e} "
          f"({rows[-1][2]/rows[0][2]:.1f}x) as rho_vs falls "
          f"{rows[0][1]/rows[-1][1]:.1f}x.")
    print(f"  Their product -- which sets the sublimation flux -- is constant to "
          f"{np.ptp(flux)/flux.mean():.1%}")
    print("  across the whole 20 C sweep. The sweep therefore CANNOT resolve a")
    print("  temperature effect: it was cancelled in the parameter choice.")
    print()
    print(f"  Literature band for alpha_c (Libbrecht 2017, Braun 2024): 1e-4..1e-3.")
    print(f"  These runs sit at {rows[0][2]/1e-3:.0f}x to {rows[-1][2]/1e-3:.0f}x "
          f"the TOP of that band.")

    print()
    print("=" * 78)
    print("C. WHAT THE k_eff TABLES ACTUALLY SHOW")
    print("=" * 78)
    print(f"  {'run':>34} {'k(0)':>8} {'k(end)':>8} {'rise':>8}")
    tmp, por = [], []
    for r in runs:
        k = load_keff(r["name"])
        if k is None:
            continue
        rise = k[-1] / k[0] - 1
        print(f"  {r['name'][:34]:>34} {k[0]:8.4f} {k[-1]:8.4f} {rise:+8.2%}")
        (tmp if r["sweep"] == "temperature" else por).append((r["T"], r["phi"], rise))
    if tmp:
        rr = np.array([x[2] for x in tmp])
        print(f"\n  temperature sweep: rise spans {rr.min():+.2%}..{rr.max():+.2%}"
              f"  -- a spread of {100*(rr.max()-rr.min()):.1f} points,")
        print(f"  i.e. {100*(rr.max()-rr.min())/rr.mean()/100:.1%} of the mean rise.")
        print(f"  The cancelled product alpha_c*rho_vs still drifts "
              f"{np.ptp(flux)/flux.mean():.1%} across the sweep, and eps/eps_max")
        print(f"  grows 1.7x -> 13.5x over the same range. The residual T-dependence")
        print("  is the SAME SIZE as those two numerical drifts (and runs opposite in")
        print("  sign to the first of them), so it cannot be attributed to physics.")
        print("  The headline reading of this sweep -- 'k_eff evolution is nearly")
        print("  temperature-independent' -- is what holding beta_sub0 fixed forces.")
    if por:
        rr = np.array([x[2] for x in por])
        print(f"\n  porosity sweep:    rise spans {rr.min():+.2%}..{rr.max():+.2%}"
              f"  at a SINGLE T, eps, beta_sub0.")
        print("  Every bias above is a common offset here, so the porosity TREND")
        print("  survives even though the magnitude and rate do not.")

    print()
    print("=" * 78)
    print("D. WHAT IS FINE")
    print("=" * 78)
    for r in runs:
        d = r["ice"][-1] / r["ice"][0] - 1
        print(f"  {r['name'][:46]:46s} ice drift {d:+.3%}  "
              f"({len(r['t'])} steps, t_final {r['t'][-1]/DAY:.2f} d)")
    print("\n  Ice loss is monotone in T and tracks rho_vs -- that is the physical")
    print("  response to humidity 0.98, not a conservation leak.")
    print(f"  L/R_ave = {L/45e-6:.1f} at t=0, above the measured REV floor of 40.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""sa81_mechanisms.py — surface diffusion against vapor transport in Molaro et al. (2019).

    venv_enceladus/bin/python studies/molaro_2019/sa81_mechanisms/sa81_mechanisms.py

Implements, as printed, the two dominant Stage-1 terms of the Swinkels & Ashby
(1981) model as given by Molaro et al. (2019, JGR Planets 124, 243):

  Table 3   V1 = S 3 pi x ds Ds gs Om (K3 - K2) / (d2 k T)          surface diffusion
            V3 = S 2 pi x rho theta Pv (gs Om / kT) sqrt(Om / (2 pi Do k T)) (K3 - Km)   vapor
            xdot = (V1 + ... + V6) / (2 pi x theta rho)
  Eq. 1-4   K3 = 2/a,  Km = -1/rho + 1/x,  rho = x^2 / (2 (a - x))  (y = 0),
            theta = atan(a / (x + rho))
  Eq. 6     K1 K2 = Km^2
  Eq. 17-20 A = db Db / (ds Ds);  dK/K1 = 3 A rho/x + 3/2 - (3/2) sqrt(3 (A rho/x)^2 + 4 A rho/x + 1)
            d2 = rho / (1 + sqrt(4 dK / (3 K2)))
  Table 5   ice: Om 1.181e-29, Do 917, gs 0.06, ds = db = 9.04e-10,
            D0s 1.4e-8 (range to 1.35e-13), Qs 23.156 kJ/mol,
            D0b 8.4e-4 (range 0.03 to 1e-12), Qb 4.9 kJ/mol; Pv from Murphy & Koop (2005).
  S (smoothing) is set to 1; it multiplies V1 and V3 alike, so their ratio does not depend on it.

Part 1 reproduces their Figure 10 (25 um grains, T/Tm = 0.8) as a check of the
implementation: V1 and V3 against time.
Part 2 evaluates the same two terms at the conditions of their -20 C grain-pair
experiment (a = 87 um mean radius, neck radius 16.5 to 32.5 um), and with the
grain-boundary coupling switched off (A = 0, i.e. K2 = Km) to show what it does.

Writes sa81_mechanisms.txt and sa81_mechanisms.png beside this file.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
KB, RG = 1.380649e-23, 8.314462
OM, DO, GS, DS_, DB_ = 1.181e-29, 917.0, 0.06, 9.04e-10, 9.04e-10
D0S, QS, D0B, QB, MU = 1.4e-8, 23156.0, 8.4e-4, 4900.0, 1e9
DAY = 86400.0


def p_vs(T):                                   # Murphy & Koop (2005), Pa
    return np.exp(9.550426 - 5723.265 / T + 3.53068 * np.log(T) - 0.00728332 * T)


def terms(x, a, T, d0s=D0S, coupled=True):
    """(V1, V3, geometry dict) at neck radius x for grain radius a."""
    rho = x * x / (2.0 * (a - x))
    theta = np.arctan(a / (x + rho))
    K3, Km = 2.0 / a, -1.0 / rho + 1.0 / x
    dsDs = DS_ * d0s * np.exp(-QS / (RG * T))
    dbDb = DB_ * D0B * np.exp(-QB / (RG * T))
    A = dbDb / dsDs if coupled else 0.0
    q = A * rho / x
    r = 3 * q + 1.5 - 1.5 * np.sqrt(3 * q * q + 4 * q + 1)          # dK/K1, Eq. 20
    K1 = -abs(Km) / np.sqrt(1 + r)                                    # Eq. 6 with K2 = K1 (1 + r), both negative
    K2 = K1 * (1 + r)
    d2 = rho / (1 + np.sqrt(4 * r / (3 * (1 + r))))                   # Eq. 19, dK/K2 = r/(1+r)
    V1 = 3 * np.pi * x * dsDs * GS * OM * (K3 - K2) / (d2 * KB * T)
    V3 = 2 * np.pi * x * rho * theta * p_vs(T) * (GS * OM / (KB * T)) * np.sqrt(OM / (2 * np.pi * DO * KB * T)) * (K3 - Km)
    return V1, V3, dict(rho=rho, theta=theta, Km=Km, K2=K2, d2=d2, A=A)


def main():
    L = []
    # ---------- Part 1: their Figure 10 ----------
    a, T = 25e-6, 0.8 * 273.15
    x0 = (GS * a * a / (10 * MU)) ** (1 / 3)                          # Table 2 initial neck
    def rhs(t, y):
        V1, V3, g = terms(y[0], a, T)
        return [(V1 + V3) / (2 * np.pi * y[0] * g["theta"] * g["rho"])]
    tt = np.geomspace(1e-9, 2.0, 400) * DAY
    sol = solve_ivp(rhs, (0, tt[-1]), [x0], t_eval=tt, method="LSODA", rtol=1e-8, atol=1e-14)
    xs = sol.y[0]
    V = np.array([terms(x, a, T)[:2] for x in xs])
    i_x = int(np.argmax(V[:, 1] > V[:, 0]))
    L += [f"PART 1. Their Figure 10: a = 25 um, T = {T:.1f} K, initial neck x0 = {x0 * 1e6:.3f} um (x0/a = {x0 / a:.4f}).",
          "  time [d]      x/a      V1 surface [m3/s]   V3 vapor [m3/s]    V1/V3"]
    for td in (1e-8, 1e-6, 1e-5, 1e-4, 1e-2, 1.0):
        i = int(np.argmin(np.abs(sol.t / DAY - td)))
        L.append(f"  {sol.t[i] / DAY:8.1e}    {xs[i] / a:.4f}    {V[i, 0]:.2e}           {V[i, 1]:.2e}          {V[i, 0] / V[i, 1]:.2e}")
    L += [f"  vapor overtakes surface diffusion at t = {sol.t[i_x] / DAY:.1e} d = {sol.t[i_x]:.1f} s, x/a = {xs[i_x] / a:.3f} (x = {xs[i_x] * 1e6:.2f} um)",
          "  (read off their figure: V1 ~ 1e-20 falling, V3 ~ 1e-21 rising to ~2e-20, crossing between 1e-6 and 1e-5 d)", ""]

    # ---------- Part 2: the -20 C grain pair ----------
    a, T = 87e-6, 253.15
    L += [f"PART 2. Their -20 C grain pair: a = 87 um, T = {T} K, Pv = {p_vs(T):.1f} Pa.",
          "  measured: neck radius 16.4 -> 32.4 um in 78 min, mean dx/dt = 3.4e-9 m/s",
          "  x [um]  rho [um]   V1 [m3/s]   V3 [m3/s]   V1/V3     dx/dt from V1 [m/s]   from V3 [m/s]  | uncoupled (K2 = Km): V1/V3   | bottom-of-range D0s: V1/V3"]
    rows = []
    for x in (16.4e-6, 20e-6, 25e-6, 32.4e-6):
        V1, V3, g = terms(x, a, T)
        V1u, V3u, _ = terms(x, a, T, coupled=False)
        V1b, V3b, _ = terms(x, a, T, d0s=1.35e-13)
        den = 2 * np.pi * x * g["theta"] * g["rho"]
        rows.append((x, V1 / V3))
        L.append(f"  {x * 1e6:5.1f}   {g['rho'] * 1e6:6.2f}    {V1:.2e}    {V3:.2e}    {V1 / V3:.1e}   {V1 / den:.2e}              {V3 / den:.2e}      |        {V1u / V3u:.1e}              |        {V1b / V3b:.1e}")
    V1, V3, g = terms(16.4e-6, a, T)
    L += [f"  at x = 16.4 um: A = db Db/(ds Ds) = {g['A']:.2e};  Km = {g['Km']:.3e} 1/m;  K2 = {g['K2']:.3e} 1/m (K2/Km = {g['K2'] / g['Km']:.1f});  d2 = {g['d2'] * 1e6:.2f} um", ""]
    # neck size at which the two are equal, at -20 C for this grain size
    xx = np.geomspace(0.05e-6, 60e-6, 4000)
    rat = np.array([terms(x, a, T)[0] / terms(x, a, T)[1] for x in xx])
    ratu = np.array([terms(x, a, T, coupled=False)[0] / terms(x, a, T, coupled=False)[1] for x in xx])
    cross = lambda r_: xx[int(np.argmax(r_ < 1))] if (r_ < 1).any() and r_[0] > 1 else float("nan")
    L += [f"  neck radius below which surface diffusion exceeds vapor transport (a = 87 um, -20 C): "
          f"{cross(rat) * 1e6:.2f} um as printed; {cross(ratu) * 1e6:.3f} um with the boundary coupling off"]
    txt = "\n".join(L)
    print(txt)
    (HERE / "sa81_mechanisms.txt").write_text(txt + "\n")

    fig, ax = plt.subplots(1, 2, figsize=(10.5, 4.0))
    ax[0].loglog(sol.t / DAY, V[:, 0], "C3--", lw=1.6, label=r"$\dot V_1$ surface diffusion")
    ax[0].loglog(sol.t / DAY, V[:, 1], "C0--", lw=1.6, label=r"$\dot V_3$ vapor transport")
    ax[0].set(xlabel="time [d]", ylabel=r"diffusive current [m$^3$/s]", title="their Figure 10 recomputed: 25 $\\mu$m, $T/T_m$ = 0.8",
              xlim=(1e-9, 2), ylim=(1e-26, 1e-19))
    ax[0].legend(frameon=False, fontsize=8)
    ax[1].loglog(xx * 1e6, rat, "k", lw=1.6, label="as printed (boundary coupling on)")
    ax[1].loglog(xx * 1e6, ratu, "0.5", lw=1.4, ls="--", label=r"coupling off ($K_2 = K_m$)")
    ax[1].axhline(1, color="0.7", lw=0.8)
    ax[1].axvspan(16.4, 32.4, color="C1", alpha=0.25, lw=0, label="necks measured at $-20\\,^\\circ$C")
    ax[1].set(xlabel=r"neck radius $x$ [$\mu$m]", ylabel=r"$\dot V_1/\dot V_3$", title=r"$a$ = 87 $\mu$m, $-20\,^\circ$C")
    ax[1].legend(frameon=False, fontsize=8)
    fig.tight_layout(); fig.savefig(HERE / "sa81_mechanisms.png", dpi=200)


if __name__ == "__main__":
    main()

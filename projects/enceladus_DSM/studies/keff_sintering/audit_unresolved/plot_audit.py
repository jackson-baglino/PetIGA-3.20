#!/usr/bin/env python3
"""Figures for the unresolved_results audit. Writes audit.png next to itself."""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
from audit import (EPS_USED, NX_USED, eps_ceiling, implied_alpha_c,  # noqa: E402
                   load_keff, load_runs, rho_vs)

TS = (-20, -25, -30, -35, -40)
runs = load_runs()

ac = np.array([implied_alpha_c(T) for T in TS])
rv = np.array([rho_vs(T) for T in TS])
emax = np.array([eps_ceiling(T, a)[0] for T, a in zip(TS, ac)])
nxr = np.array([eps_ceiling(T, a)[1] for T, a in zip(TS, ac)])

fig, ax = plt.subplots(2, 2, figsize=(11.5, 8.2))

# (a) the cancellation
a0 = ax[0, 0]
a0.semilogy(TS, ac / ac[0], "o-", label=r"implied $\alpha_c$ (rises 8.4$\times$)")
a0.semilogy(TS, rv / rv[0], "s-", label=r"$\rho_{vs}(T)$ (falls 8.0$\times$)")
a0.semilogy(TS, (ac * rv) / (ac * rv)[0], "k^-", lw=2.5,
            label=r"product $\alpha_c\,\rho_{vs}$ — flat to 4.2%")
a0.axhline(1.0, color="0.7", lw=0.8, zorder=0)
a0.set_xlabel("T [°C]"); a0.set_ylabel("normalised to T = −20 °C")
a0.set_title("(a) A fixed beta_sub0 cancels the physical T-dependence")
a0.legend(fontsize=8); a0.grid(alpha=.3)

# (b) the consequence
a1 = ax[0, 1]
rise = []
for T in TS:
    r = next(x for x in runs if x["T"] == T and x["sweep"] == "temperature")
    k = load_keff(r["name"]); rise.append(100 * (k[-1] / k[0] - 1))
a1.plot(TS, rise, "o-", color="crimson", lw=2)
a1.set_ylim(0, 60)
a1.set_xlabel("T [°C]"); a1.set_ylabel(r"$k_{eff}$ rise over 28 d [%]")
a1.set_title("(b) …so the sweep reports no T-effect: 47.8–49.5%")
a1.annotate("1.7 points across 20 °C\n= the size of the numerical drift",
            xy=(-30, rise[2]), xytext=(-37, 22), fontsize=9,
            arrowprops=dict(arrowstyle="->", color="0.4"))
a1.grid(alpha=.3)

# (c) resolution
a2 = ax[1, 0]
a2.semilogy(TS, EPS_USED / emax, "o-", color="darkorange", lw=2,
            label=r"$\epsilon_{used}/\epsilon_{max}$")
a2.semilogy(TS, nxr / NX_USED, "s-", color="purple", lw=2,
            label=r"$N_x$ required / 1142")
a2.axhline(1.0, color="k", ls="--", lw=1.2, label="requirement")
a2.set_xlabel("T [°C]"); a2.set_ylabel("violation factor")
a2.set_title("(c) One eps and one mesh for the whole sweep")
a2.legend(fontsize=8); a2.grid(alpha=.3, which="both")

# (d) porosity sweep -- the part that survives
a3 = ax[1, 1]
for r in sorted((x for x in runs if x["sweep"] == "porosity"),
                key=lambda z: z["phi"]):
    k = load_keff(r["name"])
    t = np.linspace(0, 28, len(k))
    a3.plot(t, k, lw=1.8, label=rf"$\phi$ = {r['phi']:.2f}")
a3.set_xlabel("t [days]"); a3.set_ylabel(r"$k_{eff}$ [W m$^{-1}$ K$^{-1}$]")
a3.set_title("(d) Porosity sweep — one T, one eps: the trend survives")
a3.legend(fontsize=8); a3.grid(alpha=.3)

fig.suptitle("Audit of unresolved_results (2025-09): what is and is not usable",
             fontsize=13, y=0.995)
fig.tight_layout()
fig.savefig(HERE / "audit.png", dpi=150)
print("wrote", HERE / "audit.png")

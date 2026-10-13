#!/usr/bin/env python3
"""geometric_mean.py — the volume-weighted geometric mean against the simulated k_eff.

    venv_enceladus/bin/python studies/keff_sintering/mixing_laws/geometric_mean.py

The mixing rule "product of the conductivities, each raised to its volume
fraction",

    k_geo = k_ice^(1 - phi) * k_air^phi ,   i.e.  k_geo / k_ice = exp(phi ln(k_air/k_ice)) ,

depends on porosity only. This compares it with the 125 production runs at
four sintering ages (opening frame, theta = 30, 331, and the end of the -5 C
runs, theta = 1256), alongside the arithmetic (parallel) and harmonic
(series) means and the two-dimensional Hashin-Shtrikman bounds.

Writes geometric_mean.txt and geometric_mean.png beside this file.
"""
from __future__ import annotations
import re, sys
from collections import defaultdict
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO / "postprocess"))
from plot_keff import load, read_tau_sub                              # noqa: E402

CAMP = Path.home() / "SimulationResults/HPC_results/enceladus_DSM/keff_sintering_campaign"
PAT = re.compile(r"phi([0-9.]+)_Rave50um_LR40_seed(\d+)_L2mm_eps1000nm_perxy_T(-?\d+)__")
K_I, K_A = 2.29, 0.02
AGES = (30.0, 331.0, 1256.0)


def geo(p): return K_I ** (1 - p) * K_A ** p
def arith(p): return (1 - p) * K_I + p * K_A
def harm(p): return 1.0 / ((1 - p) / K_I + p / K_A)
def hs_upper_2d(p): return K_I + p / (1.0 / (K_A - K_I) + (1 - p) / (2 * K_I))        # ice continuous
def hs_lower_2d(p): return K_A + (1 - p) / (1.0 / (K_I - K_A) + p / (2 * K_A))        # air continuous


def main():
    val = defaultdict(lambda: defaultdict(list))                      # val[phi][age] -> k per packing/temperature
    for kf in sorted(CAMP.glob("packing_*/k_eff.csv")):
        m = PAT.search(kf.parent.name)
        if not m or not 1601 <= int(m.group(2)) <= 2005:
            continue
        r = load(kf.parent)
        th = r["t"] / read_tau_sub(kf.parent)
        p = float(m.group(1))
        i0 = int(np.argmax(r["t"] >= 1.0))
        if int(m.group(3)) == -20:
            val[p]["open"].append(r["kiso"][i0])
        for a in AGES:
            if th[-1] >= 0.99 * a:
                val[p][a].append(np.interp(a, th, r["kiso"]))
    phis = sorted(val)
    L = [f"k_ice = {K_I}, k_air = {K_A} W/m/K.  k_geo/k_ice = exp({np.log(K_A / K_I):.2f} phi).",
         "Simulated k_eff: mean over packings (and temperatures where the age is reached), W/m/K; ratio to the geometric mean in brackets.", "",
         "phi     k_geo   harmonic  HS lower(2D)  HS upper(2D)  arithmetic |  opening        theta = 30     theta = 331    theta = 1256"]
    for p in phis:
        row = f"{p:.3f}   {geo(p):.3f}   {harm(p):.3f}     {hs_lower_2d(p):.3f}         {hs_upper_2d(p):.3f}         {arith(p):.3f}      | "
        for a in ("open",) + AGES:
            k = np.mean(val[p][a]); row += f" {k:.3f} ({k / geo(p):.2f})  "
        L.append(row)
    for a in AGES:
        k = np.array([np.mean(val[p][a]) for p in phis]); q = np.array(phis)
        c = np.polyfit(q, np.log(k / K_I), 1); c3 = np.polyfit(q[:3], np.log(k[:3] / K_I), 1)
        L.append(f"theta = {a:g}: fit k_eff/k_ice = {np.exp(c[1]):.2f} exp({c[0]:.2f} phi) over all five porosities;  "
                 f"{np.exp(c3[1]):.2f} exp({c3[0]:.2f} phi) over phi <= 0.375;  geometric mean: 1.00 exp({np.log(K_A / K_I):.2f} phi)")
    txt = "\n".join(L)
    print(txt)
    (HERE / "geometric_mean.txt").write_text(txt + "\n")

    pp = np.linspace(0.25, 0.50, 100)
    fig, ax = plt.subplots(figsize=(6.2, 4.4))
    ax.fill_between(pp, hs_lower_2d(pp), hs_upper_2d(pp), color="0.92", lw=0, label="Hashin-Shtrikman bounds (2D)")
    ax.plot(pp, arith(pp), color="0.6", lw=0.9, ls=":", label="arithmetic (parallel)")
    ax.plot(pp, harm(pp), color="0.6", lw=0.9, ls="-.", label="harmonic (series)")
    ax.plot(pp, geo(pp), color="k", lw=1.8, label=r"geometric mean, $k_\mathrm{ice}^{1-\varphi}k_\mathrm{air}^{\varphi}$")
    for a, c, lab in (("open", "C0", "opening frame"), (30.0, "C2", r"$\theta = 30$"), (331.0, "C1", r"$\theta = 331$"), (1256.0, "C3", r"$\theta = 1256$")):
        m_ = [np.mean(val[p][a]) for p in phis]; lo = [np.min(val[p][a]) for p in phis]; hi = [np.max(val[p][a]) for p in phis]
        ax.errorbar(phis, m_, yerr=[np.subtract(m_, lo), np.subtract(hi, m_)], fmt="o", color=c, ms=4, lw=0.9, capsize=2, label=f"simulated, {lab}")
    ax.set(yscale="log", xlabel=r"porosity $\varphi$", ylabel=r"$k_\mathrm{eff}$ [W m$^{-1}$ K$^{-1}$]", ylim=(0.03, 2.0))
    ax.legend(frameon=False, fontsize=7.5, ncol=2, loc="lower left")
    fig.tight_layout(); fig.savefig(HERE / "geometric_mean.png", dpi=200)


if __name__ == "__main__":
    main()

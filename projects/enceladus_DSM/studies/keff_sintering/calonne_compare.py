#!/usr/bin/env python3
"""Our porosity series against the density regressions of Calonne et al.

  calonne_density.png   (a) k_iso against density: every packing as a vertical
                        segment from its opening sample (t >= 1 s) to its last,
                        the Calonne et al. (2011) quadratic (solid inside its
                        fitted range, 100-550 kg/m3, dashed beyond), the
                        Calonne et al. (2019) snow-firn formula, and the 2D
                        Hashin-Shtrikman bounds that contain any 2D microstructure.
                        (b) k_iso(t), seed mean per porosity (band = seed min-max),
                        as the fraction h of the bound window it occupies,
                        h = (k - k_HS-) / (k_HS+ - k_HS-); the dashed lines are
                        the Calonne (2019) value at the same density as a
                        fraction of the 3D window. h removes the 2D/3D offset
                        in absolute k; it is a heuristic, not a mapping.
  calonne_compare.csv   per packing: density, k at start/end, Calonne values,
                        bounds, h at start/end.

Calonne formulas are fitted at -3 C (k_i 2.107, k_a 0.024 W/m/K). They are
moved to our k_i = 2.29, k_a = 0.02 with the factors of Calonne 2019 Eq. 5:
k_i k_a / (k_i,ref k_a,ref) on the snow term, k_i / k_i,ref on the firn term.

    venv_enceladus/bin/python studies/keff_sintering/calonne_compare.py \
        ~/SimulationResults/HPC_results/enceladus_DSM/keff_sintering_campaign
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
PROJ = HERE.parents[1]
sys.path.insert(0, str(PROJ / "postprocess"))
from plot_keff import load  # noqa: E402
from compare_keff import CMAP  # noqa: E402

DAY = 86400.0
MM = 1 / 25.4
PAT = re.compile(r"phi([\d.]+)_Rave50um_LR40_seed(\d+)_L2mm_eps1000nm_perxy_T(-?\d+)__")
INK, MUTED = "#1a1a1a", "#6b6b6b"

RHO_I = 917.0                      # ice density used by Calonne et al.
K_I, K_A = 2.29, 0.02              # ours (src/enceladus_main.c)
K_I_REF, K_A_REF = 2.107, 0.024    # Calonne fits, -3 C


def k_snow_2011(rho):
    """Calonne et al. (2011) Eq. 12, at -3 C."""
    return 2.5e-6 * rho**2 - 1.23e-4 * rho + 0.024


def k_calonne_2019(rho, ki=K_I, ka=K_A):
    """Calonne et al. (2019) Eq. 5, moved to (ki, ka)."""
    theta = 1.0 / (1.0 + np.exp(-2 * 0.02 * (rho - 450.0)))
    k_firn = 2.107 + 0.003618 * (rho - RHO_I)
    return ((1 - theta) * ki * ka / (K_I_REF * K_A_REF) * k_snow_2011(rho)
            + theta * ki / K_I_REF * k_firn)


def k_2011_moved(rho, ki=K_I, ka=K_A):
    return ki * ka / (K_I_REF * K_A_REF) * k_snow_2011(rho)


def hs_bounds(ice, d, ki=K_I, ka=K_A):
    """Hashin-Shtrikman bounds for an isotropic two-phase medium in d dims."""
    air = 1.0 - ice
    up = ki + air / (1.0 / (ka - ki) + ice / (d * ki))
    lo = ka + ice / (1.0 / (ki - ka) + air / (d * ka))
    return lo, up


def style(ax):
    ax.grid(True, alpha=0.25, lw=0.5)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(labelsize=8)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root", type=Path)
    ap.add_argument("--T", type=int, default=-20)
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args()
    out = a.out or a.root / "compare" / "calonne" / f"T{a.T}"
    out.mkdir(parents=True, exist_ok=True)

    G = defaultdict(list)
    for kf in sorted(a.root.glob("packing_*/k_eff.csv")):
        m = PAT.search(kf.parent.name)
        if not m or int(m.group(3)) != a.T or not 1601 <= int(m.group(2)) <= 2005:
            continue
        r = load(kf.parent)
        raw = np.genfromtxt(kf, delimiter=",", names=True)
        r["ice"] = float(np.mean(raw["phi_bar"]))
        r["seed"] = int(m.group(2))
        G[float(m.group(1))].append(r)
    phis = sorted(G)
    if not phis:
        raise SystemExit(f"no production runs at T = {a.T}")
    cm, (c0, c1) = CMAP["phi"]
    col = {p: cm(c0 + (c1 - c0) * i / max(1, len(phis) - 1)) for i, p in enumerate(phis)}

    rows = []
    for p in phis:
        for r in G[p]:
            i0 = int(np.argmax(r["t"] >= 1.0))
            rho = r["ice"] * RHO_I
            lo2, up2 = hs_bounds(r["ice"], 2)
            lo3, up3 = hs_bounds(r["ice"], 3)
            kc = float(k_calonne_2019(rho))
            k0, k1 = r["kiso"][i0], r["kiso"][-1]
            rows.append(dict(phi=p, seed=r["seed"], ice=r["ice"], rho=rho, k_start=k0, k_end=k1,
                             k_cal2011=float(k_2011_moved(rho)), k_cal2019=kc,
                             hs2_lo=lo2, hs2_up=up2, hs3_lo=lo3, hs3_up=up3,
                             h_start=(k0 - lo2) / (up2 - lo2), h_end=(k1 - lo2) / (up2 - lo2),
                             h_cal2019=(kc - lo3) / (up3 - lo3),
                             k_end_over_cal=k1 / kc))
    with open(out / "calonne_compare.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

    fig, ax = plt.subplots(1, 2, figsize=(170 * MM, 72 * MM), constrained_layout=True)

    # (a) k against density
    rr = np.linspace(100, 917, 400)
    lo2, up2 = hs_bounds(rr / RHO_I, 2)
    ax[0].fill_between(rr, lo2, up2, color="#e6e6e6", lw=0, zorder=0)
    ax[0].plot(rr, up2, color="#b0b0b0", lw=0.8, zorder=1)
    ax[0].text(300, float(np.interp(300, rr, up2)) + 0.05, "2D bounds", fontsize=8,
               color=MUTED, rotation=0)
    inside = rr <= 550
    ax[0].plot(rr[inside], k_2011_moved(rr[inside]), color=INK, lw=1.4,
               label="Calonne et al. (2011)")
    ax[0].plot(rr[~inside], k_2011_moved(rr[~inside]), color=INK, lw=1.0, ls="--")
    ax[0].plot(rr, k_calonne_2019(rr), color=INK, lw=1.0, ls=":",
               label="Calonne et al. (2019)")
    for p in phis:
        for r in G[p]:
            i0 = int(np.argmax(r["t"] >= 1.0))
            x = r["ice"] * RHO_I
            ax[0].plot([x, x], [r["kiso"][i0], r["kiso"][-1]], color=col[p], lw=2.0,
                       solid_capstyle="butt")
    for p in phis:
        ax[0].plot([], [], color=col[p], lw=2, label=rf"$\varphi$ = {p}")
    ax[0].set(xlim=(100, 917), ylim=(0, 2.4), xlabel=r"density $\rho$ (kg m$^{-3}$)",
              ylabel=r"$k_\mathrm{eff}$ (W m$^{-1}$ K$^{-1}$)")
    ax[0].xaxis.label.set_size(9); ax[0].yaxis.label.set_size(9)
    ax[0].legend(frameon=False, fontsize=7.5, loc="upper left")
    ax[0].text(0.97, 0.03, "(a)", transform=ax[0].transAxes, ha="right", fontsize=10)
    style(ax[0])

    # (b) fraction of the bound window, in time
    for p in phis:
        g = G[p]
        hi = min(r["t"][-1] for r in g)
        lo_t = max(r["t"][int(np.argmax(r["t"] >= 1.0))] for r in g)
        tt = np.linspace(lo_t, hi, 400)
        H = []
        for r in g:
            l2, u2 = hs_bounds(r["ice"], 2)
            H.append(np.interp(tt, r["t"], (r["kiso"] - l2) / (u2 - l2)))
        H = np.array(H)
        ax[1].fill_between(tt / DAY, H.min(0), H.max(0), color=col[p], alpha=0.15, lw=0)
        ax[1].plot(tt / DAY, H.mean(0), color=col[p], lw=1.6)
        hc = np.mean([x["h_cal2019"] for x in rows if x["phi"] == p])
        ax[1].axhline(hc, color=col[p], lw=1.0, ls="--")
    ax[1].set(xlim=(0, None), ylim=(0, 1), xlabel="time (d)",
              ylabel=r"$h = (k - k_\mathrm{HS}^-)/(k_\mathrm{HS}^+ - k_\mathrm{HS}^-)$")
    ax[1].xaxis.label.set_size(9); ax[1].yaxis.label.set_size(9)
    ax[1].text(0.97, 0.03, "(b)", transform=ax[1].transAxes, ha="right", fontsize=10)
    style(ax[1])
    fig.savefig(out / "calonne_density.png", dpi=300, transparent=True)
    fig.savefig(out / "calonne_density.pdf", transparent=True)

    print(f"{'phi':>6} {'rho':>6} {'k_start':>8} {'k_end':>7} {'Cal11':>6} {'Cal19':>6} "
          f"{'HS2+':>6} {'h_end':>6} {'h_C19':>6} {'end/C19':>7}")
    for p in phis:
        v = [x for x in rows if x["phi"] == p]
        f_ = lambda k: np.mean([x[k] for x in v])
        print(f"{p:6.3f} {f_('rho'):6.0f} {f_('k_start'):8.3f} {f_('k_end'):7.3f} "
              f"{f_('k_cal2011'):6.3f} {f_('k_cal2019'):6.3f} {f_('hs2_up'):6.3f} "
              f"{f_('h_end'):6.3f} {f_('h_cal2019'):6.3f} {f_('k_end_over_cal'):7.3f}")
    print(f"wrote {out}/calonne_density.png, .pdf, calonne_compare.csv")


if __name__ == "__main__":
    main()

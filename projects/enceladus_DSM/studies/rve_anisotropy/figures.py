#!/usr/bin/env python3
"""Figures for the RVE / anisotropy check. Existing data only -- no solves.

    venv_enceladus/bin/python studies/rve_anisotropy/geometry.py     # first
    venv_enceladus/bin/python studies/rve_anisotropy/figures.py

Reads geometry.csv / seam_profile.csv (this directory), the earlier t = 0
finite-volume results in studies/packing_design/rev_bias.csv, and the PetIGA
k_eff CSVs of the pilot (arith, L/R 40), the warm-end batch (tensor, L/R 40)
and rev64 (arith, L/R 64). Writes rve_anisotropy.png and k_data.csv.
"""
from __future__ import annotations

import csv
import glob
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
PROJ = HERE.parents[1]
RES = Path.home() / "SimulationResults/HPC_results/enceladus_DSM"
C = ["#0072B2", "#D55E00", "#009E73", "#CC79A7"]      # Okabe-Ito, figstyle order
INK, MUTED = "#1a1a1a", "#777777"


def geometry():
    return {(r["family"], r["name"]): r for r in csv.DictReader(open(HERE / "geometry.csv"))}


def k_rows(G):
    """One row per (packing, source, law, time): k anisotropy with its geometry."""
    out = []

    def add(src, law, L, name, g, f):
        k = np.atleast_1d(np.genfromtxt(f, delimiter=",", names=True))
        for i, when in ((0, "t0"), (-1, "30d")):
            kx, ky, kxy = k["k_00"][i], k["k_11"][i], k["k_01"][i]
            out.append(dict(src=src, law=law, L_over_R=L, packing=name, when=when,
                            t_d=k["time"][i] / 86400, kxx_kyy=kx / ky,
                            kxy_kiso=kxy / (0.5 * (kx + ky)),
                            F_ratio=float(g["F_ratio"]), F_xy=float(g["F_xy"]),
                            seam=float(g["seam_ratio"])))

    for s in (1, 2, 3, 4):
        g = G[("pilot_LR40", f"pilot_phi0.325_Rave50um_LR40_seed{s}")]
        ft = glob.glob(str(RES / f"GrainPackingSintering/keff_T_warm_phi0.325/*seed{s}*T-20*/k_eff*.csv"))
        if ft:
            add("PetIGA", "tensor", 40, f"pilot seed{s}", g, ft[0])
        elif s == 2:   # no warm-end run; its tensor replay starts at ~16 d
            fr = glob.glob(str(RES / "GrainPackingSintering/batch_2026-09-16*/*seed2*T-20/k_eff_tensor.csv"))
            if fr:
                add("PetIGA", "tensor(16d+)", 40, "pilot seed2", g, fr[0])
        fa = glob.glob(str(RES / f"GrainPackingSintering/batch_2026-09-16*/*seed{s}*T-20/k_eff.csv"))
        if fa and s != 2:
            add("PetIGA", "arith", 40, f"pilot seed{s}", g, fa[0])
    for f in glob.glob(str(RES / "rev64/*/k_eff.csv")):
        add("PetIGA", "arith", 64, "rev64 seed1",
            G[("rev_LR64", "rev_phi0.325_Rave50um_LR64_seed1")], f)
    for r in csv.DictReader(open(PROJ / "studies/packing_design/rev_bias.csv")):
        g = G[("design", r["name"])]
        out.append(dict(src="FV", law="arith", L_over_R=int(r["L_over_R"]), packing=r["name"],
                        when="t0", t_d=0.0, kxx_kyy=1.0 / float(r["k_anis"]), kxy_kiso=np.nan,
                        F_ratio=float(g["F_ratio"]), F_xy=float(g["F_xy"]),
                        seam=float(g["seam_ratio"])))
    return out


def style(ax, xl, yl, title):
    ax.set_xlabel(xl, fontsize=11)
    ax.set_ylabel(yl, fontsize=11)
    ax.set_title(title, fontsize=12, loc="left")
    ax.grid(True, alpha=0.25, lw=0.6)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


def main():
    G = geometry()
    K = k_rows(G)
    with open(HERE / "k_data.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(K[0]))
        w.writeheader()
        w.writerows(K)

    xy = [r for r in G.values() if r["periodic"] == "xy" and abs(float(r["porosity"]) - 0.325) < 0.01
          and not (r["family"] == "design" and r["name"].startswith("LR40"))]   # = prod_LR40_xy
    L = np.array([float(r["L_over_R"]) for r in xy])
    Fr = np.array([float(r["F_ratio"]) for r in xy])
    Fx = np.array([float(r["F_xy"]) for r in xy])
    noseam = [r for r in G.values() if r["family"] in ("prod_LR40_x", "prod_LR40_none")]

    fig, axes = plt.subplots(2, 3, figsize=(16, 9.5))

    # (a) fabric ratio vs domain size
    ax = axes[0, 0]
    ax.axhline(1.0, color=MUTED, lw=0.8, ls=":")
    ax.plot(L, Fr, "o", color=C[0], ms=6, label="xy-periodic packings")
    ax.plot([float(r["L_over_R"]) for r in noseam], [float(r["F_ratio"]) for r in noseam],
            "s", mfc="none", color=C[1], ms=6, label="no y seam (x-periodic / none)")
    for Lc in (10, 20, 40, 64):
        m = np.abs(L - Lc) < 0.3 * Lc
        if m.sum() > 1:
            ax.errorbar(Lc * 1.12, Fr[m].mean(), yerr=Fr[m].std(ddof=1), fmt="_", color=INK,
                        ms=14, capsize=4, lw=1.5)
    style(ax, r"$L/R_\mathrm{ave}$", r"$F_{xx}/F_{yy}$",
          "(a) contact fabric: horizontal bias at every size\n     (black = mean ± sd)")
    ax.set_xscale("log")
    ax.set_xticks([10, 20, 40, 64]); ax.set_xticklabels(["10", "20", "40", "64"])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.legend(fontsize=9, frameon=False, loc="upper right")

    # (b) fabric tilt vs domain size, with 1/L envelope
    ax = axes[0, 1]
    ax.axhline(0.0, color=MUTED, lw=0.8, ls=":")
    ax.plot(L, Fx, "o", color=C[0], ms=6)
    sd40 = Fx[np.abs(L - 40) < 12].std(ddof=1)
    Ls = np.linspace(9, 75, 100)
    ax.fill_between(Ls, -2 * sd40 * 40 / Ls, 2 * sd40 * 40 / Ls, color=C[0], alpha=0.10, lw=0,
                    label=r"$\pm 2\sigma$, scaled from L/R = 40 as $1/L$")
    style(ax, r"$L/R_\mathrm{ave}$", r"$F_{xy}$",
          "(b) fabric tilt: zero mean, shrinks as 1/L\n     (a finite-sample fluctuation)")
    ax.set_xscale("log")
    ax.set_xticks([10, 20, 40, 64]); ax.set_xticklabels(["10", "20", "40", "64"])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.legend(fontsize=9, frameon=False, loc="upper right")

    # (c) seam contact profile
    ax = axes[0, 2]
    prof = {}
    for r in csv.DictReader(open(HERE / "seam_profile.csv")):
        if r["packing"].startswith(("pilot_LR40", "rev_LR64", "prod_LR40_xy")):
            prof.setdefault(r["packing"], []).append((float(r["y_over_L"]), float(r["hline_rel"]),
                                                      float(r["vline_rel"])))
    Y = None
    H, V = [], []
    for v in prof.values():
        a = np.array(v)
        Y = a[:, 0]
        H.append(a[:, 1]); V.append(a[:, 2])
    H, V = np.array(H), np.array(V)
    Ys = np.where(Y > 0.5, Y - 1.0, Y)          # centre the seam at 0
    o = np.argsort(Ys)
    ax.fill_between(Ys[o], H.min(0)[o], H.max(0)[o], color=C[1], alpha=0.15, lw=0)
    ax.plot(Ys[o], H.mean(0)[o], "-", color=C[1], lw=1.8, label="horizontal lines (cross y-seam)")
    ax.plot(Ys[o], V.mean(0)[o], "-", color=C[0], lw=1.2, label="vertical lines (x boundary, control)")
    ax.axhline(1.0, color=MUTED, lw=0.8, ls=":")
    style(ax, "position / L   (0 = periodic boundary)", "contacts crossing a line / interior mean",
          f"(c) the y seam is a contact-poor layer\n     ({len(H)} packings, L/R 40 and 64; band = min–max)")
    ax.legend(fontsize=9, frameon=False, loc="lower right")

    # (d) k anisotropy vs fabric
    ax = axes[1, 0]
    t0 = [r for r in K if r["when"] == "t0"]
    for (src, law), col, mk in ((("FV", "arith"), C[2], "^"), (("PetIGA", "arith"), C[0], "o"),
                                (("PetIGA", "tensor"), C[1], "s")):
        s = [r for r in t0 if r["src"] == src and r["law"] == law]
        if s:
            ax.plot([r["F_ratio"] for r in s], [r["kxx_kyy"] for r in s], mk, color=col, ms=7,
                    mfc="none" if law == "tensor" else col, label=f"{src} {law}, t = 0")
    a = np.array([[r["F_ratio"], r["kxx_kyy"]] for r in t0 if not (r["src"] == "PetIGA" and r["law"] == "arith" and r["L_over_R"] == 40)])
    rr = np.corrcoef(a[:, 0], a[:, 1])[0, 1]
    big = np.array([[r["F_ratio"], r["kxx_kyy"]] for r in t0 if r["L_over_R"] >= 40
                    and not (r["src"] == "PetIGA" and r["law"] == "arith" and r["L_over_R"] == 40)])
    rr40 = np.corrcoef(big[:, 0], big[:, 1])[0, 1]
    ax.axhline(1.0, color=MUTED, lw=0.8, ls=":"); ax.axvline(1.0, color=MUTED, lw=0.8, ls=":")
    style(ax, r"contact fabric $F_{xx}/F_{yy}$", r"$k_{xx}/k_{yy}$",
          f"(d) both lean horizontal, but per seed only at small L\n     (r = {rr:.2f} all; r = {rr40:.2f} at L/R >= 40)")
    ax.legend(fontsize=9, frameon=False, loc="upper left")

    # (e) k_xy sign vs F_xy sign
    ax = axes[1, 1]
    kk = [r for r in K if r["src"] == "PetIGA" and r["when"] == "30d" and r["law"] != "arith"] + \
         [r for r in K if r["packing"] == "rev64 seed1" and r["when"] == "30d"]
    ax.axhline(0, color=MUTED, lw=0.8, ls=":"); ax.axvline(0, color=MUTED, lw=0.8, ls=":")
    for r in kk:
        ax.plot(r["F_xy"], r["kxy_kiso"], "o", color=C[0] if r["L_over_R"] == 40 else C[1], ms=8)
        ax.annotate(r["packing"], (r["F_xy"], r["kxy_kiso"]), xytext=(5, 4),
                    textcoords="offset points", fontsize=9, color="#333333")
    agree = sum(np.sign(r["F_xy"]) == np.sign(r["kxy_kiso"]) for r in kk)
    style(ax, r"contact fabric $F_{xy}$ (from grains.dat)", r"$k_{xy}/k_\mathrm{iso}$ at 30 d (PetIGA)",
          f"(e) k_xy has the sign of each packing's fabric tilt\n     ({agree} of {len(kk)} agree; blue L/R 40, orange 64)")

    # (f) seam vs k anisotropy, directly
    ax = axes[1, 2]
    pts = [r for r in t0 if not (r["src"] == "PetIGA" and r["law"] == "arith" and r["L_over_R"] == 40)]
    for r in pts:
        col = C[2] if r["src"] == "FV" else (C[1] if r["law"] == "tensor" else C[0])
        ax.plot(r["seam"], r["kxx_kyy"], "o", color=col, ms=6,
                mfc="none" if r["law"] == "tensor" else col)
    rs = np.corrcoef([r["seam"] for r in pts], [r["kxx_kyy"] for r in pts])[0, 1]
    ax.axhline(1.0, color=MUTED, lw=0.8, ls=":")
    ax.axvline(1.0, color=MUTED, lw=0.8, ls=":")
    style(ax, "seam contact density / interior", r"$k_{xx}/k_{yy}$ at t = 0",
          f"(f) seam depth barely predicts k anisotropy\n     (r = {rs:.2f}; colours as in (d))")

    fig.suptitle("Is the domain an RVE, and where does the anisotropy come from?  "
                 "(φ = 0.325, existing packings and runs only)", fontsize=14)
    fig.tight_layout()
    fig.savefig(HERE / "rve_anisotropy.png", dpi=150, bbox_inches="tight")
    print(f"  wrote {HERE / 'rve_anisotropy.png'} and k_data.csv")


if __name__ == "__main__":
    main()

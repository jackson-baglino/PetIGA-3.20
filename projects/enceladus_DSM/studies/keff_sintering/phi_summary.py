#!/usr/bin/env python3
"""Porosity-series figures for the k_eff campaign (production packings, one T).

    venv_enceladus/bin/python studies/keff_sintering/phi_summary.py <campaign dir> \\
        [--T -20] [--out <dir>]

Uses the production packings only (L/R 40, seeds 1601-2005) at one temperature
and writes four figures (default <campaign>/compare/phi_summary/T<T>/):

  phi_trends.png     the porosity trends, seed mean +- sd per phi:
                     (a) k_iso at the opening sample, 11 tau_sub and the end
                     (b) the sintering rise k_end/k_11 - 1
                     (c) k_xx/k_yy at t = 0 and the end, with the packings'
                         contact-fabric ratio F_xx/F_yy
                     (d) the SSA sensitivity d ln k / d ln SSA after 11 tau_sub
                         (least-squares slope of ln k_iso on ln SSA, per seed)
  seeds_by_phi.png   every packing's k_iso vs SSA (after 11 tau_sub), one panel
                     per porosity, so outliers and spread are visible
  anisotropy_time.png  k_xx/k_yy vs time per porosity (seed mean, band =
                     seed min..max), on common times only
  ssa_by_phi.png     SSA, SSA/SSA_0 and k/k_0 against sintering age, by porosity

Colour: porosity on the map compare_keff.py defines (black, then viridis blue to lime).
Also writes phi_trends.csv.
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
from plot_keff import load, read_tau_sub  # noqa: E402
from compare_keff import CMAP  # noqa: E402

DAY = 86400.0
PAT = re.compile(r"phi([\d.]+)_Rave50um_LR40_seed(\d+)_L2mm_eps1000nm_perxy_T(-?\d+)__")
INK, MUTED = "#1a1a1a", "#6b6b6b"


def style(ax):
    ax.grid(True, alpha=0.25)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root", type=Path)
    ap.add_argument("--T", type=int, default=-20)
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args()
    out = a.out or a.root / "compare" / "phi_summary" / f"T{a.T}"
    out.mkdir(parents=True, exist_ok=True)
    fab = {int(r["base_seed"]): float(r["F_ratio"]) for r in
           csv.DictReader(open(PROJ / "inputs/packings/keff_LR40/packings_summary.csv"))}

    G = defaultdict(list)
    for kf in sorted(a.root.glob("packing_*/k_eff.csv")):
        m = PAT.search(kf.parent.name)
        if not m or int(m.group(3)) != a.T or not 1601 <= int(m.group(2)) <= 2005:
            continue
        r = load(kf.parent)
        r["tau"] = read_tau_sub(kf.parent)
        r["seed"] = int(m.group(2))
        G[float(m.group(1))].append(r)
    phis = sorted(G)
    if not phis:
        raise SystemExit(f"no production runs at T = {a.T}")
    cm, (c0, c1) = CMAP["phi"]
    col = {p: cm(c0 + (c1 - c0) * i / max(1, len(phis) - 1)) for i, p in enumerate(phis)}

    # ---- per-seed metrics ----
    rows = []
    for p in phis:
        for r in G[p]:
            t, k = r["t"], r["kiso"]
            i0 = int(np.argmax(t >= 1.0))
            t11 = 11.05 * r["tau"]
            live = t >= t11
            slope = np.polyfit(np.log(r["ssa"][live]), np.log(k[live]), 1)[0]
            rows.append(dict(phi=p, seed=r["seed"], k_0=k[i0], k_11=float(np.interp(t11, t, k)),
                             k_end=k[-1], rise_pct=100 * (k[-1] / np.interp(t11, t, k) - 1),
                             kr_0=r["kxx"][i0] / r["kyy"][i0], kr_end=r["kxx"][-1] / r["kyy"][-1],
                             F_ratio=fab.get(r["seed"], np.nan), dlnk_dlnssa=slope))
    with open(out / "phi_trends.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

    def ms(key):
        m_, s_ = [], []
        for p in phis:
            v = np.array([x[key] for x in rows if x["phi"] == p])
            m_.append(v.mean()); s_.append(v.std(ddof=1) if len(v) > 1 else 0.0)
        return np.array(m_), np.array(s_)

    # ---- figure 1: trends ----
    fig, ax = plt.subplots(2, 2, figsize=(10.5, 7.6), constrained_layout=True)
    X = np.array(phis)
    for key, lab, ls, mk in (("k_0", "opening sample", ":", "o"),
                             ("k_11", r"11 $\tau_\mathrm{sub}$ (1 d)", "--", "s"),
                             ("k_end", "30 d", "-", "^")):
        m_, s_ = ms(key)
        ax[0, 0].errorbar(X, m_, yerr=s_, fmt=mk, ls=ls, color=INK, ms=6, capsize=3, lw=1.4, label=lab)
    ax[0, 0].set(ylabel=r"$k_\mathrm{iso}$  [W m$^{-1}$ K$^{-1}$]", title="(a) conductivity")
    m_, s_ = ms("rise_pct")
    ax[0, 1].errorbar(X, m_, yerr=s_, fmt="o-", color=INK, ms=6, capsize=3, lw=1.6)
    for p in phis:
        v = [x["rise_pct"] for x in rows if x["phi"] == p]
        ax[0, 1].scatter([p] * len(v), v, s=14, color=col[p], alpha=0.7, zorder=0)
    ax[0, 1].set(ylabel=r"rise  $k_{30\,\mathrm{d}}/k_{11\tau}-1$  [%]", title="(b) sintering rise")
    for key, lab, mk, ls in (("kr_0", r"$k_{xx}/k_{yy}$, t = 0", "o", ":"),
                             ("kr_end", r"$k_{xx}/k_{yy}$, 30 d", "^", "-"),
                             ("F_ratio", r"contact fabric $F_{xx}/F_{yy}$", "D", "--")):
        m_, s_ = ms(key)
        ax[1, 0].errorbar(X, m_, yerr=s_, fmt=mk, ls=ls, color=INK if key != "F_ratio" else MUTED,
                          ms=6, capsize=3, lw=1.4, label=lab)
    ax[1, 0].axhline(1, color="#bbbbbb", lw=0.8)
    ax[1, 0].set(ylabel="x / y ratio", title="(c) anisotropy")
    m_, s_ = ms("dlnk_dlnssa")
    ax[1, 1].errorbar(X, m_, yerr=s_, fmt="o-", color=INK, ms=6, capsize=3, lw=1.6)
    ax[1, 1].set(ylabel=r"$d\ln k_\mathrm{iso}\,/\,d\ln\mathrm{SSA}$  (after 11 $\tau_\mathrm{sub}$)",
                 title="(d) sensitivity of k to SSA")
    for x in ax.flat:
        style(x); x.set_xlabel(r"porosity $\varphi$"); x.set_xticks(phis)
    ax[0, 0].legend(frameon=False); ax[1, 0].legend(frameon=False)
    fig.suptitle(f"Porosity series, T = {a.T} °C, L/R 40, {len(rows)} packings "
                 "(seed mean ± sd)", fontsize=12)
    fig.savefig(out / "phi_trends.png", dpi=160)

    # ---- figure 2: every seed, k vs SSA ----
    fig, ax = plt.subplots(1, len(phis), figsize=(3.1 * len(phis), 3.6), sharey=False,
                           constrained_layout=True)
    for x, p in zip(np.atleast_1d(ax), phis):
        for r in G[p]:
            live = r["t"] >= 11.05 * r["tau"]
            x.plot(r["ssa"][live], r["kiso"][live], color=col[p], lw=1.4)
            x.text(r["ssa"][live][-1], r["kiso"][live][-1], f" {r['seed']}", fontsize=7,
                   color=MUTED, va="center", ha="right")
        x.invert_xaxis(); style(x)
        x.set_title(f"φ = {p}  ({len(G[p])} packings)", loc="left", fontsize=10)
        x.set_xlabel(r"SSA  [m$^{-1}$]")
    np.atleast_1d(ax)[0].set_ylabel(r"$k_\mathrm{iso}$  [W m$^{-1}$ K$^{-1}$]")
    fig.suptitle(f"Every packing, k_iso vs SSA after 11 τ_sub, T = {a.T} °C (time runs left)",
                 fontsize=11)
    fig.savefig(out / "seeds_by_phi.png", dpi=160)

    # ---- figure 3: anisotropy vs time ----
    fig, ax = plt.subplots(figsize=(8, 4.6), constrained_layout=True)
    for p in phis:
        g = G[p]
        hi = min(r["t"][-1] for r in g)
        lo = max(r["t"][int(np.argmax(r["t"] >= 1.0))] for r in g)
        tt = np.linspace(lo, hi, 300)
        R = np.array([np.interp(tt, r["t"], r["kxx"] / r["kyy"]) for r in g])
        ax.fill_between(tt / DAY, R.min(0), R.max(0), color=col[p], alpha=0.12, lw=0)
        ax.plot(tt / DAY, R.mean(0), color=col[p], lw=2, label=f"φ = {p}")
    ax.axhline(1, color="#bbbbbb", lw=0.8)
    ax.set(xlabel="time [d]", ylabel=r"$k_{xx}/k_{yy}$",
           title=f"Anisotropy vs time, T = {a.T} °C (seed mean, band = seed min–max)")
    style(ax); ax.legend(frameon=False, ncol=5, fontsize=8, loc="lower center")
    fig.savefig(out / "anisotropy_time.png", dpi=160)

    # ---- figure 4: SSA in time, by porosity ----
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6), constrained_layout=True)
    for p in phis:
        g = G[p]
        th_lo = max(r["t"][int(np.argmax(r["t"] >= 1.0))] / r["tau"] for r in g)
        th_hi = min(r["t"][-1] / r["tau"] for r in g)
        th = np.geomspace(max(th_lo, 1e-3), th_hi, 300)
        def curves(key, norm):
            out_ = []
            for r in g:
                i0 = int(np.argmax(r["t"] >= 1.0))
                v = r[key] / (r[key][i0] if norm else 1.0)
                out_.append(np.interp(th, r["t"] / r["tau"], v))
            return np.array(out_)
        for j, (key, norm) in enumerate((("ssa", False), ("ssa", True), ("kiso", True))):
            V = curves(key, norm)
            ax[j].fill_between(th, V.min(0), V.max(0), color=col[p], alpha=0.12, lw=0)
            ax[j].plot(th, V.mean(0), color=col[p], lw=2, label=f"φ = {p}")
    ax[0].set(ylabel=r"SSA  [m$^{-1}$]", title="(a) SSA")
    ax[1].set(ylabel=r"SSA / SSA$_0$", title="(b) SSA, normalized")
    ax[2].set(ylabel=r"$k_\mathrm{iso}/k_{\mathrm{iso},0}$", title="(c) k_iso, normalized, same axis")
    for x in ax:
        x.set_xscale("log"); x.set_xlim(1, None); style(x)
        x.set_xlabel(r"sintering age  $\theta=t/\tau_\mathrm{sub}$")
    ax[0].legend(frameon=False, fontsize=8)
    fig.suptitle(f"SSA and k_iso against sintering age, T = {a.T} °C (seed mean, band = seed min–max)", fontsize=12)
    fig.savefig(out / "ssa_by_phi.png", dpi=160)

    print(f"{'phi':>6} {'n':>2} {'dlnk/dlnSSA':>14}")
    m_, s_ = ms("dlnk_dlnssa")
    for p, mm, ss in zip(phis, m_, s_):
        print(f"{p:6.3f} {len(G[p]):2d} {mm:8.3f}±{ss:.3f}")
    print(f"wrote {out}/phi_trends.png, seeds_by_phi.png, anisotropy_time.png, phi_trends.csv")


if __name__ == "__main__":
    main()

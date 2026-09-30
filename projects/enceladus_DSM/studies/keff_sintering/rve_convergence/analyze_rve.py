#!/usr/bin/env python3
"""Domain-size convergence of k_eff (supplement): one curve as L grows?

    venv_enceladus/bin/python studies/keff_sintering/rve_convergence/analyze_rve.py \\
        <batch dirs or run dirs ...> [--T -20] [--out <dir>]

Finds every run of phi 0.325 at temperature --T whose name carries LR<N>.
The convergence ensemble is the UNGATED rve_phi0.325 set (seed numbers >= 1100)
at L/R 20/30/40/56/80. The production keff_LR40 packings (seeds < 1000, built
WITH the homogeneity gates) are kept apart as "production 40" and drawn as a
separate dashed curve and point, so what the production recipe does to k_eff
is visible against the unfiltered ensemble at the same size. Draws:

  (a) k_iso vs SSA: seed mean per L/R, band = seed min..max
  (b) k_iso vs t/tau_sub, the same
  (c) the gap between each size's seed mean and the largest size's, at three
      states (t = 11, 100, 330 tau_sub), with the standard error of the gap
  (d) seed-to-seed CV of k_iso at t_final vs L/R, against the 1/L reference

"Converged" is judged in (c): the production size is close enough when its
gap to the largest domain is inside that gap's standard error. SSA is length
per cell area, so it is intensive and comparable across sizes.
Writes rve_convergence.png and rve_convergence.csv (or --out).
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
sys.path.insert(0, str(HERE.parents[2] / "postprocess"))
from plot_keff import DAY, load, read_tau_sub      # noqa: E402

STATES_TAU = (11.05, 100.0, 330.0)
N_GRID = 300


def discover(roots, T):
    runs = defaultdict(list)
    for root in roots:
        for kf in Path(root).rglob("k_eff.csv"):
            d = kf.parent
            n = d.name
            m_lr = re.search(r"LR(\d+)_", n)
            m_T = re.search(r"_T(-?\d+)__", n)
            if not (m_lr and m_T and "phi0.325" in n and int(m_T.group(1)) == T):
                continue
            r = load(d)
            tau = read_tau_sub(d)
            if r is None or not tau:
                continue
            r["tau"] = tau
            r["name"] = n
            seed = int(re.search(r"seed(\d+)", n).group(1))
            key = int(m_lr.group(1)) if seed >= 1100 else "production 40"
            runs[key].append(r)
    ens = {k: v for k, v in sorted(((k, v) for k, v in runs.items() if k != "production 40"))}
    return ens, runs.get("production 40", [])


def seed_curves(g, x_of, n=N_GRID):
    """Seed curves resampled on a common grid of x (t/tau or SSA)."""
    lo = max(np.min(x_of(r)) for r in g)
    hi = min(np.max(x_of(r)) for r in g)
    xs = np.linspace(lo, hi, n)
    K = []
    for r in g:
        x = x_of(r)
        o = np.argsort(x)
        K.append(np.interp(xs, x[o], r["kiso"][o]))
    return xs, np.array(K)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("roots", nargs="+", type=Path)
    ap.add_argument("--T", type=int, default=-20)
    ap.add_argument("--out", type=Path, default=HERE)
    a = ap.parse_args()
    runs, prod = discover(a.roots, a.T)
    if not runs:
        raise SystemExit("no ungated phi 0.325 LR<N> runs (seed >= 1100) with k_eff.csv found")
    for L, g in runs.items():
        print(f"  L/R {L:3d}: {len(g)} seeds (ungated)")
    print(f"  production L/R 40 (gated): {len(prod)} seeds")
    Ls = list(runs)
    cmap = plt.get_cmap("viridis")
    col = {L: cmap(0.1 + 0.8 * i / max(1, len(Ls) - 1)) for i, L in enumerate(Ls)}

    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    rows = []
    # (a) vs SSA after the baseline, (b) vs t/tau
    for L, g in runs.items():
        live = lambda r: r["t"] >= STATES_TAU[0] * r["tau"]
        xs, K = seed_curves([{**r, "kiso": r["kiso"][live(r)], "ssa": r["ssa"][live(r)]} for r in g],
                            lambda r: r["ssa"])
        axes[0, 0].fill_between(xs, K.min(0), K.max(0), color=col[L], alpha=0.15, lw=0)
        axes[0, 0].plot(xs, K.mean(0), "-", color=col[L], lw=2, label=f"L/R = {L} ({len(g)} seeds)")
        xs, K = seed_curves(g, lambda r: r["t"] / r["tau"])
        axes[0, 1].fill_between(xs, K.min(0), K.max(0), color=col[L], alpha=0.15, lw=0)
        axes[0, 1].plot(xs, K.mean(0), "-", color=col[L], lw=2, label=f"L/R = {L}")

    if prod:
        live = lambda r: r["t"] >= STATES_TAU[0] * r["tau"]
        xs, K = seed_curves([{**r, "kiso": r["kiso"][live(r)], "ssa": r["ssa"][live(r)]} for r in prod],
                            lambda r: r["ssa"])
        axes[0, 0].plot(xs, K.mean(0), "--", color="#c0392b", lw=1.6,
                        label=f"production L/R 40, gated ({len(prod)})")
        xs, K = seed_curves(prod, lambda r: r["t"] / r["tau"])
        axes[0, 1].plot(xs, K.mean(0), "--", color="#c0392b", lw=1.6, label="production 40, gated")

    # (c) gap to the largest size, (d) CV at t_final
    Lref = Ls[-1]
    def at(g, st):
        return np.array([np.interp(st * r["tau"], r["t"], r["kiso"]) for r in g])
    for st, mk in zip(STATES_TAU, ("o", "s", "^")):
        ref = at(runs[Lref], st)
        gaps, errs = [], []
        for L in Ls:
            v = at(runs[L], st)
            gap = (v.mean() / ref.mean() - 1) * 100
            se = 100 * np.sqrt((v.std(ddof=1) ** 2 / len(v) if len(v) > 1 else 0)
                               + (ref.std(ddof=1) ** 2 / len(ref) if len(ref) > 1 else 0)) / ref.mean()
            gaps.append(gap); errs.append(se)
            rows.append(dict(L_over_R=L, state_tau=st, n=len(v), kiso_mean=v.mean(),
                             kiso_sd=v.std(ddof=1) if len(v) > 1 else np.nan,
                             gap_to_ref_pct=gap, gap_se_pct=se))
        axes[1, 0].errorbar(Ls, gaps, yerr=errs, fmt=mk + "-", capsize=4,
                            label=f"t = {st:g} τ_sub")
        if prod:
            v = at(prod, st)
            gp = (v.mean() / ref.mean() - 1) * 100
            sp = 100 * np.sqrt(v.std(ddof=1) ** 2 / len(v) + ref.std(ddof=1) ** 2 / len(ref)) / ref.mean() \
                if len(v) > 1 and len(ref) > 1 else 0.0
            axes[1, 0].errorbar([41.5], [gp], yerr=[sp], fmt=mk, mfc="none", color="#c0392b", capsize=3)
            rows.append(dict(L_over_R="production 40", state_tau=st, n=len(v), kiso_mean=v.mean(),
                             kiso_sd=v.std(ddof=1) if len(v) > 1 else np.nan,
                             gap_to_ref_pct=gp, gap_se_pct=sp))
    cv = []
    for L in Ls:
        v = np.array([r["kiso"][-1] for r in runs[L]])
        cv.append(100 * v.std(ddof=1) / v.mean() if len(v) > 1 else np.nan)
    axes[1, 1].plot(Ls, cv, "o-", color="#1a1a1a", label="measured")
    if not np.isnan(cv[Ls.index(40)] if 40 in Ls else np.nan):
        ref40 = cv[Ls.index(40)]
        Lg = np.linspace(min(Ls), max(Ls), 50)
        axes[1, 1].plot(Lg, ref40 * 40 / Lg, "--", color="#888888", label="∝ 1/L (through L/R = 40)")

    t = {(0, 0): ("(a) k_iso vs SSA, after 11 τ_sub", r"SSA  [m$^{-1}$]", r"$k_\mathrm{iso}$  [W m$^{-1}$ K$^{-1}$]"),
         (0, 1): ("(b) k_iso vs t/τ_sub", r"$t\,/\,\tau_\mathrm{sub}$", r"$k_\mathrm{iso}$  [W m$^{-1}$ K$^{-1}$]"),
         (1, 0): (f"(c) seed-mean gap to L/R = {Lref}", "L / R_ave", "gap  [%]"),
         (1, 1): ("(d) seed-to-seed CV of k_iso at t_final", "L / R_ave", "CV  [%]")}
    for (i, j), (ti, xl, yl) in t.items():
        ax = axes[i, j]
        ax.set_title(ti, loc="left", fontsize=12); ax.set_xlabel(xl); ax.set_ylabel(yl)
        ax.grid(True, alpha=0.25)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.legend(frameon=False, fontsize=9)
    axes[1, 0].axhline(0, color="#888888", lw=0.8, ls=":")
    axes[1, 0].axvline(40, color="#c0392b", lw=0.8, ls=":")
    axes[1, 1].axvline(40, color="#c0392b", lw=0.8, ls=":")
    fig.suptitle(f"k_eff domain-size convergence, φ = 0.325, T = {a.T} °C "
                 "(ungated ensembles; red: production L/R 40, gated)", fontsize=14)
    fig.tight_layout()
    a.out.mkdir(parents=True, exist_ok=True)
    fig.savefig(a.out / "rve_convergence.png", dpi=150, bbox_inches="tight")
    with open(a.out / "rve_convergence.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    for r in rows:
        if r["L_over_R"] == 40:
            ok = abs(r["gap_to_ref_pct"]) <= r["gap_se_pct"]
            print(f"  L/R 40 at {r['state_tau']:g} tau: gap {r['gap_to_ref_pct']:+.2f}% "
                  f"± {r['gap_se_pct']:.2f}%  -> {'within' if ok else 'OUTSIDE'} the standard error")
    print(f"  wrote {a.out / 'rve_convergence.png'} and rve_convergence.csv")


if __name__ == "__main__":
    main()

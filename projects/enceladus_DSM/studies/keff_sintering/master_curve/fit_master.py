#!/usr/bin/env python3
"""First pass at the closure k_iso = F(theta, connectivity).

    venv_enceladus/bin/python studies/keff_sintering/master_curve/fit_master.py <campaign dir>
        [--theta-ref 30] [--phi-max 0.375]

theta = t / tau_sub(T) is the sintering age. For every production run (L/R 40,
seeds 1601-2005, every temperature) and theta >= theta_ref:

    k_iso(theta) / k_iso(theta_ref) = 1 + a * log10(theta / theta_ref)      (log law)
    k_iso(theta) / k_iso(theta_ref) = (theta / theta_ref)^n                  (power law)

theta_ref = 30 by default: after the early transient, which is interface-width
dependent (eps_sensitivity/README.md) and happens in unresolved necks.

Per PACKING (all its temperatures pooled, since they coincide) the fit gives a
growth rate a and a level k_ref = k_iso(theta_ref). Both are then regressed on
the packing's descriptors -- porosity and contacts per grain at the band
(z_band) -- over phi <= --phi-max, where the 2D solid is well connected.

    F:  k_iso = k_ref(phi) * [ 1 + a * log10(theta / theta_ref) ]

Writes master_curve.png, master_packings.csv, master_fit.txt here.
"""
from __future__ import annotations

import argparse, csv, re, sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
PROJ = HERE.parents[2]
sys.path.insert(0, str(PROJ / "postprocess"))
from plot_keff import load, read_tau_sub  # noqa: E402
from compare_keff import CMAP  # noqa: E402

PAT = re.compile(r"phi([\d.]+)_Rave50um_LR40_seed(\d+)_L2mm_eps1000nm_perxy_T(-?\d+)__")
K_ICE = 2.29


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root", type=Path)
    ap.add_argument("--theta-ref", type=float, default=30.0)
    ap.add_argument("--phi-max", type=float, default=0.375)
    a = ap.parse_args()
    desc = {int(r["base_seed"]): r for r in
            csv.DictReader(open(PROJ / "inputs/packings/keff_LR40/packings_summary.csv"))}
    runs = defaultdict(list)
    for kf in sorted(a.root.glob("packing_*/k_eff.csv")):
        m = PAT.search(kf.parent.name)
        if not m or not 1601 <= int(m.group(2)) <= 2005:
            continue
        r = load(kf.parent)
        th = r["t"] / read_tau_sub(kf.parent)
        if th[-1] < 1.5 * a.theta_ref:
            continue
        kref = float(np.interp(a.theta_ref, th, r["kiso"]))
        w = th >= a.theta_ref
        runs[(float(m.group(1)), int(m.group(2)))].append(
            dict(T=int(m.group(3)), x=np.log10(th[w] / a.theta_ref), y=r["kiso"][w] / kref, kref=kref))

    rows = []
    for (phi, seed), g in sorted(runs.items()):
        x = np.concatenate([q["x"] for q in g]); y = np.concatenate([q["y"] for q in g])
        al = float(np.sum(x * (y - 1)) / np.sum(x * x))                 # log law through (0, 1)
        n = float(np.sum(x * np.log10(y)) / np.sum(x * x))              # power law through (0, 1)
        rows.append(dict(phi=phi, seed=seed, n_T=len(g), theta_max=float(a.theta_ref * 10 ** x.max()),
                         a_log=al, n_pow=n,
                         rms_log_pct=float(100 * np.sqrt(np.mean((y - 1 - al * x) ** 2))),
                         rms_pow_pct=float(100 * np.sqrt(np.mean((y - (10 ** x) ** n) ** 2))),
                         k_ref=float(np.mean([q["kref"] for q in g])),
                         kref_spread_pct=float(100 * np.ptp([q["kref"] for q in g]) / np.mean([q["kref"] for q in g])),
                         z_band=float(desc[seed]["z_band"])))
    with open(HERE / "master_packings.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

    A = lambda k, sel=None: np.array([r[k] for r in rows if sel is None or sel(r)])
    good = lambda r: r["phi"] <= a.phi_max
    out = []
    P = lambda s: (print(s), out.append(s))
    P(f"theta_ref = {a.theta_ref:g}; {len(rows)} packings, {sum(r['n_T'] for r in rows)} runs; "
      f"fit quality per packing (rms of k/k_ref): log law {A('rms_log_pct').mean():.2f}%, "
      f"power law {A('rms_pow_pct').mean():.2f}%")
    P(f"{'phi':>6} {'n':>2} {'a (per decade)':>16} {'n (power)':>14} {'k_ref [W/m/K]':>16} {'z_band':>7}")
    for p in sorted({r["phi"] for r in rows}):
        s = lambda r, p=p: r["phi"] == p
        P(f"{p:6.3f} {len(A('a_log', s)):2d} {A('a_log', s).mean():8.4f}±{A('a_log', s).std(ddof=1):.4f} "
          f"{A('n_pow', s).mean():7.4f}±{A('n_pow', s).std(ddof=1):.4f} "
          f"{A('k_ref', s).mean():8.3f}±{A('k_ref', s).std(ddof=1):.3f} {A('z_band', s).mean():7.2f}")
    ag = A("a_log", good)
    P(f"\nwell-connected set (phi <= {a.phi_max}, {len(ag)} packings): a = {ag.mean():.4f} ± {ag.std(ddof=1):.4f} "
      f"(se {ag.std(ddof=1) / np.sqrt(len(ag)):.4f});  n = {A('n_pow', good).mean():.4f}")
    for key in ("phi", "z_band"):
        c = np.corrcoef(A(key, good), ag)[0, 1]
        ca = np.corrcoef(A(key), A("a_log"))[0, 1]
        P(f"  corr(a, {key}): {c:+.2f} within phi <= {a.phi_max}; {ca:+.2f} over all porosities")
    # level: ln k_ref linear in phi (all), and in z
    cf = np.polyfit(A("phi", good), np.log(A("k_ref", good) / K_ICE), 1)
    res = np.log(A("k_ref", good) / K_ICE) - np.polyval(cf, A("phi", good))
    P(f"  level (phi <= {a.phi_max}): k_ref/k_ice = {np.exp(cf[1]):.3f} * exp({cf[0]:.2f} * phi); "
      f"packing scatter about it {100 * res.std(ddof=2):.1f}%")
    cz = np.corrcoef(res, A("z_band", good) - np.polyval(np.polyfit(A("phi", good), A("z_band", good), 1), A("phi", good)))[0, 1]
    P(f"  does z_band explain the scatter at fixed phi?  corr(level residual, z_band residual) = {cz:+.2f}")
    a0 = ag.mean()
    P(f"\nF (first pass, valid for theta >= {a.theta_ref:g}, phi <= {a.phi_max}):")
    P(f"   k_iso = {np.exp(cf[1]):.3f} k_ice exp({cf[0]:.2f} phi) * [1 + {a0:.3f} log10(theta/{a.theta_ref:g})]")
    # how well does F predict every run?
    err = []
    for (phi, seed), g in runs.items():
        if phi > a.phi_max:
            continue
        for q in g:
            pred = np.exp(cf[1]) * K_ICE * np.exp(cf[0] * phi) * (1 + a0 * q["x"])
            err.append(100 * (pred / (q["y"] * q["kref"]) - 1))
    err = np.concatenate(err)
    P(f"   absolute k_iso predicted for every sample of every run: rms {np.sqrt(np.mean(err ** 2)):.1f}% "
      f"(the packing-to-packing scatter); with each packing's own k_ref: "
      f"{np.mean(A('rms_log_pct', good)):.2f}% + {100 * ag.std(ddof=1) / 1:.1f}%-of-a spread")
    (HERE / "master_fit.txt").write_text("\n".join(out) + "\n")

    phis = sorted({r["phi"] for r in rows})
    cm, (c0, c1) = CMAP["phi"]
    col = {p: cm(c0 + (c1 - c0) * i / max(1, len(phis) - 1)) for i, p in enumerate(phis)}
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6), constrained_layout=True)
    for (phi, seed), g in runs.items():
        for q in g:
            ax[0].plot(a.theta_ref * 10 ** q["x"], q["y"], color=col[phi], lw=0.9, alpha=0.75)
    xx = np.linspace(0, np.log10(1300 / a.theta_ref), 50)
    ax[0].plot(a.theta_ref * 10 ** xx, 1 + a0 * xx, "k--", lw=1.6, label=f"1 + {a0:.3f} log10(θ/{a.theta_ref:g})")
    for p in phis:
        ax[0].plot([], [], color=col[p], lw=2, label=f"φ = {p}")
    ax[0].set(xscale="log", xlabel=r"sintering age  $\theta=t/\tau_\mathrm{sub}$",
              ylabel=rf"$k_\mathrm{{iso}}(\theta)\,/\,k_\mathrm{{iso}}(\theta={a.theta_ref:g})$",
              title="(a) every run, every temperature")
    ax[0].legend(frameon=False, fontsize=8)
    for p in phis:
        s = lambda r, p=p: r["phi"] == p
        ax[1].scatter(A("z_band", s), A("a_log", s), color=col[p], s=34, edgecolor="white", lw=0.6)
        ax[2].scatter(A("phi", s), A("k_ref", s) / K_ICE, color=col[p], s=34, edgecolor="white", lw=0.6)
    ax[1].axhline(a0, color="k", ls="--", lw=1)
    ax[1].set(xlabel="contacts per grain at the band, $z$", ylabel="growth per decade of age, $a$",
              title="(b) growth rate per packing")
    pp = np.linspace(min(phis), max(phis), 50)
    ax[2].plot(pp, np.exp(np.polyval(cf, pp)), "k--", lw=1, label=f"fit, φ ≤ {a.phi_max}")
    ax[2].axvspan(a.phi_max + 0.025, max(phis) + 0.02, color="#f2f2f2", zorder=0, lw=0)
    ax[2].set(xlabel=r"porosity $\varphi$", ylabel=rf"level  $k_\mathrm{{iso}}(\theta={a.theta_ref:g})/k_\mathrm{{ice}}$",
              title="(c) level per packing", yscale="log")
    ax[2].legend(frameon=False, fontsize=8)
    for x in ax:
        x.grid(True, alpha=0.25)
        for sp in ("top", "right"):
            x.spines[sp].set_visible(False)
    fig.savefig(HERE / "master_curve.png", dpi=150)
    print(f"wrote {HERE}/master_curve.png, master_packings.csv, master_fit.txt")


if __name__ == "__main__":
    main()

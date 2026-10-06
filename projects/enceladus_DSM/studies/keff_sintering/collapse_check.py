#!/usr/bin/env python3
"""Temperature collapse, packing by packing: is T only a time rescaling?

    venv_enceladus/bin/python studies/keff_sintering/collapse_check.py <campaign dir> [--ref -20]

For every production packing that ran at the reference temperature and at
another one, over the SSA range both runs cover after 11 tau_sub:

  dk      max |k_iso(T) / k_iso(ref) - 1| at matched SSA     -> same PATH?
  ratio   median of t_ref(SSA) / t_T(SSA)                     -> the speed-up
  tau     tau_sub(ref) / tau_sub(T)                           -> what the model predicts
  dev     ratio / tau - 1

Prints one line per (phi, T) with the worst packing, writes collapse_check.csv
(every packing) and collapse_check.png into <campaign>/compare/collapse/.
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
sys.path.insert(0, str(HERE.parents[1] / "postprocess"))
from plot_keff import load, read_tau_sub  # noqa: E402
from compare_keff import CMAP  # noqa: E402

PAT = re.compile(r"phi([\d.]+)_Rave50um_LR40_seed(\d+)_L2mm_eps1000nm_perxy_T(-?\d+)__")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root", type=Path)
    ap.add_argument("--ref", type=int, default=-20)
    a = ap.parse_args()
    R = {}
    for kf in sorted(a.root.glob("packing_*/k_eff.csv")):
        m = PAT.search(kf.parent.name)
        if not m or not 1601 <= int(m.group(2)) <= 2005:
            continue
        r = load(kf.parent)
        r["tau"] = read_tau_sub(kf.parent)
        R[(float(m.group(1)), int(m.group(2)), int(m.group(3)))] = r
    rows = []
    for (phi, seed, T), r in sorted(R.items()):
        if T == a.ref or (phi, seed, a.ref) not in R:
            continue
        ref = R[(phi, seed, a.ref)]
        la, lb = ref["t"] >= 11.05 * ref["tau"], r["t"] >= 11.05 * r["tau"]
        lo = max(ref["ssa"].min(), r["ssa"].min())
        hi = min(ref["ssa"][la].max(), r["ssa"][lb].max())
        if hi <= lo:
            continue
        S = np.linspace(lo, hi, 300)
        f = lambda q, key: np.interp(S, q["ssa"][::-1], q[key][::-1])
        dk = (f(r, "kiso") / f(ref, "kiso") - 1) * 100
        ratio = f(ref, "t") / f(r, "t")
        rows.append(dict(phi=phi, seed=seed, T=T, ssa_lo=lo, ssa_hi=hi,
                         dk_max_pct=float(np.abs(dk).max()), dk_mean_pct=float(dk.mean()),
                         time_ratio=float(np.median(ratio)), tau_ratio=ref["tau"] / r["tau"],
                         dev_pct=float((np.median(ratio) / (ref["tau"] / r["tau"]) - 1) * 100)))
    out = a.root / "compare" / "collapse"
    out.mkdir(parents=True, exist_ok=True)
    with open(out / "collapse_check.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    G = defaultdict(list)
    for x in rows:
        G[(x["phi"], x["T"])].append(x)
    print(f"reference T = {a.ref} C; {len(rows)} paired runs")
    print(f"{'phi':>6} {'T':>4} {'n':>2} {'max dk %':>9} {'time ratio':>11} {'tau ratio':>10} {'dev %':>12}")
    for (phi, T), g in sorted(G.items()):
        dv = [x["dev_pct"] for x in g]
        print(f"{phi:6.3f} {T:4d} {len(g):2d} {max(x['dk_max_pct'] for x in g):9.3f} "
              f"{np.mean([x['time_ratio'] for x in g]):11.3f} {g[0]['tau_ratio']:10.3f} "
              f"{np.mean(dv):+6.2f} ({min(dv):+.2f}..{max(dv):+.2f})")
    print(f"worst k mismatch at matched SSA: {max(x['dk_max_pct'] for x in rows):.3f}%   "
          f"worst rate deviation from tau_sub: {max(abs(x['dev_pct']) for x in rows):.2f}%")

    phis = sorted({x["phi"] for x in rows}); Ts = sorted({x["T"] for x in rows})
    cm, (c0, c1) = CMAP["phi"]
    col = {p: cm(c0 + (c1 - c0) * i / max(1, len(phis) - 1)) for i, p in enumerate(phis)}
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.4), constrained_layout=True)
    for i, p in enumerate(phis):
        g = [x for x in rows if x["phi"] == p]
        off = (i - (len(phis) - 1) / 2) * 0.6
        ax[0].scatter([x["T"] + off for x in g], [x["dk_max_pct"] for x in g], color=col[p], s=26, label=f"φ = {p}")
        ax[1].scatter([x["T"] + off for x in g], [x["dev_pct"] for x in g], color=col[p], s=26)
    ax[0].set(xlabel="T [°C]", ylabel="max |Δk| at matched SSA  [%]", title=f"(a) same path as {a.ref} °C?")
    ax[1].axhline(0, color="#888", lw=0.8)
    ax[1].set(xlabel="T [°C]", ylabel=r"speed-up vs $\tau_\mathrm{sub}$ ratio  [%]", title="(b) rate = τ_sub ratio?")
    for x in ax:
        x.set_xticks(Ts); x.grid(True, alpha=0.25)
        for sp in ("top", "right"):
            x.spines[sp].set_visible(False)
    ax[0].legend(frameon=False, fontsize=8)
    fig.suptitle("Temperature collapse, every packing (one point per packing)")
    fig.savefig(out / "collapse_check.png", dpi=160)
    print(f"wrote {out}/collapse_check.png and .csv")


if __name__ == "__main__":
    main()

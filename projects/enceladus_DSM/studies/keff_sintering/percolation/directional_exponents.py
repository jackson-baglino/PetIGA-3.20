#!/usr/bin/env python3
"""directional_exponents.py — the state law for k_xx and k_yy separately.

    venv_enceladus/bin/python studies/keff_sintering/percolation/directional_exponents.py

Figure 7 fits k_eff = (k_xx + k_yy)/2 against SSA. At porosity 0.425 and 0.475
the packings are anisotropic (k_xx < k_yy; y is the direction of deposition),
so this asks whether each direction follows its own power law. Same fit as
the figure (sample_figures.fig_porosity): after the reference age theta = 30,
log(k/k_r) = p log(SSA/SSA_r) through the origin, temperatures -30 to -5 C
pooled per packing.

Writes directional_exponents.txt beside this file.
"""
from __future__ import annotations
import re, sys
from collections import defaultdict
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO / "postprocess"))
from plot_keff import load, read_tau_sub                              # noqa: E402

CAMP = Path.home() / "SimulationResults/HPC_results/enceladus_DSM/keff_sintering_campaign"
PAT = re.compile(r"phi([0-9.]+)_Rave50um_LR40_seed(\d+)_L2mm_eps1000nm_perxy_T(-?\d+)__")
TH_R = 30.0


def main():
    per = defaultdict(lambda: defaultdict(list))
    for kf in sorted(CAMP.glob("packing_*/k_eff.csv")):
        m = PAT.search(kf.parent.name)
        if not m or not 1601 <= int(m.group(2)) <= 2005:
            continue
        r = load(kf.parent)
        th = r["t"] / read_tau_sub(kf.parent)
        if th[-1] < 1.5 * TH_R:
            continue
        w = th >= TH_R
        X = np.log(r["ssa"][w] / np.interp(TH_R, th, r["ssa"]))
        for key in ("kxx", "kyy", "kiso"):
            per[(float(m.group(1)), int(m.group(2)))][key].append((X, np.log(r[key][w] / np.interp(TH_R, th, r[key]))))
        per[(float(m.group(1)), int(m.group(2)))]["aniso"].append(np.interp(TH_R, th, r["kxx"] / r["kyy"]))
    by = defaultdict(list)
    for (p, seed), d in per.items():
        e = {}
        for key in ("kxx", "kyy", "kiso"):
            X = np.concatenate([a for a, _ in d[key]]); Y = np.concatenate([b for _, b in d[key]])
            e[key] = float(np.sum(X * Y) / np.sum(X * X))
        by[p].append((e["kxx"], e["kyy"], e["kiso"], float(np.mean(d["aniso"]))))
    L = ["Exponent p of k ~ SSA^p after theta = 30, per packing (temperatures pooled); mean +- sd over the five packings.", "",
         "phi      k_xx            k_yy            (k_xx+k_yy)/2     k_xx/k_yy at theta = 30    corr(p_xx - p_yy, k_xx/k_yy)"]
    allv = []
    for p, v in sorted(by.items()):
        v = np.array(v); allv.append(v)
        L.append(f"{p:.3f}   {v[:,0].mean():+.2f} +- {v[:,0].std(ddof=1):.2f}   {v[:,1].mean():+.2f} +- {v[:,1].std(ddof=1):.2f}   "
                 f"{v[:,2].mean():+.2f} +- {v[:,2].std(ddof=1):.2f}      {v[:,3].mean():.2f} +- {v[:,3].std(ddof=1):.2f}             "
                 f"{np.corrcoef(v[:,0]-v[:,1], v[:,3])[0,1]:+.2f}")
    a = np.vstack(allv)
    L += ["", f"all 25 packings: corr(p_xx - p_yy, k_xx/k_yy) = {np.corrcoef(a[:,0]-a[:,1], a[:,3])[0,1]:+.2f}",
          f"                 corr(p of the weaker direction, min/max k) = "
          f"{np.corrcoef(np.where(a[:,3] < 1, a[:,0], a[:,1]), np.minimum(a[:,3], 1/a[:,3]))[0,1]:+.2f}",
          f"                 exponent of the STRONGER direction: {np.where(a[:,3] < 1, a[:,1], a[:,0]).mean():+.2f} +- {np.where(a[:,3] < 1, a[:,1], a[:,0]).std(ddof=1):.2f}",
          f"                 exponent of the WEAKER direction:   {np.where(a[:,3] < 1, a[:,0], a[:,1]).mean():+.2f} +- {np.where(a[:,3] < 1, a[:,0], a[:,1]).std(ddof=1):.2f}"]
    for p, v in sorted(by.items()):
        v = np.array(v); s = np.where(v[:,3] < 1, v[:,1], v[:,0]); wk = np.where(v[:,3] < 1, v[:,0], v[:,1])
        L.append(f"   phi {p:.3f}: stronger {s.mean():+.2f} +- {s.std(ddof=1):.2f}   weaker {wk.mean():+.2f} +- {wk.std(ddof=1):.2f}")
    txt = "\n".join(L)
    print(txt)
    (HERE / "directional_exponents.txt").write_text(txt + "\n")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Off-diagonal k_xy vs domain size: a zero-mean fluctuation that shrinks ~1/L?

    venv_enceladus/bin/python studies/keff_sintering/rve_convergence/kxy_vs_L.py <run dirs ...>

phi 0.325, -20 C, ungated seeds (>= 1100) at L/R 20-80. For every run,
k_xy/k_iso at the opening sample and at the end; per L/R the seed mean (should
be ~0) and the RMS (should fall ~1/L). Writes kxy_vs_L.png/.csv here.
"""
import csv, re, sys
from collections import defaultdict
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2] / "postprocess"))
from plot_keff import load  # noqa: E402

G = defaultdict(list)
for d in map(Path, sys.argv[1:]):
    m = re.search(r"phi0\.325_Rave50um_LR(\d+)_seed(\d+)_.*eps1000nm.*_T-20__", d.name)
    if not m or int(m.group(2)) < 1100 or not (d / "k_eff.csv").is_file():
        continue
    r = load(d)
    i0 = int(np.argmax(r["t"] >= 1.0))
    G[int(m.group(1))].append((r["kxy"][i0] / r["kiso"][i0], r["kxy"][-1] / r["kiso"][-1]))
Ls = sorted(G)
rows = []
for L in Ls:
    a = 100 * np.array(G[L])
    rows.append(dict(L_over_R=L, n=len(a), mean_0=a[:, 0].mean(), rms_0=np.sqrt((a[:, 0] ** 2).mean()),
                     mean_end=a[:, 1].mean(), rms_end=np.sqrt((a[:, 1] ** 2).mean())))
    print(f"L/R {L:3d} n={len(a)}  k_xy/k_iso  t0: mean {rows[-1]['mean_0']:+.2f}% rms {rows[-1]['rms_0']:.2f}%   "
          f"end: mean {rows[-1]['mean_end']:+.2f}% rms {rows[-1]['rms_end']:.2f}%")
with open(HERE / "kxy_vs_L.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
fig, ax = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
for L in Ls:
    a = 100 * np.array(G[L])
    ax[0].scatter([L] * len(a), a[:, 1], color="#0072B2", s=22, alpha=0.8)
ax[0].plot(Ls, [r["mean_end"] for r in rows], "k-o", ms=5, label="seed mean")
ax[0].axhline(0, color="#888", lw=0.8)
ax[0].set(xlabel="L / R_ave", ylabel="k_xy / k_iso at 30 d  [%]", title="(a) every seed: zero-mean")
ax[1].plot(Ls, [r["rms_end"] for r in rows], "k-o", label="RMS, 30 d")
ax[1].plot(Ls, [r["rms_0"] for r in rows], "o:", color="#888", label="RMS, t = 0")
ref = rows[Ls.index(40)]["rms_end"] if 40 in Ls else rows[0]["rms_end"]
Lg = np.linspace(min(Ls), max(Ls), 50)
ax[1].plot(Lg, ref * 40 / Lg, "--", color="#D55E00", label="∝ 1/L (through L/R 40)")
ax[1].set(xlabel="L / R_ave", ylabel="RMS k_xy / k_iso  [%]", title="(b) magnitude shrinks with L")
for x in ax:
    x.grid(True, alpha=0.25)
    for s in ("top", "right"):
        x.spines[s].set_visible(False)
    x.legend(frameon=False)
fig.suptitle("Off-diagonal conductivity vs domain size, φ = 0.325, −20 °C (ungated)")
fig.savefig(HERE / "kxy_vs_L.png", dpi=160)
print("wrote", HERE / "kxy_vs_L.png")

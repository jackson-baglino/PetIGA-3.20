#!/usr/bin/env python3
"""eps x2 sensitivity: does the k_eff trajectory depend on the interface width?

    venv_enceladus/bin/python studies/keff_sintering/eps_sensitivity/compare_eps.py <campaign dir>

phi 0.325, -20 C, seeds 1701-1703, run at eps = 1 um (production) and at
eps = 2 um (batch_eps2.txt) with the SAME experiment file, so dtmax is the same
in seconds. Paired per packing, from t = 1 d on:

  k_iso(t)      eps 2 um interpolated onto the 1 um run's samples
  SSA(t)        the same
  rise          k(30 d)/k(1 d) - 1, both
  k_iso(SSA)    on the common SSA range

Writes eps_sensitivity.png and eps_sensitivity.csv next to this file.
"""
import csv, sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2] / "postprocess"))
from plot_keff import load  # noqa: E402

DAY = 86400.0
C = Path(sys.argv[1])
rows = []
fig, ax = plt.subplots(1, 3, figsize=(15, 4.6), constrained_layout=True)
cols = ["#0072B2", "#D55E00", "#009E73"]
for seed, c in zip((1701, 1702, 1703), cols):
    a = load(next(C.glob(f"packing_*seed{seed}_L2mm_eps1000nm_perxy_T-20__*")))
    b = load(next(C.glob(f"packing_*seed{seed}_L2mm_eps2.00um_perxy_T-20__*")))
    m = a["t"] >= DAY
    kb = np.interp(a["t"][m], b["t"], b["kiso"]); sb = np.interp(a["t"][m], b["t"], b["ssa"])
    dk = (kb / a["kiso"][m] - 1) * 100; ds = (sb / a["ssa"][m] - 1) * 100
    k1 = lambda r: np.interp(DAY, r["t"], r["kiso"])
    te = min(a["t"][-1], b["t"][-1])
    ke = lambda r: np.interp(te, r["t"], r["kiso"])
    ra, rb = 100 * (ke(a) / k1(a) - 1), 100 * (ke(b) / k1(b) - 1)
    i0a, i0b = int(np.argmax(a["t"] >= 1.0)), int(np.argmax(b["t"] >= 1.0))
    rows.append(dict(seed=seed, k0_eps1=a["kiso"][i0a], k0_eps2=b["kiso"][i0b],
                     k0_diff_pct=100 * (b["kiso"][i0b] / a["kiso"][i0a] - 1),
                     dk_t_mean_pct=dk.mean(), dk_t_max_pct=np.abs(dk).max(),
                     dssa_t_mean_pct=ds.mean(), rise_eps1=ra, rise_eps2=rb, rise_diff_pts=rb - ra,
                     ssa0_diff_pct=100 * (b["ssa"][i0b] / a["ssa"][i0a] - 1)))
    ax[0].plot(a["t"] / DAY, a["kiso"], color=c, lw=2, label=f"seed {seed}, ε = 1 µm")
    ax[0].plot(b["t"] / DAY, b["kiso"], color=c, lw=1.6, ls="--", label=f"seed {seed}, ε = 2 µm")
    ax[1].plot(a["ssa"], a["kiso"], color=c, lw=2); ax[1].plot(b["ssa"], b["kiso"], color=c, lw=1.6, ls="--")
    ax[2].plot(a["t"][m] / DAY, dk, color=c, lw=1.8, label=f"seed {seed}")
ax[0].set(xlabel="time [d]", ylabel=r"$k_\mathrm{iso}$ [W m$^{-1}$ K$^{-1}$]", title="(a) k_iso vs time")
ax[1].set(xlabel=r"SSA [m$^{-1}$]", ylabel=r"$k_\mathrm{iso}$", title="(b) k_iso vs SSA"); ax[1].invert_xaxis()
ax[2].axhline(0, color="#888", lw=0.8)
ax[2].set(xlabel="time [d]", ylabel="k_iso(ε = 2 µm) / k_iso(ε = 1 µm) − 1  [%]", title="(c) paired difference, t ≥ 1 d")
for x in ax:
    x.grid(True, alpha=0.25)
    for s in ("top", "right"):
        x.spines[s].set_visible(False)
ax[0].legend(frameon=False, fontsize=7, ncol=2); ax[2].legend(frameon=False)
fig.savefig(HERE / "eps_sensitivity.png", dpi=150)
with open(HERE / "eps_sensitivity.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
print(f"{'seed':>5} {'k0 diff':>8} {'SSA0 diff':>9} {'k(t) mean':>10} {'k(t) max':>9} {'SSA(t) mean':>11} {'rise 1um':>9} {'rise 2um':>9} {'diff pts':>9}")
for r in rows:
    print(f"{r['seed']:5d} {r['k0_diff_pct']:+7.2f}% {r['ssa0_diff_pct']:+8.2f}% {r['dk_t_mean_pct']:+9.2f}% {r['dk_t_max_pct']:8.2f}% "
          f"{r['dssa_t_mean_pct']:+10.2f}% {r['rise_eps1']:8.1f}% {r['rise_eps2']:8.1f}% {r['rise_diff_pts']:+9.2f}")
d = np.array([r["rise_diff_pts"] for r in rows])
print(f"paired rise difference: {d.mean():+.2f} ± {d.std(ddof=1):.2f} points (seed scatter of the rise itself ~2.3)")

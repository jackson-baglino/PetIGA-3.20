"""Time-lapse of one type-6 run for the run plan: phi = 1/2 and the ice field at
chosen days, with the pore walls from the mesh. Checked with figcheck at the
placed width (6.5 in).

    venv_cryosuction/bin/python docs/run_plan/t6_timelapse.py <run folder> <out.png> <day> [<day> ...]
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "postprocess"))
import pplib  # noqa: E402
from figcheck import check_figure  # noqa: E402

run, out, days = sys.argv[1], sys.argv[2], [float(d) for d in sys.argv[3:]]
snaps = list(pplib.iter_snapshots(run))
ts = np.array([s.t for s in snaps]) / 86400
pick = [int(np.argmin(abs(ts - d))) for d in days]
plt.rcParams.update({"font.family": "Arial", "font.size": 7.5})
fig, axs = plt.subplots(len(pick), 1, figsize=(6.5, 0.62 * len(pick) + 0.45), sharex=True,
                        constrained_layout=True)
for ax, i in zip(axs, pick):
    s = snaps[i]
    ax.contourf(s.X * 1e6, s.Y * 1e6, s.phi, levels=[0.5, 2], colors=["#bfe0f5"])
    ax.contour(s.X * 1e6, s.Y * 1e6, s.phi, levels=[0.5], colors=["#1f6fb2"], linewidths=0.8)
    k = 0 if np.ptp(s.Y[0, :]) < np.ptp(s.Y[:, 0]) else 1      # index along the walls
    for j in (0, -1):
        xw, yw = (s.X[j, :], s.Y[j, :]) if k == 0 else (s.X[:, j], s.Y[:, j])
        ax.plot(xw * 1e6, yw * 1e6, color="#8c7350", lw=1.4)
    for side in ("top", "bottom", "right"):
        ax.spines[side].set_visible(False)
    ax.set_aspect("equal")
    ax.set_xlim(0, s.X.max() * 1e6)
    ax.set_yticks([0, 150])
    ax.text(1.01, 0.5, f"{ts[i]:.0f} d", transform=ax.transAxes, va="center")
axs[-1].set_xlabel("x [µm]   (warm end x = 0)")
bad = check_figure(fig, min_pt=7.0, placed_width_in=6.5)
print("figcheck:", bad or "OK")
fig.savefig(out, dpi=250)
sys.exit(1 if bad else 0)

#!/usr/bin/env python3
"""Does the homogenization assumption survive the run? Measure the PORE.

    venv_enceladus/bin/python studies/keff_sintering/pore_scale_evolution.py <run_dir>

THE POINT. The REV requirement is L >> (every heterogeneity length in the
medium), and the design so far has sized L against the GRAIN radius only. But
a two-phase medium has two structures, and the pore is the one that coarsens
fastest: grains merge, their pores merge with them, and a pore network can
develop a feature that spans a good fraction of the domain long before the
grains do. A thin pore cutting across the cell breaks homogenization
completely, and nothing in a grain-size criterion would see it coming.

So this measures, per snapshot:

  xi_pore, xi_solid   correlation lengths from the two-point autocovariance of
                      each phase, radially averaged. The FFT is periodic, which
                      is exactly right for these cells. xi is where the
                      autocovariance first falls to 1/e of its zero-lag value.
                      THIS, not the grain radius, is the length L must dwarf.
  max pore radius     the largest circle that fits in the pore, from a periodic
                      distance transform. Catches a single large void that a
                      correlation length would average away.
  L/xi_pore           the homogenization ratio, tracked in time.
  open-pore fraction  pore further than band/2 from ice -- the part vapour can
                      actually use -- and how many pieces it is in.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import ndimage

HERE = Path(__file__).parent
PROJ = HERE.parents[1]
sys.path.insert(0, str(PROJ / "preprocess"))
sys.path.insert(0, str(PROJ / "postprocess"))
import figstyle as fs                                       # noqa: E402
import pplib                                                # noqa: E402
import packing_lib as pl                                    # noqa: E402

DAY = 86400.0


def corr_length(mask: np.ndarray, dx: float) -> float:
    """Correlation length from the radially averaged autocovariance.

    Periodic FFT, which matches the cell. Returns the lag at which the
    autocovariance first drops to 1/e of its zero-lag value -- the standard
    definition, and the length the domain must be large compared to.
    """
    f = mask.astype(float)
    f -= f.mean()
    if not np.any(f):
        return 0.0
    F = np.fft.rfft2(f)
    ac = np.fft.irfft2(F * np.conj(F), s=mask.shape).real / f.size
    ac = np.fft.fftshift(ac)
    ny, nx = ac.shape
    yy, xx = np.indices(ac.shape)
    r = np.hypot(yy - ny // 2, xx - nx // 2).astype(int)
    prof = np.bincount(r.ravel(), ac.ravel()) / np.maximum(np.bincount(r.ravel()), 1)
    zero = prof[0]
    if zero <= 0:
        return 0.0
    below = np.flatnonzero(prof < zero / np.e)
    return float(below[0] * dx) if below.size else float(len(prof) * dx)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run", type=Path)
    ap.add_argument("--n", type=int, default=24, help="snapshots to sample")
    ap.add_argument("--eps", type=float, default=1.0e-6)
    ap.add_argument("--L", type=float, default=2.0e-3)
    a = ap.parse_args()
    band = 9.2 * a.eps

    vts = sorted((a.run / "vtkOut").glob("solV_*.vts"),
                 key=lambda p: int(re.search(r"solV_(\d+)", p.name).group(1)))
    if not vts:
        raise SystemExit(f"no vtkOut/*.vts under {a.run}")
    pick = [vts[i] for i in np.unique(np.linspace(0, len(vts) - 1, a.n).astype(int))]
    ssa = pplib.load_ssa(str(a.run))
    step2t = dict(zip(ssa[:, 3].astype(int), ssa[:, 2]))

    print(f"{'t[d]':>7} {'xi_pore':>9} {'ice_w50':>9} {'pore_w90':>8} {'L/xi_p':>8} "
          f"{'maxpore':>9} {'open%':>7} {'pieces':>7}")
    rows = []
    for f in pick:
        step = int(re.search(r"solV_(\d+)", f.name).group(1))
        t = step2t.get(step, np.nan)
        fl, X, _ = pplib.read_vts(f, want=["IcePhase"])
        phi = fl["IcePhase"]
        dx = a.L / phi.shape[1]
        solid = phi >= 0.5
        pore = ~solid
        # ONE correlation length, not two. For a two-phase medium the pore and
        # the solid have the SAME autocovariance -- pore = 1 - solid, so the
        # fluctuation fields differ only in sign and the autocovariance is
        # identical. Reporting both looked like agreement between independent
        # measurements; it was the same number twice.
        xp = corr_length(pore, dx)
        dist = pl.periodic_edt(solid, True, True) * dx          # pore half-widths
        dsol = pl.periodic_edt(pore, True, True) * dx           # solid half-widths
        maxpore = float(dist.max())
        # these DO differ between the phases, which is the point: the pore and
        # the ice coarsen at different rates
        xs = float(np.percentile(dsol[solid], 50)) * 2.0        # median solid width
        p90 = float(np.percentile(dist[pore], 90)) * 2.0        # wide-pore width
        op = pore & (dist > 0.5 * band)
        _, _, _, npiece = pl.percolates(op, diagonal=True)
        openfrac = op.sum() / max(pore.sum(), 1)
        rows.append((t / DAY, xp, xs, a.L / max(xp, 1e-12), maxpore, openfrac,
                     npiece, p90))
        print(f"{t/DAY:7.2f} {xp*1e6:9.2f} {xs*1e6:9.2f} {p90*1e6:8.2f} "
              f"{a.L/max(xp,1e-12):8.1f} {maxpore*1e6:9.2f} {openfrac:7.1%} {npiece:7d}")

    r = np.array(rows)
    fig, ax = plt.subplots(1, 3, figsize=(10.0, 3.8))
    ax[0].plot(r[:, 0], r[:, 1] * 1e6, "-o", color=fs.C[1], ms=4, label="pore")
    ax[0].plot(r[:, 0], r[:, 2] * 1e6, "-s", color=fs.C[0], ms=4,
               label="median ice width")
    ax[0].plot(r[:, 0], r[:, 7] * 1e6, "-^", color=fs.C[2], ms=4,
               label="90th-pct pore width")
    fs.style(ax[0], "time [days]", r"correlation length $\xi$  [$\mu$m]",
             "(a)  both phases coarsen", logy=False)
    ax[0].legend(fontsize=fs.FS_LEG, frameon=False)

    ax[1].plot(r[:, 0], r[:, 3], "-o", color=fs.C[1], ms=4)
    ax[1].axhline(10, color=fs.MUTED, ls="--", lw=1.2)
    ax[1].text(r[0, 0], 10.6, r"$L/\xi = 10$, a common REV floor",
               fontsize=fs.FS_NOTE, color=fs.MUTED)
    fs.style(ax[1], "time [days]", r"$L/\xi_{\rm pore}$",
             "(b)  homogenization margin", logy=False)
    ax[1].set_ylim(bottom=0)

    ax[2].plot(r[:, 0], r[:, 5] * 100, "-o", color=fs.C[1], ms=4, label="open pore %")
    a2 = ax[2].twinx()
    a2.plot(r[:, 0], r[:, 6], "-s", color=fs.C[2], ms=4)
    a2.set_ylabel("disconnected pieces", fontsize=fs.FS_LABEL, color=fs.C[2])
    a2.tick_params(labelsize=fs.FS_TICK, labelcolor=fs.C[2])
    a2.spines[["top"]].set_visible(False)
    fs.style(ax[2], "time [days]", "pore vapour can use  [%]",
             "(c)  and fragments further", logy=False)
    fig.tight_layout()
    fs.save(fig, HERE, "pore_scale", dpi=190)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

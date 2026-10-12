#!/usr/bin/env python3
"""percolation.py — does the ice (and the pore) span the periodic cell?

    venv_enceladus/bin/python studies/keff_sintering/percolation/percolation.py \
        [<campaign dir>] [--T -20] [--also-final -5] [--stride 2]

For each production run at temperature --T, and each snapshot stored locally
(four per run: the opening frame and about 10, 19 and 30 d), the phase field
is thresholded at phi = 0.5 and the connected clusters are found.

A phase "wraps" in x when one of its clusters is connected to its own periodic
image one cell to the right, i.e. heat (ice) or vapor (pore) can cross the
cell in x without leaving the phase. Tested on a 3 x 3 tiling of the periodic
cell: the cluster holding pixel (i, j) of the centre tile must also hold
(i + N, j). Same in y. This is the property that matters for a periodic cell;
"touches both edges" is weaker and is not used.

Connectivity: 8-neighbour for the ice and 4-neighbour for the pore, the
complementary pair on a square lattice (with 4 and 4 a diagonal pinch would
cut both phases; with 8 and 8 it would join both).

The snapshot is a PETSc Vec: 8-byte header, then (phi, T, rho_v) per control
point, x fastest; 2831 points a side of which the first 2829 are one period.

--also-final T2 adds the last snapshot of the runs at T2 (the most sintered
state available: -5 C ends at t/tau_sub = 1256 against 331 at -20 C).

Writes percolation.csv, percolation.txt and percolation.png beside this file.
"""
from __future__ import annotations
import argparse, csv, re, sys
from collections import defaultdict
from pathlib import Path
import numpy as np
from scipy import ndimage
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO / "postprocess"))
import pplib                                                        # noqa: E402
from compare_keff import CMAP                                        # noqa: E402

TARGETS_D = (0.0, 10.0, 19.0, 30.0)                                   # the gallery's instants
S8 = np.ones((3, 3), int)                                            # 8-neighbour
S4 = ndimage.generate_binary_structure(2, 1)                         # 4-neighbour
RUN = re.compile(r"packing_2D_phi([0-9.]+)_Rave50um_LR40_seed(\d+)_L2mm_eps1000nm_perxy_T(-?\d+)__snow_T-?\d+_h1\.00_30d$")


def read_phi(f: Path, stride: int) -> np.ndarray:
    a = np.fromfile(f, dtype=">f8", offset=8)
    n = int(round((a.size / 3) ** 0.5))
    return a.reshape(n, n, 3)[: n - 2 : stride, : n - 2 : stride, 0]  # [y, x], one period


def wraps(mask: np.ndarray, structure) -> tuple[bool, bool, float]:
    """(wraps in x, wraps in y, fraction of the phase in wrapping clusters)."""
    ny, nx = mask.shape
    lab, _ = ndimage.label(np.tile(mask, (3, 3)), structure=structure)
    c = lab[ny:2 * ny, nx:2 * nx]
    on = c > 0
    wx = on & (c == lab[ny:2 * ny, 2 * nx:3 * nx])
    wy = on & (c == lab[2 * ny:3 * ny, nx:2 * nx])
    ids = np.union1d(np.unique(c[wx]), np.unique(c[wy]))
    frac = float(np.isin(c, ids).sum() / max(on.sum(), 1))
    return bool(wx.any()), bool(wy.any()), frac


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("campaign", nargs="?", type=Path,
                    default=Path.home() / "SimulationResults/HPC_results/enceladus_DSM/keff_sintering_campaign")
    ap.add_argument("--T", type=int, default=-20)
    ap.add_argument("--also-final", type=int, default=-5)
    ap.add_argument("--stride", type=int, default=2)
    a = ap.parse_args()

    rows = []
    for d in sorted(a.campaign.iterdir()):
        m = RUN.match(d.name)
        if not m:
            continue
        phi0, seed, T = float(m.group(1)), int(m.group(2)), int(m.group(3))
        if T not in (a.T, a.also_final) or seed < 1000:      # seeds < 1000: earlier, gated packings
            continue
        sols = sorted(d.glob("sol_*.dat"))
        tm = pplib.step_times(str(d))
        if T == a.also_final and T != a.T:
            sols = sols[-1:]
        elif len(sols) > len(TARGETS_D):                     # runs that kept every snapshot: the four instants only
            td = np.array([tm.get(int(re.search(r"sol_(\d+)", f.name).group(1)), np.nan) / 86400.0 for f in sols])
            sols = [sols[int(np.nanargmin(np.abs(td - x)))] for x in TARGETS_D]
        tau = None
        mm = re.search(r"^\s*tau_sub\s+([0-9.eE+-]+)\s*s", (d / "outp.txt").read_text(errors="replace"), re.M)
        if mm:
            tau = float(mm.group(1))
        for f in sols:
            step = int(re.search(r"sol_(\d+)", f.name).group(1))
            t = tm.get(step, float("nan"))
            ice = read_phi(f, a.stride) > 0.5
            ix, iy, fi = wraps(ice, S8)
            px, py, fp = wraps(~ice, S4)
            rows.append(dict(phi=phi0, seed=seed, T=T, step=step, t_d=t / 86400.0,
                             theta=(t / tau if tau else float("nan")), porosity=1 - ice.mean(),
                             ice_x=int(ix), ice_y=int(iy), ice_frac=fi, pore_x=int(px), pore_y=int(py), pore_frac=fp))
            r = rows[-1]
            print(f"phi {phi0:.3f} seed {seed} T {T:>3} t {r['t_d']:5.1f} d  ice x{ix:d} y{iy:d} ({fi:.2f})  "
                  f"pore x{px:d} y{py:d} ({fp:.2f})", flush=True)
    with open(HERE / "percolation.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

    # ---- summary: packings (of 5) in which each phase wraps, by porosity and instant ----
    L = [f"Packings (of 5 per porosity) in which a phase wraps the periodic cell; phi = 0.5 threshold, stride {a.stride}.", ""]
    phis = sorted({r["phi"] for r in rows})
    for T in sorted({r["T"] for r in rows}, reverse=True):
        sub = [r for r in rows if r["T"] == T]
        by = defaultdict(list)
        for r in sub:
            by[(r["phi"], r["seed"])].append(r)
        nsnap = max(len(v) for v in by.values())
        for k in range(nsnap):
            grp = [sorted(v, key=lambda r: r["step"])[k] for v in by.values() if len(v) > k]
            td = np.mean([r["t_d"] for r in grp]); th = np.nanmean([r["theta"] for r in grp])
            L.append(f"T = {T} C, snapshot {k + 1} (t = {td:.1f} d, t/tau_sub = {th:.0f})")
            L.append("  porosity   ice x  ice y  ice both   pore x  pore y  pore both   ice in wrapping cluster")
            for p in phis:
                g = [r for r in grp if r["phi"] == p]
                if not g:
                    continue
                L.append(f"  {p:.3f}      {sum(r['ice_x'] for r in g)}/{len(g)}    {sum(r['ice_y'] for r in g)}/{len(g)}    "
                         f"{sum(r['ice_x'] and r['ice_y'] for r in g)}/{len(g)}        "
                         f"{sum(r['pore_x'] for r in g)}/{len(g)}     {sum(r['pore_y'] for r in g)}/{len(g)}     "
                         f"{sum(r['pore_x'] and r['pore_y'] for r in g)}/{len(g)}         "
                         f"{np.mean([r['ice_frac'] for r in g]):.2f}")
            L.append("")
    txt = "\n".join(L)
    print(txt)
    (HERE / "percolation.txt").write_text(txt + "\n")

    # ---- figure: fraction of packings whose ice wraps, against time, per porosity ----
    sub = [r for r in rows if r["T"] == a.T]
    cm, (lo, hi) = CMAP["phi"]
    col = {p: cm(lo + (hi - lo) * i / max(1, len(phis) - 1)) for i, p in enumerate(phis)}
    fig, ax = plt.subplots(1, 3, figsize=(10.5, 3.4), sharey=True)
    for j, (key, ttl) in enumerate((("ice_x", "ice wraps in $x$"), ("ice_y", "ice wraps in $y$"), ("pore_any", "pore wraps in $x$ or $y$"))):
        for p in phis:
            by = defaultdict(list)
            for r in sub:
                if r["phi"] == p:
                    by[r["seed"]].append(r)
            n = min(len(v) for v in by.values())
            t = [np.mean([sorted(v, key=lambda r: r["step"])[k]["t_d"] for v in by.values()]) for k in range(n)]
            f = [np.mean([(sorted(v, key=lambda r: r["step"])[k]["pore_x"] or sorted(v, key=lambda r: r["step"])[k]["pore_y"])
                          if key == "pore_any" else sorted(v, key=lambda r: r["step"])[k][key] for v in by.values()]) for k in range(n)]
            ax[j].plot(t, f, "-o", color=col[p], ms=4, lw=1.4, label=f"{p:g}")
        ax[j].set(title=ttl, xlabel="$t$ [d]", ylim=(-0.05, 1.05))
    ax[0].set_ylabel("fraction of the 5 packings")
    ax[2].legend(title=r"$\varphi$", frameon=False, fontsize=8)
    fig.tight_layout(); fig.savefig(HERE / "percolation.png", dpi=200)
    print(f"wrote {HERE}/percolation.csv, percolation.txt, percolation.png")


if __name__ == "__main__":
    main()

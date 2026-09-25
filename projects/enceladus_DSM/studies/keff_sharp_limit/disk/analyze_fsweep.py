#!/usr/bin/env python3
"""How the two conductivity laws behave as the solid fraction rises.

    venv_enceladus/bin/python studies/keff_sharp_limit/disk/analyze_fsweep.py <batch_dir>

The single-cylinder ladder in this directory fixes the area fraction at
f = 0.196 and varies eps. That establishes the ORDER of each law's error but
says nothing about how it behaves as the solid crowds -- and the packing the
campaign runs sits at an ice fraction of 0.675, far denser than any single
isolated cylinder.

This sweeps BOTH: R = 125..400 um in a 1 mm periodic cell (f = 0.049..0.503)
x eps/R = 0.01, 0.02, 0.04 x {arithmetic, tensor}. 36 `-keff_only` solves.
Reference is Rayleigh's square array (keff_disk_analytic.k_sharp), which is
exact to O(f^4) and good to 1.8e-4 relative at f = 0.196.

WHAT IT SHOWS. The arithmetic law's error is not a fixed offset: it grows
steeply with f, because the diffuse band bridges an ever larger share of the
shrinking gap between neighbours. The tensor law's does not.

Writes fsweep.csv and fsweep.png next to this script.
"""
from __future__ import annotations

import argparse
import csv
import glob
import os
import re
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import keff_disk_analytic as ana                          # noqa: E402

L_DEFAULT = 1.0e-3
RUN_RE = re.compile(r"__(\w+?)_R(\d+)_e([\d.]+)$")


def collect(batch: Path, L: float) -> list:
    rows = []
    for d in sorted(glob.glob(str(batch / "singleice*"))):
        m = RUN_RE.search(os.path.basename(d))
        if not m:
            continue
        law, R_um, eps_over_R = m.group(1), float(m.group(2)), float(m.group(3))
        csvs = glob.glob(os.path.join(d, "k_eff*.csv"))
        if not csvs:
            continue
        r = list(csv.DictReader(open(csvs[0])))[-1]
        # (k_00+k_11)/2: the cell is isotropic by symmetry, so the average is
        # the scalar conductivity and their difference is a solver check.
        k = 0.5 * (float(r["k_00"]) + float(r["k_11"]))
        aniso = abs(float(r["k_00"]) - float(r["k_11"])) / k
        R = R_um * 1e-6
        f = ana.area_fraction(R, L)
        rows.append(dict(law=law, R_um=R_um, eps_over_R=eps_over_R, f=f,
                         k=k, k_sharp=ana.k_sharp(f), aniso=aniso,
                         phi_bar=float(r["phi_bar"]), its=int(float(r["ksp_its"]))))
    for r in rows:
        r["err_pct"] = 100.0 * (r["k"] / r["k_sharp"] - 1.0)
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("batch", type=Path)
    ap.add_argument("--L", type=float, default=L_DEFAULT)
    ap.add_argument("--out", type=Path, default=HERE)
    args = ap.parse_args()

    rows = collect(args.batch, args.L)
    if not rows:
        print(f"no runs found under {args.batch}", file=sys.stderr)
        return 1

    eps_vals = sorted({r["eps_over_R"] for r in rows})
    fs = sorted({r["f"] for r in rows})

    print(f"{'R[um]':>6} {'f':>7} {'eps/R':>6} {'k_sharp':>9} "
          f"{'arith err':>10} {'tensor err':>11} {'tensor/arith':>13} {'|k00-k11|/k':>12}")
    worst_aniso = 0.0
    for f in fs:
        for e in eps_vals:
            a = next((r for r in rows if r["f"] == f and r["eps_over_R"] == e
                      and r["law"].startswith("arith")), None)
            t = next((r for r in rows if r["f"] == f and r["eps_over_R"] == e
                      and r["law"] == "tensor"), None)
            if not a or not t:
                continue
            worst_aniso = max(worst_aniso, a["aniso"], t["aniso"])
            print(f"{a['R_um']:6.0f} {f:7.4f} {e:6.2f} {a['k_sharp']:9.5f} "
                  f"{a['err_pct']:+9.2f}% {t['err_pct']:+10.3f}% "
                  f"{abs(t['err_pct']/a['err_pct']):13.3f} {max(a['aniso'],t['aniso']):12.2e}")
    print(f"\nworst off-isotropy |k00-k11|/k over all 36 runs: {worst_aniso:.2e}")

    out_csv = args.out / "fsweep.csv"
    with out_csv.open("w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    # -- figure -----------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.3))
    colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(eps_vals)))
    for ax, law, ttl in ((axes[0], "arith", "(a) arithmetic law"),
                         (axes[1], "tensor", "(b) tensor law")):
        for c, e in zip(colors, eps_vals):
            sel = sorted([r for r in rows if r["law"].startswith(law)
                          and r["eps_over_R"] == e], key=lambda r: r["f"])
            ax.plot([r["f"] for r in sel], [r["err_pct"] for r in sel],
                    "o-", color=c, lw=1.4, ms=4, label=rf"$\epsilon/R$ = {e:g}")
        ax.axhline(0, color="0.6", lw=0.8, ls=":")
        ax.set_xlabel("solid area fraction $f$")
        ax.set_ylabel(r"$k_{\mathrm{eff}}$ error vs Rayleigh  [%]")
        ax.set_title(ttl, fontsize=10)
        ax.grid(alpha=0.25, lw=0.6)
        ax.legend(frameon=False, fontsize=9)
    # Same y-scale would hide panel (b) entirely: arith reaches +132%, tensor
    # +5.8%. Log scale on (b) instead, and the ranges are stated in the title.
    axes[1].set_yscale("log")
    axes[1].set_ylabel(r"$k_{\mathrm{eff}}$ error vs Rayleigh  [%, log]")
    fig.suptitle("Conductivity-law error against solid fraction — "
                 "single cylinder, 1 mm periodic cell", fontsize=11)
    fig.tight_layout()
    out_png = args.out / "fsweep.png"
    fig.savefig(out_png, dpi=160, bbox_inches="tight")
    print(f"\nwrote {out_csv}\n      {out_png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

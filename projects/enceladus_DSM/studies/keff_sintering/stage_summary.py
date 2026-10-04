#!/usr/bin/env python3
"""Seed means and scatter per (porosity, temperature) for the k_eff campaign.

    venv_enceladus/bin/python studies/keff_sintering/stage_summary.py <campaign dir> \\
        [--out <dir>] [--family keff_LR40]

One row per condition, over the production packings (seed numbers 1601-2005,
L/R 40; gated and convergence runs are left out):

  n          packings
  k_0        k_iso at the opening sample (first t >= 1 s)
  k_11       k_iso at 11 tau_sub (1 d at -20 C), interpolated
  k_30       k_iso at the last sample
  rise       k_30/k_11 - 1, the sintering-driven rise   (mean, sd, min..max)
  ssa_ratio  SSA_30/SSA_11
  kxx/kyy    at the last sample, seed mean and sd
  Fxx/Fyy    contact-fabric ratio of the same packings at t = 0
             (inputs/packings/keff_LR40/packings_summary.csv)

Also the per-run table. Writes stage_summary.csv and stage_summary_runs.csv
into --out (default: <campaign dir>/compare/summary/).
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PROJ = HERE.parents[1]
sys.path.insert(0, str(PROJ / "postprocess"))
from plot_keff import load, read_tau_sub  # noqa: E402

PAT = re.compile(r"phi([\d.]+)_Rave50um_LR40_seed(\d+)_L2mm_eps1000nm_perxy_T(-?\d+)__")


def fabric():
    f = PROJ / "inputs/packings/keff_LR40/packings_summary.csv"
    out = {}
    if f.is_file():
        for r in csv.DictReader(open(f)):
            out[int(r["base_seed"])] = (float(r["F_ratio"]), float(r["z_band"]))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root", type=Path)
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args()
    out = a.out or a.root / "compare" / "summary"
    F = fabric()
    runs = []
    for kf in sorted(a.root.glob("packing_*/k_eff.csv")):
        d = kf.parent
        m = PAT.search(d.name)
        if not m or not (1601 <= int(m.group(2)) <= 2005):
            continue
        r = load(d)
        tau = read_tau_sub(d)
        if r is None or not tau:
            continue
        t, k = r["t"], r["kiso"]
        i0 = int(np.argmax(t >= 1.0))
        t11 = 11.05 * tau
        seed = int(m.group(2))
        runs.append(dict(
            phi=float(m.group(1)), T=int(m.group(3)), seed=seed,
            k_0=k[i0], k_11=float(np.interp(t11, t, k)), k_30=k[-1],
            t_end_d=float(np.loadtxt(d / "SSA_evo.dat")[-1, 2]) / 86400,   # last STEP, not last sample
            ssa_ratio=r["ssa"][-1] / float(np.interp(t11, t, r["ssa"])),
            kxx_kyy=r["kxx"][-1] / r["kyy"][-1],
            kxy_kxx=r["kxy"][-1] / r["kxx"][-1],
            F_ratio=F.get(seed, (np.nan, np.nan))[0], z_band=F.get(seed, (np.nan, np.nan))[1]))
        runs[-1]["rise_pct"] = 100 * (runs[-1]["k_30"] / runs[-1]["k_11"] - 1)

    out.mkdir(parents=True, exist_ok=True)
    with open(out / "stage_summary_runs.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(runs[0])); w.writeheader(); w.writerows(runs)

    G = defaultdict(list)
    for r in runs:
        G[(r["phi"], r["T"])].append(r)
    rows = []
    sd = lambda v: float(np.std(v, ddof=1)) if len(v) > 1 else float("nan")
    print(f"{'phi':>6} {'T':>4} {'n':>2} {'k_0':>6} {'k_11':>13} {'k_30':>13} "
          f"{'rise %':>16} {'kxx/kyy':>12} {'Fxx/Fyy':>8} {'z_band':>6}")
    for (phi, T), g in sorted(G.items()):
        v = lambda k: np.array([x[k] for x in g])
        row = dict(phi=phi, T=T, n=len(g), k_0=v("k_0").mean(),
                   k_11=v("k_11").mean(), k_11_sd=sd(v("k_11")),
                   k_30=v("k_30").mean(), k_30_sd=sd(v("k_30")),
                   rise_pct=v("rise_pct").mean(), rise_sd=sd(v("rise_pct")),
                   rise_min=v("rise_pct").min(), rise_max=v("rise_pct").max(),
                   ssa_ratio=v("ssa_ratio").mean(),
                   kxx_kyy=v("kxx_kyy").mean(), kxx_kyy_sd=sd(v("kxx_kyy")),
                   F_ratio=np.nanmean(v("F_ratio")), z_band=np.nanmean(v("z_band")))
        rows.append(row)
        print(f"{phi:6.3f} {T:4d} {len(g):2d} {row['k_0']:6.3f} "
              f"{row['k_11']:6.3f}±{row['k_11_sd']:5.3f} {row['k_30']:6.3f}±{row['k_30_sd']:5.3f} "
              f"{row['rise_pct']:5.1f}±{row['rise_sd']:3.1f} ({row['rise_min']:4.1f}-{row['rise_max']:4.1f}) "
              f"{row['kxx_kyy']:5.3f}±{row['kxx_kyy_sd']:5.3f} {row['F_ratio']:8.3f} {row['z_band']:6.2f}")
    with open(out / "stage_summary.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    short = [r for r in runs if r["t_end_d"] < 29.99]
    for r in short:
        print(f"  NOTE phi {r['phi']} seed {r['seed']} T {r['T']}: last step at {r['t_end_d']:.2f} d "
              f"(final step rejected by the CFL limiter and not retried)")
    print(f"wrote {out / 'stage_summary.csv'} and stage_summary_runs.csv")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Read the DoF/core scaling test (scripts/HPC/submit_scaling_test.sh).

    venv_enceladus/bin/python studies/keff_sintering/scaling/analyze_scaling.py \\
        <downloaded batch dirs ...> [--out <dir>]

For every job: ranks, the phase-field time per step (-log_view's TSStep,
max over ranks), and the k_eff time per sample (k_eff.csv wall_s, median).
Then, per DoF/core target, the COST per step and per sample in core-seconds
(ranks x seconds) -- the quantity the bill is made of -- and the projected
core-hours of one production run at -20 C and -5 C from their measured step
counts (3a: 368 / 1216) and the SSA-trigger sample counts (179 / 308).

Writes scaling.csv and scaling.png next to this file (or --out).
Fetch per job: SSA_evo.dat, k_eff.csv, the .o log (holds -log_view), *.opts.
"""
from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
RATE = 0.012                                   # $/core-hour, Resnick tier 1
# steps, k_eff samples per 30-day run: 3a steps, SSA-trigger cadence (predict_cadence.py)
RUNS = {"-20 C": (368, 179), "-5 C": (1216, 308)}
C = ["#0072B2", "#D55E00", "#009E73"]


def one(d: Path):
    logs = sorted(d.glob("*.o[0-9]*"))
    if not logs:
        return None
    txt = logs[-1].read_text(errors="replace")
    m = re.search(r"(\d+) MPI ranks", txt)
    ranks = int(m.group(1)) if m else None
    m = re.search(r"^TSStep\s+(\d+)\s+[\d.]+\s+([\d.eE+-]+)", txt, re.M)
    if not (ranks and m):
        return None
    nstep, t_step = int(m.group(1)), float(m.group(2))
    kf = d / "k_eff.csv"
    ks = np.atleast_1d(np.genfromtxt(kf, delimiter=",", names=True)) if kf.is_file() else None
    s_samp = float(np.median(ks["wall_s"])) if ks is not None and len(ks) else np.nan
    its = float(np.median(ks["ksp_its"])) if ks is not None and len(ks) else np.nan
    return dict(dir=d.name, ranks=ranks, steps=nstep, s_step=t_step / max(nstep, 1),
                s_sample=s_samp, keff_its=its, n_samples=0 if ks is None else len(ks))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dirs", nargs="+", type=Path)
    ap.add_argument("--out", type=Path, default=HERE)
    a = ap.parse_args()
    rows = []
    for root in a.dirs:
        for so in root.rglob("SSA_evo.dat"):
            r = one(so.parent)
            if r:
                rows.append(r)
    if not rows:
        raise SystemExit("no scaling jobs found (need SSA_evo.dat + .o with -log_view)")

    by = defaultdict(list)
    for r in rows:
        by[r["ranks"]].append(r)
    out = []
    print(f"{'ranks':>5} {'DoF/core':>8} {'n':>2} {'s/step':>7} {'s/sample':>8} "
          f"{'its':>4} {'core-s/step':>11} {'core-s/sample':>13} "
          + " ".join(f"{'$/run ' + k:>12}" for k in RUNS))
    for ranks in sorted(by):
        g = by[ranks]
        ss = np.mean([r["s_step"] for r in g]); sk = np.nanmean([r["s_sample"] for r in g])
        its = np.nanmean([r["keff_its"] for r in g])
        row = dict(ranks=ranks, dof_per_core=round(24_009_723 / ranks), n=len(g),
                   s_step=ss, s_step_sd=np.std([r["s_step"] for r in g]),
                   s_sample=sk, s_sample_sd=np.nanstd([r["s_sample"] for r in g]),
                   keff_its=its, core_s_step=ss * ranks, core_s_sample=sk * ranks)
        for k, (ns, nk) in RUNS.items():
            row[f"usd_{k}"] = (ns * ss + nk * sk) * ranks / 3600 * RATE
        out.append(row)
        print(f"{ranks:5d} {row['dof_per_core']:8d} {len(g):2d} {ss:7.1f} {sk:8.1f} {its:4.0f} "
              f"{row['core_s_step']:11.0f} {row['core_s_sample']:13.0f} "
              + " ".join(f"{row['usd_' + k]:12.2f}" for k in RUNS))

    a.out.mkdir(parents=True, exist_ok=True)
    with open(a.out / "scaling.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0]))
        w.writeheader(); w.writerows(out)

    X = np.array([o["dof_per_core"] for o in out]) / 1e3      # kDoF per core
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))
    ax = axes[0]
    ax.errorbar(X, [o["s_step"] for o in out], yerr=[o["s_step_sd"] for o in out], fmt="o-",
                color=C[0], label="phase-field step")
    ax.errorbar(X, [o["s_sample"] for o in out], yerr=[o["s_sample_sd"] for o in out], fmt="s-",
                color=C[1], label="k_eff sample")
    ax.set(ylabel="wall seconds", title="(a) time per step / per sample")
    ax.legend(frameon=False)
    ax = axes[1]
    ax.plot(X, [o["core_s_step"] for o in out], "o-", color=C[0], label="phase-field step")
    ax.plot(X, [o["core_s_sample"] for o in out], "s-", color=C[1], label="k_eff sample")
    ax.set(ylabel="core-seconds", title="(b) cost per step / per sample")
    ax.set_yscale("log")
    ax.legend(frameon=False)
    ax = axes[2]
    for (k, _), col in zip(RUNS.items(), C):
        ax.plot(X, [o[f"usd_{k}"] for o in out], "o-", color=col, label=k)
    ax.set(ylabel="$ per 30-day run", title="(c) projected cost per run")
    ax.set_yscale("log")
    ax.legend(frameon=False)
    for ax in axes:
        ax.set_xscale("log")
        ax.set_xticks(X)
        ax.set_xticklabels([f"{x:.0f}k\n({o['ranks']})" for x, o in zip(X, out)], fontsize=9)
        ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.set_xlabel("DoF per core  (MPI ranks)")
        ax.grid(True, alpha=0.25)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    fig.tight_layout()
    fig.savefig(a.out / "scaling.png", dpi=150, bbox_inches="tight")
    print(f"  wrote {a.out / 'scaling.csv'} and scaling.png")


if __name__ == "__main__":
    main()

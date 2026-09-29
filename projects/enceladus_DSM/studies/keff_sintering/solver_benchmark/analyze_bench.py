#!/usr/bin/env python3
"""Read the k_eff solver benchmark: speed and accuracy per configuration.

    venv_enceladus/bin/python studies/keff_sintering/solver_benchmark/analyze_bench.py <run_dir> [<bench batch dir> ...]

<run_dir> holds the bench_<cfg>.csv files the replays wrote (the 3a phi 0.475
run). The optional batch dirs hold the replay jobs' SLURM .o logs, from which
the -log_view times of the main events are pulled (PCSetUp = building the
multigrid hierarchy, KSPSolve = the whole solve, MatMult = one operator apply).

Accuracy is judged against 'base' (CG + GAMG, rtol 1e-10, the production
setting): max |k_ij - k_ij,base| / k_iso over the replayed samples. The
campaign's seed-to-seed scatter is ~3-4%, so anything below ~1e-5 is exact for
every purpose here; a configuration is only adoptable if it clears that.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np

EVENTS = ("KSPSolve", "PCSetUp", "PCApply", "MatMult")


def log_times(batch_dirs, cfg):
    """Max-over-ranks time [s] per event from the -log_view of config `cfg`."""
    for b in batch_dirs:
        for o in Path(b).rglob(f"*__{cfg}.o*"):
            txt = o.read_text(errors="replace")
            out = {}
            for ev in EVENTS:
                m = re.search(rf"^{ev}\s+\d+\s+[\d.]+\s+([\d.eE+-]+)", txt, re.M)
                if m:
                    out[ev] = float(m.group(1))
            m = re.search(r"^Time \(sec\):\s+([\d.eE+-]+)", txt, re.M)
            if m:
                out["total"] = float(m.group(1))
            m = re.search(r"(\d+) MPI ranks", txt)
            if m:
                out["ranks"] = int(m.group(1))
            return out
    return {}


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    run = Path(sys.argv[1])
    batches = sys.argv[2:]
    files = sorted(run.glob("bench_*.csv"))
    if not files:
        sys.exit(f"no bench_*.csv in {run}")
    data = {f.stem[len("bench_"):]: np.atleast_1d(np.genfromtxt(f, delimiter=",", names=True))
            for f in files}
    base = data.get("base")
    print(f"{'config':10s} {'n':>2} {'its/sample':>10} {'s/sample':>9} {'s/it':>6} "
          f"{'max|dk|/k_iso':>13}   log_view (max over ranks): "
          + "  ".join(EVENTS))
    for cfg, k in data.items():
        its = np.median(k["ksp_its"]); w = np.median(k["wall_s"])
        if base is not None and cfg != "base":
            n = min(len(k), len(base))
            dk = max(np.max(np.abs(k[c][:n] - base[c][:n]) / base["k_iso"][:n])
                     for c in ("k_00", "k_01", "k_11"))
            acc = f"{dk:13.2e}"
        else:
            acc = f"{'(reference)':>13}"
        lt = log_times(batches, cfg)
        ev = "  ".join(f"{lt.get(e, float('nan')):8.1f}" for e in EVENTS)
        rk = f" [{lt['ranks']} ranks]" if "ranks" in lt else ""
        print(f"{cfg:10s} {len(k):2d} {its:10.0f} {w:9.1f} {w / its:6.3f} {acc}   {ev}{rk}")


if __name__ == "__main__":
    main()

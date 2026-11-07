"""How well each run of a batch ran: one row per run, from the files the
solver wrote, nothing fitted.

    completed   the job log says "Workflow completed successfully" and the
                last SSA_evo.dat time reaches -t_final (to 0.1 %)
    wall, ranks the job log's cost summary ("Wall time", "MPI ranks")
    steps       accepted steps (solver_evo.dat rows after step 0)
    newton      mean Newton iterations per step; max in brackets
    krylov      median Krylov iterations per Newton iteration
    rejected    total rejected steps (solver_evo.dat, last column, max)
    rollbacks   band rollbacks ("Phase field out of bounds") and interface-CFL
                rollbacks ("Interface-CFL violated ... rolling back") in the log
    at_dtmax    share of steps with dt >= 0.99 dtmax (-dtmax from the opts)
    phi range   min and max of phi_ice over every BOUNDS line in the log
    Pi_max      largest K&P Eq. 45 number after the start-up window
                (resolution.dat; the window is the solver's own, from its
                RESOLUTION SUMMARY line); Pi_gt_0.5: steps above the design limit
    res_face    largest phi on a vapor-reservoir face (resolution.dat)
    mass_drift  (total_mass_end - total_mass_0) / total_mass_0 from SSA_evo.dat;
                for an open pore this includes the vapor that left, so it is
                only a closure check when the run is closed. Use water_budget
                for the closure of an open run.

Usage:
    python3 postprocess/run_health.py --batch <batch folder> [--md out.md] [--csv out.csv]
    python3 postprocess/run_health.py --run <run folder> [...]

The markdown table lists runs by folder name with --label to shorten them:
--label 'channel_2D_' removes that prefix.
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

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pplib  # noqa: E402

_BOUNDS = re.compile(r"BOUNDS: phi_ice \[\s*([-+0-9.eE]+),\s*([-+0-9.eE]+)\]")
_RESSUM = re.compile(r"RESOLUTION SUMMARY \(after the ([0-9.eE+-]+) s start-up window")


def job_log(run):
    logs = sorted(f for f in glob.glob(os.path.join(run, "*.o*"))
                  if f.rsplit(".o", 1)[-1].isdigit())
    return logs[-1] if logs else None


def run_dirs(batch):
    return sorted(d for d in glob.glob(os.path.join(batch, "*"))
                  if os.path.isfile(os.path.join(d, "igasol.dat")))


def health(run):
    r = dict(run=os.path.basename(run.rstrip("/")))
    opts = pplib.read_opts(run)
    t_final = pplib.opt_float(opts, "-t_final", float("nan"))
    dtmax = pplib.opt_float(opts, "-dtmax", float("nan"))

    jl = job_log(run)
    text = Path(jl).read_text(errors="replace") if jl else ""
    r["job_ok"] = "Workflow completed successfully" in text
    m = re.search(r"Wall time:\s*(.+)", text)
    r["wall"] = m.group(1).strip() if m else ""
    m = re.search(r"MPI ranks:\s*(\d+)", text)
    r["ranks"] = int(m.group(1)) if m else None

    ssa = pplib.load_ssa(run)
    t_end = float(ssa[:, pplib.TIME].max()) if ssa is not None else float("nan")
    r["t_end_d"] = t_end / 86400
    r["completed"] = bool(r["job_ok"] and np.isfinite(t_final) and t_end >= 0.999 * t_final)
    if ssa is not None and ssa.shape[1] > 7:
        m0, m1 = ssa[0, 7], ssa[-1, 7]
        r["mass_drift"] = (m1 - m0) / m0 if m0 else float("nan")
    else:
        r["mass_drift"] = float("nan")

    se = os.path.join(run, "solver_evo.dat")
    if os.path.isfile(se):
        a = np.genfromtxt(se, comments="#", invalid_raise=False, ndmin=2)
        a = a[a[:, 0] > 0]
        r["steps"] = int(a.shape[0])
        r["newton_mean"] = float(a[:, 3].mean())
        r["newton_max"] = int(a[:, 3].max())
        r["krylov_median"] = float(np.median(a[:, 5]))
        r["rejected"] = int(a[:, 6].max())
        r["at_dtmax"] = float(np.mean(a[:, 2] >= 0.99 * dtmax)) if np.isfinite(dtmax) else float("nan")

    outp = os.path.join(run, "outp.txt")
    log = Path(outp).read_text(errors="replace") if os.path.isfile(outp) else text
    r["rb_band"] = log.count("Phase field out of bounds at step")
    r["rb_cfl"] = len(re.findall(r"Interface-CFL violated at step \d+: .*rolling back", log))
    lo = [float(x) for x, _ in _BOUNDS.findall(log)]
    hi = [float(y) for _, y in _BOUNDS.findall(log)]
    r["phi_min"] = min(lo) if lo else float("nan")
    r["phi_max"] = max(hi) if hi else float("nan")

    m = _RESSUM.search(log)
    t_skip = float(m.group(1)) if m else 0.0
    rd = os.path.join(run, "resolution.dat")
    if os.path.isfile(rd):
        res = np.genfromtxt(rd, comments="#", invalid_raise=False, ndmin=2)
        after = res[res[:, 0] > t_skip]
        r["Pi_max"] = float(after[:, 4].max()) if after.size else float("nan")
        r["Pi_gt_0.5"] = int((after[:, 4] > 0.5).sum()) if after.size else 0
        r["res_face"] = float(res[:, 5].max())
    return r


def fmt_row(r, label):
    name = r["run"]
    for p in label:
        name = name.replace(p, "")
    ok = "yes" if r["completed"] else f"**no** (t = {r['t_end_d']:.2f} d)"
    return (f"| {name} | {ok} | {r['wall']} / {r['ranks']} | {r.get('steps', '')} | "
            f"{r.get('newton_mean', float('nan')):.3f} ({r.get('newton_max', '')}) | "
            f"{r.get('krylov_median', float('nan')):.0f} | {r.get('rejected', '')} | "
            f"{r['rb_band']} / {r['rb_cfl']} | {100 * r.get('at_dtmax', float('nan')):.1f} % | "
            f"[{r['phi_min']:.1e}, {r['phi_max']:.4f}] | {r.get('Pi_max', float('nan')):.3f} "
            f"({r.get('Pi_gt_0.5', '')}) | {r.get('res_face', float('nan')):.1e} |")


HEADER = ("| Run | Completed | Wall time / ranks | Steps | Newton per step (max) | "
          "Krylov per Newton | Rejected steps | Rollbacks band / CFL | Steps at dtmax | "
          "phi_ice range | Pi max (steps > 0.5) | phi on reservoir face |\n"
          "|---|---|---|---|---|---|---|---|---|---|---|---|")


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--batch")
    g.add_argument("--run", action="append")
    p.add_argument("--md", help="write the markdown table here")
    p.add_argument("--csv", help="write the rows here")
    p.add_argument("--label", action="append", default=[], help="text to remove from run names")
    a = p.parse_args(argv)
    runs = run_dirs(a.batch) if a.batch else a.run
    if not runs:
        raise SystemExit("no run folders (with igasol.dat) found")
    rows = [health(r) for r in runs]
    table = HEADER + "\n" + "\n".join(fmt_row(r, a.label) for r in rows)
    print(table)
    if a.md:
        with open(a.md, "w") as fh:
            fh.write(table + "\n")
    if a.csv:
        keys = list(dict.fromkeys(k for r in rows for k in r))
        with open(a.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=keys)
            w.writeheader()
            w.writerows(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

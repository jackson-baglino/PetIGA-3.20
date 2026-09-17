#!/usr/bin/env python3
"""Which limiter is actually setting dt, and what will the run cost?

Run this against a LIVE or finished run directory. It reads SSA_evo.dat,
which the solver flushes every step, so it works mid-run:

    venv_enceladus/bin/python postprocess/diagnose_dt.py <run_dir>

WHY THIS MATTERS. Three different things can be holding dt down, and they
have completely different fixes:

  -dtmax           a hard ceiling from the opts file. If dt sits exactly on
                   it, the run is paying for a number somebody chose, and
                   raising it is free speed -- IF the CFL limiter agrees.
  -dtCFL           the adaptive interface-CFL limiter: no point may change
                   phi by more than -dtCFL_dphimax (0.2) in a step. If this
                   binds, dt is tracking real interface motion and raising
                   dtmax buys nothing.
  NRmin/NRmax      the Newton-iteration heuristic. Grows dt when Newton
                   converges in fewer than NRmin iterations. Rarely binds
                   once the other two are in play.

The output says which, and projects the remaining cost at the current rate.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np


def read_ssa(run: Path):
    f = run / "SSA_evo.dat"
    if not f.is_file():
        raise SystemExit(f"no SSA_evo.dat in {run}")
    a = np.loadtxt(f)
    if a.ndim == 1:
        a = a[None, :]
    # columns: ssa/eps, tot_ice, t, step, dt, tot_air, tot_rhov, tot_mass
    return {"t": a[:, 2], "step": a[:, 3].astype(int), "dt": a[:, 4],
            "ssa": a[:, 0], "ice": a[:, 1]}


def opt_from_run(run: Path, key: str):
    """Read a flag out of whatever .opts the run staged next to its output."""
    for p in list(run.glob("*.opts")) + list(run.glob("**/*.opts")):
        for line in p.read_text().splitlines():
            line = line.split("#", 1)[0].strip()
            m = re.match(rf"^-{re.escape(key)}\s+(\S+)$", line)
            if m:
                try:
                    return float(m.group(1))
                except ValueError:
                    return None
    return None


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir", type=Path)
    ap.add_argument("--dtmax", type=float, default=None,
                    help="override; otherwise read from the staged opts")
    ap.add_argument("--t-final", dest="t_final", type=float, default=None)
    ap.add_argument("--hours", type=float, default=None,
                    help="wall-clock hours elapsed so far, for the projection")
    a = ap.parse_args()

    d = read_ssa(a.run_dir)
    dtmax = a.dtmax if a.dtmax is not None else opt_from_run(a.run_dir, "dtmax")
    t_final = a.t_final if a.t_final is not None else opt_from_run(a.run_dir, "t_final")
    dt, t = d["dt"], d["t"]
    n = len(dt)

    print(f"run           {a.run_dir}")
    print(f"steps taken   {n}")
    print(f"sim time      {t[-1]:.4e} s"
          + (f"  of {t_final:.4e} s  ({t[-1]/t_final:.1%} done)" if t_final else ""))
    print(f"dt            min {dt.min():.4e}   median {np.median(dt):.4e}   "
          f"max {dt.max():.4e}")

    print("\n--- what is holding dt down? ---")
    if dtmax:
        at_cap = np.isclose(dt, dtmax, rtol=1e-6)
        frac = at_cap.mean()
        print(f"  -dtmax = {dtmax:.4e} s")
        print(f"  steps sitting exactly on it: {at_cap.sum()}/{n} ({frac:.1%})")
        # the recent history is what matters for the projection
        tail = at_cap[-min(200, n):]
        print(f"  of the last {len(tail)} steps: {tail.mean():.1%} at the cap")
        if tail.mean() > 0.8:
            print("\n  >> -dtmax IS THE BINDING CONSTRAINT.")
            print("     dt is not tracking the physics, it is sitting on a number")
            print("     from the opts file. The interface-CFL limiter is not")
            print("     asking for anything smaller, so raising -dtmax would")
            print("     speed the run up roughly in proportion -- up to the point")
            print("     where CFL starts to bind instead.")
        elif tail.mean() < 0.2:
            print("\n  >> -dtCFL IS THE BINDING CONSTRAINT (or Newton is).")
            print("     dt is below the cap almost every step, so it is tracking")
            print("     real interface motion. Raising -dtmax buys NOTHING.")
            print("     The levers are -dtCFL_dphimax (accuracy/cost trade) or a")
            print("     coarser eps (fewer steps AND cheaper steps).")
        else:
            print("\n  >> MIXED: the cap binds during quiet phases, CFL during")
            print("     events. Raising -dtmax gives a partial speedup.")
    else:
        print("  could not find -dtmax in the staged opts; pass --dtmax")

    if t_final and n > 1:
        remaining = t_final - t[-1]
        rate = (t[-1] - t[0]) / max(n - 1, 1)          # sim seconds per step
        steps_left = remaining / max(rate, 1e-30)
        print(f"\n--- projection at the current average dt ({rate:.4e} s/step) ---")
        print(f"  steps still to run   ~{steps_left:,.0f}")
        print(f"  total steps          ~{n + steps_left:,.0f}")
        if a.hours:
            per_step = a.hours * 3600.0 / n
            print(f"  wall per step        {per_step:.2f} s  (from --hours {a.hours})")
            print(f"  time remaining       ~{steps_left*per_step/3600:.1f} h")
            print(f"  total wall           ~{(n+steps_left)*per_step/3600:.1f} h")
        if dtmax:
            print(f"\n  if dt could sit at 10x -dtmax throughout:"
                  f"  ~{(t_final-t[0])/(10*dtmax):,.0f} steps total")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

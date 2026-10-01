#!/usr/bin/env python3
"""Score a dtmax ladder: does a larger timestep change the answer, and what does it cost?

Every rung is the same run (channel, theta = 60, alpha_c = 1e-3, sigma_inf = 0,
30 days) at a different -dtmax. The comparison uses the ice area the solver
writes EVERY step to SSA_evo.dat, not the snapshots, so rungs with very few
snapshots are still compared along the whole trajectory.

Reference = the rung with the smallest dtmax present. For each other rung:

  err_growth   max over its own step times of |ice(t) - ice_ref(t)|, as a
               fraction of the reference's total ice change. This is the error
               in the quantity the study measures (growth, i.e. velocity).
  d_final      the same at the final time only.

Solver cost comes from solver_evo.dat (Newton and Krylov iterations per step,
step rejections).

Usage:
    venv_lunar/bin/python3 studies/contact_angle/velocity_plan_2026-10-01/analyze_dtmax_ladder.py <batch_dir>
"""
import glob
import os
import re
import sys

import numpy as np


def load(run):
    ssa = np.loadtxt(os.path.join(run, "SSA_evo.dat"), ndmin=2)
    t, ice, dt = ssa[:, 2], ssa[:, 1], ssa[:, 4]
    sol = None
    p = os.path.join(run, "solver_evo.dat")
    if os.path.isfile(p):
        sol = np.loadtxt(p, ndmin=2)        # step t dt newton krylov k/n rejections
    return t, ice, dt, sol


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    runs = {}
    for d in sorted(glob.glob(os.path.join(sys.argv[1], "*_dt*tau"))):
        m = re.search(r"_dt([0-9.]+)tau$", d)
        if m and os.path.isfile(os.path.join(d, "SSA_evo.dat")):
            runs[float(m.group(1))] = load(d)
    if len(runs) < 2:
        sys.exit("need at least two finished rungs in %s" % sys.argv[1])

    ref_k = min(runs)
    tr, ir, _, _ = runs[ref_k]
    growth = ir[-1] - ir[0]
    print(f"reference: {ref_k:g} x tau_sub, {len(tr)-1} steps, "
          f"ice change {100*growth/ir[0]:+.4f} % over {tr[-1]/86400:.1f} d\n")
    print(f"{'dtmax/tau':>9} {'steps':>6} {'dt_max[s]':>10} {'newton':>7} {'kry/newt':>8} "
          f"{'reject':>6} {'ice chg %':>10} {'err_growth':>10} {'d_final':>9}")
    for k in sorted(runs):
        t, ice, dt, sol = runs[k]
        ref = np.interp(t, tr, ir)
        err = np.max(np.abs(ice - ref)) / abs(growth)
        tend = min(t[-1], tr[-1])
        dfin = (np.interp(tend, t, ice) - np.interp(tend, tr, ir)) / abs(growth)
        if sol is not None and len(sol):
            newt, kpn, rej = int(sol[:, 3].sum()), sol[:, 4].sum() / max(sol[:, 3].sum(), 1), int(sol[:, 6].sum())
            cost = f"{newt:7d} {kpn:8.1f} {rej:6d}"
        else:
            cost = f"{'-':>7} {'-':>8} {'-':>6}"
        print(f"{k:9g} {len(t)-1:6d} {dt.max():10.3e} {cost} "
              f"{100*(ice[-1]-ice[0])/ice[0]:+10.4f} {err:10.2e} {dfin:+9.2e}")
    print("\nerr_growth and d_final are fractions of the reference's total ice change.\n"
          "SSA_evo.dat prints 7 significant figures, so differences below ~1e-6 of the\n"
          "ice area itself are at the precision floor.")


if __name__ == "__main__":
    main()

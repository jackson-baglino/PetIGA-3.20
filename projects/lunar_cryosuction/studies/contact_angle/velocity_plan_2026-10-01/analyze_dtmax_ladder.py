#!/usr/bin/env python3
"""Score a dtmax ladder: does a larger timestep change the answer, and what does it cost?

Every rung is the same run at a different -dtmax. The comparison uses a
quantity the solver writes EVERY step to SSA_evo.dat, not the snapshots, so
rungs with very few snapshots are still compared along the whole trajectory.

Which quantity depends on the case, and is chosen automatically:

  ice area          reservoir runs (channel ladder): the ice grows or shrinks,
                    and that growth is what the velocity study measures.
  interface length  sealed runs (sintering / ripening pairs): total ice is
                    conserved, so the trajectory of the interface length is the
                    signal. It falls as the small grain is consumed.

Reference = the rung with the smallest dtmax present. For each other rung:

  err_traj   max over its own step times of |q(t) - q_ref(t)|, as a fraction
             of the reference's total change in q.
  d_final    the same at the final time only.
  t_half     time at which q has covered half of the REFERENCE's total change:
             a robust event time. For the grain pairs it tracks when the small
             grain disappears; the 2026-07-10 stress test found that 32 % late
             at a 12x larger dtmax even though the run stayed clean.

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
    t, ice, dt, interf = ssa[:, 2], ssa[:, 1], ssa[:, 4], ssa[:, 0]
    sol = None
    p = os.path.join(run, "solver_evo.dat")
    if os.path.isfile(p):
        sol = np.loadtxt(p, ndmin=2)        # step t dt newton krylov k/n rejections
    return t, ice, dt, sol, interf


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
    tr, ice_r, _, _, int_r = runs[ref_k]
    use_ice = abs(ice_r[-1] - ice_r[0]) / ice_r[0] > 1e-3
    col, label = (1, "ice area") if use_ice else (4, "interface length")
    qr = runs[ref_k][col]
    change = qr[-1] - qr[0]
    half = qr[0] + 0.5 * change

    def t_half(t, q):
        past = np.where((q - half) * np.sign(change) >= 0)[0]
        if len(past) == 0:
            return np.nan
        i = past[0]
        if i == 0:
            return t[0]
        return np.interp(half, sorted([q[i-1], q[i]]),
                         [t[i-1], t[i]] if q[i] > q[i-1] else [t[i], t[i-1]])

    th_ref = t_half(tr, qr)
    print(f"scored on {label}; reference {ref_k:g} x tau_sub, {len(tr)-1} steps, "
          f"change {100*change/qr[0]:+.3f} % over {tr[-1]/86400:.1f} d, "
          f"t_half = {th_ref/86400:.2f} d\n")
    print(f"{'dtmax/tau':>9} {'steps':>6} {'dt_max[s]':>10} {'newton':>7} {'kry/newt':>8} "
          f"{'reject':>6} {'err_traj':>9} {'d_final':>9} {'t_half[d]':>9} {'vs ref':>7}")
    for k in sorted(runs):
        t, dt, sol = runs[k][0], runs[k][2], runs[k][3]
        q = runs[k][col]
        err = np.max(np.abs(q - np.interp(t, tr, qr))) / abs(change)
        tend = min(t[-1], tr[-1])
        dfin = (np.interp(tend, t, q) - np.interp(tend, tr, qr)) / abs(change)
        th = t_half(t, q)
        if sol is not None and len(sol):
            newt, kpn, rej = int(sol[:, 3].sum()), sol[:, 4].sum() / max(sol[:, 3].sum(), 1), int(sol[:, 6].sum())
            cost = f"{newt:7d} {kpn:8.1f} {rej:6d}"
        else:
            cost = f"{'-':>7} {'-':>8} {'-':>6}"
        print(f"{k:9g} {len(t)-1:6d} {dt.max():10.3e} {cost} {err:9.2e} {dfin:+9.2e} "
              f"{th/86400:9.2f} {100*(th/th_ref-1):+6.1f}%")
    print("\nerr_traj and d_final are fractions of the reference's total change.\n"
          "SSA_evo.dat prints 7 significant figures, so differences below ~1e-6 of the\n"
          "quantity itself are at the precision floor.")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Predict where the SSA-triggered k_eff cadence samples, from a run's SSA_evo.dat.

    venv_enceladus/bin/python studies/keff_sintering/predict_cadence.py <run_dir> ...
        [--dlnssa 0.001] [--t0-tau 11.05] [--freq 5] [--max-gap-tau 20]

A line-by-line mirror of KeffDue in src/keff_sample.c (-keff_dlnssa), driven by
the per-step SSA every run already records. Used to choose the production
setting (2026-09-29) and to check a run's k_eff.csv against it: after a run
with the trigger on, the steps listed here must be exactly the steps in its
k_eff.csv (--check).
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np


def schedule(steps, t, ssa, tau, dlnssa, t0_tau, freq, max_gap_tau):
    have = False
    last = tl = None
    out = []
    for st, tt, s in zip(steps, t, ssa):
        if dlnssa <= 0:                     # trigger off: the plain step cadence
            if st == 0 or (freq > 0 and st % freq == 0):
                out.append(int(st))
            continue
        ln = np.log(s)
        if st == 0:
            due = True
        elif not have:                      # first step after a restart
            have, last, tl = True, ln, tt
            continue
        elif tau <= 0 or tt < t0_tau * tau:
            due = (freq > 0 and st % freq == 0)
        else:
            due = (last - ln >= dlnssa) or (max_gap_tau > 0 and tt - tl >= max_gap_tau * tau)
        if due or st == 0:
            have, last, tl = True, ln, tt
        if due:
            out.append(int(st))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+", type=Path)
    ap.add_argument("--dlnssa", type=float, default=0.001)
    ap.add_argument("--t0-tau", type=float, default=86400.0 / 7822.3)
    ap.add_argument("--freq", type=int, default=5)
    ap.add_argument("--max-gap-tau", type=float, default=20.0)
    ap.add_argument("--check", action="store_true",
                    help="compare with the run's k_eff.csv steps")
    a = ap.parse_args()
    for d in a.runs:
        s = np.loadtxt(d / "SSA_evo.dat")
        m = re.search(r"tau_sub\s+([0-9.eE+-]+)", (d / "outp.txt").read_text(errors="replace"))
        tau = float(m.group(1)) if m else 0.0
        sch = schedule(s[:, 3].astype(int), s[:, 2], s[:, 0], tau,
                       a.dlnssa, a.t0_tau, a.freq, a.max_gap_tau)
        ts = s[np.isin(s[:, 3].astype(int), sch), 2]
        pre = int(np.sum(ts < a.t0_tau * tau))
        msg = f"{d.name[:70]:70s} {len(s):5d} steps -> {len(sch):4d} samples ({pre} before t0)"
        if a.check:
            k = np.atleast_1d(np.genfromtxt(d / "k_eff.csv", delimiter=",", names=True))
            got = [int(x) for x in k["step"]]
            msg += "  MATCH" if got == sch else f"  MISMATCH ({len(got)} in k_eff.csv)"
        print(msg)


if __name__ == "__main__":
    main()

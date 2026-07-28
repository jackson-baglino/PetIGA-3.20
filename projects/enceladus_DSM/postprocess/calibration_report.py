#!/usr/bin/env python3
"""
calibration_report.py — turn the Phase-4 sweep into an accept/reject table.

For each arm, compares its observables against the FINEST arm of the same test
and reports the relative deviation. The recommendation is the cheapest arm whose
every observable stays inside tolerance.

Observables
  ripen  final large/small grain ice split -> Ostwald transfer;
         ice-mass drift over the run (conservation)
  pack   SSA(t) endpoint and trajectory; ice-mass drift.
         k_eff must be added separately by running effective_thermal_cond over
         the snapshots -- pass --keff <dir> to fold in its k_eff.csv.
  molaro neck width, measured by postprocess/neck_width.py (run separately;
         this script only reports what it finds in neck_width.txt)

SSA_evo.dat columns (monitoring.c:280):
    sub_interf/eps   tot_ice   t   step   dt   tot_air   tot_rhov   tot_mass
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import sys
from pathlib import Path


def read_ssa(run: Path):
    f = run / "SSA_evo.dat"
    if not f.is_file():
        return None
    rows = []
    for line in f.read_text().splitlines():
        p = line.split()
        if len(p) < 8:
            continue
        try:
            rows.append([float(x) for x in p])
        except ValueError:
            continue
    if not rows:
        return None
    return {
        "ssa": [r[0] for r in rows],
        "ice": [r[1] for r in rows],
        "t": [r[2] for r in rows],
        "mass": [r[7] for r in rows],
    }


def summarize(run: Path):
    d = read_ssa(run)
    if d is None:
        return None
    ice0, ice1 = d["ice"][0], d["ice"][-1]
    m0, m1 = d["mass"][0], d["mass"][-1]
    return {
        "t_end": d["t"][-1],
        "n": len(d["t"]),
        "ssa_end": d["ssa"][-1],
        "ssa_drop": (d["ssa"][0] - d["ssa"][-1]) / d["ssa"][0] if d["ssa"][0] else float("nan"),
        "ice_drift": abs(ice1 - ice0) / ice0 if ice0 else float("nan"),
        "mass_drift": abs(m1 - m0) / abs(m0) if m0 else float("nan"),
    }


def keff_last(path: Path):
    """Last row of a k_eff.csv -> mean diagonal."""
    if not path.is_file():
        return None
    rows = list(csv.DictReader(path.open()))
    if not rows:
        return None
    r = rows[-1]
    try:
        return 0.5 * (float(r["k_00"]) + float(r["k_11"]))
    except (KeyError, ValueError):
        return None


ARM = re.compile(r"(calib_(?:ripen|pack))_(s\d{3})|_(xiv[0-9.e-]+)")


def classify(name: str):
    """(test, arm) from a run-directory name."""
    test = None
    for t in ("calib_ripen", "calib_pack", "molaro"):
        if t in name:
            test = t
            break
    m = re.search(r"_(s\d{3})(?:_|$)", name)
    if m:
        return test, m.group(1)
    m = re.search(r"_(xiv[0-9.e+-]+)", name)
    if m:
        return test, m.group(1)
    m = re.search(r"(eps(?:loose|mid|strict))", name)
    if m:
        return test, m.group(1)
    return test, None


# finest -> coarsest, so index 0 is the reference
ORDER = {
    "safety": ["s025", "s050", "s075", "s100"],
    "xiv": ["xiv1e-1", "xiv1e-2", "xiv1e-3", "xiv1e-4"],
    "molaro": ["epsstrict", "epsmid", "epsloose"],
}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("results", type=Path,
                    help="directory holding the calibration run folders")
    ap.add_argument("--tol", type=float, default=0.02,
                    help="relative tolerance on the observables")
    ap.add_argument("--mass-tol", type=float, default=0.01,
                    help="tolerance on ice-mass drift")
    ap.add_argument("--keff-root", type=Path, default=None,
                    help="effective_thermal_cond output root, for k_eff.csv")
    args = ap.parse_args(argv)

    runs = {}
    for d in sorted(args.results.rglob("SSA_evo.dat")):
        run = d.parent
        test, arm = classify(run.name)
        if test is None or arm is None:
            continue
        s = summarize(run)
        if s is None:
            continue
        if args.keff_root:
            s["keff"] = keff_last(args.keff_root / run.name / "k_eff.csv")
        runs.setdefault(test, {})[arm] = (run, s)

    if not runs:
        print(f"No calibration runs with SSA_evo.dat under {args.results}",
              file=sys.stderr)
        return 1

    verdict = {}
    for test, arms in sorted(runs.items()):
        order = (ORDER["molaro"] if test == "molaro"
                 else ORDER["safety"] if any(a.startswith("s") for a in arms)
                 else ORDER["xiv"])
        present = [a for a in order if a in arms]
        if not present:
            continue
        ref_arm = present[0]
        _, ref = arms[ref_arm]

        print(f"\n=== {test}   (reference = finest arm: {ref_arm}) ===")
        hdr = f"{'arm':>10} {'t_end[d]':>9} {'SSA_end':>10} {'dSSA':>8} " \
              f"{'ice drift':>10} {'mass drift':>11}"
        if any("keff" in s for _, s in arms.values()):
            hdr += f" {'k_eff':>9} {'dk':>8}"
        print(hdr)
        print("  " + "-" * (len(hdr) - 2))

        for a in present:
            _, s = arms[a]
            dssa = abs(s["ssa_end"] - ref["ssa_end"]) / abs(ref["ssa_end"]) \
                if ref["ssa_end"] else float("nan")
            line = (f"{a:>10} {s['t_end']/86400:>9.2f} {s['ssa_end']:>10.4e} "
                    f"{dssa*100:>7.2f}% {s['ice_drift']*100:>9.3f}% "
                    f"{s['mass_drift']*100:>10.3f}%")
            ok = (dssa <= args.tol and s["mass_drift"] <= args.mass_tol)
            if "keff" in s and s["keff"] is not None and ref.get("keff"):
                dk = abs(s["keff"] - ref["keff"]) / ref["keff"]
                line += f" {s['keff']:>9.4f} {dk*100:>7.2f}%"
                ok = ok and dk <= args.tol
            print(line + ("   OK" if ok else "   FAIL"))
            verdict.setdefault(test, []).append((a, ok))

    print("\n=== recommendation ===")
    print(f"  tolerance: {args.tol*100:g}% on observables, "
          f"{args.mass_tol*100:g}% on ice-mass drift")
    for test, res in sorted(verdict.items()):
        passing = [a for a, ok in res if ok]
        # ORDER runs finest -> coarsest, so the LAST passing arm is the cheapest
        best = passing[-1] if passing else None
        print(f"  {test:>12}: cheapest passing arm = {best or 'NONE — all arms differ'}")
    print("\n  Adopt the coarsest setting that passes in EVERY test, not just one.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""Merge the two halves of a restarted run into one directory.

A run stopped and resumed via -initial_cond leaves two directories whose step
numbers both start at 0 and whose time ranges OVERLAP -- leg 1 kept running for
a few steps past the snapshot the restart was taken from. Analysing them
separately is awkward and analysing them naively double-counts the overlap.

This produces one directory per seed that looks like an uninterrupted run:

    SSA_evo.dat     leg 1 truncated at leg 2's start time, then leg 2 with its
                    steps renumbered to continue. Strictly increasing in t.
    k_eff.csv       same treatment, same schema.
    vtkOut/         snapshots from both legs, leg 2's renumbered, no collisions.
    outp.txt        concatenated with a banner marking the join.
    <name>.opts     copied from leg 1; MERGE_INFO.json records both legs,
                    the join time, and what was dropped.

THE OVERLAP IS DROPPED, NOT AVERAGED. Leg 1's rows after the restart point are
from a trajectory that was abandoned; leg 2 recomputed that interval from the
snapshot. Keeping both would put two different states at the same t.

Idempotent and re-runnable: safe to call while the remaining legs are still
downloading, and safe to call again once they arrive (pass --force to rebuild).

    venv_enceladus/bin/python postprocess/merge_restart_legs.py <batch_dir>
"""
from __future__ import annotations

import argparse
import json
import re
import shutil
import sys
from pathlib import Path

import numpy as np

T_COL, STEP_COL = 2, 3          # SSA_evo.dat columns


def find_pairs(batch: Path):
    """(name, leg1_dir, [leg2_dirs...]) for every run in the batch."""
    pairs = []
    for l1 in sorted(batch.glob("*__*")):
        if not l1.is_dir() or not (l1 / "SSA_evo.dat").is_file():
            continue
        geom = l1.name.split("__")[0]
        legs = []
        gdir = batch / geom
        if gdir.is_dir():
            # every timestamped subdir that actually ran; sorted by name, which
            # is sorted by timestamp because the name leads with it
            legs = sorted(d for d in gdir.iterdir()
                          if d.is_dir() and (d / "SSA_evo.dat").is_file())
        pairs.append((geom, l1, legs))
    return pairs


def merge_table(rows1, rows2, t_col, join_t):
    """Leg 1 truncated below join_t, then leg 2. Returns the stacked array."""
    keep1 = rows1[rows1[:, t_col] < join_t]
    return np.vstack([keep1, rows2]) if len(rows2) else keep1


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("batch", type=Path)
    ap.add_argument("--out", type=Path, default=None,
                    help="output dir (default: <batch>/merged)")
    ap.add_argument("--force", action="store_true",
                    help="rebuild a merged run that already exists")
    ap.add_argument("--no-vtk", action="store_true",
                    help="skip copying vtkOut (much faster; tables only)")
    a = ap.parse_args()

    out_root = a.out or (a.batch / "merged")
    out_root.mkdir(parents=True, exist_ok=True)
    n_done = n_skip = 0

    for geom, l1, legs in find_pairs(a.batch):
        short = re.sub(r"^packing_2D_pilot_", "", geom)
        dest = out_root / geom
        if dest.is_dir() and not a.force:
            print(f"  EXISTS   {short}  (--force to rebuild)")
            n_skip += 1
            continue
        if not legs:
            print(f"  WAITING  {short}  — no resume leg downloaded yet")
            n_skip += 1
            continue

        ssa = [np.loadtxt(l1 / "SSA_evo.dat")]
        keff = [np.genfromtxt(l1 / "k_eff.csv", delimiter=",", names=True)]
        joins, step_off = [], int(ssa[0][:, STEP_COL].max())

        for leg in legs:
            s2 = np.loadtxt(leg / "SSA_evo.dat")
            if s2.ndim == 1:
                s2 = s2[None, :]
            join_t = float(s2[0, T_COL])
            # renumber this leg's steps to continue the merged record
            n_dropped = int((ssa[-1][:, T_COL] >= join_t).sum())
            prev_keep = ssa[-1][ssa[-1][:, T_COL] < join_t]
            base = int(prev_keep[:, STEP_COL].max()) if len(prev_keep) else 0
            s2 = s2.copy()
            s2[:, STEP_COL] += base
            ssa[-1] = prev_keep
            ssa.append(s2)

            k2 = np.genfromtxt(leg / "k_eff.csv", delimiter=",", names=True)
            k2 = np.atleast_1d(k2)
            kprev = keff[-1][keff[-1]["time"] < join_t]
            k2 = k2.copy()
            k2["step"] += base
            keff[-1] = kprev
            keff.append(k2)

            # counted BEFORE the truncation; counting after always gives 0,
            # which is what the first version did
            joins.append({"leg": leg.name, "join_time_s": join_t,
                          "rows_superseded": n_dropped})
            step_off = base

        ssa_all = np.vstack(ssa)
        keff_all = np.concatenate(keff)
        # strictly increasing in t is the whole point; assert rather than hope
        if not np.all(np.diff(ssa_all[:, T_COL]) > 0):
            bad = int(np.argmin(np.diff(ssa_all[:, T_COL]) > 0))
            print(f"  ERROR    {short}  — merged time is not increasing at row {bad}")
            n_skip += 1
            continue

        dest.mkdir(parents=True, exist_ok=True)
        np.savetxt(dest / "SSA_evo.dat", ssa_all,
                   fmt="%e %e %e %d %e %e %e %e")
        hdr = ",".join(keff_all.dtype.names)
        np.savetxt(dest / "k_eff.csv", keff_all.view(np.float64).reshape(len(keff_all), -1),
                   delimiter=",", header=hdr, comments="",
                   fmt=["%d", "%.12e"] + ["%.12e"] * 6 + ["%d", "%d", "%.3f"])

        for opt in l1.glob("*.opts"):
            shutil.copy2(opt, dest / opt.name)
        for sub in ("inputs",):
            if (l1 / sub).is_dir() and not (dest / sub).exists():
                shutil.copytree(l1 / sub, dest / sub)

        with (dest / "outp.txt").open("w") as fh:
            for lab, d in [("LEG 1", l1)] + [(f"LEG {i+2}", x) for i, x in enumerate(legs)]:
                if (d / "outp.txt").is_file():
                    fh.write(f"\n{'='*78}\n=== {lab}: {d.name}\n{'='*78}\n")
                    fh.write((d / "outp.txt").read_text())

        n_vtk = 0
        if not a.no_vtk:
            vdest = dest / "vtkOut"
            vdest.mkdir(exist_ok=True)
            cut = joins[0]["join_time_s"] if joins else np.inf
            # leg 1 keeps only snapshots at or before the join
            keep_steps = set(ssa[0][:, STEP_COL].astype(int).tolist())
            for f in sorted((l1 / "vtkOut").glob("solV_*.vts")):
                if int(f.stem.split("_")[1]) in keep_steps:
                    shutil.copy2(f, vdest / f.name); n_vtk += 1
            base = int(ssa[0][:, STEP_COL].max()) if len(ssa[0]) else 0
            for leg in legs:
                for f in sorted((leg / "vtkOut").glob("solV_*.vts")):
                    new = base + int(f.stem.split("_")[1])
                    shutil.copy2(f, vdest / f"solV_{new:05d}.vts"); n_vtk += 1
                base += int(np.loadtxt(leg / "SSA_evo.dat")[:, STEP_COL].max())

        (dest / "MERGE_INFO.json").write_text(json.dumps({
            "geometry": geom, "leg1": str(l1),
            "legs": [str(x) for x in legs], "joins": joins,
            "n_ssa_rows": int(len(ssa_all)), "n_keff_rows": int(len(keff_all)),
            "t_start_s": float(ssa_all[0, T_COL]), "t_end_s": float(ssa_all[-1, T_COL]),
            "n_vtk": n_vtk,
        }, indent=2))

        print(f"  MERGED   {short}")
        print(f"             {len(ssa_all)} steps, {len(keff_all)} k_eff samples, "
              f"t = {ssa_all[0, T_COL]:.3e} -> {ssa_all[-1, T_COL]:.3e} s"
              + (f", {n_vtk} vts" if n_vtk else ""))
        for j in joins:
            print(f"             join at t = {j['join_time_s']:.4e} s "
                  f"({j['rows_superseded']} superseded rows dropped)")
        n_done += 1

    print(f"\n{n_done} merged, {n_skip} skipped -> {out_root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

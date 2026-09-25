#!/usr/bin/env python3
"""Work out WHICH sol_*.dat are worth downloading, and emit an rsync list.

    # 1. pull the tables only (a few MB)
    rsync -av --include='*/' --include='k_eff.csv' --include='SSA_evo.dat' \
          --include='outp.txt' --include='igasol.dat' --include='*.opts' \
          --exclude='*' <host>:<remote>/ ./local/

    # 2. decide what else is needed
    venv_enceladus/bin/python postprocess/plan_download.py ./local --out files.txt

    # 3. pull exactly that
    rsync -av --files-from=files.txt <host>:<remote>/ ./local/

WHY NOT JUST PULL EVERYTHING. At L/R_ave = 64 the mesh is 4526^2, so a
snapshot is ~460 MB and a run writing one per step is ~170 GB. Four of those
is 680 GB. The redundancy is in the SAMPLING: k_eff is recorded ~54 times
while snapshots are written every step, so one snapshot per k_eff sample is
~55 files and ~2% of the bytes, with nothing analytical lost.

WHY ANY SNAPSHOTS AT ALL. The tables alone answer the REV question -- the
L/R_ave = 40 and 64 ensembles carry the same k_eff bias, so comparing them is
valid even though both are biased. They do NOT let k_eff be recomputed with
the sharp or tensorial coefficient, which needs the field. Since that
recomputation moved the headline by a factor of three, the fields are worth
having for every k_eff sample.

Matching is on TIME, not step: k_eff samples land between snapshots, so the
nearest snapshot is taken and the worst gap is reported.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pplib                                                # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("local", type=Path, help="dir holding the downloaded tables")
    ap.add_argument("--out", type=Path, default=Path("files.txt"))
    ap.add_argument("--extra", type=int, default=0,
                    help="additionally keep every Nth snapshot (0 = none)")
    ap.add_argument("--snapshot-mb", type=float, default=460.0,
                    help="approx size of one sol_*.dat, for the estimate")
    a = ap.parse_args()

    lines, total = [], 0
    for kf in sorted(a.local.glob("**/k_eff.csv")):
        d = kf.parent
        ssa = pplib.load_ssa(str(d))
        if ssa is None:
            print(f"  SKIP {d.name}: no SSA_evo.dat"); continue
        k = np.atleast_1d(np.genfromtxt(kf, delimiter=",", names=True))
        steps = ssa[:, 3].astype(int)
        want, worst = set(), 0
        for tk in np.atleast_1d(k["time"]):
            i = int(np.argmin(np.abs(ssa[:, 2] - tk)))
            want.add(int(steps[i]))
            worst = max(worst, abs(float(ssa[i, 2] - tk)))
        want |= {int(steps[0]), int(steps[-1])}
        # Step 1 is the first SOLVED state, and what movies start from: the IC
        # (step 0) has a uniform vapour field, so opening a movie on it flashes
        # from flat to structured. OutputMonitor has written it unconditionally
        # since 2026-08-13; a run older than that simply won't have the file.
        if 1 in set(steps.tolist()):
            want.add(1)
        if a.extra > 1:
            want |= {int(s) for s in steps[::a.extra]}
        rel = d.relative_to(a.local)
        for s in sorted(want):
            lines.append(f"{rel}/sol_{s:05d}.dat")
        total += len(want)
        print(f"  {str(rel)[:52]:52s} {len(k):3d} k_eff samples -> "
              f"{len(want):3d} snapshots  (worst time gap {worst:.3e} s)")

    a.out.write_text("\n".join(lines) + "\n")
    print(f"\n  {total} files, ~{total*a.snapshot_mb/1024:.0f} GiB at "
          f"{a.snapshot_mb:.0f} MB each")
    print(f"  wrote {a.out}")
    print(f"\n  rsync -av --files-from={a.out} <host>:<remote>/ {a.local}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

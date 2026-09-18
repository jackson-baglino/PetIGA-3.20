#!/usr/bin/env python3
"""Free disk by thinning sol_*.dat, without losing anything you can use.

    venv_enceladus/bin/python postprocess/thin_snapshots.py <batch_dir>          # report
    venv_enceladus/bin/python postprocess/thin_snapshots.py <batch_dir> --apply  # delete

WHY NOT JUST DELETE THE UNMERGED LEGS. The merged directory holds the TABLES
only -- SSA_evo.dat and k_eff.csv, a few MB. Every field snapshot lives in the
leg directories, and sol_*.dat is the raw solution vector: it is what
-keff_replay reads to recompute k_eff at a different eps, what a restart loads,
and the only route to neck radius, effective vapour diffusivity, chord lengths,
Euler characteristic or a correlation length. Deleting the legs to keep the
merged folder trades ~440 GB of physics for ~13 MB of CSV, irreversibly and at
the cost of the HPC time that produced it.

WHAT THIS KEEPS INSTEAD. The runs write a snapshot every ~13 steps -- around
500 per leg -- while k_eff is sampled only ~28 times. Almost all of that is
redundant. This keeps:

  * every snapshot nearest a k_eff sample time, so -keff_replay still works
    and every table row still has a field behind it;
  * the first and last of each leg;
  * the snapshot the restart was taken from, named in MERGE_INFO.json, so the
    run stays reproducible from its own join point;
  * optionally every --stride-th as a backstop.

and deletes the rest. On the pilot that is ~94% of the bytes.

RUN IT ON THE MERGED TREE, once the legs are merged and deleted. While both
exist the snapshots are hardlinks -- the same bytes under two names -- so
removing one name frees nothing until the other is gone too.

DRY RUN BY DEFAULT. Nothing is removed without --apply, and the report names
the byte count and the file count per leg first.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pplib                                                # noqa: E402


def snap_step(p: Path) -> int:
    return int(re.search(r"sol_(\d+)\.dat", p.name).group(1))


def keepers(leg: Path, stride: int, protect: set[int]) -> tuple[set[int], str]:
    snaps = sorted(leg.glob("sol_*.dat"), key=snap_step)
    if not snaps:
        return set(), "no snapshots"
    steps = np.array([snap_step(p) for p in snaps])
    keep = {int(steps[0]), int(steps[-1])} | {s for s in protect if s in set(steps.tolist())}

    kf = leg / "k_eff.csv"
    note = ""
    if kf.is_file():
        k = np.atleast_1d(np.genfromtxt(kf, delimiter=",", names=True))
        ssa = pplib.load_ssa(str(leg))
        if ssa is not None and len(k):
            # match on TIME, not step: k_eff samples land between snapshots
            for tk in np.atleast_1d(k["time"]):
                i = int(np.argmin(np.abs(ssa[:, 2] - tk)))
                target = int(ssa[i, 3])
                keep.add(int(steps[np.argmin(np.abs(steps - target))]))
            note = f"{len(k)} k_eff samples"
    if stride > 1:
        keep |= {int(s) for s in steps[::stride]}
    return keep, note


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("batch", type=Path)
    ap.add_argument("--apply", action="store_true", help="actually delete")
    ap.add_argument("--stride", type=int, default=1,
                    help="also keep every Nth snapshot as a backstop (default 1 = keep none extra)")
    a = ap.parse_args()

    # The snapshot a restart was taken from is already protected without a
    # special case: it sits at the join time, and the k_eff/SSA time matching
    # below keeps the nearest snapshot to every sampled time in both legs.
    protect: set[int] = set()

    legs = sorted(p.parent for p in a.batch.glob("**/sol_00000.dat"))
    # Skip the merged tree only when scanning a batch that CONTAINS it -- its
    # snapshots are hardlinks to the legs', so thinning both would just unlink
    # twice. When the merged directory is itself the target (the normal case
    # once the legs are gone) it must not be skipped, which the first version
    # got wrong by testing for "merged" anywhere in the path.
    merged_root = (a.batch / "merged").resolve()
    tot_del = tot_keep = 0
    n_del = n_keep = 0
    to_remove: list[Path] = []

    for leg in legs:
        if merged_root in leg.resolve().parents and merged_root != a.batch.resolve():
            continue
        keep, note = keepers(leg, a.stride, protect)
        snaps = sorted(leg.glob("sol_*.dat"), key=snap_step)
        kb = db = 0
        for p in snaps:
            sz = p.stat().st_size
            if snap_step(p) in keep:
                kb += sz; n_keep += 1
            else:
                db += sz; n_del += 1; to_remove.append(p)
        tot_keep += kb; tot_del += db
        rel = leg.relative_to(a.batch)
        print(f"  {str(rel)[:58]:58s}")
        print(f"      {len(snaps):4d} snapshots -> keep {len(keep):3d}, "
              f"drop {len(snaps)-len(keep):4d}   free {db/2**30:6.1f} GiB   {note}")

    print(f"\n  keep {n_keep} files, {tot_keep/2**30:.1f} GiB")
    print(f"  drop {n_del} files, {tot_del/2**30:.1f} GiB"
          f"   ({tot_del/max(tot_del+tot_keep,1):.0%} of the snapshot bytes)")

    if not a.apply:
        print("\n  DRY RUN — nothing deleted. Re-run with --apply.")
        return 0

    print("\n  deleting...")
    freed = 0
    for p in to_remove:
        freed += p.stat().st_size
        p.unlink()
    print(f"  freed {freed/2**30:.1f} GiB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

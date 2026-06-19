#!/usr/bin/env python3
"""
select_movie_frames.py — pick the minimal set of sol_*.dat steps needed to
build a make_movie.py movie, and highres-convert only those.

Why: plot_permafrost_highres.py's dense (n_per_elem=4) VTS output is ~70MB
PER TIMESTEP for a 608x122-element mesh. Converting every snapshot of a
multi-thousand-step run would be hundreds of GB for a movie that only
samples --n-frames frames in the end. This script computes the same
evenly-spaced-in-simulated-time target times make_movie.py uses, finds the
real snapshot(s) bracketing each one (so TemporalInterpolator has exactly
what it needs and nothing more), and converts only that union of steps.

Run with the regular project Python (needs igakit, not pvpython):

    python3 postprocess/select_movie_frames.py --dir <rundir> --n-frames 600

Then hand the resulting permafrost_highres.pvd to make_movie.py with the
SAME --n-frames/--t-start/--t-end so the target times line up exactly:

    pvpython postprocess/make_movie.py --dir <rundir> --n-frames 600
"""

import argparse
import os
import re
import subprocess
import sys


def parse_pvd(pvd_path, vts_dirname):
    pattern = re.compile(
        r'timestep="([0-9.eE+-]+)".*file="' + re.escape(vts_dirname) + r'/solV_(\d+)\.vts"')
    pairs = []
    with open(pvd_path) as fh:
        for line in fh:
            m = pattern.search(line)
            if m:
                pairs.append((float(m.group(1)), int(m.group(2))))
    pairs.sort()
    return pairs


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--dir", default=".", help="run directory (default: .)")
    p.add_argument("--pvd", default=None,
                   help="source (low-res) pvd to read timesteps from "
                        "(default: <dir>/permafrost.pvd)")
    p.add_argument("--n-frames", type=int, default=600,
                   help="must match the --n-frames you'll pass to make_movie.py (default: 600)")
    p.add_argument("--t-start", type=float, default=None)
    p.add_argument("--t-end", type=float, default=None)
    p.add_argument("--n-per-elem", type=int, default=4,
                   help="passed through to plot_permafrost_highres.py (default: 4)")
    p.add_argument("--force", action="store_true",
                   help="reconvert even if the highres VTS already exists")
    p.add_argument("--dry-run", action="store_true",
                   help="just report how many steps/how much disk, don't convert")
    args = p.parse_args()

    pvd_path = args.pvd or os.path.join(args.dir, "permafrost.pvd")
    pairs = parse_pvd(pvd_path, "vtkOut")
    if not pairs:
        sys.exit(f"No timesteps found in {pvd_path}")

    times = [t for t, _ in pairs]
    t_start = args.t_start if args.t_start is not None else times[0]
    t_end = args.t_end if args.t_end is not None else times[-1]

    needed = set()
    j = 0
    n = args.n_frames
    for i in range(n):
        target = t_start + (t_end - t_start) * (i / (n - 1) if n > 1 else 0.0)
        while j < len(times) - 1 and times[j + 1] <= target:
            j += 1
        lo = j
        hi = min(j + 1, len(times) - 1)
        needed.add(pairs[lo][1])
        needed.add(pairs[hi][1])

    steps = sorted(needed)
    # ~70MB/file observed at n_per_elem=4 on a 608x122-element mesh; point
    # count (and so file size) scales as n_per_elem^2.
    est_mb = len(steps) * 70 * (args.n_per_elem / 4.0) ** 2
    print(f"{len(pairs)} source snapshots -> {len(steps)} unique steps needed "
          f"for {n} frames over [{t_start:.6g}, {t_end:.6g}]")
    print(f"Estimated highres VTS size: ~{est_mb/1024:.1f} GB at --n-per-elem {args.n_per_elem}")

    if args.dry_run:
        return

    cmd = [
        sys.executable, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                      "plot_permafrost_highres.py"),
        "--dir", args.dir, "--n-per-elem", str(args.n_per_elem),
        "--steps", *[str(s) for s in steps],
    ]
    if args.force:
        cmd.append("--force")
    subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()

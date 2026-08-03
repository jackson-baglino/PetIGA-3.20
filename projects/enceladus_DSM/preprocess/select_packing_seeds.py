#!/usr/bin/env python3
"""Generate many packing seeds with the grains IN CONTACT, and rank them on how
much open pore they retain.

WHY
---
With `--contact-gap <= 0` the grains actually touch, which is what k_eff needs:
conductivity is largely a function of grain connectivity, and the rise of k_eff
through a run IS the growth of that connectivity. Starting from a pack whose
grains are held apart makes the early rise partly a numerical gap closing.

The cost is that in 2D a spanning contact network cuts the void into cells, so
the pore no longer percolates. That cannot be avoided by geometry -- but it can
be MITIGATED by choosing which random configuration to use. Seeds at the same
target porosity differ substantially in how much of the void stays in one
connected cluster and how wide the surviving throats are. This script generates
a batch and ranks them, so the study runs on the most open configurations rather
than on whichever seed happened to come first.

RANKING -- UNIFORMITY, not pore connectivity  (changed 2026-08-03)
------------------------------------------------------------------
Seeds are ranked by how EVENLY the grains are spread:

  max_void_radius_per_mean_r  radius of the largest circle fitting entirely in
                              the pore space, in units of the mean grain radius.
                              The direct measure of "one overly large void".
  local_density_cv            coefficient of variation of the solid fraction
                              over an 8x8 coarse grid; clumping raises it.

This REPLACES ranking by pore_largest_cluster_frac, which was actively harmful.
Maximising the fraction of void held in one connected cluster selects for
exactly the configurations that concentrate all the empty space in a single
hole -- so the winners were clumped grains beside a large void, which is the
opposite of a representative snow microstructure. The 0.362/seed12 packing
picked under the old rule had a void running the full height of the domain and
was chosen BECAUSE of it (pore 69%, the highest in its group).

The trade is explicit: uniform packings hold less of the void in one connected
cluster, so vapour transport is more obstructed at t=0. Pore connectivity is
still measured and reported, just not optimised for. In 2D it cannot be had
together with contacting grains anyway, and on Enceladus the vapour budget is
small.

Usage:
    python3 select_packing_seeds.py --porosity 0.325 --Lx 2.0e-3 \\
        --mean-r 100e-6 --seeds 1-20 --out inputs/packings_contact --keep 5
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
GEN = HERE / "generate_packing.py"


def parse_seeds(spec: str) -> list[int]:
    out: list[int] = []
    for part in spec.split(","):
        part = part.strip()
        if "-" in part:
            a, b = part.split("-")
            out.extend(range(int(a), int(b) + 1))
        elif part:
            out.append(int(part))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--porosity", type=float, required=True)
    ap.add_argument("--Lx", type=float, required=True)
    ap.add_argument("--Ly", type=float, default=None)
    ap.add_argument("--mean-r", type=float, default=100e-6,
                    help="mean grain radius [m]; 100 um matches the Molaro "
                         "et al. (2019) grain experiments")
    ap.add_argument("--sigma-ln", type=float, default=0.5)
    ap.add_argument("--contact-gap", type=float, default=0.0,
                    help="0 = exact tangency: grains touch where each one's phi "
                         "crosses 0.5, which is eps-independent and needs no "
                         "calibration. See generate_packing.py's docstring "
                         "before using anything else")
    ap.add_argument("--seeds", default="1-20", help="e.g. '1-20' or '1,3,7'")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--keep", type=int, default=5,
                    help="how many of the best seeds to retain")
    ap.add_argument("--non-periodic", action="store_true",
                    help="closed-box packing with every grain wholly inside the "
                         "walls; use for -periodic 0 runs so nothing is clipped")
    ap.add_argument("--raster", type=int, default=1024)
    ap.add_argument("--python", default=sys.executable)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args(argv)

    Ly = args.Ly if args.Ly is not None else args.Lx
    seeds = parse_seeds(args.seeds)
    args.out.mkdir(parents=True, exist_ok=True)

    print(f"generating {len(seeds)} seeds at porosity {args.porosity:.3f}, "
          f"{args.Lx*1e3:.2f} x {Ly*1e3:.2f} mm, mean r = {args.mean_r*1e6:.0f} um, "
          f"contact gap {args.contact_gap*1e6:+.2f} um")
    if args.dry_run:
        return 0

    rows = []
    for s in seeds:
        # Namespaced by porosity. With a bare "cand_seed{N}" a candidate that
        # FAILED for one porosity was left on disk, and the next porosity's run
        # skipped regeneration (the "metadata.json exists" test) and adopted it
        # -- filing a packing built for a different target under the new one.
        d = args.out / f"cand_phi{args.porosity:.3f}_seed{s}"
        if not (d / "metadata.json").is_file():
            # `--contact-gap=-4e-06`, not `--contact-gap -4e-06`: argparse reads
            # a leading '-' on the value as the start of another option.
            cmd = [args.python, str(GEN),
                   "--Lx", repr(args.Lx), "--porosity", repr(args.porosity),
                   "--seed", str(s), "--mean_r_m", repr(args.mean_r),
                   "--sigma_ln", repr(args.sigma_ln),
                   f"--contact-gap={args.contact_gap!r}",
                   "--out", str(d), "--raster", str(args.raster)]
            if args.non_periodic:
                cmd.append("--non-periodic")
            if args.Ly is not None:
                cmd += ["--Ly", repr(args.Ly)]
            r = subprocess.run(cmd, capture_output=True, text=True)
            # A non-zero exit means a real `problems` entry. Pore
            # non-percolation at gap <= 0 is a `warnings` entry and does NOT
            # fail, which is the whole point of the reclassification.
            if r.returncode != 0:
                print(f"  seed {s}: FAILED\n{r.stderr.strip()[-300:]}")
                continue
        m = json.loads((d / "metadata.json").read_text())
        if m.get("problems"):
            print(f"  seed {s}: skipped ({m['problems'][0][:60]})")
            continue
        rows.append({
            "seed": s, "dir": d,
            "phi": m["porosity_achieved"],
            "n": m["n_grains"],
            "Z": m["coordination_number"],
            "pore_frac": m["pore_largest_cluster_frac"],
            "pore_n": m["pore_n_clusters"],
            "p50": m["throat_gap_m"]["p50"],
            "p75": m["throat_gap_m"]["p75"],
            "void_r": m.get("max_void_radius_per_mean_r", float("inf")),
            "cv": m.get("local_density_cv", float("inf")),
        })

    if not rows:
        print("no usable packings generated", file=sys.stderr)
        return 1

    # Most uniform first: smallest worst-void, then least clumped.
    rows.sort(key=lambda r: (r["void_r"], r["cv"]))

    print(f"\n{'rank':>4} {'seed':>5} {'phi':>7} {'N':>5} {'Z':>5} "
          f"{'maxvoid/r':>10} {'densCV':>7} {'pore%':>6} {'p50[um]':>8}")
    for i, r in enumerate(rows):
        mark = "  <- keep" if i < args.keep else ""
        print(f"{i+1:>4} {r['seed']:>5} {r['phi']:>7.4f} {r['n']:>5} {r['Z']:>5.2f} "
              f"{r['void_r']:>10.2f} {r['cv']:>7.3f} "
              f"{r['pore_frac']*100:>5.1f}% {r['p50']*1e6:>8.1f}{mark}")

    # Promote the winners to stable names; drop the rest so the directory holds
    # only what the study will actually run.
    tag = f"phi{args.porosity:.3f}"
    for i, r in enumerate(rows):
        if i < args.keep:
            dest = args.out / f"{tag}_seed{r['seed']}"
            if dest.resolve() != r["dir"].resolve():
                if dest.exists():
                    shutil.rmtree(dest)
                r["dir"].rename(dest)
        else:
            shutil.rmtree(r["dir"], ignore_errors=True)

    # Anything still named cand_* for THIS porosity failed or was rejected;
    # leaving it behind is what caused the reuse bug above.
    for leftover in args.out.glob(f"cand_phi{args.porosity:.3f}_seed*"):
        shutil.rmtree(leftover, ignore_errors=True)

    print(f"\nkept {min(args.keep, len(rows))} packing(s) in {args.out}")
    print("Ranked by UNIFORMITY: smallest largest-void first, then least clumped.\n"
          "Pore percolation is not expected -- the grains are in contact by "
          "design -- and is no longer\noptimised for, because maximising it "
          "selects for a single large hole.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

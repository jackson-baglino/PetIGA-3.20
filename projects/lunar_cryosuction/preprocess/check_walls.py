#!/usr/bin/env python3
"""check_walls.py: are the walls of an igakit mesh smooth? (user, 2026-11-13: no sharp valleys, cusps or
other corners in the pore walls; they cause numerical trouble).

Evaluates the two wall curves (v = 0 and v = 1) at N and 4N points and reports, per wall, the smallest
radius of curvature in valleys (wall concave toward the pore) and at tips (convex), and the largest change
of tangent direction between neighboring samples. On a smooth curve that turn shrinks about 4x when the
sampling is 4x finer; at a cusp or corner it does not. Exit status 1 if any radius is below --min-eps * eps
or the turn does not shrink (ratio < 3), else 0.

    venv_cryosuction/bin/python preprocess/check_walls.py <mesh.dat> [--eps 8.584e-7] [--min-eps 20]
"""
from __future__ import annotations

import argparse
import sys

import numpy as np
from igakit.io import PetIGA


def wall_stats(nrb, v, n):
    u = np.linspace(0.0, 1.0, n)
    p = nrb(u, [v])[:, 0, :2]
    d1 = np.gradient(p, u, axis=0)
    d2 = np.gradient(d1, u, axis=0)
    sp = np.hypot(d1[:, 0], d1[:, 1])
    k = (d1[:, 0] * d2[:, 1] - d1[:, 1] * d2[:, 0]) / sp ** 3
    k = k[5:-5]
    turn = np.abs(np.diff(np.unwrap(np.arctan2(d1[:, 1], d1[:, 0])))).max()
    kv = k if v == 0.0 else -k                       # > 0: valley (concave toward the pore)
    r_valley = 1.0 / kv[kv > 0].max() if (kv > 0).any() else np.inf
    r_tip = 1.0 / (-kv[kv < 0]).max() if (kv < 0).any() else np.inf
    return r_valley, r_tip, turn


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("mesh")
    ap.add_argument("--eps", type=float, default=8.584e-7)
    ap.add_argument("--min-eps", type=float, default=20.0, help="smallest allowed radius, in eps")
    ap.add_argument("-n", type=int, default=40001)
    a = ap.parse_args(argv)
    nrb = PetIGA().read(a.mesh)
    bad = False
    for v, side in ((0.0, "bottom"), (1.0, "top")):
        rv, rt, t1 = wall_stats(nrb, v, a.n)
        _, _, t4 = wall_stats(nrb, v, 4 * a.n - 3)
        ratio = t1 / t4 if t4 > 0 else np.inf
        ok = min(rv, rt) >= a.min_eps * a.eps and (ratio >= 3.0 or t1 < 1e-6)
        bad |= not ok
        print(f"  {'OK  ' if ok else 'FAIL'} {side} wall: smallest valley radius {rv * 1e6:.1f} um "
              f"({rv / a.eps:.1f} eps), smallest tip radius {rt * 1e6:.1f} um ({rt / a.eps:.1f} eps), "
              f"tangent-turn ratio at 4x sampling {ratio:.2f} (smooth ~4)")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())

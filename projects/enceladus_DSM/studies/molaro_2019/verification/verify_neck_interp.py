#!/usr/bin/env python3
"""
verify_neck_interp.py — gate for neck_width.py's sub-cell crossing.

WHY THIS EXISTS
---------------
The neck-growth curves carried a second, oscillatory mode riding on the
power law. It was not physics: it was the sub-cell interpolation in
`chord_width()`.

The `.vts` snapshots are written on a COARSER grid than the solve -- 1080 x 541
against a 5394 x 2697 mesh -- so the sample spacing is dy = 4.17e-7 m = 3.5 eps.
Interpolating phi LINEARLY across 3.5 eps of a sigmoid is not a small
correction: the two samples bracketing phi = 0.5 can sit at phi = 0.03 and 0.97,
and the crossing error then depends on where the true interface falls inside the
cell. That error is periodic in the sub-cell offset, and a neck growing steadily
outward sweeps through offsets at a steady rate -- so the artefact appears as a
sinusoid on top of the growth curve. On the round-2 arms it correlated with the
radial sub-cell phase at |R| = 0.66-0.88.

The fix is to interpolate in logit(phi). The model's own 1D equilibrium profile
is phi = 1/(1 + exp(-s/eps)), so logit(phi) = s/eps is LINEAR in distance and a
straight line through two samples is EXACT at any spacing.

WHAT THIS CHECKS
----------------
A synthetic axisymmetric neck with an analytically known width, sampled the way
a .vts is sampled, at a spacing swept across the range the real runs use. For
each sub-cell offset the measured width is compared against the exact answer.

  1. logit interpolation is exact for the equilibrium profile (to ~1e-12 um);
  2. it is UNBIASED in the sub-cell offset -- which is the property that kills
     the oscillation, and is separate from being small on average;
  3. linear-in-phi is neither, and the gate records how large its error is so
     the regression cannot silently come back.

Usage:  python studies/molaro_2019/verification/verify_neck_interp.py
Writes: neck_interp.csv next to this script; non-zero exit on failure.
"""

import csv
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent.parent / "postprocess"))

from neck_width import chord_width, _cross, _logit  # noqa: E402

EPS = 1.18e-7          # production interface parameter
R_NECK = 27.0e-6       # a representative late-time neck RADIUS
TOL_LOGIT_UM = 1.0e-6  # exactness tolerance on the width [um]
TOL_BIAS_UM = 1.0e-3   # spread across sub-cell offsets [um]


def linear_cross(ya, yb, pa, pb, level):
    """The OLD estimator: linear in phi. Kept here, not imported, precisely so
    this gate still measures it after it is gone from the shipping code."""
    if pa == pb:
        return ya
    return ya + (pa - level) / (pa - pb) * (yb - ya)


def synthetic(dy, offset):
    """(y, phi) for one column of a synthetic axisymmetric neck, in metres.

    The GRID is fixed and starts at y = 0 -- that is what an axisymmetric run
    looks like, and chord_width() relies on it: at the axis the ice is already
    above the contour level, so index 0 is taken as the inner end with no
    interpolation. Moving the grid instead (the obvious first way to write this)
    puts y[0] > 0 and silently loses that offset from every measured width.

    It is the INTERFACE that sweeps through the cell: the profile is centred at
    R_NECK + offset*dy, so `offset` scans one full sub-cell period exactly as a
    growing neck does.
    """
    y = np.arange(0, 200) * dy
    r_true = R_NECK + offset * dy
    # clipped so the far tail does not overflow exp; phi there is 0 either way
    phi = 1.0 / (1.0 + np.exp(np.clip((y - r_true) / EPS, -500.0, 500.0)))
    return y, phi, r_true


def measure(dy, offset, cross):
    """Measured neck WIDTH [um] for the synthetic column, using `cross`.

    Axisymmetric, so width = 2 x radius -- matching neck_width.py's --axisym
    doubling. Returns micrometres, like everything else this gate prints.
    """
    y, phi, _ = synthetic(dy, offset)
    hi = np.flatnonzero(phi >= 0.5)[-1]
    r = cross(y[hi], y[hi + 1], phi[hi], phi[hi + 1], 0.5)
    return 2.0 * r * 1e6


def main():
    offsets = np.linspace(0.0, 1.0, 101, endpoint=False)
    spacings = {"production .vts (3.5 eps)": 2.25e-4 / 540,
                "half that (1.8 eps)": 2.25e-4 / 1080,
                "coarse (7.1 eps)": 2.25e-4 / 270}

    rows = [("spacing", "dy_over_eps", "estimator", "max_abs_err_um",
             "ptp_over_offset_um", "verdict")]
    ok = True
    print(f"nominal neck width = {2*R_NECK*1e6:.4f} um, swept over one sub-cell "
          f"period   (eps = {EPS:.3e} m)\n")
    print(f"{'spacing':26s} {'dy/eps':>7s} {'estimator':16s} "
          f"{'max|err|':>11s} {'ptp':>11s}  verdict")
    print("-" * 88)

    for name, dy in spacings.items():
        for est, cross, tol in (("linear-in-phi", linear_cross, None),
                                ("logit", _cross, TOL_LOGIT_UM)):
            # the exact answer MOVES with the offset -- comparing against a fixed
            # width would report the sweep itself as error
            exact = np.array([2.0 * (R_NECK + o * dy) * 1e6 for o in offsets])
            errs = np.array([measure(dy, o, cross) for o in offsets]) - exact
            mx, ptp = np.abs(errs).max(), errs.ptp()
            if tol is None:
                verdict = "reference"          # not gated: this is the old one
            else:
                good = mx < tol and ptp < TOL_BIAS_UM
                verdict = "PASS" if good else "FAIL"
                ok &= good
            print(f"{name:26s} {dy/EPS:7.2f} {est:16s} "
                  f"{mx:11.3e} {ptp:11.3e}  {verdict}")
            rows.append((name, f"{dy/EPS:.4f}", est, f"{mx:.6e}",
                         f"{ptp:.6e}", verdict))
        print()

    with (HERE / "neck_interp.csv").open("w", newline="") as fh:
        csv.writer(fh, lineterminator="\n").writerows(rows)

    # The end-to-end path must agree with the direct call, or the gate is
    # testing something chord_width() does not actually do.
    dy = 2.25e-4 / 540
    y, phi, r_true = synthetic(dy, 0.37)
    w_chord = 2.0 * chord_width(phi, y, 0.5) * 1e6
    w_exact = 2.0 * r_true * 1e6
    agree = abs(w_chord - w_exact) < TOL_LOGIT_UM
    print(f"chord_width() end-to-end, offset 0.37: "
          f"{'PASS' if agree else 'FAIL'}  "
          f"({w_chord:.9f} vs exact {w_exact:.9f} um)")
    ok &= agree

    print()
    print("ALL CHECKS PASSED" if ok else "SOME CHECKS FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

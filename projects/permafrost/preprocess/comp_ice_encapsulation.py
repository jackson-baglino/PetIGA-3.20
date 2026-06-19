#!/usr/bin/env python3
"""
comp_ice_encapsulation.py — fit a circle to a SedimentBump() floor bump and
compute the enlarged concentric circle that encapsulates it in ice.

Background: sediment grains intruding into the domain are modeled as
boundary bumps (SedimentBump() in src/initial_conditions.c), not as a
separate phase field -- direct phase-field modeling of the sediment grain
caused numerical issues. To make the ICE around such a bump look like the
grain is actually encapsulated (constant-thickness coating on every exposed
side), approximate the bump's curve with a circle through its two edges and
its peak (the standard sagitta/circular-segment formula), then use the SAME
center with a LARGER radius (+ the desired ice thickness) as a plain
-ice_grain_cx/cy/ax/ay circle (ax=ay=R_ice). Because it's concentric with
the sediment-fit circle, the ice wraps the grain at uniform radial thickness
with no separate lateral window -- no notches, no overhang, no new C code.

The fitted circle's center sits BELOW y=0 (cy_sed is negative) for any bump
that's wider than it is tall (height_ratio < 1), which is the case for
every bump in this project so far -- that's expected, not a bug: a shallow,
wide bump corresponds to a large, nearly-flat circle whose center is far
below the domain.

Usage:
    python3 preprocess/comp_ice_encapsulation.py --R 0.4e-5 --H 0.2e-5 --thickness 0.2e-5

Prints the -ice_grain_cx/cy/ax/ay values to paste into an opts file (cx is
left as a placeholder matching -sed_grain_x for that bump; only cy and the
radius are derived here).
"""

import argparse
import math


def fit_sediment_circle(R: float, H: float):
    """Circle through (-R,0), (R,0), (0,H) -- center directly below/above
    the bump's peak at x=0 (relative to the bump's own center)."""
    R_sed = (H**2 + R**2) / (2.0 * H)
    cy_sed = H - R_sed
    return R_sed, cy_sed


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--cx", type=float, default=0.0,
                   help="bump center x [m] (matches -sed_grain_x for this bump; "
                        "default 0.0, i.e. relative)")
    p.add_argument("--R", type=float, required=True, help="bump half-width [m] (-sed_grain_R)")
    p.add_argument("--H", type=float, required=True, help="bump peak height [m] (-sed_grain_h)")
    p.add_argument("--thickness", type=float, required=True,
                   help="desired ice thickness [m], measured radially from the "
                        "sediment-fit circle")
    args = p.parse_args()

    R_sed, cy_sed = fit_sediment_circle(args.R, args.H)
    R_ice = R_sed + args.thickness

    footprint_ice = math.sqrt(max(R_ice**2 - cy_sed**2, 0.0))
    footprint_sed = math.sqrt(max(R_sed**2 - cy_sed**2, 0.0))

    print(f"R_sed (sediment-fit circle radius): {R_sed:.6e}")
    print(f"cy_sed (circle center, y, below 0): {cy_sed:.6e}")
    print(f"R_ice (ice circle radius):          {R_ice:.6e}")
    print(f"ice footprint half-width at y=0:    {footprint_ice:.6e}")
    print(f"sed-fit footprint half-width at y=0:{footprint_sed:.6e}  "
          f"(should equal --R={args.R:.4e})")
    print()
    print("Paste into the opts file (append to existing -ice_grain_* lists):")
    print(f"  cx={args.cx + 0.0:.6e}  cy={cy_sed:.6e}  ax=ay={R_ice:.6e}")


if __name__ == "__main__":
    main()

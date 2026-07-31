#!/usr/bin/env python3
"""Colour-blindness separation check for categorical plot palettes.

A Python port of the essential checks from the data-viz palette validator,
because this machine has no JS runtime and the project's plots are matplotlib.

What it checks, for every pair of series colours:
  * separation under normal vision
  * separation under deuteranopia and protanopia (the two common forms)

Separation is Euclidean distance in OKLab x100, which is near-perceptually
uniform, so one threshold means the same thing across the whole gamut.

  dE >= 8   pass
  6 <= dE < 8   floor: legal only with a second encoding (marker shape, dash)
  dE < 6    fail
  normal-vision dE < 15 is a hard fail regardless of CVD performance

Usage:
    python3 validate_palette.py "#2a78d6,#eb6834,#1baf7a"
"""
import sys
import numpy as np

# Viénot, Brettel & Mollon (1999) dichromat simulation, applied in linear RGB.
LMS_FROM_LRGB = np.array([[0.31399022, 0.63951294, 0.04649755],
                          [0.15537241, 0.75789446, 0.08670142],
                          [0.01775239, 0.10944209, 0.87256922]])
LRGB_FROM_LMS = np.linalg.inv(LMS_FROM_LRGB)
CVD_LMS = {
    "deuteranopia": np.array([[1.0, 0.0, 0.0],
                              [0.9513092, 0.0, 0.04866992],
                              [0.0, 0.0, 1.0]]),
    "protanopia":   np.array([[0.0, 1.05118294, -0.05116099],
                              [0.0, 1.0, 0.0],
                              [0.0, 0.0, 1.0]]),
}

# OKLab (Björn Ottosson)
LMS_FROM_LRGB_OK = np.array([[0.4122214708, 0.5363325363, 0.0514459929],
                             [0.2119034982, 0.6806995451, 0.1073969566],
                             [0.0883024619, 0.2817188376, 0.6299787005]])
OKLAB_FROM_LMSP = np.array([[0.2104542553, 0.7936177850, -0.0040720468],
                            [1.9779984951, -2.4285922050, 0.4505937099],
                            [0.0259040371, 0.7827717662, -0.8086757660]])


def srgb_to_linear(c):
    c = np.asarray(c, dtype=float)
    return np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)


def hex_to_linear(h):
    h = h.strip().lstrip("#")
    return srgb_to_linear(np.array([int(h[i:i + 2], 16) for i in (0, 2, 4)]) / 255.0)


def simulate(lrgb, kind):
    if kind == "normal":
        return lrgb
    lms = LMS_FROM_LRGB @ lrgb
    return LRGB_FROM_LMS @ (CVD_LMS[kind] @ lms)


def oklab(lrgb):
    lms = LMS_FROM_LRGB_OK @ np.clip(lrgb, 0.0, 1.0)
    return OKLAB_FROM_LMSP @ np.cbrt(lms)


def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return 2
    hexes = [h for h in argv[1].split(",") if h.strip()]
    lin = [hex_to_linear(h) for h in hexes]

    worst_overall, ok = None, True
    for kind in ("normal", "deuteranopia", "protanopia"):
        labs = [oklab(simulate(c, kind)) for c in lin]
        print(f"\n{kind}")
        for i in range(len(hexes)):
            for j in range(i + 1, len(hexes)):
                dE = float(np.linalg.norm(labs[i] - labs[j])) * 100.0
                floor = 15.0 if kind == "normal" else 8.0
                if dE >= floor:
                    verdict = "PASS"
                elif kind != "normal" and dE >= 6.0:
                    verdict = "FLOOR (needs 2nd encoding)"
                    ok = False
                else:
                    verdict = "FAIL"
                    ok = False
                print(f"  {hexes[i]} vs {hexes[j]}   dE = {dE:6.2f}   {verdict}")
                if worst_overall is None or dE < worst_overall[0]:
                    worst_overall = (dE, kind, hexes[i], hexes[j])

    d, kind, a, b = worst_overall
    print(f"\nworst pair overall: {a} vs {b} under {kind}, dE = {d:.2f}")
    print("RESULT:", "PASS" if ok else "NOT CLEAN")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))

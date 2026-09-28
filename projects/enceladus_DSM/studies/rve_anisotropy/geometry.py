#!/usr/bin/env python3
"""Geometry of every packing on disk: contact fabric and the y-seam.

    venv_enceladus/bin/python studies/rve_anisotropy/geometry.py

No solver, no simulation: grains.dat only. Two questions.

1. CONTACT FABRIC. In a granular solid heat crosses from grain to grain through
   contacts, and the conductivity tensor follows the orientation distribution
   of the contacts (Batchelor & O'Brien 1977, one of the classical results for
   a packed bed). So

        F = < n n^T >   over contacts,   n = unit branch vector i -> j

   predicts, from geometry alone and independently of any solver, which way
   k_eff is anisotropic: F_xx > F_yy means more contacts carry heat
   horizontally, i.e. k_xx > k_yy; F_xy != 0 tilts the principal axes, i.e.
   k_xy != 0 with the same sign. "Contact" = gap below the diffuse band
   9.2*eps, eps = R_ave/50, because the solver joins every such pair with ice
   (studies/packing_design/README.md section 1).

2. THE Y SEAM. x is periodic because deposition itself wraps sideways, so
   contacts across x = 0 are real deposition contacts. y is periodic by
   construction: the window is cut from a taller bed and y = 0 is identified
   with y = Ly, overlaps along that line are relaxed apart and gaps are left
   for the void filler. If the seam carries fewer contacts than the interior,
   it is a weak layer in SERIES with vertical heat flow and lowers k_yy -- a
   generator artifact that would masquerade as anisotropy. Measured as the
   number of contacts whose branch vector crosses a horizontal line y = c, per
   unit length, for 256 lines; the seam is c = 0. Vertical lines x = c are the
   control, since the x boundary is not a seam.

Writes geometry.csv and seam_profile.csv next to this file.
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PROJ = HERE.parents[1]
sys.path.insert(0, str(PROJ / "preprocess"))
import packing_lib as pl                                    # noqa: E402

EPS_OVER_R = 1.0 / 50.0      # the campaign's R_feat/R_ave = 1/25 at safety 0.5
BAND = 9.2 * EPS_OVER_R      # in mean radii
N_LINES = 256


def load(d: Path):
    m = json.loads((d / "metadata.json").read_text())
    rows = [l.split() for l in (d / "grains.dat").read_text().splitlines()
            if l.strip() and not l.startswith("#")]
    raw = np.array([[float(v) for v in r] for r in rows[1:]])
    c, r = raw[:, :2], raw[:, 2]
    keep = ((c[:, 0] >= 0) & (c[:, 0] < m["Lx"]) & (c[:, 1] >= 0) & (c[:, 1] < m["Ly"]))
    return m, c[keep].copy(), r[keep].copy()


def contacts(c, r, Lx, Ly, px, py, band_m):
    """Bonds (i, j, ox, oy) with gap <= band_m, and their branch vectors."""
    b = pl.delaunay_bonds(c, Lx, Ly, px, py)
    if len(b) == 0:
        return np.empty((0, 4), int), np.empty((0, 2))
    i, j = b[:, 0], b[:, 1]
    cj = c[j] + b[:, 2:4] * np.array([Lx, Ly])
    v = cj - c[i]
    gap = np.hypot(v[:, 0], v[:, 1]) - (r[i] + r[j])
    k = gap <= band_m
    return b[k], v[k]


def crossings(c, bonds, v, L, axis, lines):
    """Contacts per unit length crossing each line coordinate = `lines` along `axis`."""
    a = c[bonds[:, 0], axis]
    bb = a + v[:, axis]
    lo, hi = np.minimum(a, bb), np.maximum(a, bb)
    out = np.zeros(len(lines))
    for s in (-L, 0.0, L):                       # a branch may straddle the wrap
        out += ((lo[None, :] < lines[:, None] + s) &
                (hi[None, :] > lines[:, None] + s)).sum(1)
    return out / L


def analyse(d: Path, family: str):
    m, c, r = load(d)
    Lx, Ly = m["Lx"], m["Ly"]
    px, py = bool(m.get("periodic_x", True)), bool(m.get("periodic_y", True))
    R = float(np.mean(r))
    bonds, v = contacts(c, r, Lx, Ly, px, py, BAND * R)
    n = v / np.hypot(v[:, 0], v[:, 1])[:, None]
    F = (n[:, :, None] * n[:, None, :]).mean(0)
    w, e = np.linalg.eigh(F)
    ang = float(np.degrees(np.arctan2(e[1, 1], e[0, 1])) % 180.0)

    row = dict(family=family, name=d.name, L_over_R=round(Lx / R, 1),
               porosity=m.get("porosity_achieved", m.get("porosity_raster")),
               periodic=("x" if px else "") + ("y" if py else "") or "none",
               n_grains=len(r), n_contacts=len(bonds),
               z_band=2 * len(bonds) / len(r),
               F_xx=F[0, 0], F_yy=F[1, 1], F_xy=F[0, 1],
               F_ratio=F[0, 0] / F[1, 1], F_major_deg=ang)
    prof = None
    if px and py:
        lines = (np.arange(N_LINES) + 0.5) * Ly / N_LINES
        hy = crossings(c, bonds, v, Ly, 1, lines)            # horizontal lines
        hx = crossings(c, bonds, v, Lx, 0, lines * Lx / Ly)  # vertical lines
        seam = (lines < R) | (lines > Ly - R)                 # within 1 R of y = 0
        interior = ~seam
        row.update(
            seam_contacts_per_L=hy[seam].mean(), interior_contacts_per_L=hy[interior].mean(),
            seam_ratio=hy[seam].mean() / hy[interior].mean(),
            seam_z=(hy[seam].mean() - hy[interior].mean()) / hy[interior].std(),
            xline_seam_ratio=hx[seam].mean() / hx[interior].mean(),
            seam_max_shift_R=m.get("seam_max_shift_m", 0.0) / R)
        prof = (lines / Ly, hy / hy[interior].mean(), hx / hx[interior].mean())
    return row, prof


def main():
    sets = [
        ("pilot_LR40",   PROJ / "inputs/packings/pilot_LR40"),
        ("rev_LR64",     PROJ / "inputs/packings/rev_LR64"),
        ("prod_LR40_xy", PROJ / "inputs/packings/periodic_xy"),
        ("prod_LR40_x",  PROJ / "inputs/packings/periodic_x"),
        ("prod_LR40_y",  PROJ / "inputs/packings/periodic_y"),
        ("prod_LR40_none", PROJ / "inputs/packings/periodic_none"),
        ("design",       PROJ / "studies/packing_design/packings"),
    ]
    rows, profs = [], {}
    for fam, root in sets:
        for d in sorted(root.glob("*/")):
            if not (d / "grains.dat").is_file():
                continue
            row, prof = analyse(d, fam)
            rows.append(row)
            if prof is not None:
                profs[f"{fam}/{d.name}"] = prof
            print(f"  {fam:15s} {d.name:42s} L/R {row['L_over_R']:5.1f}  "
                  f"F_xx/F_yy {row['F_ratio']:.3f}  F_xy {row['F_xy']:+.4f}"
                  + (f"  seam {row['seam_ratio']:.3f} (z {row['seam_z']:+.2f})"
                     if "seam_ratio" in row else ""))

    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(HERE / "geometry.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    with open(HERE / "seam_profile.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["packing", "y_over_L", "hline_rel", "vline_rel"])
        for name, (y, hy, hx) in profs.items():
            for a, b, cc in zip(y, hy, hx):
                w.writerow([name, f"{a:.5f}", f"{b:.5f}", f"{cc:.5f}"])
    print(f"  wrote {HERE / 'geometry.csv'} and seam_profile.csv")


if __name__ == "__main__":
    main()

"""
cone2.py — 3D cone volume via revolution of a triangular cross-section.

Builds a solid cone by:
  1. Defining two line segments from the base (y=0) to the apex (y=H):
       l1 — inner edge: (Rin, 0) to (0, H)
       l2 — outer edge: (Rout, 0) to (0, H)
  2. Ruling between l1 and l2 to get a triangular patch (the cross-section).
  3. Revolving the patch 360° around the Y axis to produce the 3D cone volume.
  4. Uniformly refining in all three parametric directions.
  5. Degree elevation to cubic in x/y and quadratic in z.

Parameters (hardcoded):
  H    = 1.0  — cone height
  Rin  = 0.5  — inner radius (creates a hollow cone / frustum if Rin > 0)
  Rout = 1.0  — outer radius

Outputs:
  cone2.dat — PetIGA binary mesh file

Usage:
  python3 cone2.py
"""

from igakit.cad import *
from numpy import linspace

# ---------------------------------------------------------------------------
# Geometry parameters
# ---------------------------------------------------------------------------
H    = 1.0   # cone height
Rin  = 0.5   # inner radius (> 0 gives a hollow cone)
Rout = 1.0   # outer radius

# ---------------------------------------------------------------------------
# Build 2D cross-section by ruling between inner and outer slant edges
# ---------------------------------------------------------------------------
l1  = line((Rin,  0), (0, H))   # inner slant: base-inner to apex
l2  = line((Rout, 0), (0, H))   # outer slant: base-outer to apex
sqr = ruled(l1, l2)             # triangular patch (cross-section)

# Knot insertions for uniform mesh after revolution
insert_x0 = linspace(0, 1, 10)[1:-1]
insert_y0 = linspace(0, 1, 10)[1:-1]
insert_z0 = linspace(0, 1, 10)[1:-1]

# ---------------------------------------------------------------------------
# Visualize the 2D cross-section
# ---------------------------------------------------------------------------
from igakit.plot import plt
plt.kplot(sqr)
plt.show()

# ---------------------------------------------------------------------------
# Revolve around Y axis to produce the 3D cone volume
# ---------------------------------------------------------------------------
vol = revolve(sqr, point=(0, H, 0), axis=(0, 1, 0), angle=2 * Pi)

plt.kplot(vol)
plt.show()

# ---------------------------------------------------------------------------
# Uniform refinement in all three parametric directions
# ---------------------------------------------------------------------------
vol.refine(0, insert_x0)
vol.refine(1, insert_y0)
vol.refine(2, insert_z0)

# Degree elevation: cubic in x and y, quadratic in z
vol.elevate(0, 2)
vol.elevate(1, 2)
vol.elevate(2, 1)

plt.kplot(vol)
plt.show()
print(vol.degree)

# ---------------------------------------------------------------------------
# Export
# ---------------------------------------------------------------------------
from igakit.io import PetIGA
PetIGA().write("cone2.dat", vol)

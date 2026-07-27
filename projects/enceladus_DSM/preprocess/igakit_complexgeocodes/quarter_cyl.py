"""
quarter_cyl.py — 3D quarter-cylinder volume mesh.

Builds a solid quarter-cylinder (90° wedge of a full cylinder) by:
  1. Ruling between a degenerate circle at r=0 (a point) and a circle
     of radius R at angle [0, pi/2] to get a quarter-disc cross-section.
  2. Extruding the disc in the z-direction by length Lz.
  3. Elevating the degree in x and z to quadratic.
  4. Uniformly refining all three parametric directions.

The result is a single-patch NURBS volume suitable for axisymmetric or
3D IGA simulations in a cylindrical domain (e.g. flow in a pipe sector).

Parameters (hardcoded):
  R  = 0.5   — cylinder radius
  Lz = 1.0   — cylinder length (z-extent)
  Refinement: 8 knots in r/theta, 100 knots in z

Outputs:
  quarter_cylinder.dat — PetIGA binary mesh file
  quarter_cylinder.vtk — VTK structured grid for visualization

Usage:
  python3 quarter_cyl.py
"""

from igakit.io import PetIGA, VTK
from igakit.plot import plt
from numpy import *
import matplotlib.pyplot as plt2
import glob
from igakit.nurbs import NURBS
from igakit.transform import transform
from igakit.cad import *

# ---------------------------------------------------------------------------
# Geometry parameters
# ---------------------------------------------------------------------------
R  = 0.5   # cylinder radius
Lz = 1.0   # cylinder length

# ---------------------------------------------------------------------------
# Build the quarter-disc cross-section by ruling between r=0 and r=R arcs.
# circle(radius=0) degenerates to a point at the origin.
# ---------------------------------------------------------------------------
from numpy import pi
c1  = circle(radius=0., angle=(0, pi / 2.))   # degenerate (point at origin)
c2  = circle(radius=R,  angle=(0, pi / 2.))   # quarter-circle arc at radius R
srf = ruled(c1, c2).transpose()

plt.plot(srf)
plt.show()

# ---------------------------------------------------------------------------
# Extrude disc in z to form the quarter-cylinder volume
# ---------------------------------------------------------------------------
vol = extrude(srf, displ=Lz, axis=2)

# Elevate degree: quadratic in r/theta (axis 0) and z (axis 2)
vol.elevate(0, 1)
vol.elevate(2, 1)

plt.plot(vol)
plt.show()
print(vol.degree)
print(vol.knots)

# ---------------------------------------------------------------------------
# Uniform refinement
#   r/theta (axes 0 and 1): 8 interior knots each
#   z (axis 2): 100 interior knots
# ---------------------------------------------------------------------------
to_insert  = linspace(0, 1, 8)[1:-1]         # 6 interior knots in r/theta
to_insertZ = linspace(0, 1, 100 + 1)[1:-1]   # 99 interior knots in z

vol.refine(0, to_insert)
vol.refine(1, to_insert)
vol.refine(2, to_insertZ)

print(vol.degree)
print(vol.knots)

# ---------------------------------------------------------------------------
# Export
# ---------------------------------------------------------------------------
PetIGA().write("quarter_cylinder.dat", vol)
VTK().write("quarter_cylinder.vtk", vol)

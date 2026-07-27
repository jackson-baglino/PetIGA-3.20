"""
nsch_axisym_solidtopextra_cosine.py — 2D axisymmetric tube with cosine neck.

Generates a 2D NURBS mesh for a fluid-solid domain consisting of a tube
with a smooth necking constriction. The neck profile is defined by a cosine
function, giving a C^1-smooth transition between the full-radius tube and
the minimum-radius neck.

Domain layout (y-direction, bottom to top):
  [0, r_min + solid_thickness]           — fluid + inner solid at the neck
  [r_min + solid_thickness, H]           — outer solid

Three control-point rows define the surface:
  C2 : bottom axis (y = 0), flat
  C3 : fluid-solid interface (y varies with the cosine neck profile)
  C1 : outer wall (y = H in the tube, varies at the neck)

The knot vector in y contains a repeated knot at the fluid-solid interface
(y = r_max/H in parametric space), creating a C0 join there so that
different material properties can be applied across the interface.

Outputs:
  nsch_axisym_solidtopextra.dat — PetIGA binary mesh file
  nsch_axisym_solidtopextra.vtk — VTK structured grid for visualization

Usage:
  python3 nsch_axisym_solidtopextra_cosine.py <Lx> <L> <Nx> <r_min> <solid_thickness>

  Lx              : total domain length [m]
  L               : length of the constriction region [m]
  Nx              : number of elements in x (actual Nx + p control points)
  r_min           : minimum (neck) fluid radius [m]
  solid_thickness : solid wall thickness [m]

Example:
  python3 nsch_axisym_solidtopextra_cosine.py 12 5.0 1200 0.25 0.1
"""

import numpy as np
import copy
import math
import sys
from scipy.io import savemat

from numpy import linspace
from igakit.cad import *
from igakit.io import PetIGA, VTK
from igakit.plot import plt
from igakit.nurbs import NURBS
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Parameters — read from command line
# ---------------------------------------------------------------------------
Lx              = float(sys.argv[1])          # total domain length
L               = float(sys.argv[2])          # constriction length
p               = 2                           # NURBS degree
p_basis         = 2                           # target degree after elevation
Nx              = int(sys.argv[3]) + p        # total x control points
r_min           = float(sys.argv[4])          # minimum (neck) fluid radius
r_max           = 1.0                         # maximum (tube) fluid radius
solid_thickness = float(sys.argv[5])
H               = r_max + solid_thickness     # total domain height
Hd              = r_min + solid_thickness     # domain height at the neck
Ny              = int(float(H / Lx) * Nx)    # y elements (aspect-ratio matched)
print(Lx, L, H, Ny)

# ---------------------------------------------------------------------------
# Build x-coordinate array and compute cosine neck profiles
# ---------------------------------------------------------------------------
x          = np.linspace(0, Lx, Nx)
neck_start = 0.5 * (Lx - L)   # x-position where the neck begins

# Cosine blending factor: +1 in the tube, -1 at the neck centre
cosine_term = np.cos(np.pi * (x - neck_start - 0.5 * L) / (0.5 * L))

# Fluid-solid interface radius (bottom of the solid layer) along x
tube_neck_min = 0.5 * r_min * ((1. + (r_max / r_min)) + (1. - (r_max / r_min)) * cosine_term)

# Outer wall height along x
tube_neck_max = 0.5 * Hd * ((1. + (H / Hd)) + (1. - (H / Hd)) * cosine_term)

print(neck_start)

# Outside the constriction region, revert to straight tube profiles
for i in range(Nx):
    if (x[i] < neck_start) or (x[i] > (neck_start + L)):
        tube_neck_min[i] = r_max
        tube_neck_max[i] = H

# ---------------------------------------------------------------------------
# Initialise three boundary curves and apply the neck profiles
# ---------------------------------------------------------------------------
C1 = np.linspace((0, H), (Lx, H), Nx)                              # outer wall
C3 = np.linspace((0, H - solid_thickness), (Lx, H - solid_thickness), Nx)  # fluid-solid interface
C2 = np.linspace((0, 0), (Lx, 0), Nx)                              # bottom axis (flat)

for i in range(Nx):
    C3[i, 1] = tube_neck_min[i]   # fluid-solid interface follows neck
    C1[i, 1] = tube_neck_max[i]   # outer wall follows neck

# ---------------------------------------------------------------------------
# Stack curves and build NURBS surface
# ---------------------------------------------------------------------------
C = np.stack((C2, C3, C1))   # bottom → interface → top

# Open-uniform knot vector in x
Ux = np.linspace(0, 1, Nx + 1 - p)
for i in range(p):
    Ux = np.append(0, Ux)
    Ux = np.append(Ux, 1)

# y knot vector: repeated interior knot at r_max/H creates a C0 join at the
# fluid-solid interface, allowing separate material properties on each side.
Uy = [0, 0, float(r_max / H), 1, 1]

srf = NURBS([Uy, Ux], C).transpose()
print("Initial degree:", srf.degree)
srf.elevate(1, p_basis - 1)    # elevate y to p_basis
srf.elevate(0, p_basis - p)    # elevate x to p_basis
print("Elevated degree:", srf.degree)
print("Knots:", srf.knots)

# ---------------------------------------------------------------------------
# Uniform y-refinement (remove any y-knots that coincide with existing ones)
# ---------------------------------------------------------------------------
insert_y = np.linspace(0, 1, Ny + 1)[1:-1]
for idx, j in enumerate(Uy):
    index    = np.argwhere(insert_y == j)
    insert_y = np.delete(insert_y, index)

srf.refine(1, insert_y)

# ---------------------------------------------------------------------------
# Export
# ---------------------------------------------------------------------------
VTK().write("nsch_axisym_solidtopextra.vtk", srf)
PetIGA().write("nsch_axisym_solidtopextra.dat", srf)

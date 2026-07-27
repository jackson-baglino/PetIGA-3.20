"""
MeshGenerate_Wedge.py — 2D mesh with a triangular wedge on the bottom boundary.

Generates a 2D NURBS patch for a rectangular channel [0,L] x [0,H] whose
bottom boundary contains a symmetric triangular wedge (linear ramp up, then
ramp down). The mesh uses anisotropic refinement in the wall-normal (y)
direction to resolve boundary layer physics near the bottom.

Geometry parameters (hardcoded):
  L    = 4.5 m   — domain length
  H    = 0.6 m   — domain height
  w_h  = 0.2 m   — wedge peak height
  x_ws = 0.2 m   — wedge start x-position
  x_wp = x_ws + 3*w_h  — wedge peak x-position
  x_we = x_wp + 4.5*w_h — wedge end x-position
  N    = 716 + p  — total control points along x
  Ny   = 72       — baseline y-refinement elements
  Ny_BL = 20      — number of boundary layer elements near the wall
  Scale = 0.3     — boundary layer thickness (fraction of H)
  TransGrid = 60  — transition elements from BL to uniform spacing

The bottom boundary C1 is modified in-place: for nodes inside [x_ws, x_wp]
the y-coordinate ramps linearly up to w_h; for nodes in [x_wp, x_we] it
ramps back down to zero (piecewise-linear wedge). Outside this range the
bottom is flat (y=0).

Outputs (provided by Tianyi):
  Geometry.dat — PetIGA binary mesh file
  Geometry.vtk — VTK structured grid for visualization
  Geometry.mat — MATLAB .mat file with control points and knot vectors

Usage:
  python3 MeshGenerate_Wedge.py
"""

import numpy as np
import math
import copy
from scipy.io import savemat

from numpy import linspace
from igakit.cad import *
from igakit.io import PetIGA, VTK
from igakit.plot import plt
from igakit.nurbs import NURBS

# ---------------------------------------------------------------------------
# Geometry parameters
# ---------------------------------------------------------------------------
p        = 2      # geometry polynomial degree
p_basis  = 2      # target basis polynomial degree after elevation
N        = 716 + p
Ny       = 72     # number of uniform y-refinement intervals
Ny_BL    = 20     # number of boundary layer (BL) elements near y=0
Scale    = 0.3    # BL thickness as a fraction of H
TransGrid = 60    # number of transition intervals from BL to uniform spacing

L    = 4.5        # domain length [m]
H    = 0.6        # domain height [m]
w_h  = 0.2        # wedge peak height [m]
x_ws = 0.2        # wedge start x [m]
x_wp = x_ws + 3. * w_h   # wedge peak x [m]
x_we = x_wp + 4.5 * w_h  # wedge end x [m]
print(x_we)

# ---------------------------------------------------------------------------
# Build bottom (C1) and top (C2) boundary control-point arrays.
# C1 starts flat at y=0; the wedge shape is applied in a loop below.
# C2 is flat at y=H throughout.
# ---------------------------------------------------------------------------
C1 = np.linspace((0, 0), (L, 0), N)
C2 = np.linspace((0, H), (L, H), N)

for i in range(N):
    if   C1[i, 0] >= x_ws and C1[i, 0] <= x_wp:
        # Rising ramp: y goes from 0 at x_ws to w_h at x_wp
        C1[i, 1] = (C1[i, 0] - x_ws) * w_h / (x_wp - x_ws)
    elif C1[i, 0] >= x_wp and C1[i, 0] <= x_we:
        # Falling ramp: y goes from w_h at x_wp to 0 at x_we
        C1[i, 1] = w_h - w_h / (x_we - x_wp) * (C1[i, 0] - x_wp)

# ---------------------------------------------------------------------------
# Stack the two boundary curves and build the NURBS surface
# ---------------------------------------------------------------------------
C  = np.stack((C1, C2))

# Open-uniform knot vector in x (N control points, degree p)
Ux = np.linspace(0, 1, N + 1 - p)
for i in range(p):
    Ux = np.append(0, Ux)
    Ux = np.append(Ux, 1)

Uy = [0, 0, 1, 1]   # linear (p=1) in y before elevation

srf = NURBS([Uy, Ux], C).transpose()
srf.elevate(1, p_basis - 1)     # elevate y to degree p_basis
srf.elevate(0, p_basis - p)     # elevate x to degree p_basis
print(srf.degree)

# ---------------------------------------------------------------------------
# Anisotropic y-refinement: three-zone strategy
#   Zone 1 (insert_y1): geometrically stretched BL elements near y=0
#   Zone 2 (insert_y2): transition from BL spacing to uniform spacing
#   Zone 3 (insert_y3): uniform spacing above the transition (currently unused)
# ---------------------------------------------------------------------------
insert_x = np.linspace(0, 1, 1000 + 1)[1:-1]   # dense x-insertions (unused)
insert_y  = np.linspace(0, 1, Ny + 1)[1:-1]     # uniform y-insertions (unused)

# Zone 1: geometrically stretched from 2.5e-4 to Scale/Ny over Ny_BL elements
grid_space  = np.linspace(2.5e-4, Scale / Ny, Ny_BL + 1)
insert_y1   = np.linspace(0, 0, Ny_BL + 1)[1:]
for i in range(Ny_BL):
    if i == 0:
        insert_y1[i] = grid_space[i]
    else:
        insert_y1[i] = insert_y1[i - 1] + grid_space[i]
print(insert_y1[-1])

# Zone 2: transition from BL spacing (Scale/Ny) to uniform (1/Ny)
grid_space2 = np.linspace(Scale / Ny, 1.0 / Ny, TransGrid + 1)
insert_y2 = [0]
for i in range(Ny + 100):
    gridIndex = i
    if i > TransGrid:
        gridIndex = -1   # hold at the final (uniform) spacing
    if i == 0:
        insert_y2[i] = insert_y1[-1] + grid_space2[gridIndex]
    else:
        if (insert_y2[i - 1] + grid_space2[gridIndex]) >= 1.0:
            break
        insert_y2.append(insert_y2[i - 1] + grid_space2[gridIndex])

# Remove any duplicate knot values that already exist in Ux
for idx, j in enumerate(Ux):
    index   = np.argwhere(insert_x == j)
    insert_x = np.delete(insert_x, index)

# Zone 3: uniform spacing from end of transition to y=1 (unused — kept for reference)
insert_y3 = np.linspace(insert_y2[-1], 1, Ny)[1:-1]

# Apply refinement in y only (x-refinement commented out for now)
srf.refine(1, insert_y1)
srf.refine(1, insert_y2)
print(srf.control.shape)

# ---------------------------------------------------------------------------
# Export
# ---------------------------------------------------------------------------
savemat("Geometry.mat", {'coefs': srf.points.transpose(), 'knots': srf.knots})
PetIGA().write("Geometry.dat", srf)

uniform = lambda U: linspace(U[0], U[-1], N * 1 * U[-1] + 1)
VTK().write("Geometry.vtk", srf)

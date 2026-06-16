"""
MeshGenerate.py — 2D curved-arc mesh with boundary layer refinement.

Generates a 2D NURBS patch swept along an arc of a given angle. The two
boundary curves (at r=0 and r=1) are discretized as sequences of N points
along the arc, forming a ruled surface that maps [0,1]^2 to the physical
domain. Polynomial degree is p=9 throughout (high-order IGA).

The mesh is refined non-uniformly:
  - Along the arc (direction 1): dense near the curved wall via two
    linspace passes (coarse interior, fine near r=1 boundary).
  - Radially (direction 0): two interior knots inserted for boundary
    layer resolution near the wall.

Outputs:
  Geometry.dat  — PetIGA binary mesh file
  Geometry.vtk  — VTK structured grid for visualization
  Geometry.mat  — MATLAB .mat file with control points and knot vectors

Usage:
  python3 MeshGenerate.py <angle_deg> <N_base>

  angle_deg : arc angle in degrees (e.g. 90 for a quarter-circle)
  N_base    : base number of points along the arc (actual N = N_base + p)
"""

import numpy as np
import glob, sys, os
import math
import copy
from scipy.io import savemat

from numpy import linspace
from igakit.cad import *
from igakit.io import PetIGA, VTK
from igakit.plot import plt
from igakit.nurbs import NURBS

# ---------------------------------------------------------------------------
# Parameters — read from command line
# ---------------------------------------------------------------------------
p     = 9                              # polynomial degree (high-order IGA)
Angle = float(sys.argv[1]) / 180.0    # arc angle, converted from degrees to [0,1] fraction of pi
N     = int(sys.argv[2]) + p          # total number of control points along arc

# Two radii: r=1 (outer/wall) and r=0 (inner/axis)
R = [1.0, 0.0]

# ---------------------------------------------------------------------------
# Build control point arrays for both boundary curves
# ---------------------------------------------------------------------------
C = []
for j in range(np.size(R)):
    Position = []
    for i in range(N):
        # Parametrize arc: theta runs from 0 to -pi*Angle (clockwise)
        theta = -float(i) / (N - 1) * math.pi * Angle
        Position.append([R[j] * math.sin(theta), R[j] * math.cos(theta)])
    C.append(np.array(Position))
    if j == 0:
        Position1 = Position
    else:
        Position2 = Position

# ---------------------------------------------------------------------------
# Construct open-uniform knot vector in x (along arc) and trivial knot in y
# ---------------------------------------------------------------------------
# p+1 zeros, then N-p-1 interior knots, then p+1 ones
Ux = np.linspace(0, 1, N + 1 - p)[1:-1]
for i in range(p + 1):
    Ux = np.append(0, Ux)
    Ux = np.append(Ux, 1)

Uy = [0, 0, 1, 1]   # linear (p=1) in the radial direction

# ---------------------------------------------------------------------------
# Build NURBS surface and elevate to target degree p in both directions
# ---------------------------------------------------------------------------
srf = NURBS([Uy, Ux], C).transpose()
srf.elevate(0, p - srf.degree[0])   # elevate radial direction to degree p
srf.elevate(1, p - srf.degree[1])   # elevate arc direction to degree p
print(srf.degree)

# ---------------------------------------------------------------------------
# k-refinement: insert knots for boundary layer resolution near r=1 (wall)
# ---------------------------------------------------------------------------
scale = 3.0

# Two interior knots near each end of the radial direction for boundary layer
insert = [(Ux[p + 1]) / scale, 1.0 - (1.0 - Ux[-2 - p]) / scale]

# Along the arc: dense interior (coarse) + fine cluster near the wall end
insert_y  = np.linspace(0, 1.0 - 5e-2, N * Angle)[1:]      # coarse interior
insert_y2 = np.linspace(1.0 - 5e-2, 1.0, 512)[1:-1]         # fine near wall

srf.refine(0, insert)
srf.refine(1, insert_y)
srf.refine(1, insert_y2)

# ---------------------------------------------------------------------------
# Export
# ---------------------------------------------------------------------------
PetIGA().write("Geometry.dat", srf)
print(srf.shape)

savemat("Geometry.mat", {'coefs': srf.points.transpose(), 'knots': srf.knots})

uniform = lambda U: linspace(U[0], U[-1], N * 1 * U[-1] + 1)
VTK().write("Geometry.vtk", srf)

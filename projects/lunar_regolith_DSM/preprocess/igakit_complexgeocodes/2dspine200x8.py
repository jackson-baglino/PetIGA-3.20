"""
2dspine200x8.py — 2D B-spline surface for a symmetric fin/spine cross-section.

Identical to 2dspine100x8.py except the u-direction uses 200 knot insertions
instead of 100, giving a finer mesh along the spine for higher-resolution runs.

Defines a 6-control-point surface whose shape resembles a tapered fin:
  upper half: (0,5) -> (0,4) -> (45,0)
  lower half: (44,0) -> (0,-4) -> (0,-5)
The surface is built with geomdl (degree_u=2, degree_v=1), then imported
into igakit for k-refinement and export.

Mesh refinement:
  u-direction: 200 uniformly-spaced knot insertions (fine along the spine)
  v-direction:   8 uniformly-spaced knot insertions (coarse across the width)

Outputs:
  2dspline.dat  — PetIGA binary mesh file (overwrites 2dspine100x8 output)
  2dspline.vtk  — VTK structured grid for visualization

See also: 2dspine100x8.py (identical except 100 u-direction insertions).

Usage:
  python3 2dspine200x8.py
"""

import numpy as np
from numpy import linspace
import matplotlib.pyplot as plt
from matplotlib import cm
import os
from tqdm import tqdm
import GPy
import pandas as pd
import seaborn as sns
from geomdl import BSpline
from geomdl.visualization import VisMPL
from geomdl import utilities
from geomdl import operations

# ---------------------------------------------------------------------------
# Control point coordinates for the fin cross-section.
# Shape: tapered symmetric profile; upper (+y) and lower (-y) arms meet at
# the tip on the right side (~x=45).
# ---------------------------------------------------------------------------
x = [0.0, 0.0, 45.0, 44.0, 0.0, 0.0]
y = [5.0, 4.0,  0.0,  0.0, -5.0, -4.0]
z = [0.0, 0.0,  0.0,  0.0,  0.0,  0.0]

x = np.ravel(x)
y = np.ravel(y)
z = np.ravel(z)

# Pack into (3 rows) x (2 columns) x (3 coords) array, then convert to list
xyz = []
for i in range(len(x)):
    xyz.append([x[i], y[i], z[i]])
xyz = np.ravel(xyz)
xyz = xyz.reshape(3, 2, 3)
ctrlpts = xyz.tolist()

# ---------------------------------------------------------------------------
# Build geomdl B-spline surface
# ---------------------------------------------------------------------------
surf = BSpline.Surface()
surf.degree_u = 2   # quadratic along the spine
surf.degree_v = 1   # linear across the width
surf.ctrlpts2d = ctrlpts

# Auto-generate open-uniform knot vectors from the control point grid size
surf.knotvector_u = utilities.generate_knot_vector(surf.degree_u, surf.ctrlpts_size_u)
surf.knotvector_v = utilities.generate_knot_vector(surf.degree_v, surf.ctrlpts_size_v)

# ---------------------------------------------------------------------------
# Visualize with geomdl (renders an interactive matplotlib window)
# ---------------------------------------------------------------------------
surf.vis = VisMPL.VisSurface(ctrlpts=False, axes_equal=True)
surf.render()

# ---------------------------------------------------------------------------
# Convert to igakit NURBS and apply k-refinement
# ---------------------------------------------------------------------------
C = np.array(ctrlpts)

from igakit.cad import *
from igakit.io import VTK

igasurf = NURBS([surf.knotvector_u, surf.knotvector_v], C)

# Insert knots uniformly: 200 in u (along spine), 8 in v (across width).
# Interior points only — linspace endpoints [0] and [-1] are excluded.
insert_x0 = linspace(0, 1, 200)[1:-1]
insert_y0 = linspace(0, 1, 8)[1:-1]
igasurf.refine(0, insert_x0)
igasurf.refine(1, insert_y0)

# ---------------------------------------------------------------------------
# Export
# ---------------------------------------------------------------------------
from igakit.io import PetIGA

PetIGA().write("2dspline.dat", igasurf)
VTK().write("2dspline.vtk", igasurf)

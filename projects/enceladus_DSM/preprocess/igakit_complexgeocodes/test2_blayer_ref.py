"""
test2_blayer_ref.py — 2D rectangular mesh with boundary layer refinement.

Generates a 2D NURBS patch for a rectangular domain [0,Lx] x [0,Ly] split
into a fluid region and a solid region separated at y = fluid_thickness.
The mesh uses non-uniform refinement in the y-direction with three zones:

  Zone 1 (boundary layer): dense cluster of knots near y=0, starting at
      y = 0.015 * fluid_thickness and using 3 finely-spaced intervals.
  Zone 2 (fluid interior): uniform spacing from the BL edge to y = fluid_thickness.
  Zone 3 (solid): uniform spacing from y = fluid_thickness to y = 1.

The x-direction uses uniform refinement with elements_x intervals.

Note: avoid using kplot / cplot for visualization — the mesh is too dense
and the plot becomes congested.

Outputs:
  test.dat — PetIGA binary mesh file
  test.vtk — VTK structured grid for visualization

Usage:
  python3 test2_blayer_ref.py <Lx> <Ly> <solid_thickness> <elements_x>

  Lx              : domain length in x
  Ly              : domain height in y
  solid_thickness : non-dimensional solid layer thickness (fraction of Ly)
  elements_x      : number of elements in x (determines elements_y by aspect ratio)

Example:
  python3 test2_blayer_ref.py 4.0 1.0 0.1 100

Note on scaling: solid_thickness is non-dimensional. The split between fluid
and solid occurs at y = fluid_thickness = Ly - solid_thickness * Ly (i.e. the
solid occupies the top solid_thickness fraction of the domain).
"""

from igakit.io import PetIGA, VTK
from igakit.plot import plt
from numpy import *
import matplotlib.pyplot as plt2
import glob
from igakit.nurbs import NURBS
from igakit.transform import transform
from igakit.cad import *
import sys

# ---------------------------------------------------------------------------
# Parameters — read from command line
# ---------------------------------------------------------------------------
Lx              = float(sys.argv[1])
Ly              = float(sys.argv[2])
solid_thickness = float(sys.argv[3])   # non-dimensional (fraction of Ly)
elements_x      = int(sys.argv[4])

# y-elements scaled to match the x aspect ratio
elements_y = int(float(Ly / Lx) * elements_x)

# ---------------------------------------------------------------------------
# Build the rectangular patch by ruling two horizontal lines
# ---------------------------------------------------------------------------
l1  = line((0., 0.), (Lx, 0.))    # bottom edge
l2  = line((0., Ly), (Lx, Ly))    # top edge
nrb = ruled(l1, l2)

# Elevate to degree 2 in both directions (quadratic basis)
nrb.elevate(1, 1)
nrb.elevate(0, 1)

# ---------------------------------------------------------------------------
# Compute fluid/solid split in parametric space
# ---------------------------------------------------------------------------
# Solid occupies [fluid_thickness, 1.0] in normalised y.
fluid_thickness = 1. - float(solid_thickness)
solid_elements  = int(elements_y * solid_thickness)
fluid_elements  = elements_y - solid_elements + 1
print(fluid_thickness, solid_elements)

# ---------------------------------------------------------------------------
# Three-zone y-refinement
# ---------------------------------------------------------------------------
# Zone 1: boundary layer — 3 elements tightly clustered near y=0, up to
#         y = 0.015 * fluid_thickness in parametric space.
insert_blayer = linspace(0., 0.015 * fluid_thickness, 3 + 1, endpoint=True)[1:]
print(insert_blayer)
nrb.refine(1, insert_blayer)

# Zone 2: fluid interior — uniform from BL edge to fluid-solid interface.
#         fluid_elements - 3 knots (3 already placed by BL zone above).
insert_fluid = concatenate((
    linspace(0.015 * fluid_thickness, fluid_thickness, fluid_elements - 3, endpoint=True)[1:],
    [fluid_thickness]
))
nrb.refine(1, insert_fluid)
print(insert_fluid)

# Zone 3: solid — uniform from fluid-solid interface to y=1.
insert_solid = linspace(fluid_thickness, 1.0, solid_elements + 1)[1:-1]
nrb.refine(1, insert_solid)

# ---------------------------------------------------------------------------
# Uniform x-refinement
# ---------------------------------------------------------------------------
insert_y = linspace(0., 1.0, elements_x + 1)[1:-1]
nrb.refine(0, insert_y)
print(nrb.degree)

# ---------------------------------------------------------------------------
# Export
# ---------------------------------------------------------------------------
PetIGA().write("test.dat", nrb)
VTK().write("test.vtk", nrb)

# Diagnostic output: control array size and knot vector lengths
print(size(nrb.control))
print(size(nrb.knots[0]), size(nrb.knots[1]))
print(nrb.knots[0])
print(nrb.knots[1])

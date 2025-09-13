from dolfinx import mesh, fem
import ufl
import numpy as np
from mpi4py import MPI
from petsc4py import PETSc
from dolfinx.fem.petsc import LinearProblem

# Unit square mesh
domain = mesh.create_unit_square(MPI.COMM_WORLD, 16, 16, mesh.CellType.triangle)

# API has changed: FunctionSpace is created from domain using fem.functionspace
V = fem.functionspace(domain, ("Lagrange", 1))

u = ufl.TrialFunction(V)
v = ufl.TestFunction(V)
f = fem.Constant(domain, PETSc.ScalarType(1.0))
a = ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx
L = f * v * ufl.dx

# Dirichlet u=0 on boundary
facets = mesh.locate_entities_boundary(domain, domain.topology.dim - 1,
                                       lambda x: np.full(x.shape[1], True))
bc_dofs = fem.locate_dofs_topological(V, domain.topology.dim - 1, facets)
zero = fem.Constant(domain, PETSc.ScalarType(0.0))
bc = fem.dirichletbc(zero, bc_dofs, V)

uh = fem.Function(V)
problem = LinearProblem(a, L, bcs=[bc], u=uh)
uh = problem.solve()

# Compute L2 norm safely via assembly (works across dolfinx versions)
uh_L2 = np.sqrt(fem.assemble_scalar(fem.form(uh * uh * ufl.dx)))
if MPI.COMM_WORLD.rank == 0:
    print("Solved Poisson on 1x1 with u=0 boundary; ||u||_L2 =", uh_L2)
#include "assembly.h"
#include "material.h"

/*-----------------------------------------------------------------------------
  AssembleStiffnessMatrix
  IGA point-wise integrand for the periodic homogenization cell problem.

  For each direction m the correction temperature t_m satisfies:
    ∫ k(x) ∇N_i · ∇t_m dV  =  -∫ k(x) ∂N_i/∂x_m dV    ∀ test functions N_i

  Assembled as a vector-valued system with dim DOFs per node:
    K[i*dim+m, j*dim+m] += k * ∑_g  ∂N_i/∂x_g * ∂N_j/∂x_g
    F[i*dim+m]           += -k * ∂N_i/∂x_m

  NOTE ON QUADRATURE WEIGHTING
  ─────────────────────────────
  PetIGA calls this function once per Gauss point and then accumulates the
  returned K and F via IGAPointAddMat / IGAPointAddVec, which internally
  multiply by JW = detJac * weight.  The user function must therefore return
  bare (un-weighted) integrand values — do NOT multiply by dV here.
-----------------------------------------------------------------------------*/
PetscErrorCode AssembleStiffnessMatrix(IGAPoint pnt,
                                       PetscScalar *K,
                                       PetscScalar *F,
                                       void *ctx)
{
  PetscErrorCode        ierr;
  AppCtx               *user = (AppCtx *)ctx;
  PetscInt              nen  = pnt->nen, dim = user->dim;
  PetscInt              ndof = nen * dim;
  const PetscScalar    *N0;
  const PetscReal     (*N1)[dim];
  PetscInt              i, j, m, g;

  PetscFunctionBegin;

  ierr = IGAPointGetShapeFuns(pnt, 0, &N0);                    CHKERRQ(ierr);
  ierr = IGAPointGetShapeFuns(pnt, 1, (const PetscReal**)&N1); CHKERRQ(ierr);

  PetscInt  indGP = pnt->index + pnt->count * pnt->parent->index;
  PetscReal thcond;
  ThermalCond(user, user->ice[indGP], &thcond, NULL);

  for (PetscInt q = 0; q < ndof * ndof; ++q) K[q] = 0.0;
  for (PetscInt q = 0; q < ndof;        ++q) F[q] = 0.0;

  for (i = 0; i < nen; i++) {
    for (j = 0; j < nen; j++) {
      for (m = 0; m < dim; m++) {
        PetscScalar stab = 0.0;
        for (g = 0; g < dim; g++)
          stab += N1[i][g] * N1[j][g] * thcond;
        K[(i*dim + m) * ndof + (j*dim + m)] += stab;
      }
    }
    for (m = 0; m < dim; m++)
      F[i*dim + m] += -thcond * N1[i][m];
  }

  PetscFunctionReturn(0);
}

/*-----------------------------------------------------------------------------
  ApplyBoundaryConditions
  No-op for fully periodic problems.  Retained for a self-documenting call site.
-----------------------------------------------------------------------------*/
PetscErrorCode ApplyBoundaryConditions(IGA iga, AppCtx *user)
{
  PetscFunctionBegin;
  /* All axes are set periodic in SetupIGA — nothing to do here. */
  (void)iga; (void)user;
  PetscFunctionReturn(0);
}

/*-----------------------------------------------------------------------------
  ComputeKeffIntegrand
  IGA scalar integrand that accumulates the effective thermal conductivity
  tensor entry-by-entry:
    S[i*dim + j] += k(x) * (∂t_i/∂x_j + δ_{ij})

  IGAComputeScalar calls this function once per Gauss point and accumulates
  the result via IGAPointAddArray, which multiplies by JW = detJac * weight.
  The user function returns the bare integrand value — no dV multiplication.

  Here t_i is component i of the correction-temperature vector T_sol.
-----------------------------------------------------------------------------*/
PetscErrorCode ComputeKeffIntegrand(IGAPoint p, const PetscScalar *U,
                                    PetscInt n, PetscScalar *S, void *ctx)
{
  AppCtx   *user = (AppCtx *)ctx;
  PetscInt  dim  = user->dim;

  PetscFunctionBegin;

  if (!p->atboundary) {
    /* Flat buffer: grad_t[i*dim+j] = ∂T_i/∂x_j
     * PetIGA fills this as grad[dof*dim + dir], so row-stride = dim.
     * A fixed-size [3][3] array would have stride 3 and corrupt 2-D results. */
    PetscReal grad_t[9] = {0.0};
    IGAPointFormGrad(p, U, &grad_t[0]);

    PetscInt  indGP = p->index + p->count * p->parent->index;
    PetscReal thcond;
    ThermalCond(user, user->ice[indGP], &thcond, NULL);

    /* No dV here — IGAPointAddArray (called by IGAComputeScalar internally)
     * multiplies by JW = detJac * weight automatically. */
    for (PetscInt i = 0; i < dim; i++)
      for (PetscInt j = 0; j < dim; j++)
        S[i * dim + j] += thcond * (grad_t[i * dim + j] + (i == j ? 1.0 : 0.0));
  }

  PetscFunctionReturn(0);
}

/*-----------------------------------------------------------------------------
  ComputeKeffective
  Integrate ComputeKeffIntegrand over the domain, then divide by volume to
  obtain the dim×dim effective thermal conductivity tensor (row-major).
-----------------------------------------------------------------------------*/
PetscErrorCode ComputeKeffective(IGA iga, Vec t_vec,
                                 PetscReal *keff, AppCtx *user)
{
  PetscErrorCode ierr;
  PetscInt       dim = user->dim;
  PetscReal      volume;

  PetscFunctionBegin;

  ierr = PetscMemzero(keff, sizeof(PetscReal) * dim * dim); CHKERRQ(ierr);
  ierr = IGAComputeScalar(iga, t_vec, dim * dim, keff,
                          ComputeKeffIntegrand, (void *)user); CHKERRQ(ierr);

  volume = (dim == 2) ? user->Lx * user->Ly
                      : user->Lx * user->Ly * user->Lz;

  for (PetscInt k = 0; k < dim * dim; k++)
    keff[k] /= volume;

  PetscFunctionReturn(0);
}

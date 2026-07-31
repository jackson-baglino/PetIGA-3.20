#include "keff.h"
#include "material_properties.h"

/* ---------------------------------------------------------------------------
 * keff_cell.c -- the periodic homogenization cell problem.
 *
 * For each macroscopic direction m we solve for the corrector t_m on the cell Y
 *
 *     -div( k(x) grad t_m )  =  div( k(x) e_m )      t_m  Y-periodic
 *
 * whose weak form is, for all periodic test functions v,
 *
 *     INT_Y k grad v . grad t_m  =  - INT_Y k  dv/dx_m
 *
 * and then average the microscopic flux to get the effective tensor
 *
 *     k_eff[i][j] = (1/|Y|) INT_Y k(x) ( d t_i/d x_j + delta_ij ) dV.
 *
 * ONE OPERATOR, dim RIGHT-HAND SIDES. The bilinear form does not reference m at
 * all -- only the load does. The original blocked implementation built a
 * dim-dof system and copied the SAME scalar Laplacian into each of its dim
 * diagonal blocks, so three quarters of the matrix at dim=2 was structural
 * zero. Here the operator is assembled once as a scalar matrix and reused for
 * every direction, which at Nx=Ny=2829 is the difference between ~9.6 GB and
 * ~2.4 GB, and hands the preconditioner a plain scalar Laplacian.
 *
 * PetIGA applies the quadrature weight and Jacobian (JW = detJac*weight) to
 * whatever the form callbacks accumulate, so no dV factor appears below.
 * ------------------------------------------------------------------------- */

/* Thermal conductivity at the current quadrature point, from the projected
 * phase field. The flat index convention is documented in keff_field.c. */
static inline PetscReal KeffPointCond(KeffCtx *kc, IGAPoint point)
{
  PetscInt    idx = point->index + point->count * point->parent->index;
  PetscScalar cond;
  ThermalCond(kc->app, (PetscScalar)kc->ice[idx], &cond, NULL);
  return PetscRealPart(cond);
}

/* ---------------------------------------------------------------------------
 * KeffFormMatrix:  K[i][j] += k(x) * grad N_i . grad N_j
 * ------------------------------------------------------------------------- */
PetscErrorCode KeffFormMatrix(IGAPoint point, PetscScalar K[], void *ctx)
{
  KeffCtx        *kc  = (KeffCtx *)ctx;
  const PetscInt  nen = point->nen;
  const PetscInt  dim = point->dim;
  PetscReal     (*N1)[dim];
  PetscReal       cond;
  PetscErrorCode  ierr;

  PetscFunctionBegin;
  if (point->atboundary) PetscFunctionReturn(0);

  ierr = IGAPointGetShapeFuns(point, 1, (const PetscReal **)&N1); CHKERRQ(ierr);
  cond = KeffPointCond(kc, point);

  for (PetscInt i = 0; i < nen; i++)
    for (PetscInt j = 0; j < nen; j++) {
      PetscReal g = 0.0;
      for (PetscInt d = 0; d < dim; d++) g += N1[i][d] * N1[j][d];
      K[i * nen + j] += cond * g;
    }

  PetscFunctionReturn(0);
}

/* ---------------------------------------------------------------------------
 * KeffFormVector:  F[i] += -k(x) * dN_i/dx_m   for the current direction m
 * ------------------------------------------------------------------------- */
PetscErrorCode KeffFormVector(IGAPoint point, PetscScalar F[], void *ctx)
{
  KeffCtx        *kc  = (KeffCtx *)ctx;
  const PetscInt  nen = point->nen;
  const PetscInt  dim = point->dim;
  PetscReal     (*N1)[dim];
  PetscReal       cond;
  PetscErrorCode  ierr;

  PetscFunctionBegin;
  if (point->atboundary) PetscFunctionReturn(0);

  ierr = IGAPointGetShapeFuns(point, 1, (const PetscReal **)&N1); CHKERRQ(ierr);
  cond = KeffPointCond(kc, point);

  for (PetscInt i = 0; i < nen; i++)
    F[i] += -cond * N1[i][kc->cur_dir];

  PetscFunctionReturn(0);
}

/* ---------------------------------------------------------------------------
 * KeffScalarIntegrand -- one row of the tensor, for direction i = kc->cur_dir:
 *
 *     S[j]   += k * ( d t_i/d x_j + delta_ij )      j = 0 .. dim-1
 *     S[dim] += phi                                 (cell-mean ice, free)
 *
 * U is the corrector t_i on the scalar corrector mesh, so IGAPointFormGrad
 * yields the dim components of grad t_i directly.
 * ------------------------------------------------------------------------- */
PetscErrorCode KeffScalarIntegrand(IGAPoint point, const PetscScalar U[],
                                   PetscInt n, PetscScalar S[], void *ctx)
{
  KeffCtx        *kc  = (KeffCtx *)ctx;
  const PetscInt  dim = point->dim;
  PetscScalar     grad_t[3];
  PetscReal       cond;
  PetscInt        idx;

  PetscFunctionBegin;
  (void)n;
  if (point->atboundary) PetscFunctionReturn(0);

  IGAPointFormGrad(point, U, &grad_t[0]);

  idx  = point->index + point->count * point->parent->index;
  cond = KeffPointCond(kc, point);

  for (PetscInt j = 0; j < dim; j++)
    S[j] += cond * (PetscRealPart(grad_t[j]) + ((j == kc->cur_dir) ? 1.0 : 0.0));

  S[dim] += kc->ice[idx];

  PetscFunctionReturn(0);
}

/* ---------------------------------------------------------------------------
 * KeffComputeTensor
 *
 * Assemble the tensor from the dim solved correctors. keff is row-major,
 * keff[i*dim + j] = k_eff[i][j]. phi_bar is the cell-mean ice fraction, taken
 * from the first pass (every pass computes it; they agree).
 *
 * Divides by the MEASURED cell volume from KeffCreate, not by Lx*Ly[*Lz].
 * ------------------------------------------------------------------------- */
PetscErrorCode KeffComputeTensor(AppCtx *app, PetscReal keff[9], PetscReal *phi_bar)
{
  PetscErrorCode ierr;
  KeffCtx       *kc  = app->keff;
  const PetscInt dim = kc->dim;

  PetscFunctionBegin;

  for (PetscInt i = 0; i < dim; i++) {
    PetscScalar S[4] = {0.0, 0.0, 0.0, 0.0};

    kc->cur_dir = i;
    ierr = IGAComputeScalar(kc->iga, kc->T[i], dim + 1, &S[0],
                            KeffScalarIntegrand, kc); CHKERRQ(ierr);

    for (PetscInt j = 0; j < dim; j++)
      keff[i * dim + j] = PetscRealPart(S[j]) / kc->vol;

    if (i == 0 && phi_bar) *phi_bar = PetscRealPart(S[dim]) / kc->vol;
  }

  PetscFunctionReturn(0);
}

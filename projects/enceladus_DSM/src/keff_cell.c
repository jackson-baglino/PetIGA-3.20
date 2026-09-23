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
 * THE COEFFICIENT IS A TENSOR. Under -keff_interp tensor, K is anisotropic
 * inside the band, so every "k * grad" below is a matrix-vector product. The
 * isotropic laws fill K = k I and give exactly the scalar assembly. The weak
 * form, the load and the flux average carry over unchanged with k -> K:
 *
 *     INT_Y grad v . K grad t_m  =  - INT_Y grad v . K e_m
 *     k_eff[i][j] = (1/|Y|) INT_Y [ K ( grad t_i + e_i ) ]_j dV
 *
 * K is symmetric positive definite (its eigenvalues are k_arith and k_harm,
 * both in [k_a, k_i]), so the operator stays SPD and CG + GAMG still apply.
 *
 * PetIGA applies the quadrature weight and Jacobian (JW = detJac*weight) to
 * whatever the form callbacks accumulate, so no dV factor appears below.
 * ------------------------------------------------------------------------- */

/* Conductivity tensor at the current quadrature point, from the projected phase
 * field, row-major K[a*dim + b]. The flat index convention is documented in
 * keff_field.c.
 *
 * TENSOR: K = k_arith (I - n n) + k_harm n n, n = grad phi / |grad phi|. Where
 * grad phi vanishes n is undefined; K falls back to k_arith I. That happens in
 * the bulk, where phi is 0 or 1 and k_arith = k_harm anyway, so the choice is
 * immaterial there. phi is clamped to [0,1] first so that both branches, and
 * k_harm's denominator, see the same admissible value. */
static inline void KeffPointCond(KeffCtx *kc, IGAPoint point, PetscReal K[9])
{
  const PetscInt dim = kc->dim;
  PetscInt       idx = point->index + point->count * point->parent->index;
  PetscReal      phi = kc->ice[idx];
  PetscScalar    cond;

  for (PetscInt a = 0; a < dim * dim; a++) K[a] = 0.0;

  if (kc->interp == KEFF_INTERP_SHARP) phi = (phi >= 0.5) ? 1.0 : 0.0;
  if (kc->interp == KEFF_INTERP_TENSOR) phi = PetscMin(PetscMax(phi, 0.0), 1.0);

  ThermalCond(kc->app, (PetscScalar)phi, &cond, NULL);
  for (PetscInt a = 0; a < dim; a++) K[a * dim + a] = PetscRealPart(cond);

  if (kc->interp == KEFF_INTERP_TENSOR) {
    const PetscReal *g      = &kc->grad_ice[idx * dim];
    const PetscReal  k_arith = PetscRealPart(cond);
    const PetscReal  k_harm  = 1.0 / (phi / kc->app->thcond_ice
                                    + (1.0 - phi) / kc->app->thcond_air);
    PetscReal        g2 = 0.0;

    for (PetscInt a = 0; a < dim; a++) g2 += g[a] * g[a];
    if (g2 > 0.0) {
      for (PetscInt a = 0; a < dim; a++)
        for (PetscInt b = 0; b < dim; b++)
          K[a * dim + b] += (k_harm - k_arith) * g[a] * g[b] / g2;
    }
  }
}

/* ---------------------------------------------------------------------------
 * KeffFormMatrix:  K[i][j] += grad N_i . K(x) grad N_j
 * ------------------------------------------------------------------------- */
PetscErrorCode KeffFormMatrix(IGAPoint point, PetscScalar K[], void *ctx)
{
  KeffCtx        *kc  = (KeffCtx *)ctx;
  const PetscInt  nen = point->nen;
  const PetscInt  dim = point->dim;
  PetscReal     (*N1)[dim];
  PetscReal       Kc[9];
  PetscErrorCode  ierr;

  PetscFunctionBegin;
  if (point->atboundary) PetscFunctionReturn(0);

  ierr = IGAPointGetShapeFuns(point, 1, (const PetscReal **)&N1); CHKERRQ(ierr);
  KeffPointCond(kc, point, Kc);

  for (PetscInt j = 0; j < nen; j++) {
    PetscReal KgradNj[3] = {0.0, 0.0, 0.0};
    for (PetscInt a = 0; a < dim; a++)
      for (PetscInt b = 0; b < dim; b++) KgradNj[a] += Kc[a * dim + b] * N1[j][b];
    for (PetscInt i = 0; i < nen; i++) {
      PetscReal g = 0.0;
      for (PetscInt d = 0; d < dim; d++) g += N1[i][d] * KgradNj[d];
      K[i * nen + j] += g;
    }
  }

  PetscFunctionReturn(0);
}

/* ---------------------------------------------------------------------------
 * KeffFormVector:  F[i] += -grad N_i . K(x) e_m   for the current direction m
 * ------------------------------------------------------------------------- */
PetscErrorCode KeffFormVector(IGAPoint point, PetscScalar F[], void *ctx)
{
  KeffCtx        *kc  = (KeffCtx *)ctx;
  const PetscInt  nen = point->nen;
  const PetscInt  dim = point->dim;
  const PetscInt  m   = kc->cur_dir;
  PetscReal     (*N1)[dim];
  PetscReal       Kc[9];
  PetscErrorCode  ierr;

  PetscFunctionBegin;
  if (point->atboundary) PetscFunctionReturn(0);

  ierr = IGAPointGetShapeFuns(point, 1, (const PetscReal **)&N1); CHKERRQ(ierr);
  KeffPointCond(kc, point, Kc);

  /* K e_m is column m of K. */
  for (PetscInt i = 0; i < nen; i++) {
    PetscReal g = 0.0;
    for (PetscInt d = 0; d < dim; d++) g += N1[i][d] * Kc[d * dim + m];
    F[i] += -g;
  }

  PetscFunctionReturn(0);
}

/* ---------------------------------------------------------------------------
 * KeffScalarIntegrand -- one row of the tensor, for direction i = kc->cur_dir:
 *
 *     S[j]   += [ K ( grad t_i + e_i ) ]_j           j = 0 .. dim-1
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
  PetscReal       Kc[9], G[3];
  PetscInt        idx;

  PetscFunctionBegin;
  (void)n;
  if (point->atboundary) PetscFunctionReturn(0);

  IGAPointFormGrad(point, U, &grad_t[0]);

  idx = point->index + point->count * point->parent->index;
  KeffPointCond(kc, point, Kc);

  for (PetscInt b = 0; b < dim; b++)
    G[b] = PetscRealPart(grad_t[b]) + ((b == kc->cur_dir) ? 1.0 : 0.0);
  for (PetscInt j = 0; j < dim; j++)
    for (PetscInt b = 0; b < dim; b++) S[j] += Kc[j * dim + b] * G[b];

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

#ifndef MATERIAL_PROPERTIES_H
#define MATERIAL_PROPERTIES_H

#include "NASA_types.h"

/* Computes effective thermal conductivity and its derivative with respect to ice */
void ThermalCond(AppCtx *user, PetscScalar ice, PetscScalar *cond, PetscScalar *dcond_ice);

/* Computes effective heat capacity and its derivative with respect to ice */
void HeatCap(AppCtx *user, PetscScalar ice, PetscScalar *cp, PetscScalar *dcp_ice);

/* Computes effective density and its derivative with respect to ice */
void Density(AppCtx *user, PetscScalar ice, PetscScalar *rho, PetscScalar *drho_ice);

/* Computes vapor diffusivity and its temperature derivative */
void VaporDiffus(AppCtx *user, PetscScalar tem, PetscScalar *difvap, PetscScalar *d_difvap);

/* Computes the saturation vapor density and its derivative */
void RhoVS_I(AppCtx *user, PetscScalar tem, PetscScalar *rho_vs, PetscScalar *d_rhovs);

/* computes the phase dependent mobility*/
void Mobility(AppCtx *user, PetscScalar ice, PetscScalar *mob);

/* Computes sigma0 using logarithmic interpolation based on temperature */
void Sigma0(PetscScalar temp, PetscScalar *sigm0);

/* Curvature of an isosurface of phi, computed from its gradient and Hessian:
 *   kappa = -div(grad_phi / |grad_phi|)
 *         = -L/G + (g.H.g)/G^3
 * where g = grad_phi, H = Hessian, L = trace(H) (Laplacian), G = |g|_reg.
 * G is regularized: G^2 = |g|^2 + eps_reg^2, to keep kappa finite where
 * |g| -> 0 (bulk regions). The ice^2*air^2 localizer in S_sub kills the
 * (numerically uninteresting) bulk contribution.
 *
 * dim must match the spatial dimension; dim==1 returns kappa=0 (no curvature
 * in 1D, but the regularized formula would otherwise produce a small artifact).
 *
 * Pass NULL for any of kappa / dkappa_dg / dkappa_dH that you don't need.
 *   dkappa_dg[l]              = d kappa / d g_l          (length dim)
 *   dkappa_dH[k*dim + l]      = d kappa / d H_{kl}       (length dim*dim)
 */
void Curvature(PetscInt dim,
               const PetscScalar grad_phi[],
               const PetscScalar hess_phi[],
               PetscReal eps_reg,
               PetscScalar *kappa,
               PetscScalar dkappa_dg[],
               PetscScalar dkappa_dH[]);

#endif // MATERIAL_PROPERTIES_H

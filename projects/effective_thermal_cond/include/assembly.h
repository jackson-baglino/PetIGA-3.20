#ifndef ASSEMBLY_H
#define ASSEMBLY_H

#include "app_ctx.h"

/*
 * AssembleStiffnessMatrix — IGA point-wise integrand for the periodic
 * homogenization cell problem:
 *   K[row,col] += k(x) * ∇N_i · ∇N_j   (per component)
 *   F[row]     += -k(x) * ∂N_i/∂x_m     (unit-gradient forcing)
 */
PetscErrorCode AssembleStiffnessMatrix(IGAPoint pnt,
                                       PetscScalar *K,
                                       PetscScalar *F,
                                       void *ctx);

/*
 * ApplyBoundaryConditions — no-op for fully-periodic problems.
 * Retained so the call-site in main() is self-documenting.
 */
PetscErrorCode ApplyBoundaryConditions(IGA iga, AppCtx *user);

/*
 * ComputeKeffIntegrand — IGA scalar integrand that accumulates the
 * effective thermal conductivity tensor:
 *   S[i,j] += k(x) * (∂t_i/∂x_j + δ_{ij}) * dV
 */
PetscErrorCode ComputeKeffIntegrand(IGAPoint p, const PetscScalar *U,
                                    PetscInt n, PetscScalar *S, void *ctx);

/*
 * ComputeKeffective — integrate ComputeKeffIntegrand over the domain and
 * divide by volume to get the dim×dim keff tensor (row-major).
 */
PetscErrorCode ComputeKeffective(IGA iga, Vec t_vec,
                                 PetscReal *keff, AppCtx *user);

#endif /* ASSEMBLY_H */

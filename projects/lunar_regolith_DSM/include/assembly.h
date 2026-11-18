#ifndef ASSEMBLY_H
#define ASSEMBLY_H

#include "NASA_types.h"

/* Residual implementation for the 2-phase (ice/temperature/vapor) system */
PetscErrorCode Residual_A1(IGAPoint pnt, PetscReal shift, const PetscScalar *V,
                            PetscReal t, const PetscScalar *U, PetscScalar *Re,
                            void *ctx);

/* The derived model (src/assembly_derived.c), selected with -residual_model derived */
PetscErrorCode Residual_D1(IGAPoint pnt, PetscReal shift, const PetscScalar *V,
                           PetscReal t, const PetscScalar *U, PetscScalar *Re,
                           void *ctx);
PetscErrorCode Jacobian_D1(IGAPoint pnt, PetscReal shift, const PetscScalar *V,
                           PetscReal t, const PetscScalar *U, PetscScalar *Je,
                           void *ctx);

/* Dispatchers: Residual_A1 / Jacobian_A1, or the D1 pair when user->residual_model == 1 */
PetscErrorCode Residual(IGAPoint pnt, PetscReal shift, const PetscScalar *V,
                        PetscReal t, const PetscScalar *U, PetscScalar *Re,
                        void *ctx);

/* Jacobian evaluation */
PetscErrorCode Jacobian(IGAPoint pnt, PetscReal shift, const PetscScalar *V,
                        PetscReal t, const PetscScalar *U, PetscScalar *Je,
                        void *ctx);

/* Computes integrated scalar quantities over the domain */
PetscErrorCode Integration(IGAPoint pnt, const PetscScalar *U, PetscInt n,
                           PetscScalar *S, void *ctx);

#endif // ASSEMBLY_H

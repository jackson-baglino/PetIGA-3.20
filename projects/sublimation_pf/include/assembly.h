#ifndef ASSEMBLY_H
#define ASSEMBLY_H

#include "NASA_types.h"

/* Residual implementation for the 2-phase (ice/temperature/vapor) system */
PetscErrorCode Residual_A1(IGAPoint pnt, PetscReal shift, const PetscScalar *V,
                            PetscReal t, const PetscScalar *U, PetscScalar *Re,
                            void *ctx);

/* Residual for the 3-phase (ice/temperature/vapor/sediment) system, dof=4
 * (defined in assembly_sed.c) */
PetscErrorCode Residual_A2(IGAPoint pnt, PetscReal shift, const PetscScalar *V,
                            PetscReal t, const PetscScalar *U, PetscScalar *Re,
                            void *ctx);

/* Dispatcher: forwards to Residual_A1 (dof=3) or Residual_A2 (dof=4) */
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

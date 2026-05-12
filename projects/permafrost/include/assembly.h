#ifndef ASSEMBLY_H
#define ASSEMBLY_H

#include "NASA_types.h"

/* Avenue 1: Allen-Cahn + penalty vapor; sediment RHS zeroed after t_sed_freeze */
PetscErrorCode Residual_A1(IGAPoint pnt, PetscReal shift, const PetscScalar *V,
                            PetscReal t, const PetscScalar *U, PetscScalar *Re,
                            void *ctx);

/* Avenue 2: Allen-Cahn + penalty vapor; sediment penalty after t_sed_freeze */
PetscErrorCode Residual_A2(IGAPoint pnt, PetscReal shift, const PetscScalar *V,
                            PetscReal t, const PetscScalar *U, PetscScalar *Re,
                            void *ctx);

/* Avenue 3: Cahn-Hilliard phase-field evolution; no penalty parameters */
PetscErrorCode Residual_A3(IGAPoint pnt, PetscReal shift, const PetscScalar *V,
                            PetscReal t, const PetscScalar *U, PetscScalar *Re,
                            void *ctx);

/* Dispatcher: calls Residual_A1/A2/A3 based on user->flag_avenue */
PetscErrorCode Residual(IGAPoint pnt, PetscReal shift, const PetscScalar *V,
                        PetscReal t, const PetscScalar *U, PetscScalar *Re,
                        void *ctx);

/* Jacobian evaluation (FD approximation used; stub returns 0) */
PetscErrorCode Jacobian(IGAPoint pnt, PetscReal shift, const PetscScalar *V,
                        PetscReal t, const PetscScalar *U, PetscScalar *Je,
                        void *ctx);

/* Computes integrated scalar quantities over the domain */
PetscErrorCode Integration(IGAPoint pnt, const PetscScalar *U, PetscInt n,
                           PetscScalar *S, void *ctx);

#endif // ASSEMBLY_H

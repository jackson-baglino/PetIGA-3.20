#ifndef INITIAL_CONDITIONS_H
#define INITIAL_CONDITIONS_H

#include "NASA_types.h"

/* Form initial soil distribution for 2D and 3D problems */
PetscErrorCode FormInitialLayeredPermafrost2D(IGA iga, Vec U, AppCtx *user);
PetscErrorCode FormInitialEnclosedPermafrost2D(IGA iga, Vec U, AppCtx *user);
PetscErrorCode FormInitialFlatSedIceCap2D(IGA iga, Vec U, AppCtx *user);
PetscErrorCode FormInitialRandomEnclosedPermafrost2D(IGA iga, Vec U, AppCtx *user);
PetscErrorCode FormInitialRandomPackedPermafrost2D(IGA iga, Vec U, AppCtx *user);
PetscErrorCode FormInitialSoil2D(IGA iga, Vec U, AppCtx *user);
PetscErrorCode FormInitialSoil3D(IGA iga, Vec U, AppCtx *user);

/* Form initial condition for the primary field variables in 2D and 3D */
PetscErrorCode FormInitialCondition2D(IGA iga, PetscReal t, Vec U, AppCtx *user,
                                      const char datafile[], const char dataPF[]);
PetscErrorCode FormLayeredInitialCondition2D(IGA iga, PetscReal t, Vec U,
                                            AppCtx *user, const char datafile[],
                                            const char dataPF[]);
PetscErrorCode FormInitialCondition3D(IGA iga, PetscReal t, Vec U, AppCtx *user,
                                      const char datafile[], const char dataPF[]);
PetscErrorCode LoadInputSolutionVec(const char *filename, Vec *U_in_seq);
PetscErrorCode InitializeFromInputSolution(IGA iga, Vec U, AppCtx *user);
PetscErrorCode FormIC_grain_ana(IGA iga, Vec U, AppCtx *user);

/* 1D initial conditions */
PetscErrorCode FormInitialCondition1D(IGA iga, Vec U, AppCtx *user);

#endif // INITIAL_CONDITIONS_H

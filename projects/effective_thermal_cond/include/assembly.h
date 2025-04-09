#ifndef ASSEMBLY_H
#define ASSEMBLY_H

#include "user_context.h"

// ========================= Function Declarations =========================
PetscErrorCode AssembleStiffnessMatrix(IGAPoint pnt, PetscScalar *K, PetscScalar *F, void *ctx);
PetscErrorCode ComputeInitialCondition(Vec T, AppCtx *user);
PetscErrorCode ApplyBoundaryConditions(IGA iga, AppCtx *user);

#endif // ASSEMBLY_H
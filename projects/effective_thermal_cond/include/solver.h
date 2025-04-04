#ifndef SOLVER_H
#define SOLVER_H

#include "user_context.h"

// ========================= Function Declarations =========================
PetscErrorCode AssembleStiffnessMatrix(IGAPoint pnt, PetscScalar *K, PetscScalar *F, void *ctx);
PetscErrorCode ComputeInitialCondition(Vec T, AppCtx *user);

#endif // SOLVER_H
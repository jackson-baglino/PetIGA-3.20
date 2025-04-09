#ifndef SETUP_THERMAL_H
#define SETUP_THERMAL_H

#include "user_context.h"

// ========================= Function Declarations =========================
PetscErrorCode FormInitialCondition(AppCtx *user);
PetscErrorCode InitializeFields(AppCtx *user, IGA iga);
PetscErrorCode SetupIGA(AppCtx *user, IGA *iga);
PetscErrorCode SetupAndSolve(AppCtx *user, IGA iga);

#endif // ENV_CONFIG_H
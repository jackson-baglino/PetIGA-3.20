#ifndef SOLVER_H
#define SOLVER_H

#include "user_context.h"
#include "assembly.h"

/* -------------------------- Function Declarations ------------------------- */
PetscErrorCode SetupAndSolve(AppCtx *user, IGA iga);

#endif // SOLVER_H

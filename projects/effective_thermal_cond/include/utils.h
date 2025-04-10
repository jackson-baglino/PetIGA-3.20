#ifndef UTILS_H
#define UTILS_H

#include "user_context.h"

PetscErrorCode AllocateAppCtxFields(IGA iga, AppCtx *user, PetscScalar **field);
PetscErrorCode FreeAppCtxFields(PetscScalar **field);

#endif // UTILS_H
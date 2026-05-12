// field_helpers.c
#include "utils.h"

/*--------------------------------------------------------------------------------------------------
  Function: AllocateAppCtxFields
  Allocates a PetscScalar field with size based on user dimensions and polynomial order.
  This function supports 2D and 3D cases and assumes tensor-product basis structure.
--------------------------------------------------------------------------------------------------*/
PetscErrorCode AllocateAppCtxFields(IGA iga, AppCtx *user, PetscScalar **field) {
  PetscFunctionBegin;

  PetscInt nmb = iga->elem_width[0] * iga->elem_width[1] * SQ(user->p + 1);
  if (user->dim == 3) {
    nmb = iga->elem_width[0] * iga->elem_width[1] * iga->elem_width[2] * CU(user->p + 1);
  }

  PetscCall(PetscMalloc1(nmb, field));
  PetscCall(PetscMemzero(*field, sizeof(PetscScalar) * nmb));

  PetscFunctionReturn(0);
}

/*--------------------------------------------------------------------------------------------------
  Function: FreeAppCtxFields
  Frees a PetscScalar field and sets its pointer to NULL to avoid dangling references.
--------------------------------------------------------------------------------------------------*/
PetscErrorCode FreeAppCtxFields(PetscScalar **field) {
  PetscFunctionBegin;

  if (*field) {
    PetscCall(PetscFree(*field));
    *field = NULL;
  }

  PetscFunctionReturn(0);
}
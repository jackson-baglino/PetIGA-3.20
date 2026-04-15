#include "field_init.h"
#include "material.h"

/*-----------------------------------------------------------------------------
  AllocateAppCtxFields
  Allocate a PetscScalar array sized for the LOCAL Gauss-point count on this
  rank.  Supports both 2-D and 3-D tensor-product element structures.
-----------------------------------------------------------------------------*/
PetscErrorCode AllocateAppCtxFields(IGA iga, AppCtx *user, PetscScalar **field)
{
  PetscFunctionBegin;

  PetscInt nmb = iga->elem_width[0] * iga->elem_width[1] * SQ(user->p + 1);
  if (user->dim == 3)
    nmb = iga->elem_width[0] * iga->elem_width[1] * iga->elem_width[2]
          * CU(user->p + 1);

  PetscCall(PetscMalloc1(nmb, field));
  PetscCall(PetscMemzero(*field, sizeof(PetscScalar) * nmb));

  PetscFunctionReturn(0);
}

/*-----------------------------------------------------------------------------
  EvaluateFieldAtGaussPoints
  Project the first DOF of vec_phase onto the local user->ice array by
  iterating over all local IGA elements and Gauss points.
-----------------------------------------------------------------------------*/
PetscErrorCode EvaluateFieldAtGaussPoints(AppCtx *user, IGA iga, Vec vec_phase)
{
  PetscErrorCode    ierr;
  IGAElement        element;
  IGAPoint          point;
  Vec               localU;
  const PetscScalar *arrayU;
  PetscScalar       *U;
  PetscInt          idx = 0;

  PetscFunctionBegin;

  ierr = IGAGetLocalVecArray(iga, vec_phase, &localU, &arrayU); CHKERRQ(ierr);
  ierr = IGABeginElement(iga, &element); CHKERRQ(ierr);
  while (IGANextElement(iga, element)) {
    ierr = IGAElementGetValues(element, arrayU, &U); CHKERRQ(ierr);
    ierr = IGAElementBeginPoint(element, &point); CHKERRQ(ierr);
    while (IGAElementNextPoint(element, point)) {
      PetscScalar sol[3];
      IGAPointFormValue(point, U, &sol[0]);
      user->ice[idx++] = sol[0];
    }
    ierr = IGAElementEndPoint(element, &point); CHKERRQ(ierr);
  }
  ierr = IGAEndElement(iga, &element); CHKERRQ(ierr);
  ierr = IGARestoreLocalVecArray(iga, vec_phase, &localU, &arrayU); CHKERRQ(ierr);

  /* Sanity check — local count must match allocation */
  PetscInt nmb_expected = iga->elem_width[0] * iga->elem_width[1] * SQ(user->p + 1);
  if (user->dim == 3)
    nmb_expected = iga->elem_width[0] * iga->elem_width[1] * iga->elem_width[2]
                   * CU(user->p + 1);
  if (idx != nmb_expected)
    PetscPrintf(PETSC_COMM_WORLD,
                "Warning: assigned %d Gauss-point values, expected %d\n",
                idx, nmb_expected);

  PetscFunctionReturn(0);
}

/*-----------------------------------------------------------------------------
  ReadSolutionVec
  Read an IGA descriptor and a solution vector from file, then project the
  first DOF onto user->ice via EvaluateFieldAtGaussPoints.
-----------------------------------------------------------------------------*/
PetscErrorCode ReadSolutionVec(const char *iga_file, const char *vec_file,
                               IGA *iga_out, AppCtx *user)
{
  PetscErrorCode ierr;
  IGA            iga_input;
  Vec            ice_phase;
  IGAAxis        axis0, axis1, axis2;

  PetscFunctionBegin;

  ierr = IGACreate(PETSC_COMM_WORLD, &iga_input); CHKERRQ(ierr);
  ierr = IGASetDim(iga_input, user->dim); CHKERRQ(ierr);
  ierr = IGASetDof(iga_input, 3); CHKERRQ(ierr);

  ierr = IGAGetAxis(iga_input, 0, &axis0); CHKERRQ(ierr);
  ierr = IGAAxisSetDegree(axis0, user->p); CHKERRQ(ierr);
  ierr = IGAAxisInitUniform(axis0, user->Nx, 0.0, user->Lx, user->C); CHKERRQ(ierr);

  ierr = IGAGetAxis(iga_input, 1, &axis1); CHKERRQ(ierr);
  ierr = IGAAxisSetDegree(axis1, user->p); CHKERRQ(ierr);
  ierr = IGAAxisInitUniform(axis1, user->Ny, 0.0, user->Ly, user->C); CHKERRQ(ierr);

  if (user->dim == 3) {
    ierr = IGAGetAxis(iga_input, 2, &axis2); CHKERRQ(ierr);
    ierr = IGAAxisSetDegree(axis2, user->p); CHKERRQ(ierr);
    ierr = IGAAxisInitUniform(axis2, user->Nz, 0.0, user->Lz, user->C); CHKERRQ(ierr);
  }

  ierr = IGASetFromOptions(iga_input); CHKERRQ(ierr);
  ierr = IGASetUp(iga_input); CHKERRQ(ierr);

  user->iga_input = iga_input;

  ierr = AllocateAppCtxFields(iga_input, user, &user->ice); CHKERRQ(ierr);

  ierr = IGACreateVec(iga_input, &ice_phase); CHKERRQ(ierr);
  ierr = IGAReadVec(iga_input, ice_phase, vec_file); CHKERRQ(ierr);
  ierr = EvaluateFieldAtGaussPoints(user, iga_input, ice_phase); CHKERRQ(ierr);

  ierr = VecDestroy(&ice_phase); CHKERRQ(ierr);
  ierr = IGADestroy(&iga_input); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*-----------------------------------------------------------------------------
  ComputeCircleIceField
  Populate user->ice with a single circular grain centred in the domain.
  The ice fraction uses a hyperbolic-tangent interface of width user->eps.
  Radius = min(Lx, Ly) / 16.
-----------------------------------------------------------------------------*/
PetscErrorCode ComputeCircleIceField(AppCtx *user)
{
  PetscErrorCode ierr;
  IGAElement     element;
  IGAPoint       point;
  PetscInt       idx, grainID;

  const PetscReal centX = user->Lx / 2.0;
  const PetscReal centY = user->Ly / 2.0;
  const PetscReal radius = PetscMin(user->Lx, user->Ly) / 16.0;

  PetscFunctionBegin;

  ierr = AllocateAppCtxFields(user->iga, user, &user->ice); CHKERRQ(ierr);

  ierr = IGABeginElement(user->iga, &element); CHKERRQ(ierr);
  while (IGANextElement(user->iga, element)) {
    ierr = IGAElementBeginPoint(element, &point); CHKERRQ(ierr);
    while (IGAElementNextPoint(element, point)) {
      idx = point->index + point->count * point->parent->index;
      PetscReal dist = PetscSqrtReal(SQ(point->mapX[0][0] - centX) +
                                     SQ(point->mapX[0][1] - centY)) - radius;
      PetscReal ice = 0.5 - 0.5 * PetscTanhReal(0.5 / user->eps * dist);
      user->ice[idx] = PetscMax(0.0, PetscMin(1.0, ice));
      (void)grainID; /* suppress unused-variable warning */
    }
    ierr = IGAElementEndPoint(element, &point); CHKERRQ(ierr);
  }
  ierr = IGAEndElement(user->iga, &element); CHKERRQ(ierr);

  PetscPrintf(PETSC_COMM_WORLD, "Ice field: circle mode (radius = %g m).\n", radius);
  PetscFunctionReturn(0);
}

/*-----------------------------------------------------------------------------
  ComputeLayeredIceField
  Populate user->ice with a horizontal layer: ice (phi=1) below y = Ly/2,
  air (phi=0) above.  Interface width controlled by user->eps.
-----------------------------------------------------------------------------*/
PetscErrorCode ComputeLayeredIceField(AppCtx *user)
{
  PetscErrorCode ierr;
  IGAElement     element;
  IGAPoint       point;
  PetscInt       indGP;

  PetscFunctionBegin;

  ierr = AllocateAppCtxFields(user->iga, user, &user->ice); CHKERRQ(ierr);

  ierr = IGABeginElement(user->iga, &element); CHKERRQ(ierr);
  while (IGANextElement(user->iga, element)) {
    ierr = IGAElementBeginPoint(element, &point); CHKERRQ(ierr);
    while (IGAElementNextPoint(element, point)) {
      indGP = point->index + point->count * point->parent->index;
      PetscReal dist     = point->mapX[0][1] - user->Ly / 2.0;
      PetscReal ice      = 0.5 - 0.5 * PetscTanhReal(0.5 / user->eps * dist);
      user->ice[indGP]   = PetscMax(0.0, PetscMin(1.0, ice));
    }
    ierr = IGAElementEndPoint(element, &point); CHKERRQ(ierr);
  }
  ierr = IGAEndElement(user->iga, &element); CHKERRQ(ierr);

  PetscPrintf(PETSC_COMM_WORLD, "Ice field: layered mode.\n");
  PetscFunctionReturn(0);
}

/*-----------------------------------------------------------------------------
  FormInitialCondition
  Dispatch to the appropriate ice-field initialiser based on user->init_mode.
-----------------------------------------------------------------------------*/
PetscErrorCode FormInitialCondition(AppCtx *user)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;

  if (strcmp(user->init_mode, "circle") == 0) {
    ierr = ComputeCircleIceField(user); CHKERRQ(ierr);
  } else if (strcmp(user->init_mode, "layered") == 0) {
    ierr = ComputeLayeredIceField(user); CHKERRQ(ierr);
  } else {
    /* File-based: read IGA + solution vector from user->init_dir */
    char iga_file[PETSC_MAX_PATH_LEN];
    char sol_file[PETSC_MAX_PATH_LEN];
    snprintf(iga_file, sizeof(iga_file), "%s/igasol.dat",         user->init_dir);
    snprintf(sol_file, sizeof(sol_file), "%s/sol_%05d.dat",
             user->init_dir, user->sol_index);

    PetscPrintf(PETSC_COMM_WORLD, "Reading IGA :        %s\n", iga_file);
    PetscPrintf(PETSC_COMM_WORLD, "Reading ice field:   %s\n\n", sol_file);
    ierr = ReadSolutionVec(iga_file, sol_file, &user->iga_input, user); CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

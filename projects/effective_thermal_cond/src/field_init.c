#include "field_init.h"
#include "material.h"

/*-----------------------------------------------------------------------------
  CheckVecMatchesIGA
  Compare the global length of the IGA-shaped vector we are about to fill
  against the length actually stored in the PETSc binary file, and fail with
  the discretisation parameters spelled out if they disagree.

  IGAReadVec would catch the mismatch on its own, but only as a bare
  "Vec sizes differ" -- and the cause is always one of a handful of parameters
  disagreeing with the run that produced the file (-p, -C, -periodic, -Nx/-Ny,
  or a snapshot with a different DOF count). Printing them turns a puzzling
  failure into an obvious one.

  Compare against the NATURAL vector, not the global one. IGAWriteVec stores
  the natural (unwrapped tensor-grid) layout, which for a periodic axis has
  the p duplicated wrap nodes present: measured on a 32x32, p=2, C=1 mesh, the
  global vector is 3072 while the file holds 3468 either way. A consequence
  worth knowing is that the file length is IDENTICAL periodic or not, so this
  check cannot detect a -periodic mismatch -- it catches -dim/-Nx/-Ny/-Nz/-p/-C
  and DOF-count errors, which is what actually bit us (this solver hardcoded
  p=2, C=1 while the producer was configured p=1, C=0).

  PETSc's binary Vec format is [VEC_FILE_CLASSID, global_rows, values...],
  written big-endian, so the length is the second int in the file.
-----------------------------------------------------------------------------*/
PetscErrorCode CheckVecMatchesIGA(IGA iga, const char *vec_file, AppCtx *user)
{
  PetscErrorCode ierr;
  PetscViewer    viewer;
  Vec            natural;
  PetscInt       header[2] = {0, 0};
  PetscInt       expected;

  PetscFunctionBegin;

  ierr = IGAGetNaturalVec(iga, &natural); CHKERRQ(ierr);
  ierr = VecGetSize(natural, &expected);  CHKERRQ(ierr);

  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, vec_file,
                               FILE_MODE_READ, &viewer); CHKERRQ(ierr);
  ierr = PetscViewerBinarySkipInfo(viewer); CHKERRQ(ierr);
  ierr = PetscViewerBinaryRead(viewer, header, 2, NULL, PETSC_INT); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

  if (header[1] != expected) {
    PetscInt nodes = expected / PF_INPUT_DOF;
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_UNEXPECTED,
            "\n*** INPUT MESH MISMATCH ***\n"
            "  %s holds %d values, but the configured mesh wants %d\n"
            "  (%d nodes x %d dof).\n"
            "  Configured: -dim %d -Nx %d -Ny %d -Nz %d -p %d -C %d -periodic %d\n"
            "  These must match the run that produced the file. Point this\n"
            "  solver at that run's own staged solver.opts + geometry .opts\n"
            "  rather than re-specifying them by hand.\n"
            "  (-periodic does not change the on-disk length, so a periodicity\n"
            "   mismatch is NOT detected here -- check it yourself.)\n",
            vec_file, (int)header[1], (int)expected, (int)nodes,
            (int)PF_INPUT_DOF, (int)user->dim, (int)user->Nx, (int)user->Ny,
            (int)user->Nz, (int)user->p, (int)user->C, (int)user->periodic);
  }

  PetscFunctionReturn(0);
}

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
  /* The producing solver writes 3 DOF per node: (ice phi, temperature,
     vapour density). Only DOF 0 is used here -- see
     EvaluateFieldAtGaussPoints, which takes sol[0]. */
  ierr = IGASetDof(iga_input, PF_INPUT_DOF); CHKERRQ(ierr);

  /* Periodicity must match the producing run. It does NOT change the on-disk
     length (IGAWriteVec stores the unwrapped natural grid either way), but it
     does change how those values scatter into the global vector, and hence
     what field gets evaluated at the Gauss points. IGARead() cannot supply
     it -- IGALoad stores only degree and knots, and IGAReset clears the
     periodic flag -- so it comes from the options database instead. */
  ierr = IGAGetAxis(iga_input, 0, &axis0); CHKERRQ(ierr);
  ierr = IGAAxisSetDegree(axis0, user->p); CHKERRQ(ierr);
  ierr = IGAAxisSetPeriodic(axis0, user->periodic); CHKERRQ(ierr);
  ierr = IGAAxisInitUniform(axis0, user->Nx, 0.0, user->Lx, user->C); CHKERRQ(ierr);

  ierr = IGAGetAxis(iga_input, 1, &axis1); CHKERRQ(ierr);
  ierr = IGAAxisSetDegree(axis1, user->p); CHKERRQ(ierr);
  ierr = IGAAxisSetPeriodic(axis1, user->periodic); CHKERRQ(ierr);
  ierr = IGAAxisInitUniform(axis1, user->Ny, 0.0, user->Ly, user->C); CHKERRQ(ierr);

  if (user->dim == 3) {
    ierr = IGAGetAxis(iga_input, 2, &axis2); CHKERRQ(ierr);
    ierr = IGAAxisSetDegree(axis2, user->p); CHKERRQ(ierr);
    ierr = IGAAxisSetPeriodic(axis2, user->periodic); CHKERRQ(ierr);
    ierr = IGAAxisInitUniform(axis2, user->Nz, 0.0, user->Lz, user->C); CHKERRQ(ierr);
  }

  ierr = IGASetFromOptions(iga_input); CHKERRQ(ierr);
  ierr = IGASetUp(iga_input); CHKERRQ(ierr);

  user->iga_input = iga_input;

  ierr = AllocateAppCtxFields(iga_input, user, &user->ice); CHKERRQ(ierr);

  ierr = IGACreateVec(iga_input, &ice_phase); CHKERRQ(ierr);
  ierr = CheckVecMatchesIGA(iga_input, vec_file, user); CHKERRQ(ierr);
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
  Populate user->ice with a PERIODIC horizontal ice stripe.

  Two tanh interfaces at y = Ly/4 and y = 3*Ly/4 create an ice band in the
  centre of the domain [Ly/4, 3*Ly/4] with air at the periodic boundaries:
    phi(0) = phi(Ly) = 0  (air)
    phi(Ly/2)           = 1  (ice)

  This gives exactly 50 % ice fraction and is consistent with the periodic
  boundary conditions imposed by the homogenisation cell problem.

  A single-interface tanh (phi(0)=1, phi(Ly)=0) is NOT periodic and yields
  incorrect effective conductivities.

  Interface width is controlled by user->eps.
-----------------------------------------------------------------------------*/
PetscErrorCode ComputeLayeredIceField(AppCtx *user)
{
  PetscErrorCode ierr;
  IGAElement     element;
  IGAPoint       point;
  PetscInt       indGP;

  PetscFunctionBegin;

  ierr = AllocateAppCtxFields(user->iga, user, &user->ice); CHKERRQ(ierr);

  const PetscReal a = 0.5 / user->eps;

  ierr = IGABeginElement(user->iga, &element); CHKERRQ(ierr);
  while (IGANextElement(user->iga, element)) {
    ierr = IGAElementBeginPoint(element, &point); CHKERRQ(ierr);
    while (IGAElementNextPoint(element, point)) {
      indGP = point->index + point->count * point->parent->index;
      PetscReal y   = point->mapX[0][1];
      /* Two rising/falling tanh transitions at Ly/4 and 3*Ly/4:
         phi = 0 outside [Ly/4, 3*Ly/4], phi = 1 in the interior */
      PetscReal ice = 0.5 * PetscTanhReal(a * (y - user->Ly / 4.0))
                   - 0.5 * PetscTanhReal(a * (y - 3.0 * user->Ly / 4.0));
      user->ice[indGP] = PetscMax(0.0, PetscMin(1.0, ice));
    }
    ierr = IGAElementEndPoint(element, &point); CHKERRQ(ierr);
  }
  ierr = IGAEndElement(user->iga, &element); CHKERRQ(ierr);

  PetscPrintf(PETSC_COMM_WORLD,
              "Ice field: layered mode (periodic stripe, two interfaces at Ly/4 and 3Ly/4).\n");
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

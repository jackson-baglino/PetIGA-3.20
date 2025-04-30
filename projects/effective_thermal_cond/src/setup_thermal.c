#include "setup_thermal.h"
#include "assembly.h"
#include "utils.h"
#include "io_thermal.h"

/*-----------------------------------------------------------
   Function: ComputeCircleIceField
   Purpose:  Computes the initial ice field using a circular profile.
             The field is computed at each Gauss (integration) point based on
             the distance from the center of the domain, using a hyperbolic
             tangent function.
   Parameters:
     - user: Pointer to the AppCtx structure containing grid and simulation parameters.
   Returns:
     - PetscErrorCode (0 on success).
-----------------------------------------------------------*/
static PetscErrorCode ComputeCircleIceField(AppCtx *user) {
    PetscErrorCode ierr;
    PetscFunctionBegin;

    IGAElement element;
    IGAPoint point;
    PetscReal dist;
    PetscInt idx = 0;

    // Allocate memory for user->ice (if not already allocated)
    ierr = AllocateAppCtxFields(user->iga, user, &user->ice); CHKERRQ(ierr);

    ierr = IGABeginElement(user->iga, &element); CHKERRQ(ierr);
    while (IGANextElement(user->iga, element)) {
        ierr = IGAElementBeginPoint(element, &point); CHKERRQ(ierr);
        while (IGAElementNextPoint(element, point)) {
          idx = point->index + point->count * point->parent->index;

          PetscReal radius = PetscMin(user->Lx / 6.0, user->Ly / 6.0);
          dist = PetscSqrtReal(SQ(point->mapX[0][0] - user->Lx / 2.0) +
                                SQ(point->mapX[0][1] - user->Ly / 2.0)) - radius;
          user->ice[idx] = 0.5 - 0.5 * PetscTanhReal(0.5 / user->eps * dist);
          user->ice[idx] = PetscMax(0.0, PetscMin(1.0, user->ice[idx]));
        }
        ierr = IGAElementEndPoint(element, &point); CHKERRQ(ierr);
    }
    ierr = IGAEndElement(user->iga, &element); CHKERRQ(ierr);

    PetscPrintf(PETSC_COMM_WORLD, "Circle mode: ice field computed with %d points.\n", idx);
    PetscFunctionReturn(0);
}

/*-----------------------------------------------------------
   Function: ComputeLayeredIceField
   Purpose:  Computes the initial ice field using a layered (vertical)
             profile. The field is computed at each Gauss (integration) point,
             with variation assumed only in the y-direction.
   Parameters:
     - user: Pointer to the AppCtx structure containing grid parameters.
   Returns:
     - PetscErrorCode (0 on success).
-----------------------------------------------------------*/
static PetscErrorCode ComputeLayeredIceField(AppCtx *user) {
    PetscErrorCode ierr;
    PetscFunctionBegin;

    IGAElement element;
    IGAPoint point;
    PetscReal dist;
    PetscInt indGP;

    // Allocate memory for user->ice (if not already allocated)
    ierr = AllocateAppCtxFields(user->iga, user, &user->ice); CHKERRQ(ierr);

    ierr = IGABeginElement(user->iga, &element); CHKERRQ(ierr);
    while (IGANextElement(user->iga, element)) {
        ierr = IGAElementBeginPoint(element, &point); CHKERRQ(ierr);
        while (IGAElementNextPoint(element, point)) {
          indGP = point->index + point->count * point->parent->index;

          // Find distance from the center of the domain--above midline is air, below is ice
          dist = point->mapX[0][1] - user->Ly / 2.0;

          user->ice[indGP] = 0.5 - 0.5 * PetscTanhReal(0.5 / user->eps * dist);
          user->ice[indGP] = PetscMax(0.0, PetscMin(1.0, user->ice[indGP]));   

          // user->ice[indGP] = 1.0;
        }
        ierr = IGAElementEndPoint(element, &point); CHKERRQ(ierr);
    }
    ierr = IGAEndElement(user->iga, &element); CHKERRQ(ierr);

    PetscPrintf(PETSC_COMM_WORLD, "--- Layered mode. ---\n\n");
    PetscFunctionReturn(0);
}

/*-----------------------------------------------------------
   Function: FormInitialCondition
   Purpose:  Determines the initialization mode and calls the
             appropriate helper function:
               - "circle": Calls ComputeCircleIceField.
               - "layered": Calls ComputeLayeredIceField.
               - Otherwise: Assumes user->init_mode is a file path (nodal data),
                 reads raw nodal data using ReadInputFile, and then interpolates
                 these nodal values onto the integration (Gauss) points using
                 the element’s shape functions.
   Parameters:
     - user: Pointer to the AppCtx structure.
   Returns:
     - PetscErrorCode (0 on success).
-----------------------------------------------------------*/
PetscErrorCode FormInitialCondition(AppCtx *user)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;

  if (strcmp(user->init_mode, "circle") == 0) {
    ierr = ComputeCircleIceField(user); CHKERRQ(ierr);
  } else if (strcmp(user->init_mode, "layered") == 0) {
    ierr = ComputeLayeredIceField(user); CHKERRQ(ierr);
  } else {
    // NEEDS TO BE CORRECTED!
    PetscPrintf(PETSC_COMM_WORLD, "Reading ice field from file: %s\n", user->init_mode);
    const char *iga_file = "/Users/jacksonbaglino/PetIGA-3.20/projects/effective_thermal_cond/inputs/igasol.dat";
    PetscPrintf(PETSC_COMM_WORLD, "[WARNING] May need to update path to IGA object.\n");
    PetscPrintf(PETSC_COMM_WORLD, "Reading solution vector from file: %s\n", iga_file);

    // const char *vec_file = "/Users/jacksonbaglino/PetIGA-3.20/projects/effective_thermal_cond/inputs/sol_00000.dat";
    ierr = ReadSolutionVec(iga_file, user->init_mode, &user->iga_input, user); CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

/*-----------------------------------------------------------
   Function: InitializeFields
   Purpose:  Initializes the ice field by calling FormInitialCondition.
             (The Vec for the ice field is allocated within FormInitialCondition.)
   Parameters:
     - user: Pointer to the AppCtx structure.
     - iga: The IGA object.
   Returns:
     - PetscErrorCode (0 on success).
-----------------------------------------------------------*/
PetscErrorCode InitializeFields(AppCtx *user, IGA iga) {
    PetscErrorCode ierr;
    PetscFunctionBegin;

    ierr = FormInitialCondition(user); CHKERRQ(ierr);

    PetscPrintf(PETSC_COMM_WORLD, "Field variables initialized.\n\n");
    PetscFunctionReturn(0);
}
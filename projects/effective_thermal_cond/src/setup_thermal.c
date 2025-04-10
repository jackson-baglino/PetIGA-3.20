#include "setup_thermal.h"
#include "assembly.h"
#include "utils.h"

// /*-----------------------------------------------------------
//    Function: ReadInputFile
//    Purpose:  Opens a text file specified by user->init_mode, reads raw
//              ice field values into a dynamically allocated array, and
//              returns the number of points read.
//    Parameters:
//      - user: Pointer to the AppCtx structure containing grid parameters.
//      - iceField: Address of a pointer that will hold the raw ice field data.
//      - nPoints: Pointer to a PetscInt to store the number of points read.
//    Returns:
//      - PetscErrorCode (0 on success).
// -----------------------------------------------------------*/
// static PetscErrorCode ReadInputFile(AppCtx *user, PetscReal **iceField, IGA *iga_DSM)
// {
//     PetscErrorCode ierr;
//     Vec DSM_sol;
//     FILE *file = fopen(user->init_mode, "r");
//     if (!file) {
//         SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_OPEN, "Error opening ice field file: %s", user->init_mode);
//     }


    
//     IGAReadVec(iga_DSM, DSM_sol, user->init_mode);

//     // /* Expected number of points */
//     // PetscInt total = user->Nx * user->Ny * (user->dim == 3 ? user->Nz : 1);
//     // ierr = PetscMalloc1(total, iceField); CHKERRQ(ierr);

//     // for (PetscInt i = 0; i < total; i++) {
//     //     if (fscanf(file, "%lf", &((*iceField)[i])) != 1) {
//     //         fclose(file);
//     //         SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_READ, "Error reading ice field file at index %D", i);
//     //     }
//     // }
//     // fclose(file);
//     // *nPoints = total;
//     return 0;
// }

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
            PetscReal radius = PetscMin(user->Lx / 6.0, user->Ly / 6.0);
            dist = PetscSqrtReal(SQ(point->mapX[0][0] - user->Lx / 2.0) +
                                 SQ(point->mapX[0][1] - user->Ly / 2.0)) - radius;
            user->ice[idx] = 0.5 - 0.5 * PetscTanhReal(0.5 / user->eps * dist);
            user->ice[idx] = PetscMax(0.0, PetscMin(1.0, user->ice[idx]));
            idx++;
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
    PetscInt idx = 0;

    // Allocate memory for user->ice (if not already allocated)
    ierr = AllocateAppCtxFields(user->iga, user, &user->ice); CHKERRQ(ierr);

    ierr = IGABeginElement(user->iga, &element); CHKERRQ(ierr);
    while (IGANextElement(user->iga, element)) {
        ierr = IGAElementBeginPoint(element, &point); CHKERRQ(ierr);
        while (IGAElementNextPoint(element, point)) {
            dist = point->mapX[0][1] - user->Ly / 2.0;
            user->ice[idx] = 0.5 - 0.5 * PetscTanhReal(0.5 / user->eps * dist);
            user->ice[idx] = PetscMax(0.0, PetscMin(1.0, user->ice[idx]));
            idx++;
        }
        ierr = IGAElementEndPoint(element, &point); CHKERRQ(ierr);
    }
    ierr = IGAEndElement(user->iga, &element); CHKERRQ(ierr);

    PetscPrintf(PETSC_COMM_WORLD, "Layered mode: ice field computed with %d points.\n", idx);
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
    ierr = ReadSolutionVec(user->init_mode, user->init_mode, &user->iga_input, &user->ice, user); CHKERRQ(ierr);
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

    PetscPrintf(PETSC_COMM_WORLD, "Field variables initialized.\n");
    PetscFunctionReturn(0);
}
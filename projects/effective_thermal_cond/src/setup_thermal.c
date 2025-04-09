#include "setup_thermal.h"
#include "assembly.h"

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
static PetscErrorCode ComputeCircleIceField(AppCtx *user)
{
    PetscErrorCode ierr;
    IGAElement element;
    IGAPoint point;
    PetscReal dist;
    PetscInt idx = 0;
    PetscInt total = user->Nx * user->Ny * (user->dim == 3 ? user->Nz : 1);

    /* Create a PETSc Vec to hold the ice field at integration (Gauss) points */
    Vec ice;
    ierr = VecCreate(PETSC_COMM_WORLD, &ice); CHKERRQ(ierr);
    ierr = VecSetSizes(ice, PETSC_DECIDE, total); CHKERRQ(ierr);
    ierr = VecSetFromOptions(ice); CHKERRQ(ierr);

    /* Get array pointer for direct access */
    PetscScalar *iceArray;
    ierr = VecGetArray(ice, &iceArray); CHKERRQ(ierr);

    ierr = IGABeginElement(user->iga, &element); CHKERRQ(ierr);
    while (IGANextElement(user->iga, element)) {
        ierr = IGAElementBeginPoint(element, &point); CHKERRQ(ierr);
        while (IGAElementNextPoint(element, point)) {
            PetscReal radius = PetscMin(user->Lx / 6.0, user->Ly / 6.0);
            /* Compute distance from the circle center */
            dist = PetscSqrtReal(SQ(point->mapX[0][0] - user->Lx / 2.0) +
                                 SQ(point->mapX[0][1] - user->Ly / 2.0)) - radius;
            /* Compute ice phase using a hyperbolic tangent transition */
            iceArray[idx] = 0.5 - 0.5 * PetscTanhReal(0.5 / user->eps * dist);
            iceArray[idx] = PetscMax(0.0, PetscMin(1.0, iceArray[idx])); /* Clamp [0,1] */
            idx++;
            /* For debugging */
            if (idx > total) {
                PetscPrintf(PETSC_COMM_WORLD, "There are %d points in ice Vec. Originally set to %d.\n", idx, total);
            }
        }
        ierr = IGAElementEndPoint(element, &point); CHKERRQ(ierr);
    }
    ierr = IGAEndElement(user->iga, &element); CHKERRQ(ierr);
    ierr = VecRestoreArray(ice, &iceArray); CHKERRQ(ierr);

    PetscPrintf(PETSC_COMM_WORLD, "Circle mode: ice field computed with %d points.\n", idx);
    PetscPrintf(PETSC_COMM_WORLD, "Ice field is %g times larger than expected.\n", (PetscReal)idx / (PetscReal)total);

    /* Assign the computed Vec to user->ice */
    user->ice = ice;
    return 0;
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
static PetscErrorCode ComputeLayeredIceField(AppCtx *user)
{
    PetscErrorCode ierr;
    IGAElement element;
    IGAPoint point;
    PetscReal dist;
    PetscInt idx = 0;
    PetscInt total = user->Nx * user->Ny * (user->dim == 3 ? user->Nz : 1);

    /* Create a PETSc Vec to hold the ice field */
    Vec ice;
    ierr = VecCreate(PETSC_COMM_WORLD, &ice); CHKERRQ(ierr);
    ierr = VecSetSizes(ice, PETSC_DECIDE, total); CHKERRQ(ierr);
    ierr = VecSetFromOptions(ice); CHKERRQ(ierr);

    /* Get array pointer for direct access */
    PetscScalar *iceArray;
    ierr = VecGetArray(ice, &iceArray); CHKERRQ(ierr);

    ierr = IGABeginElement(user->iga, &element); CHKERRQ(ierr);
    while (IGANextElement(user->iga, element)) {
        ierr = IGAElementBeginPoint(element, &point); CHKERRQ(ierr);
        while (IGAElementNextPoint(element, point)) {
            /* For layered mode, assume variation in the y-direction only */
            dist = point->mapX[0][1] - user->Ly / 2.0;
            iceArray[idx] = 0.5 - 0.5 * PetscTanhReal(0.5 / user->eps * dist);
            iceArray[idx] = PetscMax(0.0, PetscMin(1.0, iceArray[idx])); /* Clamp [0,1] */
            idx++;
        }
        ierr = IGAElementEndPoint(element, &point); CHKERRQ(ierr);
    }
    ierr = IGAEndElement(user->iga, &element); CHKERRQ(ierr);
    ierr = VecRestoreArray(ice, &iceArray); CHKERRQ(ierr);

    PetscPrintf(PETSC_COMM_WORLD, "Layered mode: ice field computed with %d points.\n", idx);

    /* Assign the computed Vec to user->ice */
    user->ice = ice;
    return 0;
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
        // /* File mode: read nodal data then interpolate to integration (Gauss) points */
        // PetscReal *rawField = NULL;
        // PetscInt nRaw = 0;

        // PetscPrintf(PETSC_COMM_WORLD, "Assuming init_mode is a file path.\n");
        // PetscPrintf(PETSC_COMM_WORLD, "Reading nodal ice field from %s\n", user->init_mode);
        // ierr = ReadInputFile(user, &rawField, &nRaw); CHKERRQ(ierr);

        // /* Determine total number of integration (Gauss) points by looping over all elements */
        // PetscInt total_gp = 0;
        // IGAElement element;
        // IGAPoint point;
        // ierr = IGABeginElement(user->iga, &element); CHKERRQ(ierr);
        // while (IGANextElement(user->iga, element)) {
        //     ierr = IGAElementBeginPoint(element, &point); CHKERRQ(ierr);
        //     while (IGAElementNextPoint(element, point)) {
        //         total_gp++;
        //     }
        //     ierr = IGAElementEndPoint(element, &point); CHKERRQ(ierr);
        // }
        // ierr = IGAEndElement(user->iga, &element); CHKERRQ(ierr);

        // /* Create a Vec to hold the ice field at integration points */
        // Vec ice;
        // ierr = VecCreate(PETSC_COMM_WORLD, &ice); CHKERRQ(ierr);
        // ierr = VecSetSizes(ice, PETSC_DECIDE, total_gp); CHKERRQ(ierr);
        // ierr = VecSetFromOptions(ice); CHKERRQ(ierr);

        // PetscInt idx = 0;
        // ierr = IGABeginElement(user->iga, &element); CHKERRQ(ierr);
        // while (IGANextElement(user->iga, element)) {
        //     /* Get nodal connectivity for the element using the five-argument version */
        //     PetscInt imin, imax, nen;
        //     const PetscInt *rowmap, *colmap;
        //     ierr = IGAElementGetIndices(element, &imin, &imax, &rowmap, &nen); CHKERRQ(ierr);
 
        //     /* Allocate temporary array to hold nodal values for this element */
        //     PetscReal *nodalVals;
        //     ierr = PetscMalloc1(nen, &nodalVals); CHKERRQ(ierr);
        //     for (PetscInt a = 0; a < nen; a++) {
        //         if (rowmap[a] < nRaw) {
        //             nodalVals[a] = rawField[rowmap[a]];
        //         } else {
        //             nodalVals[a] = 0.0;
        //         }
        //     }
 
        //     ierr = IGAElementBeginPoint(element, &point); CHKERRQ(ierr);
        //     while (IGAElementNextPoint(element, point)) {
        //         const PetscScalar *N0;
        //         ierr = IGAPointGetShapeFuns(point, 0, &N0); CHKERRQ(ierr);
        //         PetscReal value = 0.0;
        //         for (PetscInt a = 0; a < nen; a++) {
        //             value += PetscRealPart(N0[a]) * nodalVals[a];
        //         }
        //         ierr = VecSetValue(ice, idx, value, INSERT_VALUES); CHKERRQ(ierr);
        //         idx++;
        //     }
        //     ierr = IGAElementEndPoint(element, &point); CHKERRQ(ierr);
        //     ierr = PetscFree(nodalVals); CHKERRQ(ierr);
        // }
        // ierr = IGAEndElement(user->iga, &element); CHKERRQ(ierr);
 
        // ierr = VecAssemblyBegin(ice); CHKERRQ(ierr);
        // ierr = VecAssemblyEnd(ice); CHKERRQ(ierr);
        // ierr = PetscFree(rawField); CHKERRQ(ierr);
 
        // PetscPrintf(PETSC_COMM_WORLD, "File mode: ice field interpolated to integration points (%d points).\n", total_gp);
        // user->ice = ice;
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
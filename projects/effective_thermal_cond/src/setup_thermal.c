#include "setup_thermal.h"
#include "solver.h"


/*-----------------------------------------------------------
   Function: InitializeUserContext
   Purpose:  Initialize the application context (AppCtx) with default
             material properties and discretization parameters.
   Parameters:
     - user: Pointer to the AppCtx structure to be initialized.
   Returns:
     - void.
-----------------------------------------------------------*/
void InitializeUserContext(AppCtx *user) {
    /* Set material properties */
    user->thcond_ice = 2.29;  
    user->thcond_air = 0.02;  
    user->cp_ice     = 1.96e3;   
    user->cp_air     = 1.044e3;  
    user->rho_ice    = 919.0;  
    user->rho_air    = 1.341;  

    /* Set default polynomial order and continuity */
    user->p = 1;
    user->C = 0;

    PetscPrintf(PETSC_COMM_WORLD, "User context initialized.\n");
}

/*-----------------------------------------------------------
   Function: ReadInputFile
   Purpose:  Opens a text file specified by user->init_mode, reads raw
             ice field values into a dynamically allocated array, and
             returns the number of points read.
   Parameters:
     - user: Pointer to the AppCtx structure containing grid parameters.
     - iceField: Address of a pointer that will hold the raw ice field data.
     - nPoints: Pointer to a PetscInt to store the number of points read.
   Returns:
     - PetscErrorCode (0 on success).
-----------------------------------------------------------*/
static PetscErrorCode ReadInputFile(AppCtx *user, PetscReal **iceField, IGA *iga_DSM)
{
    PetscErrorCode ierr;
    Vec DSM_sol;
    FILE *file = fopen(user->init_mode, "r");
    if (!file) {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_OPEN, "Error opening ice field file: %s", user->init_mode);
    }


    
    IGAReadVec(iga_DSM, DSM_sol, user->init_mode);

    // /* Expected number of points */
    // PetscInt total = user->Nx * user->Ny * (user->dim == 3 ? user->Nz : 1);
    // ierr = PetscMalloc1(total, iceField); CHKERRQ(ierr);

    // for (PetscInt i = 0; i < total; i++) {
    //     if (fscanf(file, "%lf", &((*iceField)[i])) != 1) {
    //         fclose(file);
    //         SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_READ, "Error reading ice field file at index %D", i);
    //     }
    // }
    // fclose(file);
    // *nPoints = total;
    return 0;
}

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
        /* File mode: read nodal data then interpolate to integration (Gauss) points */
        PetscReal *rawField = NULL;
        PetscInt nRaw = 0;

        PetscPrintf(PETSC_COMM_WORLD, "Assuming init_mode is a file path.\n");
        PetscPrintf(PETSC_COMM_WORLD, "Reading nodal ice field from %s\n", user->init_mode);
        ierr = ReadInputFile(user, &rawField, &nRaw); CHKERRQ(ierr);

        /* Determine total number of integration (Gauss) points by looping over all elements */
        PetscInt total_gp = 0;
        IGAElement element;
        IGAPoint point;
        ierr = IGABeginElement(user->iga, &element); CHKERRQ(ierr);
        while (IGANextElement(user->iga, element)) {
            ierr = IGAElementBeginPoint(element, &point); CHKERRQ(ierr);
            while (IGAElementNextPoint(element, point)) {
                total_gp++;
            }
            ierr = IGAElementEndPoint(element, &point); CHKERRQ(ierr);
        }
        ierr = IGAEndElement(user->iga, &element); CHKERRQ(ierr);

        /* Create a Vec to hold the ice field at integration points */
        Vec ice;
        ierr = VecCreate(PETSC_COMM_WORLD, &ice); CHKERRQ(ierr);
        ierr = VecSetSizes(ice, PETSC_DECIDE, total_gp); CHKERRQ(ierr);
        ierr = VecSetFromOptions(ice); CHKERRQ(ierr);

        PetscInt idx = 0;
        ierr = IGABeginElement(user->iga, &element); CHKERRQ(ierr);
        while (IGANextElement(user->iga, element)) {
            /* Get nodal connectivity for the element using the five-argument version */
            PetscInt imin, imax, nen;
            const PetscInt *rowmap, *colmap;
            ierr = IGAElementGetIndices(element, &imin, &imax, &rowmap, &nen); CHKERRQ(ierr);
 
            /* Allocate temporary array to hold nodal values for this element */
            PetscReal *nodalVals;
            ierr = PetscMalloc1(nen, &nodalVals); CHKERRQ(ierr);
            for (PetscInt a = 0; a < nen; a++) {
                if (rowmap[a] < nRaw) {
                    nodalVals[a] = rawField[rowmap[a]];
                } else {
                    nodalVals[a] = 0.0;
                }
            }
 
            ierr = IGAElementBeginPoint(element, &point); CHKERRQ(ierr);
            while (IGAElementNextPoint(element, point)) {
                const PetscScalar *N0;
                ierr = IGAPointGetShapeFuns(point, 0, &N0); CHKERRQ(ierr);
                PetscReal value = 0.0;
                for (PetscInt a = 0; a < nen; a++) {
                    value += PetscRealPart(N0[a]) * nodalVals[a];
                }
                ierr = VecSetValue(ice, idx, value, INSERT_VALUES); CHKERRQ(ierr);
                idx++;
            }
            ierr = IGAElementEndPoint(element, &point); CHKERRQ(ierr);
            ierr = PetscFree(nodalVals); CHKERRQ(ierr);
        }
        ierr = IGAEndElement(user->iga, &element); CHKERRQ(ierr);
 
        ierr = VecAssemblyBegin(ice); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(ice); CHKERRQ(ierr);
        ierr = PetscFree(rawField); CHKERRQ(ierr);
 
        PetscPrintf(PETSC_COMM_WORLD, "File mode: ice field interpolated to integration points (%d points).\n", total_gp);
        user->ice = ice;
    }

    PetscFunctionReturn(0);
}

/*-----------------------------------------------------------
   Function: ApplyBoundaryConditions
   Purpose:  Applies boundary conditions for the thermal model.
             It sets a fixed temperature at the top, prescribed flux at
             the bottom, and, for 3D problems, zero-flux conditions on the z-axis.
   Parameters:
     - iga: The IGA object representing the discretization.
     - user: Pointer to the AppCtx structure containing boundary parameters.
   Returns:
     - PetscErrorCode (0 on success).
-----------------------------------------------------------*/
PetscErrorCode ApplyBoundaryConditions(IGA iga, AppCtx *user) {
    PetscErrorCode ierr;
    PetscFunctionBegin;

    PetscPrintf(PETSC_COMM_WORLD, "Applying boundary conditions...\n");

    /* Fixed temperature at the top */
    ierr = IGASetBoundaryValue(iga, 1, 1, 0, user->T_top); CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD, "  - Fixed temperature applied at top: T = %g K\n", user->T_top);

    /* Prescribed flux at the bottom */
    ierr = IGASetBoundaryForm(iga, 1, 0, PETSC_TRUE); CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD, "  - Prescribed flux at bottom: q = %g W/m²\n", user->q_bottom);

    /* Zero-flux (Neumann) conditions on the z-axis (if 3D) */
    if (user->dim == 3) {
        ierr = IGASetBoundaryForm(iga, 2, 0, PETSC_TRUE); CHKERRQ(ierr);
        ierr = IGASetBoundaryForm(iga, 2, 1, PETSC_TRUE); CHKERRQ(ierr);
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

/*-----------------------------------------------------------
   Function: SetupIGA
   Purpose:  Sets up the IGA object by defining the problem domain,
             degrees of freedom, field names, and the uniform grid.
   Parameters:
     - user: Pointer to the AppCtx structure containing grid and simulation parameters.
     - iga: Address of the IGA object pointer to be created and configured.
   Returns:
     - PetscErrorCode (0 on success).
-----------------------------------------------------------*/
PetscErrorCode SetupIGA(AppCtx *user, IGA *iga) {
    PetscErrorCode ierr;
    PetscFunctionBegin;

    ierr = IGACreate(PETSC_COMM_WORLD, iga); CHKERRQ(ierr);
    ierr = IGASetDim(*iga, user->dim); CHKERRQ(ierr);
    ierr = IGASetDof(*iga, 1); CHKERRQ(ierr);
    ierr = IGASetFieldName(*iga, 0, "temperature"); CHKERRQ(ierr);

    /* Set up uniform grid along each axis */
    IGAAxis axisX, axisY, axisZ;
    ierr = IGAGetAxis(*iga, 0, &axisX); CHKERRQ(ierr);
    ierr = IGAAxisSetDegree(axisX, user->p); CHKERRQ(ierr);
    ierr = IGAGetAxis(*iga, 1, &axisY); CHKERRQ(ierr);
    ierr = IGAAxisSetDegree(axisY, user->p); CHKERRQ(ierr);
    if (user->dim == 3) {
        ierr = IGAGetAxis(*iga, 2, &axisZ); CHKERRQ(ierr);
    }
    
    /* Initialize each axis */
    ierr = IGAAxisInitUniform(axisX, user->Nx, 0.0, user->Lx, user->C); CHKERRQ(ierr);
    ierr = IGAAxisInitUniform(axisY, user->Ny, 0.0, user->Ly, user->C); CHKERRQ(ierr);
    if (user->dim == 3) {
        ierr = IGAAxisInitUniform(axisZ, user->Nz, 0.0, user->Lz, user->C); CHKERRQ(ierr);
    }

    /* Finalize IGA setup */
    ierr = IGASetFromOptions(*iga); CHKERRQ(ierr);
    ierr = IGASetUp(*iga); CHKERRQ(ierr);

    PetscPrintf(PETSC_COMM_WORLD, "IGA setup complete.\n");
    PetscFunctionReturn(0);
}

/*-----------------------------------------------------------
   Function: SetupAndSolve
   Purpose:  Assembles the system matrix and right-hand side vector,
             computes the initial condition for the temperature field,
             and solves the resulting linear system using KSP.
   Parameters:
     - user: Pointer to the AppCtx structure.
     - iga: The IGA object representing the discretization.
   Returns:
     - PetscErrorCode (0 on success).
-----------------------------------------------------------*/
PetscErrorCode SetupAndSolve(AppCtx *user, IGA iga) {
    PetscErrorCode ierr;
    Mat A;
    Vec b;
    KSP ksp;

    PetscFunctionBegin;

    ierr = IGACreateMat(iga, &A); CHKERRQ(ierr);
    ierr = IGACreateVec(iga, &user->T_sol); CHKERRQ(ierr);
    ierr = IGACreateVec(iga, &b); CHKERRQ(ierr);

    /* Assemble the system */
    ierr = IGASetFormSystem(iga, AssembleStiffnessMatrix, user); CHKERRQ(ierr);
    ierr = IGAComputeSystem(iga, A, b); CHKERRQ(ierr);

    /* Compute initial condition for the temperature field */
    ierr = ComputeInitialCondition(user->T_sol, user); CHKERRQ(ierr);

    /* Create and configure the KSP solver */
    ierr = IGACreateKSP(iga, &ksp); CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp, A, A); CHKERRQ(ierr);
    ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);
    ierr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE); CHKERRQ(ierr);
    ierr = KSPSolve(ksp, b, user->T_sol); CHKERRQ(ierr);

    PetscPrintf(PETSC_COMM_WORLD, "KSP solve complete.\n");

    ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
    ierr = MatDestroy(&A); CHKERRQ(ierr);
    ierr = VecDestroy(&b); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
#include "thermal_solver.h"

void InitializeUserContext(AppCtx *user) {
    // Set material properties
    user->thcond_ice = 2.29;  
    user->thcond_air = 0.02;  
    user->cp_ice = 1.96e3;   
    user->cp_air = 1.044e3;  
    user->rho_ice = 919.0;  
    user->rho_air = 1.341;  

    // Set default polynomial order and continuity
    user->p = 1;
    user->C = 0;

    PetscPrintf(PETSC_COMM_WORLD, "User context initialized.\n");
}

// Initialize the ice field based on a phase-field approach
PetscErrorCode FormInitialCondition(AppCtx *user) {
    PetscFunctionBegin;
    
    PetscErrorCode ierr;
    IGAElement element;
    IGAPoint point;
    PetscReal ice = 0.0, dist;
    PetscInt idx = 0;
    
    // =============================
    // 1. COMPUTE GAUSS POINTS & INITIALIZE ICE FIELD
    // =============================
    ierr = IGABeginElement(user->iga, &element); CHKERRQ(ierr);
    while (IGANextElement(user->iga, element)) {
        ierr = IGAElementBeginPoint(element, &point); CHKERRQ(ierr);
        while (IGAElementNextPoint(element, point)) {
            // Check initialization mode (init_mode)
            if (strcmp(user->init_mode, "circle") == 0) {
                PetscReal radius = PetscMin(user->Lx / 6.0, user->Ly / 6.0);

                // Initialize ice field in a circular region
                dist = PetscSqrtReal(SQ(point->mapX[0][0] - user->Lx / 2.0) + 
                                     SQ(point->mapX[0][1] - user->Ly / 2.0)) - radius;

                // PetscPrintf(PETSC_COMM_SELF, "Circle radius: %g, dist: %g\n", radius, dist);
            } else if (strcmp(user->init_mode, "layered") == 0) {
                // Initialize ice field in a layered manner
                dist = point->mapX[0][1] - user->Ly / 2.0;
            } else {
                PetscPrintf(PETSC_COMM_SELF, "Assiming init_mode is file path.\n");
                PetscPrintf(PETSC_COMM_SELF, "Reading ice field from %s\n", user->init_mode);
            }

            // Apply phase field function to compute ice phase
            ice = 0.5 - 0.5 * PetscTanhReal(2.0 / user->eps * dist);
            ice = PetscMax(0.0, PetscMin(1.0, ice)); // Clamp between [0,1]

            // Store computed ice field values 
            user->ice[idx] = ice;
            idx++;
        }
        ierr = IGAElementEndPoint(element, &point); CHKERRQ(ierr);
    }
    ierr = IGAEndElement(user->iga, &element); CHKERRQ(ierr);
    
    PetscPrintf(PETSC_COMM_WORLD, "Ice field initialized with %d Gauss points.\n", idx);

    PetscFunctionReturn(0);
}

PetscErrorCode ApplyBoundaryConditions(IGA iga, AppCtx *user) {
    PetscErrorCode ierr;
    PetscFunctionBegin;

    PetscPrintf(PETSC_COMM_WORLD, "Applying boundary conditions...\n");

    /* ===========================
       FIXED TEMPERATURE AT THE TOP
       =========================== */
    ierr = IGASetBoundaryValue(iga, 1, 1, 0, user->T_top); CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD, "  - Fixed temperature applied at top: T = %g K\n", user->T_top);

    /* ===========================
       PRESCRIBED FLUX AT THE BOTTOM
       =========================== */
    ierr = IGASetBoundaryForm(iga, 1, 0, PETSC_TRUE); CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD, "  - Prescribed flux at bottom: q = %g W/m²\n", user->q_bottom);

    /* ===========================
       ZERO-FLUX (NEUMANN) BCs on SIDES
       =========================== */
    // ierr = IGASetBoundaryForm(iga, 0, 0, PETSC_TRUE); CHKERRQ(ierr);  // x = 0
    // ierr = IGASetBoundaryForm(iga, 0, 1, PETSC_TRUE); CHKERRQ(ierr);  // x = Lx
    // PetscPrintf(PETSC_COMM_WORLD, "  - Zero-flux (Neumann) BCs applied at x = 0 and x = Lx\n");

    /* ===========================
       ZERO-FLUX (NEUMANN) BCs on Z-AXIS (ONLY IF 3D)
       =========================== */
    if (user->dim == 3) {
        ierr = IGASetBoundaryForm(iga, 2, 0, PETSC_TRUE); CHKERRQ(ierr);  // z = 0
        ierr = IGASetBoundaryForm(iga, 2, 1, PETSC_TRUE); CHKERRQ(ierr);  // z = Lz
    }

    PetscFunctionReturn(0);
}

PetscErrorCode InitializeFields(AppCtx *user, IGA iga) {
    PetscErrorCode ierr;
    PetscInt nmb = iga->elem_width[0] * iga->elem_width[1] * SQ(user->p + 1);

    PetscFunctionBegin;

    // Allocate memory for ice field
    ierr = PetscMalloc1(nmb, &user->ice); CHKERRQ(ierr);
    ierr = PetscMemzero(user->ice, nmb * sizeof(PetscReal)); CHKERRQ(ierr);

    // Call function to initialize ice field values
    ierr = FormInitialCondition(user); CHKERRQ(ierr);

    PetscPrintf(PETSC_COMM_WORLD, "Field variables initialized.\n");
    PetscFunctionReturn(0);
}

PetscErrorCode SetupIGA(AppCtx *user, IGA *iga) {
    PetscErrorCode ierr;
    PetscFunctionBegin;

    ierr = IGACreate(PETSC_COMM_WORLD, iga); CHKERRQ(ierr);
    ierr = IGASetDim(*iga, user->dim); CHKERRQ(ierr);
    ierr = IGASetDof(*iga, 1); CHKERRQ(ierr);
    ierr = IGASetFieldName(*iga, 0, "temperature"); CHKERRQ(ierr);

    // Set up uniform grid along each axis
    IGAAxis axisX, axisY, axisZ;
    ierr = IGAGetAxis(*iga, 0, &axisX); CHKERRQ(ierr);
    ierr = IGAAxisSetDegree(axisX, user->p); CHKERRQ(ierr);
    ierr = IGAGetAxis(*iga, 1, &axisY); CHKERRQ(ierr);
    ierr = IGAAxisSetDegree(axisY, user->p); CHKERRQ(ierr);
    if (user->dim == 3) {
        ierr = IGAGetAxis(*iga, 2, &axisZ); CHKERRQ(ierr);
    }
    
    // Initialize each axis
    ierr = IGAAxisInitUniform(axisX, user->Nx, 0.0, user->Lx, user->C); CHKERRQ(ierr);
    ierr = IGAAxisInitUniform(axisY, user->Ny, 0.0, user->Ly, user->C); CHKERRQ(ierr);
    if (user->dim == 3) {
        ierr = IGAAxisInitUniform(axisZ, user->Nz, 0.0, user->Lz, user->C); CHKERRQ(ierr);
    }

    // Finalize IGA setup
    ierr = IGASetFromOptions(*iga); CHKERRQ(ierr);
    ierr = IGASetUp(*iga); CHKERRQ(ierr);

    PetscPrintf(PETSC_COMM_WORLD, "IGA setup complete.\n");
    PetscFunctionReturn(0);
}

PetscErrorCode SetupAndSolve(AppCtx *user, IGA iga) {
    PetscErrorCode ierr;
    Mat A;
    Vec b;
    KSP ksp;

    PetscFunctionBegin;

    ierr = IGACreateMat(iga, &A); CHKERRQ(ierr);
    ierr = IGACreateVec(iga, &user->T_sol); CHKERRQ(ierr);
    ierr = IGACreateVec(iga, &b); CHKERRQ(ierr);

    ierr = IGASetFormSystem(iga, AssembleStiffnessMatrix, user); CHKERRQ(ierr);
    ierr = IGAComputeSystem(iga, A, b); CHKERRQ(ierr);

    ierr = ComputeInitialCondition(user->T_sol, user); CHKERRQ(ierr);

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
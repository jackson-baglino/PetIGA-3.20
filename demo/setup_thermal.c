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

PetscErrorCode FormInitialCondition(AppCtx *user) {
    PetscFunctionBegin;
    
    PetscErrorCode ierr;
    IGAElement element;
    IGAPoint point;
    PetscReal dist;
    PetscReal *ice_temp = NULL;
    PetscInt idx = 0;
    PetscBool loadFromFile = PETSC_FALSE;
    FILE *file = NULL;

    // =============================
    // 1. CHECK IF READING FROM FILE
    // =============================
    if (strcmp(user->init_mode, "circle") != 0 && strcmp(user->init_mode, "layered") != 0) {
        PetscPrintf(PETSC_COMM_WORLD, "Assuming init_mode is a file path.\n");
        PetscPrintf(PETSC_COMM_WORLD, "Reading ice field from %s\n", user->init_mode);
        loadFromFile = PETSC_TRUE;
    }

    // =============================
    // 2. LOAD ICE FIELD FROM .dat FILE (IF NEEDED)
    // =============================
    if (loadFromFile) {
        PetscPrintf(PETSC_COMM_WORLD, "Loading ice field from file...\n");

        file = fopen(user->init_mode, "r");
        if (!file) {
            SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_OPEN, "Error opening ice field file: %s", user->init_mode);
        }

        // Determine expected size
        PetscInt num_points = user->Nx * user->Ny * (user->dim == 3 ? user->Nz : 1);
        
        // Allocate memory if needed
        if (!ice_temp) {
            PetscPrintf(PETSC_COMM_WORLD, "Allocating memory for ice field with %d points.\n", num_points);
            ierr = PetscMalloc1(num_points, &ice_temp); CHKERRQ(ierr);
        }

        // Read file values into user->ice
        for (PetscInt i = 0; i < num_points; i++) {
            if (fscanf(file, "%lf", &ice_temp[i]) != 1) {
                fclose(file);
                SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_READ, 
                         "Error reading ice field file at index %D. Expected %D values.", i, num_points);
            }
        }

        // Close the file
        fclose(file);

        // Interpolate ice field to match Gauss points
        PetscInt num_gauss_points = (user->p + 1) * user->dim * num_points;
        PetscReal *ice_interpolated;
        ierr = PetscMalloc1(num_gauss_points, &ice_interpolated); CHKERRQ(ierr);

        PetscInt i;
        PetscReal alpha;
        PetscInt num_points_per_dir = (PetscInt)PetscSqrtReal((PetscReal)num_gauss_points); // Assuming 2D structured grid

        for (i = 0; i < num_gauss_points; i++) {
            PetscInt row = i / num_points_per_dir;
            PetscInt col = i % num_points_per_dir;

            // Compute interpolation weights
            PetscReal x_frac = (PetscReal)col / (num_points_per_dir - 1);
            PetscReal y_frac = (PetscReal)row / (num_points_per_dir - 1);

            // Bilinear interpolation
            PetscInt x0 = PetscMin(col, num_points_per_dir - 2);
            PetscInt x1 = x0 + 1;
            PetscInt y0 = PetscMin(row, num_points_per_dir - 2);
            PetscInt y1 = y0 + 1;

            PetscReal q11 = ice_temp[y0 * num_points_per_dir + x0];
            PetscReal q12 = ice_temp[y1 * num_points_per_dir + x0];
            PetscReal q21 = ice_temp[y0 * num_points_per_dir + x1];
            PetscReal q22 = ice_temp[y1 * num_points_per_dir + x1];

            alpha = (1.0 - x_frac) * (1.0 - y_frac) * q11 +
                    (1.0 - x_frac) * y_frac * q12 +
                    x_frac * (1.0 - y_frac) * q21 +
                    x_frac * y_frac * q22;

            ice_interpolated[i] = alpha;

        }

        // Assign interpolated values to user->ice
        ierr = PetscMalloc1(num_gauss_points, &user->ice); CHKERRQ(ierr);
        ierr = PetscMemcpy(user->ice, ice_interpolated, num_gauss_points * sizeof(PetscReal)); CHKERRQ(ierr);

        // Free temporary memory
        ierr = PetscFree(ice_interpolated); CHKERRQ(ierr);

        // Print the interpolated ice field as a 2D array
        PetscPrintf(PETSC_COMM_WORLD, "Interpolated Ice Field (size: %d x %d):\n", num_points_per_dir, num_points_per_dir);

        for (PetscInt row = 0; row < num_points_per_dir; row++) {
            for (PetscInt col = 0; col < num_points_per_dir; col++) {
                PetscPrintf(PETSC_COMM_WORLD, "%.4f ", user->ice[row * num_points_per_dir + col]);
            }
            PetscPrintf(PETSC_COMM_WORLD, "\n");
        }

        PetscPrintf(PETSC_COMM_WORLD, "✅ Ice field successfully loaded from file.\n");
        PetscPrintf(PETSC_COMM_WORLD, "----- Ice field has %d points.\n", num_gauss_points);
    }

    // =============================
    // 3. COMPUTE ICE FIELD AT NODES IF NOT LOADING FROM FILE
    // =============================
    PetscInt num_points = user->Nx * user->Ny * (user->dim == 3 ? user->Nz : 1);
    if (!loadFromFile) {

        // Allocate memory if needed
        if (!user->ice) {
            ierr = PetscMalloc1(num_points, &user->ice); CHKERRQ(ierr);
        }

        ierr = IGABeginElement(user->iga, &element); CHKERRQ(ierr);
        while (IGANextElement(user->iga, element)) {
            ierr = IGAElementBeginPoint(element, &point); CHKERRQ(ierr);
            while (IGAElementNextPoint(element, point)) {
                // ===========================
                // COMPUTE INITIAL ICE FIELD AT NODES
                // ===========================
                if (strcmp(user->init_mode, "circle") == 0) {
                    PetscReal radius = PetscMin(user->Lx / 6.0, user->Ly / 6.0);

                    // Compute distance from the circle center
                    dist = PetscSqrtReal(SQ(point->mapX[0][0] - user->Lx / 2.0) + 
                                         SQ(point->mapX[0][1] - user->Ly / 2.0)) - radius;
                } 
                else if (strcmp(user->init_mode, "layered") == 0) {
                    // Initialize ice field in a layered manner
                    dist = point->mapX[0][1] - user->Ly / 2.0;
                }

                // Apply phase field function to compute ice phase
                user->ice[idx] = 0.5 - 0.5 * PetscTanhReal(2.0 / user->eps * dist);
                user->ice[idx] = PetscMax(0.0, PetscMin(1.0, user->ice[idx])); // Clamp between [0,1]

                idx++;
            }
            ierr = IGAElementEndPoint(element, &point); CHKERRQ(ierr);
        }
        ierr = IGAEndElement(user->iga, &element); CHKERRQ(ierr);

        PetscPrintf(PETSC_COMM_WORLD, "Ice field has %d points.\n", idx);
        PetscPrintf(PETSC_COMM_WORLD, "There are %d elements total.\n", num_points);
    }

    // =============================
    // 4. SAVE ICE FIELD TO FILE (IF CIRCLE MODE)
    // =============================
    if (!loadFromFile && strcmp(user->init_mode, "circle") == 0) {
        char circle_output_file[256];
        snprintf(circle_output_file, sizeof(circle_output_file), "/Users/jacksonbaglino/PetIGA-3.20/demo/input/Thermal_IO/circle_ice_field.dat");

        FILE *circle_file = fopen(circle_output_file, "w");
        if (!circle_file) {
            SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_OPEN, "Error opening output file for writing ice field.");
        }

        for (PetscInt i = 0; i < idx; i++) {
            fprintf(circle_file, "%.8e\n", user->ice[i]);
        }

        fclose(circle_file);
        PetscPrintf(PETSC_COMM_WORLD, "✅ Ice field written to %s\n", circle_output_file);
    }

    PetscPrintf(PETSC_COMM_WORLD, "✅ Ice field initialized with %d nodal points.\n", idx);

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
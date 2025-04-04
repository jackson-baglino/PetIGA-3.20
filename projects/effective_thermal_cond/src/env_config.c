#include "thermal_solver.h"

PetscErrorCode GetEnvironment(AppCtx *user) {
    PetscFunctionBegin;

    PetscPrintf(PETSC_COMM_WORLD, "Reading simulation parameters...\n\n");

    /* Retrieve environment variables */
    const char *Nx_str = getenv("Nx"), *Ny_str = getenv("Ny"), *Nz_str = getenv("Nz");
    const char *Lx_str = getenv("Lx"), *Ly_str = getenv("Ly"), *Lz_str = getenv("Lz");
    const char *temp_str = getenv("temp"), *grad_temp0X_str = getenv("grad_temp0X");
    const char *grad_temp0Y_str = getenv("grad_temp0Y"), *grad_temp0Z_str = getenv("grad_temp0Z");
    const char *dim_str = getenv("dim"), *outputBinary_str = getenv("OUTPUT_BINARY");
    const char *eps_str = getenv("eps");

    /* New boundary condition variables */
    const char *temp_top_str = getenv("TEMP_TOP");
    const char *flux_bottom_str = getenv("FLUX_BOTTOM");

    /* Check for missing required variables */
    if (!Nx_str || !Ny_str || !Nz_str || !Lx_str || !Ly_str || !Lz_str || 
        !temp_str || !grad_temp0X_str || !grad_temp0Y_str || !grad_temp0Z_str ||
        !dim_str || !outputBinary_str || !eps_str) {
        PetscPrintf(PETSC_COMM_WORLD, "Error: Missing required environment variables.\n");
        PetscFinalize();
        return EXIT_FAILURE;
    }

    /* Convert environment variables */
    char *endptr;
    user->dim = (PetscInt)strtol(dim_str, &endptr, 10);
    user->Nx = (PetscInt)strtol(Nx_str, &endptr, 10);
    user->Ny = (PetscInt)strtol(Ny_str, &endptr, 10);
    
    /* Set Nz based on dimension */
    user->Nz = (user->dim == 3) ? (PetscInt)strtol(Nz_str, &endptr, 10) : 1;

    /* Domain size */
    user->Lx = strtod(Lx_str, &endptr);
    user->Ly = strtod(Ly_str, &endptr);
    user->Lz = strtod(Lz_str, &endptr);

    PetscPrintf(PETSC_COMM_WORLD, "Domain size: Lx = %g, Ly = %g, Lz = %g\n", user->Lx, user->Ly, user->Lz);

    /* Temperature & gradient */
    user->temp0 = strtod(temp_str, &endptr);
    user->grad_temp0[0] = strtod(grad_temp0X_str, &endptr);
    user->grad_temp0[1] = strtod(grad_temp0Y_str, &endptr);
    user->grad_temp0[2] = strtod(grad_temp0Z_str, &endptr);

    /* Interface width */
    user->eps = strtod(eps_str, &endptr);

    /* Output settings */
    user->outputBinary = (PetscBool)strtol(outputBinary_str, &endptr, 10);


    /* ==============================
       Boundary Condition Handling
       ============================== */

    /* Enforce fixed temperature at the top */
    if (temp_top_str) {
        user->T_top = strtod(temp_top_str, &endptr);
        PetscPrintf(PETSC_COMM_WORLD, "Using fixed temperature at top: %g K\n", user->T_top);
    } else {
        PetscPrintf(PETSC_COMM_WORLD, "Error: A fixed temperature at the top boundary (TEMP_TOP) is required!\n");
        PetscFinalize();
        return EXIT_FAILURE;
    }

    /* Enforce prescribed flux at the bottom */
    if (flux_bottom_str) {
        user->q_bottom = strtod(flux_bottom_str, &endptr);
        PetscPrintf(PETSC_COMM_WORLD, "Using prescribed flux at bottom: %g W/m²•K\n", user->q_bottom);
    } else {
        PetscPrintf(PETSC_COMM_WORLD, "Error: A prescribed flux at the bottom boundary (FLUX_BOTTOM) is required!\n");
        PetscFinalize();
        return EXIT_FAILURE;
    }

    PetscFunctionReturn(0);
}
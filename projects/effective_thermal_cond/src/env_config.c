#include "env_config.h"

// Function declarations (static = private to this file)
static PetscErrorCode CheckRequiredEnvVars(void);
static void ParseDomainAndMesh(AppCtx *user, char *endptr);
static void ParseTemperatureSettings(AppCtx *user, char *endptr);
static void ParseOutputSettings(AppCtx *user, char *endptr);
static PetscErrorCode ParseBoundaryConditions(AppCtx *user, char *endptr);

/*------------------------------------------------------------------------------
  Function: InitializeUserContext
  Purpose : Initialize AppCtx with default material properties and settings.
------------------------------------------------------------------------------*/
void InitializeUserContext(AppCtx *user) {
    // Material properties
    user->thcond_ice = 2.29;       // W/m·K
    user->thcond_air = 0.02;       // W/m·K
    user->cp_ice     = 1.96e3;     // J/kg·K
    user->cp_air     = 1.044e3;    // J/kg·K
    user->rho_ice    = 919.0;      // kg/m³
    user->rho_air    = 1.341;      // kg/m³

    // Discretization order and continuity
    user->p = 1; // Polynomial order
    user->C = 0; // Inter-element continuity

    PetscPrintf(PETSC_COMM_WORLD, "User context initialized.\n\n");
}

/*------------------------------------------------------------------------------
  Function: GetEnvironment
  Purpose : Read and parse environment variables into AppCtx.
------------------------------------------------------------------------------*/
PetscErrorCode GetEnvironment(AppCtx *user) {
    PetscFunctionBegin;

    PetscPrintf(PETSC_COMM_WORLD, "Reading simulation parameters from environment variables...\n\n");

    // Step 1: Check for required environment variables
    PetscErrorCode ierr = CheckRequiredEnvVars();
    if (ierr) PetscFunctionReturn(ierr);

    // Step 2: Parse all values into AppCtx
    char *endptr = NULL;
    ParseDomainAndMesh(user, endptr);
    ParseTemperatureSettings(user, endptr);
    ParseOutputSettings(user, endptr);

    // Step 3: Parse boundary conditions
    ierr = ParseBoundaryConditions(user, endptr); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

/*------------------------------------------------------------------------------
  Function: CheckRequiredEnvVars
  Purpose : Ensure all required environment variables are defined.
------------------------------------------------------------------------------*/
static PetscErrorCode CheckRequiredEnvVars(void) {
    const char *required_vars[] = {
        "Nx", "Ny", "Nz", "Lx", "Ly", "Lz", "temp",
        "grad_temp0X", "grad_temp0Y", "grad_temp0Z",
        "dim", "OUTPUT_BINARY", "eps", "TEMP_TOP", "FLUX_BOTTOM"
    };

    int num_vars = sizeof(required_vars) / sizeof(required_vars[0]);
    for (int i = 0; i < num_vars; i++) {
        if (!getenv(required_vars[i])) {
            PetscPrintf(PETSC_COMM_WORLD, "❌ Error: Missing required environment variable: %s\n", required_vars[i]);
            PetscFinalize();
            return EXIT_FAILURE;
        }
    }

    return 0;
}

/*------------------------------------------------------------------------------
  Function: ParseDomainAndMesh
  Purpose : Set domain size, mesh resolution, and problem dimension.
------------------------------------------------------------------------------*/
static void ParseDomainAndMesh(AppCtx *user, char *endptr) {
    user->dim = (PetscInt)strtol(getenv("dim"), &endptr, 10);
    user->Nx  = (PetscInt)strtol(getenv("Nx"),  &endptr, 10);
    user->Ny  = (PetscInt)strtol(getenv("Ny"),  &endptr, 10);
    user->Nz  = (user->dim == 3) ? (PetscInt)strtol(getenv("Nz"), &endptr, 10) : 1;

    user->Lx  = strtod(getenv("Lx"), &endptr);
    user->Ly  = strtod(getenv("Ly"), &endptr);
    user->Lz  = strtod(getenv("Lz"), &endptr);

    PetscPrintf(PETSC_COMM_WORLD, "Domain: Lx = %g, Ly = %g, Lz = %g\n\n", user->Lx, user->Ly, user->Lz);
}

/*------------------------------------------------------------------------------
  Function: ParseTemperatureSettings
  Purpose : Set initial temperature and its gradient field.
------------------------------------------------------------------------------*/
static void ParseTemperatureSettings(AppCtx *user, char *endptr) {
    user->temp0         = strtod(getenv("temp"), &endptr);
    user->grad_temp0[0] = strtod(getenv("grad_temp0X"), &endptr);
    user->grad_temp0[1] = strtod(getenv("grad_temp0Y"), &endptr);
    user->grad_temp0[2] = strtod(getenv("grad_temp0Z"), &endptr);
    user->eps           = strtod(getenv("eps"), &endptr);
}

/*------------------------------------------------------------------------------
  Function: ParseOutputSettings
  Purpose : Set output format and settings.
------------------------------------------------------------------------------*/
static void ParseOutputSettings(AppCtx *user, char *endptr) {
    user->outputBinary = (PetscBool)strtol(getenv("OUTPUT_BINARY"), &endptr, 10);
}

/*------------------------------------------------------------------------------
  Function: ParseBoundaryConditions
  Purpose : Set boundary condition values for top and bottom of the domain.
------------------------------------------------------------------------------*/
static PetscErrorCode ParseBoundaryConditions(AppCtx *user, char *endptr) {
    const char *temp_top_str    = getenv("TEMP_TOP");
    const char *flux_bottom_str = getenv("FLUX_BOTTOM");

    if (temp_top_str) {
        user->T_top = strtod(temp_top_str, &endptr);
        PetscPrintf(PETSC_COMM_WORLD, "Top boundary temperature: %g K\n\n", user->T_top);
    } else {
        PetscPrintf(PETSC_COMM_WORLD, "❌ Error: TEMP_TOP is required.\n");
        PetscFinalize();
        return EXIT_FAILURE;
    }

    if (flux_bottom_str) {
        user->q_bottom = strtod(flux_bottom_str, &endptr);
        PetscPrintf(PETSC_COMM_WORLD, "Bottom boundary flux: %g W/m²·K\n\n", user->q_bottom);
    } else {
        PetscPrintf(PETSC_COMM_WORLD, "❌ Error: FLUX_BOTTOM is required.\n");
        PetscFinalize();
        return EXIT_FAILURE;
    }

    return 0;
}
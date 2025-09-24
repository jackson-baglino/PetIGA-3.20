//------------------------------------------------------------------------------
// 0) Includes, macros, and system headers
//------------------------------------------------------------------------------
#include "petiga.h"
#include <petsc/private/tsimpl.h>
#include <petsc/private/snesimpl.h>
#include <petsctime.h>
#include <petscsys.h>
#include <petscviewer.h>
#include <math.h>
#include <mpi.h>
#include <sys/stat.h>
#include <regex.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>

#define SQ(x) ((x)*(x))
#define CU(x) ((x)*(x)*(x))

//------------------------------------------------------------------------------
// 1) AppCtx definition
//------------------------------------------------------------------------------
typedef struct {
  IGA iga;  // Isogeometric Analysis (IGA) object for managing the finite element discretization
  IGA iga_input;  // IGA object for reading input data

  // **Material properties**
  PetscReal thcond_ice;   // Thermal conductivity of ice (W/mK)
  PetscReal thcond_air;   // Thermal conductivity of air (W/mK)
  PetscReal cp_ice;       // Specific heat capacity of ice (J/kgK)
  PetscReal cp_air;      // Specific heat capacity of air (J/kgK)
  PetscReal rho_ice;     // Density of ice (kg/m³)
  PetscReal rho_air;     // Density of air (kg/m³)

  // **Initial conditions**
  PetscReal temp0;            // Initial temperature (K)
  PetscReal grad_temp0[3];    // Initial temperature gradient (x, y, z components)
  PetscReal *ice;             // Ice phase variable
  Vec T_sol;                  // Solution vector for temperature


  // **Domain size and mesh resolution**
  PetscInt dim;         // 2D or 3D
  PetscReal Lx, Ly, Lz;  // Domain size
  PetscInt Nx, Ny, Nz;  // Mesh resolution
  PetscReal eps;        // Interface width parameter for phase field method
  PetscInt  p;         // Polynomial orders for basis functions
  PetscInt  C;         // Global continuity order

  // **Dirichlet BC (Fixed Temperature at Top & Bottom)**
  PetscReal T_top;    // Temperature at y = Ly

  // **Neumann BC (Constant flux at top or bottom)**
  PetscReal q_bottom; // Heat flux at y = 0
  PetscBool useFluxBottom; // If true, apply flux at bottom boundary

  // **Input options**
  char init_mode[256]; // Mode for initializing the ice field 
  char init_dir[256];  // Directory for initial condition files

  // **Output options**
  PetscBool outputBinary; // Flag for binary output
  PetscInt sol_index; // -1 means loop over all, otherwise specific index
  const char *output_dir; // Directory for output files

} AppCtx;

//------------------------------------------------------------------------------
// 2) Low-level helpers & memory management
//------------------------------------------------------------------------------
/*--------------------------------------------------------------------------------------------------
  Function: AllocateAppCtxFields
  Allocates a PetscScalar field with size based on user dimensions and polynomial order.
  This function supports 2D and 3D cases and assumes tensor-product basis structure.
--------------------------------------------------------------------------------------------------*/
PetscErrorCode AllocateAppCtxFields(IGA iga, AppCtx *user, PetscScalar **field) {
  PetscFunctionBegin;

  PetscInt nmb = iga->elem_width[0] * iga->elem_width[1] * SQ(user->p + 1);
  if (user->dim == 3) {
    nmb = iga->elem_width[0] * iga->elem_width[1] * iga->elem_width[2] * CU(user->p + 1);
  }

  PetscCall(PetscMalloc1(nmb, field));
  PetscCall(PetscMemzero(*field, sizeof(PetscScalar) * nmb));

  PetscFunctionReturn(0);
}

//------------------------------------------------------------------------------
// 3) Material property functions
//------------------------------------------------------------------------------
/* Function: ThermalCond  
   Computes thermal conductivity and its derivative w.r.t. ice content. */
void ThermalCond(AppCtx *user, PetscScalar ice, PetscScalar *cond, PetscScalar *dcond_ice)
{
    PetscReal dice = 1.0, dair = 1.0;
    PetscReal air = 1.0 - ice;
    if (ice < 0.0) { ice = 0.0; dice = 0.0; }
    if (air < 0.0) { air = 0.0; dair = 0.0; }

    PetscReal cond_ice = user->thcond_ice;
    PetscReal cond_air = user->thcond_air;

    if (cond)       (*cond)      = ice * cond_ice + air * cond_air;
    if (dcond_ice)  (*dcond_ice) = cond_ice * dice - cond_air * dair;

    // if (*cond > cond_ice - 0.1) {
    //   PetscPrintf(PETSC_COMM_WORLD, "ice = %g, thcond = %g \n", ice, *cond);
    // }

    return;
}
   
//------------------------------------------------------------------------------
// 4) Field evaluation & I/O helpers
//------------------------------------------------------------------------------
/*--------------------------------------------------------------------------------------------------
  Function: EvaluateFieldAtGaussPoints
  Purpose: Evaluates the scalar field from a PETSc vector at Gauss points and stores it in 
           user->ice array.
  Inputs: 
    - user: Application context containing simulation parameters and output field pointer.
    - iga:  IGA object describing geometry and discretization.
    - vec_phase: Vector containing the solution to be evaluated.
--------------------------------------------------------------------------------------------------*/
PetscErrorCode EvaluateFieldAtGaussPoints(AppCtx *user, IGA iga, Vec vec_phase) {
  PetscErrorCode ierr;
  PetscFunctionBegin;

  IGAElement            element;
  IGAPoint              point;

  PetscInt              idx = 0;
  Vec                   localU;
  const PetscScalar     *arrayU;
  PetscScalar           *U;

  ierr = IGAGetLocalVecArray(iga, vec_phase, &localU, &arrayU); CHKERRQ(ierr);

  ierr = IGABeginElement(iga, &element); CHKERRQ(ierr);
  while (IGANextElement(iga, element)) {
    ierr = IGAElementGetValues(element, arrayU, &U); CHKERRQ(ierr);
    ierr = IGAElementBeginPoint(element, &point); CHKERRQ(ierr);

    while (IGAElementNextPoint(element, point)) {
      PetscScalar sol[3];
      IGAPointFormValue(point, U, &sol[0]);
      PetscScalar phi = sol[0];
      user->ice[idx++] = phi;
    }

    ierr = IGAElementEndPoint(element, &point); CHKERRQ(ierr);
  }
  
  PetscInt nmb_expected = iga->elem_width[0] * iga->elem_width[1] * SQ(user->p + 1);
  if (idx != nmb_expected) {
    PetscPrintf(PETSC_COMM_WORLD, "⚠️ Warning: Assigned %d Gauss point values, expected %d\n", idx, nmb_expected);
  }

  ierr = IGAEndElement(iga, &element); CHKERRQ(ierr);

  ierr = IGARestoreLocalVecArray(iga, vec_phase, &localU, &arrayU); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*--------------------------------------------------------------------------------------------------
  Function: ReadSolutionVec
  Purpose: Reads an IGA object and solution vector from file, evaluates the scalar field 
           at Gauss points, and stores it in user->ice.
  Inputs:
    - iga_file: Path to the IGA file (.iga).
    - vec_file: Path to the vector file (.vec).
    - user: Application context (to store loaded field and IGA).
  Outputs:
    - iga_out: Pointer to loaded IGA object.
    - vec_out: Pointer to loaded solution vector (before destruction).
--------------------------------------------------------------------------------------------------*/
PetscErrorCode ReadSolutionVec(const char *iga_file, const char *vec_file, 
                               IGA *iga_out, AppCtx *user) {
  PetscErrorCode ierr;
  PetscFunctionBegin;

  IGA iga_input;  // IGA object read from file (external to simulation)

  // Create the IGA object (NOTE: WE CAN PROBABLY REPLACE ALL THIS WITH A PREDEFINED FUNCTION)
  ierr = IGACreate(PETSC_COMM_WORLD, &iga_input); CHKERRQ(ierr);
  ierr = IGASetDim(iga_input, user->dim); CHKERRQ(ierr);
  ierr = IGASetDof(iga_input, 3); CHKERRQ(ierr);

  // Set up the axes
  IGAAxis axis0, axis1, axis2;
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

  // Set the IGA from options and set up
  ierr = IGASetFromOptions(iga_input); CHKERRQ(ierr);
  ierr = IGASetUp(iga_input); CHKERRQ(ierr);

  // Set the input IGA object
  user->iga_input = iga_input;

  // Allocate memory fro user->ice
  ierr = AllocateAppCtxFields(iga_input, user, &user->ice); CHKERRQ(ierr);

  // Create the solution vector
  Vec ice_phase;
  ierr = IGACreateVec(iga_input, &ice_phase); CHKERRQ(ierr);
  ierr = IGAReadVec(iga_input, ice_phase, vec_file); CHKERRQ(ierr);

  // Assign Gauss points via a function (to be created)
  ierr = EvaluateFieldAtGaussPoints(user, iga_input, ice_phase); CHKERRQ(ierr);

  // Clean up
  ierr = VecDestroy(&ice_phase); CHKERRQ(ierr);
  ierr = IGADestroy(&iga_input); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

//------------------------------------------------------------------------------
// 5) Initial-condition generators
//------------------------------------------------------------------------------
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
PetscErrorCode ComputeCircleIceField(AppCtx *user) {
  PetscErrorCode ierr;
  PetscFunctionBegin;

  IGAElement element;
  IGAPoint point;
  PetscReal dist;
  PetscInt idx = 0;

  /* Four grains, one in each quadrant */
  enum {NUM_GRAINS = 1};
  PetscInt  grainID;
  // const PetscReal centX[NUM_GRAINS] = {user->Lx/2, user->Lx/2};
  // const PetscReal centY[NUM_GRAINS] = {0.0,        user->Ly};
  // const PetscReal centX[NUM_GRAINS] = {0.0,        0.0,        user->Lx, user->Lx};
  // const PetscReal centY[NUM_GRAINS] = {0.0,        user->Ly,   0.0,      user->Ly};
  const PetscReal centX[NUM_GRAINS] = {user->Lx / 2.0};
  const PetscReal centY[NUM_GRAINS] = {user->Ly / 2.0};
  PetscReal radius = PetscMin(user->Lx / 16.0, user->Ly / 16.0);
  // PetscReal radius = user->Ly*0.6;/

  // Allocate memory for user->ice (if not already allocated)
  ierr = AllocateAppCtxFields(user->iga, user, &user->ice); CHKERRQ(ierr);

  ierr = IGABeginElement(user->iga, &element); CHKERRQ(ierr);
  while (IGANextElement(user->iga, element)) {
      ierr = IGAElementBeginPoint(element, &point); CHKERRQ(ierr);
      while (IGAElementNextPoint(element, point)) {
        idx = point->index + point->count * point->parent->index;
        PetscReal ice = 0.0;
        for (grainID = 0; grainID < NUM_GRAINS; grainID++) {
          dist = PetscSqrtReal(SQ(point->mapX[0][0] - centX[grainID]) +
                                SQ(point->mapX[0][1] - centY[grainID])) - radius;
          ice += 0.5 - 0.5 * PetscTanhReal(0.5 / user->eps * dist);
        }
        user->ice[idx] = PetscMax(0.0, PetscMin(1.0, ice));
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
PetscErrorCode ComputeLayeredIceField(AppCtx *user) {
  PetscErrorCode ierr;
  PetscFunctionBegin;

  IGAElement element;
  IGAPoint point;
  PetscReal dist;
  PetscInt indGP;
  PetscInt counts = 0;

  // Allocate memory for user->ice (if not already allocated)
  ierr = AllocateAppCtxFields(user->iga, user, &user->ice); CHKERRQ(ierr);

  PetscInt nmb = user->iga->elem_width[0] * user->iga->elem_width[1] * SQ(user->p + 1);

  ierr = IGABeginElement(user->iga, &element); CHKERRQ(ierr);
  while (IGANextElement(user->iga, element)) {
      ierr = IGAElementBeginPoint(element, &point); CHKERRQ(ierr);
      while (IGAElementNextPoint(element, point)) {
        indGP = point->index + point->count * point->parent->index;

        // Find distance from the center of the domain--above midline is air, below is ice
        dist = point->mapX[0][1] - user->Ly / 2.0;

        user->ice[indGP] = 0.5 - 0.5 * PetscTanhReal(0.5 / user->eps * dist);
        user->ice[indGP] = PetscMax(0.0, PetscMin(1.0, user->ice[indGP]));   

        counts++;

      //   user->ice[indGP] = 1.0;
      }
      ierr = IGAElementEndPoint(element, &point); CHKERRQ(ierr);
  }
  ierr = IGAEndElement(user->iga, &element); CHKERRQ(ierr);

  if (counts != nmb) {
    PetscPrintf(PETSC_COMM_WORLD, "⚠️ Warning: Assigned %d Gauss point values, expected %d\n", counts, nmb);
  }

  PetscPrintf(PETSC_COMM_WORLD, "--- Layered mode. ---\n\n");
  PetscFunctionReturn(0);
}

/*-----------------------------------------------------------
  Function: FormInitialCondition
  Purpose : Initializes the ice field in user->ice according to the mode specified in user->init_mode.
         - "circle"  : Calls ComputeCircleIceField to generate a circular grain.
         - "layered" : Calls ComputeLayeredIceField for a vertical layer profile.
         - Otherwise : Assumes file-based initialization, reads solution from files in user->init_dir.
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
   // Generate a circular grain ice field
   ierr = ComputeCircleIceField(user); CHKERRQ(ierr);
  } else if (strcmp(user->init_mode, "layered") == 0) {
   // Generate a layered (vertical) ice field
   ierr = ComputeLayeredIceField(user); CHKERRQ(ierr);
  } else {
   // File-based initialization: read IGA and solution vector from files
   char iga_file[PETSC_MAX_PATH_LEN];
   char sol_file[PETSC_MAX_PATH_LEN];

   // Construct file paths for IGA and solution vector
   snprintf(iga_file, sizeof(iga_file), "%s/igasol.dat", user->init_dir);
   snprintf(sol_file, sizeof(sol_file), "%s/sol_%05d.dat", user->init_dir, user->sol_index);

   PetscPrintf(PETSC_COMM_WORLD, "Reading IGA from file: %s\n", iga_file);
   PetscPrintf(PETSC_COMM_WORLD, "Reading ice field from file: %s\n", sol_file);

   // Read the solution vector and evaluate the field at Gauss points
   ierr = ReadSolutionVec(iga_file, sol_file, &user->iga_input, user); CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

//------------------------------------------------------------------------------
// 6) Environment & option parsing
//------------------------------------------------------------------------------
/*------------------------------------------------------------------------------
  Function: CheckRequiredEnvVars
  Purpose : Ensure all required environment variables are defined.
------------------------------------------------------------------------------*/
PetscErrorCode CheckRequiredEnvVars(void) {
  const char *required_vars[] = {
      "Nx", "Ny", "Nz", "Lx", "Ly", "Lz",
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
void ParseDomainAndMesh(AppCtx *user, char *endptr) {
  user->dim = (PetscInt)strtol(getenv("dim"), &endptr, 10);
  user->Nx  = (PetscInt)strtol(getenv("Nx"),  &endptr, 10);
  user->Ny  = (PetscInt)strtol(getenv("Ny"),  &endptr, 10);
  user->Nz  = (user->dim == 3) ? (PetscInt)strtol(getenv("Nz"), &endptr, 10) : 1;

  user->Lx  = strtod(getenv("Lx"), &endptr);
  user->Ly  = strtod(getenv("Ly"), &endptr);
  user->Lz  = strtod(getenv("Lz"), &endptr);
}

/*------------------------------------------------------------------------------
Function: ParseTemperatureSettings
Purpose : Set initial temperature and its gradient field.
------------------------------------------------------------------------------*/
void ParseTemperatureSettings(AppCtx *user, char *endptr) {
  user->eps           = strtod(getenv("eps"), &endptr);
}

/*------------------------------------------------------------------------------
Function: ParseOutputSettings
Purpose : Set output format and settings.
------------------------------------------------------------------------------*/
void ParseOutputSettings(AppCtx *user, char *endptr) {
  user->outputBinary = (PetscBool)strtol(getenv("OUTPUT_BINARY"), &endptr, 10);
}

/*------------------------------------------------------------------------------
Function: ParseBoundaryConditions
Purpose : Set boundary condition values for top and bottom of the domain.
------------------------------------------------------------------------------*/
PetscErrorCode ParseBoundaryConditions(AppCtx *user, char *endptr) {
  const char *temp_top_str    = getenv("TEMP_TOP");
  const char *flux_bottom_str = getenv("FLUX_BOTTOM");

  if (temp_top_str) {
      user->T_top = strtod(temp_top_str, &endptr);
  } else {
      PetscPrintf(PETSC_COMM_WORLD, "❌ Error: TEMP_TOP is required.\n");
      PetscFinalize();
      return EXIT_FAILURE;
  }

  if (flux_bottom_str) {
      user->q_bottom = strtod(flux_bottom_str, &endptr);
  } else {
      PetscPrintf(PETSC_COMM_WORLD, "❌ Error: FLUX_BOTTOM is required.\n");
      PetscFinalize();
      return EXIT_FAILURE;
  }

  return 0;
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
  // ParseDomainAndMesh(user, endptr);
  // ParseTemperatureSettings(user, endptr);
  // ParseOutputSettings(user, endptr);

  user->dim = (PetscInt)strtol(getenv("dim"), &endptr, 10);
  user->Nx  = (PetscInt)strtol(getenv("Nx"),  &endptr, 10);
  user->Ny  = (PetscInt)strtol(getenv("Ny"),  &endptr, 10);
  user->Nz  = (user->dim == 3) ? (PetscInt)strtol(getenv("Nz"), &endptr, 10) : 1;

  user->Lx  = strtod(getenv("Lx"), &endptr);
  user->Ly  = strtod(getenv("Ly"), &endptr);
  user->Lz  = strtod(getenv("Lz"), &endptr);

  user->eps           = strtod(getenv("eps"), &endptr);

  user->outputBinary = (PetscBool)strtol(getenv("OUTPUT_BINARY"), &endptr, 10);
  user->sol_index = (PetscInt)strtol(getenv("SOL_INDEX"), &endptr, 10);

  const char *temp_top_str    = getenv("TEMP_TOP");
  const char *flux_bottom_str = getenv("FLUX_BOTTOM");

  const char *output_dir_str = getenv("OUTPUT_DIR");
  user->output_dir = output_dir_str;

  PetscPrintf(PETSC_COMM_WORLD, "Output directory: %s\n", user->output_dir);

  if (temp_top_str) {
      user->T_top = strtod(temp_top_str, &endptr);
  } else {
      PetscPrintf(PETSC_COMM_WORLD, "❌ Error: TEMP_TOP is required.\n");
      PetscFinalize();
      return EXIT_FAILURE;
  }

  if (flux_bottom_str) {
      user->q_bottom = strtod(flux_bottom_str, &endptr);
  } else {
      PetscPrintf(PETSC_COMM_WORLD, "❌ Error: FLUX_BOTTOM is required.\n");
      PetscFinalize();
      return EXIT_FAILURE;
  }


  // // Step 3: Parse boundary conditions
  // ierr = ParseBoundaryConditions(user, endptr); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

//------------------------------------------------------------------------------
// 7) PETSc/IGA setup & solver
//------------------------------------------------------------------------------
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
}

/*------------------------------------------------------------------------------
  Function: SetupIGA
  Purpose:  Sets up the IGA object by defining the problem domain, degrees of 
            freedom, field names, and the uniform grid.
------------------------------------------------------------------------------*/
PetscErrorCode SetupIGA(AppCtx *user, IGA *iga) {
  PetscErrorCode ierr;
  PetscFunctionBegin;

  ierr = IGACreate(PETSC_COMM_WORLD, iga); CHKERRQ(ierr);
  ierr = IGASetDim(*iga, user->dim); CHKERRQ(ierr);
  ierr = IGASetDof(*iga, user->dim); CHKERRQ(ierr); /* To solve for vector field */
  ierr = IGASetFieldName(*iga, 0, "temperature_gradients"); CHKERRQ(ierr);

  /* Set up uniform grid along each axis */
  IGAAxis axisX, axisY, axisZ;
  ierr = IGAGetAxis(*iga, 0, &axisX); CHKERRQ(ierr);
  ierr = IGAAxisSetDegree(axisX, user->p); CHKERRQ(ierr);
  ierr = IGAAxisSetPeriodic(axisX, PETSC_TRUE); CHKERRQ(ierr);

  ierr = IGAGetAxis(*iga, 1, &axisY); CHKERRQ(ierr);
  ierr = IGAAxisSetDegree(axisY, user->p); CHKERRQ(ierr);
  ierr = IGAAxisSetPeriodic(axisY, PETSC_TRUE); CHKERRQ(ierr);

  if (user->dim == 3) {
      ierr = IGAGetAxis(*iga, 2, &axisZ); CHKERRQ(ierr);
      ierr = IGAAxisSetDegree(axisZ, user->p); CHKERRQ(ierr);
      ierr = IGAAxisSetPeriodic(axisZ, PETSC_TRUE); CHKERRQ(ierr);
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

  PetscPrintf(PETSC_COMM_WORLD, "IGA setup complete.\n\n");

  /* Write IGA object to file */
  char filename[256];
  sprintf(filename, "%s/igaice.dat", user->output_dir);
  ierr=IGAWrite(*iga,filename);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*------------------------------------------------------------------------------
  Function: AssembleStiffnessMatrix
  Purpose: Assemble FEM stiffness matrix and apply bottom flux BC.
------------------------------------------------------------------------------*/
PetscErrorCode AssembleStiffnessMatrix(IGAPoint pnt,
                                       PetscScalar *K,   /* size (nen*dim)^2 */
                                       PetscScalar *F,   /* size  nen*dim    */
                                       void *ctx)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;

  // Retrieve user context and define variables
  AppCtx *user = (AppCtx *)ctx;
  PetscInt nen = pnt->nen, dim = user->dim;
  PetscInt ndof = nen * dim;
  PetscInt i, j, m, g;

  const PetscScalar *N0;
  const PetscReal (*N1)[dim];

  // Retrieve basis functions
  ierr = IGAPointGetShapeFuns(pnt, 0, &N0);                       CHKERRQ(ierr);
  ierr = IGAPointGetShapeFuns(pnt, 1, (const PetscReal**)&N1);    CHKERRQ(ierr);

  // Get the index of the Gauss point and thermal conductivity
  PetscInt indGP = pnt->index + pnt->count * pnt->parent->index;
  PetscReal thcond;
  ThermalCond(user, user->ice[indGP], &thcond, NULL);

  /* zero local matrix/vector */
  for (PetscInt q = 0; q < ndof * ndof; ++q) K[q] = 0.0;
  for (PetscInt q = 0; q < ndof; ++q) F[q] = 0.0;

  // Initialize stiffness matrix and force vector
  for (i = 0; i < nen; i++) {
    for (j = 0; j < nen; j++) {
      for (m = 0; m < dim; m++) {
        PetscScalar stab = 0.0;
        for (g = 0; g < dim; g++) {
          stab += N1[i][g] * N1[j][g] * thcond;
        }
        PetscInt row = i*dim + m;
        PetscInt col = j*dim + m;
        K[row*ndof + col] += stab;
      }
    }
  }
  
  /* Initialize laod vector */
  for (i = 0; i < nen; i++) {
    for (m = 0; m < dim; m++) {
      F[i*dim + m] += -thcond * N1[i][m];
    }
  }
  PetscFunctionReturn(0);
}

/*------------------------------------------------------------------------------
   Function: ApplyBoundaryConditions
   Purpose : For periodic problems, there are no Dirichlet/Neumann terms to set.
------------------------------------------------------------------------------*/
PetscErrorCode ApplyBoundaryConditions(IGA iga, AppCtx *user)
{
  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD,
      "Periodic boundary conditions activated in all %d directions — no "
      "Dirichlet or Neumann BCs applied.\n\n", user->dim);
  PetscFunctionReturn(0);
}

/*------------------------------------------------------------------------------
  Function: InitialCondition
  Purpose: Set initial temperature using linear profile from BCs.
------------------------------------------------------------------------------*/
PetscErrorCode ComputeInitialCondition(Vec T, AppCtx *user)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = VecSet(T, 0.0); CHKERRQ(ierr);    /* set all dofs to zero */
  PetscFunctionReturn(0);
}

/*-----------------------------------------------------------
  Function: SetupAndSolve
  Assembles the system matrix and RHS vector, sets the initial
  condition for temperature, and solves the linear system using KSP.
  Returns: 0 on success, or a PETSc error code.
-----------------------------------------------------------*/
PetscErrorCode SetupAndSolve(AppCtx *user, IGA iga) {
  PetscErrorCode ierr;
  Mat A;
  Vec b;
  KSP ksp;
  PC pc;

  PetscInt dim = user->dim;
  PetscMPIInt size;

  PetscFunctionBegin;

  MPI_Comm_size(PETSC_COMM_WORLD, &size);

  // Create system matrix and vectors
  ierr = IGACreateMat(iga, &A); CHKERRQ(ierr);
  ierr = IGACreateVec(iga, &user->T_sol); CHKERRQ(ierr);
  ierr = IGACreateVec(iga, &b); CHKERRQ(ierr);

  // Assemble system matrix and RHS vector
  ierr = IGASetFormSystem(iga, AssembleStiffnessMatrix, user); CHKERRQ(ierr);
  ierr = IGAComputeSystem(iga, A, b); CHKERRQ(ierr);

  /* ---------- constant-mode null-space (one per component) ---------- */
  {
    PetscInt rows;
    Vec      basis[3];                  /* up to 3 constant vectors     */
    MatNullSpace nsp;

    MatGetSize(A,&rows,NULL);
    PetscInt nnodes = rows / dim;       /* total nodes (= rows/dim)     */

    /* build orthonormal constant vectors */
    for (PetscInt m=0; m<dim; ++m) {
      VecCreate(PETSC_COMM_WORLD,&basis[m]);
      VecSetSizes(basis[m],PETSC_DECIDE,rows);
      VecSetFromOptions(basis[m]);
      VecSet(basis[m],0.0);

      PetscScalar *y;
      VecGetArray(basis[m],&y);
      for (PetscInt n=0; n<nnodes; ++n) y[n*dim + m] = 1.0;
      VecRestoreArray(basis[m],&y);
      VecAssemblyBegin(basis[m]); VecAssemblyEnd(basis[m]);
      VecNormalize(basis[m],NULL);
    }

    MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_FALSE,dim,basis,&nsp);
    MatSetNullSpace(A,nsp);
    MatNullSpaceDestroy(&nsp);
    for (PetscInt m=0; m<dim; ++m) VecDestroy(&basis[m]);
  }
  /* ---------- end null-space block ---------- */

  /* Set initial condition for temperature field */
  ierr = ComputeInitialCondition(user->T_sol, user); CHKERRQ(ierr);

  // Set initial condition for temperature field
  ierr = ComputeInitialCondition(user->T_sol, user); CHKERRQ(ierr);
  // VecZeroEntries(user->T_sol);

  // Solve the linear system using KSP
  ierr = IGACreateKSP(iga, &ksp); CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp, A, A); CHKERRQ(ierr);

  // Configure preconditioner
  ierr = KSPGetPC(ksp, &pc); CHKERRQ(ierr);
  ierr = PCSetType(pc, PCLU); CHKERRQ(ierr); // Use LU preconditioner

  #if defined(PETSC_HAVE_MUMPS)
    if (size > 1) PetscCall(PCFactorSetMatSolverType(pc, MATSOLVERMUMPS));
  #endif

  ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);
  // Set solver tolerances
  ierr = KSPSetTolerances(ksp, PETSC_SMALL, PETSC_SMALL, PETSC_DEFAULT, 4000); CHKERRQ(ierr); // rtol, abstol, dtol, maxits
  ierr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE); CHKERRQ(ierr);

  // Start timing
  PetscLogDouble t1, t2;
  ierr = PetscTime(&t1); CHKERRQ(ierr);

  // Solve the system
  ierr = KSPSolve(ksp, b, user->T_sol); CHKERRQ(ierr);

  // End timing
  ierr = PetscTime(&t2); CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD, "KSP solve completed in %.2f seconds.\n\n", t2 - t1);

  /* ----- enforce zero–mean per component ----- */
  {
    PetscInt    dim = user->dim;
    PetscInt    gsize;
    const PetscScalar *x;
    PetscScalar       *y;
    PetscScalar        mean[3] = {0,0,0};

    /* global size and #nodes */
    ierr = VecGetSize(user->T_sol,&gsize); CHKERRQ(ierr);
    PetscInt nnodes = gsize / dim;

    /* read-only pass: accumulate sums */
    ierr = VecGetArrayRead(user->T_sol,&x); CHKERRQ(ierr);
    for (PetscInt n=0; n<nnodes; ++n)
      for (PetscInt m=0; m<dim; ++m)
        mean[m] += x[n*dim + m];
    ierr = VecRestoreArrayRead(user->T_sol,&x); CHKERRQ(ierr);

    for (PetscInt m=0; m<dim; ++m) mean[m] /= (PetscScalar)nnodes;

    /* writable pass: subtract means */
    ierr = VecGetArray(user->T_sol,&y); CHKERRQ(ierr);
    for (PetscInt n=0; n<nnodes; ++n)
      for (PetscInt m=0; m<dim; ++m)
        y[n*dim + m] -= mean[m];
    ierr = VecRestoreArray(user->T_sol,&y); CHKERRQ(ierr);
  }

  PetscPrintf(PETSC_COMM_WORLD, "Zero-mean enforced per component.\n\n");

  // Clean up
  ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
  ierr = MatDestroy(&A); CHKERRQ(ierr);
  ierr = VecDestroy(&b); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}



//------------------------------------------------------------------------------
// 8) Output routines
//------------------------------------------------------------------------------
/*--------------------------------------------------------------------------------------------------
  Function: WriteBinaryFile
  Purpose: Writes a PETSc vector to a binary file using a PetscViewer.
  Inputs:
    - field: PETSc vector to write.
    - filename: Output file name.
--------------------------------------------------------------------------------------------------*/
PetscErrorCode WriteBinaryFile(Vec field, const char *filename) {
  PetscErrorCode ierr;
  PetscViewer viewer;
  PetscBool fileExists;

  PetscFunctionBegin;
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename, FILE_MODE_WRITE, &viewer); CHKERRQ(ierr);

  ierr = VecView(field, viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

  ierr = PetscTestFile(filename, 'r', &fileExists); CHKERRQ(ierr);
  if (!fileExists) {
    PetscPrintf(PETSC_COMM_WORLD, "Warning: Binary output file %s was not created!\n", filename);
  }
  
  PetscFunctionReturn(0);
}

/*--------------------------------------------------------------------------------------------------
  Function: WriteIceFieldToFile
  Purpose: Writes the contents of user->ice to a human-readable .dat text file.
  Inputs:
    - filename: Path to output .dat file.
    - user: Application context containing the ice field.
--------------------------------------------------------------------------------------------------*/
PetscErrorCode WriteIceFieldToFile(const char *filename, AppCtx *user) {
  FILE *file;
  PetscFunctionBegin;
  PetscErrorCode ierr;
  IGAElement element;
  IGAPoint point;
  PetscInt indGP; 
  PetscReal ice_val;

  file = fopen(filename, "w");
  if (!file) {
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_FILE_OPEN, "Error opening file for writing: %s", filename);
  }

  ierr = IGABeginElement(user->iga, &element); CHKERRQ(ierr);
  while (IGANextElement(user->iga, element)) {
      ierr = IGAElementBeginPoint(element, &point); CHKERRQ(ierr);
      while (IGAElementNextPoint(element, point)) {
        indGP = point->index + point->count * point->parent->index;

        ice_val = user->ice[indGP];
        fprintf(file, "%g %g %g %g\n", point->mapX[0][0], point->mapX[0][1], point->mapX[0][2], ice_val);

        if (point->atboundary) {
          PetscPrintf(PETSC_COMM_WORLD, "Warning: Gauss point %d is at the boundary!\n", indGP);
        }

        // user->ice[indGP] = 1.0;
      }
      ierr = IGAElementEndPoint(element, &point); CHKERRQ(ierr);
  }
  ierr = IGAEndElement(user->iga, &element); CHKERRQ(ierr);

  fclose(file);
  PetscPrintf(PETSC_COMM_WORLD, "Ice field successfully written to file: %s\n", filename);
  PetscFunctionReturn(0);
}

/*--------------------------------------------------------------------------------------------------
  Function: WriteOutput
  Purpose: Writes the given PETSc vector to file in binary format if user->outputBinary is set.
  Inputs:
    - user: Application context containing the binary flag.
    - x: Vector to be written.
    - filename: Name of output file.
--------------------------------------------------------------------------------------------------*/
PetscErrorCode WriteOutput(AppCtx *user, Vec x, const char *filename) {
  PetscErrorCode ierr;
  PetscFunctionBegin;

  if (user->outputBinary) {
    ierr = VecAssemblyBegin(x); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(x); CHKERRQ(ierr);
    ierr = WriteBinaryFile(x, filename); CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode ComputeKeffIntegrand(IGAPoint p, const PetscScalar *U,
                                    PetscInt n, PetscScalar *S, void *ctx)
{
  PetscFunctionBegin;

  // Define user context and retrieve necessary variables
  AppCtx *user = (AppCtx *)ctx;

  // Initialize variables
  PetscInt       i, j;             // Loop indices
  PetscReal      Iij;              // Identity delta
  PetscReal      dV;               // Differential volume

  PetscInt       dim = user->dim;  // Problem dimension (2 or 3)
  PetscReal      grad_t[dim][dim]; // ∇t: dim x dim

  // PetscInt       indGP = p->index + p->count * p->parent->index; // Gauss point index
  PetscReal      thcond; // Thermal conductivity

  // Safety chck
  if (!p->atboundary) {
    // U has dim components at each point
    // Compute the gradient of t vector field
    IGAPointFormGrad(p, U, &grad_t[0][0]);

    // Get the index of the Gauss point and thermal conductivity
    PetscInt indGP = p->index + p->count * p->parent->index;
    PetscReal ice = user->ice[indGP];
    ThermalCond(user, ice, &thcond, NULL);

    // Compute the differential volume
    dV = *(p->weight) * (*(p->detX));

    // Compute the integrand for the effective thermal conductivity
    for (i = 0; i < dim; i++) {
      for (j = 0; j < dim; j++) {
        Iij = (i == j) ? 1.0 : 0.0; // Identity delta
        S[i * dim + j] += (thcond * (grad_t[i][j] + Iij) * dV);
      }
    }
  }

  PetscFunctionReturn(0);
}

PetscErrorCode ComputeKeffective(IGA iga, Vec t_vec, PetscReal *keff, AppCtx *user)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;
  PetscInt dim = user->dim;
  
  // Initialize keff to zero
  ierr = PetscMemzero(keff, sizeof(PetscReal) * dim * dim); CHKERRQ(ierr);

  // Call IGAComputeScalar to integrate
  PetscPrintf(PETSC_COMM_WORLD, "Computing effective thermal conductivity...\n");
  ierr = IGAComputeScalar(iga, t_vec, dim * dim, keff, ComputeKeffIntegrand, (void*)user); CHKERRQ(ierr);


  // Divide by total domain volume
  PetscReal volume = (user->dim == 2) 
                     ? user->Lx * user->Ly 
                     : user->Lx * user->Ly * user->Lz;

  for (PetscInt i = 0; i < dim * dim; i++) {
    keff[i] /= volume;
  }

  PetscFunctionReturn(0);
}

// Helper to check if a file exists
PetscBool file_exists(const char *filename)
{
  struct stat buffer;
  return (stat(filename, &buffer) == 0) ? PETSC_TRUE : PETSC_FALSE;
}

// Write keff to CSV, add header if file doesn't exist
PetscErrorCode WriteKeffToCSV(AppCtx *user, const char *filename, PetscInt dim, const PetscReal *keff)
{
  PetscFunctionBegin;
  FILE *fp;
  PetscBool exists = file_exists(filename);

  // Open file in append mode
  fp = fopen(filename, "a");
  if (!fp) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Cannot open file: %s", filename);

  // If file didn't exist, write header
  if (!exists) {
    for (PetscInt i = 0; i < dim; i++) {
      for (PetscInt j = 0; j < dim; j++) {
        if (i > 0 || j > 0) {
          fprintf(fp, ","); // Commas between columns
        } else {
          fprintf(fp, "sol_index,"); // First column header
        }
        fprintf(fp, "k_%d%d", i, j); // Column names: k_00, k_01, ..., k_22
      }
    }
    fprintf(fp, "\n");
  }

  // Write keff matrix entries row-major
  for (PetscInt i = 0; i < dim * dim; i++) {
    if (i > 0) {
      fprintf(fp, ","); // Commas between columns
    } else {
      fprintf(fp, "%d,", user->sol_index); // First column header
    }
    fprintf(fp, "%.12e", keff[i]);
  }

  fprintf(fp, "\n"); // End of line
  fclose(fp);
  PetscFunctionReturn(0);
}

//------------------------------------------------------------------------------
// 9) Main
//------------------------------------------------------------------------------
int main(int argc, char *argv[]) {
  AppCtx              user;
  PetscErrorCode      ierr;

  /* ------------------ Initialize PETSc ------------------ */
  ierr = PetscInitialize(&argc, &argv, NULL, NULL); CHKERRQ(ierr);

  PetscMPIInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  // Set up timer
  PetscLogDouble t_start, t_1, t_end;
  ierr = PetscTime(&t_start); CHKERRQ(ierr);

  /* ------------------ Set options ------------------ */
  PetscBool print_error = PETSC_TRUE;
  PetscBool check_error = PETSC_FALSE;
  PetscBool save = PETSC_FALSE;
  PetscBool draw = PETSC_FALSE;
  PetscOptionsBegin(PETSC_COMM_WORLD,"","Laplace Options","IGA");CHKERRQ(ierr);
  ierr = PetscOptionsBool("-print_error","Prints the error of the solution",__FILE__,print_error,&print_error,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool("-check_error","Checks the error of the solution",__FILE__,check_error,&check_error,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool("-save","Save the solution to file",__FILE__,save,&save,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool("-draw","If dim <= 2, then draw the solution to the screen",__FILE__,draw,&draw,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsString("-init_mode", "Set initial ice field mode (circle, layered, file)", __FILE__, "circle", user.init_mode, sizeof(user.init_mode), NULL); CHKERRQ(ierr);
  ierr = PetscOptionsString("-init_dir", "Directory for initial condition files", __FILE__, ".", user.init_dir, sizeof(user.init_dir), NULL); CHKERRQ(ierr);
  PetscOptionsEnd();CHKERRQ(ierr);

  /* ------------------ Read Environment Variables ------------------ */
  ierr = GetEnvironment(&user); CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD, "Domain: Lx = %g, Ly = %g, Lz = %g\n", user.Lx, user.Ly, user.Lz);
  PetscPrintf(PETSC_COMM_WORLD, "eps = %g\n", user.eps);
  PetscPrintf(PETSC_COMM_WORLD, "dim = %d\n", user.dim);
  PetscPrintf(PETSC_COMM_WORLD, "TEMP_TOP = %g\n", user.T_top);
  PetscPrintf(PETSC_COMM_WORLD, "FLUX_BOTTOM = %g\n", user.q_bottom);

  // Initialize keff array
  PetscReal keff[user.dim * user.dim];
  ierr = PetscMemzero(keff, sizeof(PetscReal) * user.dim * user.dim); CHKERRQ(ierr);

  // Check if sol_indx > 0
  PetscPrintf(PETSC_COMM_WORLD, "SOL_INDEX = %d\n\n", user.sol_index);
  if (user.sol_index > 0) {
    /* ------------------ Define user context ------------------ */
    InitializeUserContext(&user);
    PetscPrintf(PETSC_COMM_WORLD, "thermal conductivity ice = %g W/m·K\n", user.thcond_ice);
    PetscPrintf(PETSC_COMM_WORLD, "thermal conductivity air = %g W/m·K\n", user.thcond_air);
    PetscPrintf(PETSC_COMM_WORLD, "Specific heat ice = %g J/kg·K\n", user.cp_ice);
    PetscPrintf(PETSC_COMM_WORLD, "Specific heat air = %g J/kg·K\n", user.cp_air);
    PetscPrintf(PETSC_COMM_WORLD, "Density ice = %g kg/m³\n", user.rho_ice);
    PetscPrintf(PETSC_COMM_WORLD, "Density air = %g kg/m³\n", user.rho_air);
    PetscPrintf(PETSC_COMM_WORLD, "User context initialized.\n\n");

    /* ------------------ Initialize IGA ------------------ */
    IGA iga;
    ierr = SetupIGA(&user, &iga); CHKERRQ(ierr); // Create and set up the IGA object
    user.iga = iga;

    /* ------------------ Initialize Field Variables ------------------ */
    ierr = FormInitialCondition(&user); CHKERRQ(ierr); // Initialize the ice field

    /* ------------------ Define Boundary Conditions ------------------ */
    // Apply the boundary conditions
    ierr = ApplyBoundaryConditions(iga, &user); CHKERRQ(ierr);

    /* ------------------ Set Up KSP Solver ------------------ */
    // Creat KSP solver
    ierr = SetupAndSolve(&user, iga); CHKERRQ(ierr);
    PetscTime(&t_1); // Record time after solving
    PetscPrintf(PETSC_COMM_WORLD, "KSP solve completed in %.2f seconds.\n\n", t_1 - t_start);

    /* ------------------ Compute Effective Thermal Conductivity ------------------ */
    ierr = ComputeKeffective(iga, user.T_sol, keff, &user); CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD, "Effective thermal conductivity computed:\n");
    for (PetscInt i = 0; i < user.dim; i++) {
      for (PetscInt j = 0; j < user.dim; j++) {
        PetscPrintf(PETSC_COMM_WORLD, "k_eff[%d][%d] = %.12e\n", i, j, keff[i * user.dim + j]);
      }
    }

    /* ------------------ Write Output ------------------ */
    if (rank == 0) {
      // Write the solution vector to a file
      char t_vec_file[PETSC_MAX_PATH_LEN];
      sprintf(t_vec_file, "%s/t_vec.dat", user.output_dir);
      ierr = WriteOutput(&user, user.T_sol, t_vec_file); CHKERRQ(ierr);
      PetscPrintf(PETSC_COMM_WORLD, "Successfully wrote solution vector to t_vec.dat\n");

      // Write ice field to file
      char iceFile[PETSC_MAX_PATH_LEN];
      sprintf(iceFile, "%s/ice_data.dat", user.output_dir);
      ierr = WriteIceFieldToFile(iceFile, &user); CHKERRQ(ierr);
      PetscPrintf(PETSC_COMM_WORLD, "Successfully wrote ice field to ice_data.dat\n");

      // Write the effective thermal conductivity to a CSV file
      char csvFile[PETSC_MAX_PATH_LEN];
      sprintf(csvFile, "%s/k_eff.csv", user.output_dir);
      ierr = WriteKeffToCSV(&user, csvFile, user.dim, keff); CHKERRQ(ierr);
      PetscPrintf(PETSC_COMM_WORLD, "Effective thermal conductivity written to k_eff.csv\n");

      PetscPrintf(PETSC_COMM_WORLD, "Wrote output files.\n\n");
    }

    /* ------------------ Clean Up ------------------ */
    ierr = IGADestroy(&iga); CHKERRQ(ierr);
    ierr = PetscFree(user.ice); CHKERRQ(ierr);
  } else {
    PetscPrintf(PETSC_COMM_WORLD, "SOL_INDEX < 0. Looping over entire solution folder.\n");

    // Open the directory
    DIR *dir;
    struct dirent *entry;
    dir = opendir(user.init_dir);
    if (!dir) {
      PetscPrintf(PETSC_COMM_WORLD, "Error opening directory: %s\n", user.init_dir);
      PetscFinalize();
      return EXIT_FAILURE;
    }

    // Create a dynamic array to store all file indices
    PetscInt file_indices[10000]; // Adjust size as needed
    PetscInt num_files = 0;

    // Regex to match file names
    regex_t regex;
    regcomp(&regex, "^sol_([0-9]+)\\.dat$", REG_EXTENDED);

    while ((entry = readdir(dir)) != NULL) {
      // Check if the file name matches the regex
      if (regexec(&regex, entry->d_name, 0, NULL, 0) == 0) {
        // Extract the index from the file name
        int index;
        sscanf(entry->d_name, "sol_%05d.dat", &index);
        file_indices[num_files++] = index; // Store the index
      }
    }
    closedir(dir);
    regfree(&regex); // Free the regex resources

    // Sort the file indices
    PetscSortInt(num_files, file_indices);
    PetscPrintf(PETSC_COMM_WORLD, "Found %d solution files in directory %s\n\n", num_files, user.init_dir);

    /* ------------------ Define user context ------------------ */
    InitializeUserContext(&user);
    PetscPrintf(PETSC_COMM_WORLD, "thermal conductivity ice = %g W/m·K\n", user.thcond_ice);
    PetscPrintf(PETSC_COMM_WORLD, "thermal conductivity air = %g W/m·K\n", user.thcond_air);
    PetscPrintf(PETSC_COMM_WORLD, "Specific heat ice = %g J/kg·K\n", user.cp_ice);
    PetscPrintf(PETSC_COMM_WORLD, "Specific heat air = %g J/kg·K\n", user.cp_air);
    PetscPrintf(PETSC_COMM_WORLD, "Density ice = %g kg/m³\n", user.rho_ice);
    PetscPrintf(PETSC_COMM_WORLD, "Density air = %g kg/m³\n\n", user.rho_air);

    /* ------------------ Initialize IGA ------------------ */
    IGA iga;
    ierr = SetupIGA(&user, &iga); CHKERRQ(ierr); // Create and set up the IGA object
    user.iga = iga;

    PetscInt num_files_analyze = 40;

    // Create a list of indices from 0 to num_files-1, evenly spaced with num_files_analyze entries
    PetscInt analyze_indices[num_files_analyze];
    for (PetscInt i = 0; i < num_files_analyze; i++) {
      PetscReal pos = ((PetscReal)i) * (num_files - 1) / (num_files_analyze - 1);
      analyze_indices[i] = file_indices[(PetscInt)(pos + 0.5)]; // round to nearest integer
    }

    for (PetscInt i = 0; i < num_files_analyze; i++)
    {
      // Update sol_index
    //   user.sol_index = file_indices[i];
      user.sol_index = analyze_indices[i];

      PetscPrintf(PETSC_COMM_WORLD, "Processing sol_index = %05d\n\n", user.sol_index);

      /* ------------------ Initialize Field Variables ------------------ */
      ierr = FormInitialCondition(&user); CHKERRQ(ierr); // Initialize the ice field

      /* ------------------ Define Boundary Conditions ------------------ */
      // Apply the boundary conditions
      ierr = ApplyBoundaryConditions(iga, &user); CHKERRQ(ierr);

      /* ------------------ Set Up KSP Solver ------------------ */
      // Creat KSP solver
      ierr = SetupAndSolve(&user, iga); CHKERRQ(ierr);

      /* ------------------ Compute Effective Thermal Conductivity ------------------ */
      ierr = ComputeKeffective(iga, user.T_sol, keff, &user); CHKERRQ(ierr);
      PetscPrintf(PETSC_COMM_WORLD, "Effective thermal conductivity computed:\n");
      for (PetscInt i = 0; i < user.dim; i++) {
        for (PetscInt j = 0; j < user.dim; j++) {
          PetscPrintf(PETSC_COMM_WORLD, "k_eff[%d][%d] = %.12e\n", i, j, keff[i * user.dim + j]);
        }
      }

      /* ------------------ Write Output ------------------ */
      if (rank == 0) {
        PetscPrintf(PETSC_COMM_WORLD, "Writing output for sol_index = %05d\n", user.sol_index);
        // Write the solution vector to a file
        char t_vec_file[PETSC_MAX_PATH_LEN];
        sprintf(t_vec_file, "%s/t_vec.dat", user.output_dir);
        ierr = WriteOutput(&user, user.T_sol, t_vec_file); CHKERRQ(ierr);
        PetscPrintf(PETSC_COMM_WORLD, "Successfully wrote solution vector to t_vec.dat\n");

        // Write ice field to file
        char iceFile[PETSC_MAX_PATH_LEN];
        sprintf(iceFile, "%s/ice_data.dat", user.output_dir);
        ierr = WriteIceFieldToFile(iceFile, &user); CHKERRQ(ierr);
        PetscPrintf(PETSC_COMM_WORLD, "Successfully wrote ice field to ice_data.dat\n");

        // Write the effective thermal conductivity to a CSV file
        char csvFile[PETSC_MAX_PATH_LEN];
        sprintf(csvFile, "%s/k_eff.csv", user.output_dir);
        ierr = WriteKeffToCSV(&user, csvFile, user.dim, keff); CHKERRQ(ierr);
        PetscPrintf(PETSC_COMM_WORLD, "Effective thermal conductivity written to k_eff.csv\n");

        PetscPrintf(PETSC_COMM_WORLD, "Wrote output files.\n\n");
      }

    } // End of for loop

    /* ------------------ Clean Up ------------------ */
    ierr = IGADestroy(&iga); CHKERRQ(ierr);
    ierr = PetscFree(user.ice); CHKERRQ(ierr);

  }

  // Calculate total elapsed time
  ierr = PetscTime(&t_end); CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD, "Total elapsed time: %.2f seconds\n", t_end - t_start);
  PetscPrintf(PETSC_COMM_WORLD, "Program completed successfully.\n");

  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;
}
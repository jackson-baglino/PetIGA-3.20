#include <petsc/private/tsimpl.h>
#include <petsc/private/snesimpl.h>
#include "petiga.h"
#include <math.h>
#include <mpi.h>

#include <petscsys.h>
#include <petscviewer.h>

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

  // **Output options**
  PetscBool outputBinary; // Flag for binary output

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

  PetscPrintf(PETSC_COMM_WORLD, "✅ Gauss point ice field computed with %d values.\n", idx);

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

        user->ice[indGP] = 0.5 - 0.5 * PetscTanhReal(1.0 / user->eps * dist);
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

  PetscPrintf(PETSC_COMM_WORLD, "Domain: Lx = %g, Ly = %g, Lz = %g\n\n", user->Lx, user->Ly, user->Lz);
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

  PetscPrintf(PETSC_COMM_WORLD, "Domain: Lx = %g, Ly = %g, Lz = %g\n\n", user->Lx, user->Ly, user->Lz);

  user->eps           = strtod(getenv("eps"), &endptr);

  user->outputBinary = (PetscBool)strtol(getenv("OUTPUT_BINARY"), &endptr, 10);

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

  PetscPrintf(PETSC_COMM_WORLD, "User context initialized.\n\n");
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

  PetscPrintf(PETSC_COMM_WORLD, "IGA setup complete.\n\n");
  PetscFunctionReturn(0);
}

/*------------------------------------------------------------------------------
  Function: AssembleStiffnessMatrix
  Purpose: Assemble FEM stiffness matrix and apply bottom flux BC.
------------------------------------------------------------------------------*/
PetscErrorCode AssembleStiffnessMatrix(IGAPoint pnt, PetscScalar *K, PetscScalar *F, void *ctx) {
  PetscErrorCode ierr;
  PetscFunctionBegin;

  // Retrieve user context and defiene variables
  AppCtx *user = (AppCtx *)ctx;
  PetscInt a, b, l;
  PetscInt nen = pnt->nen, dim = user->dim;
  const PetscScalar *N0;
  const PetscReal (*N1)[dim];

  // Retrieve basis functions
  ierr = IGAPointGetShapeFuns(pnt, 0, &N0); CHKERRQ(ierr);
  ierr = IGAPointGetShapeFuns(pnt, 1, (const PetscReal**)&N1); CHKERRQ(ierr);

  // Get ice phase at Gauss point
  PetscInt indGP = pnt->index + pnt->count * pnt->parent->index;
  PetscScalar ice = user->ice[indGP];

  // Evaluate thermal conductivity
  PetscReal thcond;
  ThermalCond(user, ice, &thcond, NULL);

  // Get the weight of the Gauss point
  // PetscReal W = *pnt->weight;

  for (a = 0; a < nen; a++) {
    F[a] += 0;
  }

  // Initialize stiffness matrix and force vector
  for (a = 0; a < nen; a++) {
    for (b = 0; b < nen; b++) {
      for (l = 0; l < dim; l++) {
        // Assemble stiffness matrix
        K[a * nen + b] += thcond * N1[a][l] * N1[b][l];
      }
    }
    // Apply Neumann boundary condition for bottom flux
    if (pnt->boundary_id == 2) {
      F[a] += N0[a] * user->q_bottom;
    }
  }
    
  PetscFunctionReturn(0);
}

/*------------------------------------------------------------------------------
Function: ApplyBoundaryConditions
Purpose: Set top temperature, bottom flux, and side zero-flux BCs.
------------------------------------------------------------------------------*/
PetscErrorCode ApplyBoundaryConditions(IGA iga, AppCtx *user) {
PetscErrorCode ierr;
PetscFunctionBegin;

PetscPrintf(PETSC_COMM_WORLD, "Applying boundary conditions...\n");

/* Fixed temperature at the top */
ierr = IGASetBoundaryValue(iga, 1, 1, 0, user->T_top); CHKERRQ(ierr);
PetscPrintf(PETSC_COMM_WORLD, "  - Fixed temperature applied at top: T = %g K\n\n", user->T_top);

// /* Fixed temperature at the top */
// ierr = IGASetBoundaryValue(iga, 1, 0, 0, user->T_top-3.0); CHKERRQ(ierr);
// PetscPrintf(PETSC_COMM_WORLD, "  - Fixed temperature applied at top: T = %g K\n\n", user->T_top);

/* Prescribed flux at the bottom */
ierr = IGASetBoundaryForm(iga, 1, 0, PETSC_TRUE); CHKERRQ(ierr);
// ierr = IGASetBoundaryLoad(iga, 1, 0, 0, -user->q_bottom); CHKERRQ(ierr);
PetscPrintf(PETSC_COMM_WORLD, "  - Prescribed flux at bottom: q = %g W/m²\n", user->q_bottom);

// /* Insulated side boundaries */
// ierr = IGASetBoundaryForm(iga, 0, 0, PETSC_TRUE); CHKERRQ(ierr);
// ierr = IGASetBoundaryForm(iga, 0, 1, PETSC_TRUE); CHKERRQ(ierr);
// PetscPrintf(PETSC_COMM_WORLD, "  - Insulated side boundaries applied\n");

PetscFunctionReturn(0);
}

/*------------------------------------------------------------------------------
  Function: InitialCondition
  Purpose: Set initial temperature using linear profile from BCs.
------------------------------------------------------------------------------*/
PetscErrorCode ComputeInitialCondition(Vec T, AppCtx *user) {
  PetscErrorCode ierr;
  PetscFunctionBegin;

  PetscScalar *T_array;
  ierr = VecGetArray(T, &T_array); CHKERRQ(ierr);

  PetscInt Nx = user->Nx, Ny = user->Ny;
  PetscReal Ly = user->Ly;
  PetscReal T_top = user->T_top;
  PetscReal q_bottom = user->q_bottom;
  PetscReal k_eff = 0.5 * (user->thcond_ice + user->thcond_air);

  for (PetscInt j = 0; j <= Ny; j++) {
    PetscReal y = (PetscReal)j * Ly / (PetscReal)Ny;
    PetscReal T_init = T_top - (q_bottom / k_eff) * (Ly - y);
    for (PetscInt i = 0; i <= Nx; i++) {
      PetscInt idx = j * (Nx + 1) + i;
      // T_array[idx] = T_init;
      T_array[idx] = 240.15;
    }
  }

  ierr = VecRestoreArray(T, &T_array); CHKERRQ(ierr);
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
  Vec b, u;
  KSP ksp;
  PC pc;

  PetscFunctionBegin;

  // Create system matrix and vectors
  ierr = IGACreateMat(iga, &A); CHKERRQ(ierr);
  ierr = IGACreateVec(iga, &user->T_sol); CHKERRQ(ierr);
  ierr = IGACreateVec(iga, &b); CHKERRQ(ierr);

  // Assemble system matrix and RHS vector
  ierr = IGASetFormSystem(iga, AssembleStiffnessMatrix, user); CHKERRQ(ierr);
  ierr = IGAComputeSystem(iga, A, b); CHKERRQ(ierr);

  // Set initial condition for temperature field
  ierr = ComputeInitialCondition(user->T_sol, user); CHKERRQ(ierr);
  // VecZeroEntries(user->T_sol);

  /* Change solver type to direct solver */

  // Solve the linear system using KSP
  ierr = IGACreateKSP(iga, &ksp); CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp, A, A); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);
  // Set solver tolerances
  // ierr = KSPSetTolerances(ksp, 1.0e-8, 1.0e-05, PETSC_DEFAULT, 4000); CHKERRQ(ierr); // rtol, abstol, dtol, maxits
  ierr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE); CHKERRQ(ierr);

  // Specify solver type (e.g., GMRES, CG, etc.)
  ierr = KSPSetType(ksp, KSPGMRES); CHKERRQ(ierr); // Replace KSPGMRES with desired solver type

  // Configure preconditioner
  ierr = KSPGetPC(ksp, &pc); CHKERRQ(ierr);
  ierr = PCSetType(pc, PCLU); CHKERRQ(ierr); // Use LU preconditioner

#if defined(PETSC_HAVE_MUMPS)
  if (size > 1) PetscCall(PCFactorSetMatSolverType(pc, MATSOLVERMUMPS));
#endif
  // ierr = PCFactorSetZeroPivot(pc, 1.0e-50); CHKERRQ(ierr);

  // Solve the system
  ierr = KSPSolve(ksp, b, user->T_sol); CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD, "KSP solve complete.\n\n");

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
    PetscPrintf(PETSC_COMM_WORLD, "Writing binary output...\n\n");
    ierr = VecAssemblyBegin(x); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(x); CHKERRQ(ierr);
    ierr = WriteBinaryFile(x, filename); CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

//------------------------------------------------------------------------------
// 9) Main
//------------------------------------------------------------------------------
int main (int argc, char *argv[]) {
  AppCtx              user;
  PetscErrorCode      ierr;

  /* ------------------ Initialize PETSc ------------------ */
  ierr = PetscInitialize(&argc, &argv, NULL, NULL); CHKERRQ(ierr);

  /* ------------------ Set options ------------------ */
  // Look into all of these options. See if there are other options that can be set.
  // These options are set in the command line when running the program.
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
  PetscOptionsEnd();CHKERRQ(ierr);

  /* ------------------ Read Environment Variables ------------------ */
  ierr = GetEnvironment(&user); CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD, "Domain: Lx = %g, Ly = %g, Lz = %g\n\n", user.Lx, user.Ly, user.Lz);
  PetscPrintf(PETSC_COMM_WORLD, "eps = %g\n", user.eps);
  PetscPrintf(PETSC_COMM_WORLD, "dim = %d\n", user.dim);
  PetscPrintf(PETSC_COMM_WORLD, "TEMP_TOP = %g\n", user.T_top);
  PetscPrintf(PETSC_COMM_WORLD, "FLUX_BOTTOM = %g\n", user.q_bottom);

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

  /* ------------------ Write Output ------------------ */
  ierr = WriteOutput(&user, user.T_sol, "temperature.bin"); CHKERRQ(ierr); // Write the solution to file
  ierr = WriteIceFieldToFile("ice_data.dat", &user); CHKERRQ(ierr); // Write the ice field to a .dat file

  /* ------------------ Clean Up ------------------ */
  ierr = IGADestroy(&iga); CHKERRQ(ierr);
  ierr = PetscFree(user.ice); CHKERRQ(ierr);
  ierr = PetscFinalize(); CHKERRQ(ierr);

  return 0;
}
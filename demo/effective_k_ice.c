#include <petsc/private/tsimpl.h>
#include <petsc/private/snesimpl.h>
#include "petiga.h"
#include <math.h>
#include <mpi.h>

#define SQ(x) ((x)*(x))
#define CU(x) ((x)*(x)*(x))

// Structure to hold application-specific parameters and data
typedef struct {
    IGA iga;  // Isogeometric Analysis (IGA) object for managing the finite element discretization

    // **Material properties**
    PetscReal thcond_ice;  // Thermal conductivity of ice (W/mK)
    PetscReal thcond_air;  // Thermal conductivity of air (W/mK)
    PetscReal cp_ice;      // Specific heat capacity of ice (J/kgK)
    PetscReal cp_air;      // Specific heat capacity of air (J/kgK)
    PetscReal rho_ice;     // Density of ice (kg/m³)
    PetscReal rho_air;     // Density of air (kg/m³)

    // **Initial conditions**
    PetscReal temp0;       // Initial temperature (K)
    PetscReal grad_temp0[3]; // Initial temperature gradient (x, y, z components)
    PetscReal *ice;        // Ice phase variable

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
                
                
                SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Invalid init_mode specified.");
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

void ThermalCond(AppCtx *user, PetscScalar ice, PetscScalar *cond, 
    PetscScalar *dcond_ice)
{
    PetscReal dice=1.0, dair=1.0;
    PetscReal air = 1.0-ice;
    if(ice<0.0) {ice=0.0;dice=0.0;}
    if(air<0.0) {air=0.0;dair=0.0;}
    PetscReal cond_ice = user->thcond_ice;
    PetscReal cond_air = user->thcond_air;  
    if(cond)      (*cond)  = ice*cond_ice + air*cond_air;
    if(dcond_ice)    (*dcond_ice) = cond_ice*dice-cond_air*dair;
  
    return;
}

void HeatCap(AppCtx *user, PetscScalar ice, PetscScalar *cp, 
    PetscScalar *dcp_ice)
{
    PetscReal dice=1.0, dair=1.0;
    PetscReal air = 1.0-ice;
    if(ice<0.0) {ice=0.0;dice=0.0;}
    if(air<0.0) {air=0.0;dair=0.0;}
    PetscReal cp_ice = user->cp_ice;
    PetscReal cp_air = user->cp_air;
    if(cp)     (*cp)  = ice*cp_ice + air*cp_air;
    if(dcp_ice)    (*dcp_ice) = cp_ice*dice-cp_air*dair;
  
    return;
}

void Density(AppCtx *user, PetscScalar ice, PetscScalar *rho, 
                PetscScalar *drho_ice)
  {
    PetscReal dice=1.0, dair=1.0;
    PetscReal air = 1.0-ice;
    if(ice<0.0) {ice=0.0;dice=0.0;}
    if(air<0.0) {air=0.0;dair=0.0;}
    PetscReal rho_ice = user->rho_ice;
    PetscReal rho_air = user->rho_air;
    if(rho)     (*rho)  = ice*rho_ice + air*rho_air;
    if(drho_ice)    (*drho_ice) = rho_ice*dice-rho_air*dair;
  
    
    return;
}

PetscErrorCode AssembleStiffnessMatrix(IGAPoint pnt, PetscScalar *K, PetscScalar *F, void *ctx) {
    PetscErrorCode ierr;
    PetscFunctionBegin;

    // Retrieve user context
    AppCtx *user = (AppCtx *)ctx;

    PetscInt a, b, l;
    PetscInt nen = pnt->nen, dim = user->dim;
    const PetscScalar *N0, (*N1)[dim];

    // Get the ice phase value and compute thermal conductivity
    PetscReal ice, thcond;
    PetscInt indGP = pnt->index + pnt->count * pnt->parent->index;
    ice = user->ice[indGP];
    ThermalCond(user, ice, &thcond, NULL);

    // Get basis functions
    ierr = IGAPointGetShapeFuns(pnt, 0, &N0); CHKERRQ(ierr);
    ierr = IGAPointGetShapeFuns(pnt, 1, (const PetscReal**)&N1); CHKERRQ(ierr);

    /* ===================================================
       APPLY FLUX BC AT THE BOTTOM BOUNDARY (y=0)
       =================================================== */
    if (pnt->atboundary && pnt->boundary_id == 2 ) {  // Check if the Gauss point is on a boundary
        // Boundary in y-direction (dir = 1) and bottom side (side = 0)
        for (a = 0; a < nen; a++) {
            F[a] -= N0[a] * user->q_bottom / user->Lx; // Apply flux as a Neumann condition
        }
    } else {
        // Loop over test functions (a) and trial functions (b)
        for (a = 0; a < nen; a++) {
            for (b = 0; b < nen; b++) {
                for (l = 0; l < dim; l++) {
                    K[a * nen + b] += thcond * N1[a][l] * N1[b][l]; // Diffusion term
                }
            }
        }
    }

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
    ierr = IGASetBoundaryForm(iga, 0, 0, PETSC_TRUE); CHKERRQ(ierr);  // x = 0
    ierr = IGASetBoundaryForm(iga, 0, 1, PETSC_TRUE); CHKERRQ(ierr);  // x = Lx
    PetscPrintf(PETSC_COMM_WORLD, "  - Zero-flux (Neumann) BCs applied at x = 0 and x = Lx\n");

    /* ===========================
       ZERO-FLUX (NEUMANN) BCs on Z-AXIS (ONLY IF 3D)
       =========================== */
    if (user->dim == 3) {
        ierr = IGASetBoundaryForm(iga, 2, 0, PETSC_TRUE); CHKERRQ(ierr);  // z = 0
        ierr = IGASetBoundaryForm(iga, 2, 1, PETSC_TRUE); CHKERRQ(ierr);  // z = Lz
        PetscPrintf(PETSC_COMM_WORLD, "  - Zero-flux (Neumann) BCs applied at z = 0 and z = Lz (3D case)\n");
    }

    PetscFunctionReturn(0);
}

PetscErrorCode ComputeInitialCondition(Vec T, AppCtx *user) {
    PetscErrorCode ierr;
    PetscFunctionBegin;

    // Get array to modify Vec directly
    PetscScalar *T_array;
    ierr = VecGetArray(T, &T_array); CHKERRQ(ierr);

    // Get mesh dimensions
    PetscInt Nx = user->Nx;
    PetscInt Ny = user->Ny;
    // PetscReal Lx = user->Lx;  // CHECK WHY WE DON'T USE THIS!
    PetscReal Ly = user->Ly;

    // Get boundary conditions
    PetscReal T_top = user->T_top;
    PetscReal q_bottom = user->q_bottom;  // Prescribed heat flux at the bottom

    // Approximate an effective thermal conductivity (use average of min and max)
    PetscReal k_eff = 0.5 * (user->thcond_ice + user->thcond_air);  // Can be refined

    // Loop over the mesh points and set initial temperature
    for (PetscInt j = 0; j <= Ny; j++) {
        PetscReal y = (PetscReal)j * Ly / (PetscReal)Ny;  // Physical y-position
        PetscReal T_init = T_top - (q_bottom / k_eff) * (Ly - y); // Linear profile

        for (PetscInt i = 0; i <= Nx; i++) {
            PetscInt idx = j * (Nx + 1) + i;  // Flattened 2D index
            T_array[idx] = T_init;  // Assign temperature
        }
    }

    // Restore array back to PETSc Vec
    ierr = VecRestoreArray(T, &T_array); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode ComputeAndStoreThermalConductivity(AppCtx *user, Vec K) {
    PetscInt i, N;
    PetscScalar ice, cond, dcond_ice;
    PetscErrorCode ierr;

    N = user->Nx * user->Ny; // Total number of grid points

    for (i = 0; i < N; i++) {
        ice = user->ice[i];
        ThermalCond(user, ice, &cond, &dcond_ice);

        // Store thermal conductivity value in K
        ierr = VecSetValue(K, i, cond, INSERT_VALUES); CHKERRQ(ierr);
    }

    ierr = VecAssemblyBegin(K); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(K); CHKERRQ(ierr);

    PetscPrintf(PETSC_COMM_WORLD, "Thermal conductivity field computed and stored.\n");

    PetscFunctionReturn(0);
}

PetscErrorCode WriteBinaryFile(Vec field, const char *filename) {
    PetscErrorCode ierr;
    PetscViewer viewer;
    PetscBool fileExists;
    
    PetscFunctionBegin;

    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename, FILE_MODE_WRITE, &viewer);CHKERRQ(ierr);

    // Debug: Check if the vector is valid before writing
    PetscScalar norm;
    ierr = VecNorm(field, NORM_2, &norm);CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD, "Norm of vector %s before writing: %g\n", filename, (double)norm);

    ierr = VecView(field, viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

    ierr = PetscTestFile(filename, 'r', &fileExists);CHKERRQ(ierr);
    if (!fileExists) {
        PetscPrintf(PETSC_COMM_WORLD, "Warning: Binary output file %s was not created!\n", filename);
    } else {
        PetscPrintf(PETSC_COMM_WORLD, "Binary output successfully written to %s\n", filename);
    }

    PetscFunctionReturn(0);
}

int main (int argc, char *argv[]) {
    AppCtx              user;
    PetscErrorCode      ierr;

    /* ------------------ Initialize PETSc ------------------ */
    ierr = PetscInitialize(&argc, &argv, NULL, NULL); CHKERRQ(ierr);

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
    PetscOptionsEnd();CHKERRQ(ierr);

    /* ------------------ Read Environment Variables ------------------ */
    ierr = GetEnvironment(&user); CHKERRQ(ierr);

    /* ------------------ Define user context ------------------ */
    // Material properties
    user.thcond_ice             = 2.29;     // Thermal conductivity of ice
    user.thcond_air             = 0.02;     // Thermal conductivity of air
    user.cp_ice                 = 1.96e3;   // Specific heat capacity of ice
    user.cp_air                 = 1.044e3;  // Specific heat capacity of air
    user.rho_ice                = 919.0;    // Density of ice
    user.rho_air                = 1.341;    // Density of air

    // Basis functions order and continuity
    PetscInt p = 1, C = 0;
    user.p                      = p;        // Polynomial order for basis functions      
    user.C                      = C;        // Continuity of basis functions

    // Set boundary conditions for a simple test case
    // user.T_bottom = 265.15;  // -8°C at y = 0  // NEED TO FIX INITAL GUESS!!!!å
    // user.T_top    = 270.15;  // -3°C at y = Ly

    /* Define user context */
    PetscPrintf(PETSC_COMM_WORLD, "Initializing thermal diffusion model...\n\n");

    /* ------------------ Initialize IGA ------------------ */
    IGA iga;
    ierr = IGACreate(PETSC_COMM_WORLD, &iga); CHKERRQ(ierr);
    ierr = IGASetDim(iga, user.dim); CHKERRQ(ierr);
    ierr = IGASetDof(iga, 1); CHKERRQ(ierr);
    ierr = IGASetFieldName(iga, 0, "temperature"); CHKERRQ(ierr);

    /* Set up IGA mesh */
    IGAAxis axisX, axisY, axisZ;
    ierr = IGAGetAxis(iga, 0, &axisX); CHKERRQ(ierr);
    ierr = IGAAxisSetDegree(axisX, p); CHKERRQ(ierr);
    ierr = IGAAxisInitUniform(axisX, user.Nx, 0.0, user.Lx, C); CHKERRQ(ierr);
    ierr = IGAGetAxis(iga, 1, &axisY); CHKERRQ(ierr);
    ierr = IGAAxisSetDegree(axisY, p); CHKERRQ(ierr);
    ierr = IGAAxisInitUniform(axisY, user.Ny, 0.0, user.Ly, C); CHKERRQ(ierr);
    if (user.dim == 3) {
        ierr = IGAGetAxis(iga, 2, &axisZ); CHKERRQ(ierr);
        ierr = IGAAxisSetDegree(axisZ, p); CHKERRQ(ierr);
        ierr = IGAAxisInitUniform(axisZ, user.Nz, 0.0, user.Lz, C); CHKERRQ(ierr);
    }

    ierr = IGASetFromOptions(iga); CHKERRQ(ierr);
    ierr = IGASetUp(iga); CHKERRQ(ierr);

    // Check the IGA object has been created correctly
    if (!iga) {
        PetscPrintf(PETSC_COMM_WORLD, "Error: IGA object is NULL before creating KSP.\n\n");
        return PETSC_ERR_ARG_WRONG;
    } else {
        PetscPrintf(PETSC_COMM_WORLD, "IGA object created.\n\n");
    }
    

    user.iga = iga;

    /* ------------------ Initialize Field Variables ------------------ */
    PetscInt nmb = iga->elem_width[0] * iga->elem_width[1] * SQ(p + 1);
    ierr = PetscMalloc1(nmb, &user.ice); CHKERRQ(ierr);
    ierr = PetscMemzero(user.ice, nmb * sizeof(PetscReal)); CHKERRQ(ierr);
    ierr = FormInitialCondition(&user); CHKERRQ(ierr);

    /* ------------------ Define Boundary Conditions ------------------ */
    // Apply the boundary conditions
    ierr = ApplyBoundaryConditions(iga, &user); CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD, "Boundary conditions applied.\n\n");

    /* ------------------ Set Up KSP Solver ------------------ */
    // Creat KSP solver
    Mat A;
    Vec x, b;

    ierr = IGACreateMat(iga, &A); CHKERRQ(ierr);
    ierr = IGACreateVec(iga, &x); CHKERRQ(ierr);
    ierr = IGACreateVec(iga, &b); CHKERRQ(ierr);
    ierr = IGASetFormSystem(iga, AssembleStiffnessMatrix, &user); CHKERRQ(ierr);
    ierr = IGAComputeSystem(iga, A, b); CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD, "System assembled.\n\n");

    // Print The Matrix size
    PetscInt m, n;
    ierr = MatGetSize(A, &m, &n); CHKERRQ(ierr);

    // Compute initial condition before solving
    ierr = ComputeInitialCondition(x, &user); CHKERRQ(ierr);

    KSP ksp;
    ierr = IGACreateKSP(iga, &ksp); CHKERRQ(ierr);
    // PetscPrintf(PETSC_COMM_WORLD, "KSP solver created.\n");

    ierr = KSPSetType(ksp, KSPCG); CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD, "KSP type set to GMRES.\n\n");

    ierr = KSPSetOperators(ksp, A, A); CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD, "KSP operators set.\n\n");    

    ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD, "KSP options set.\n\n");

    ierr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE); CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD, "KSP initial guess set.\n\n");

    ierr = KSPSolve(ksp, b, x); CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD, "KSP solve complete.\n\n");

    ierr = VecView(x,PETSC_VIEWER_DRAW_WORLD);CHKERRQ(ierr);

    /* ------------------ Write Output ------------------ */
    if (user.outputBinary) {
        PetscPrintf(PETSC_COMM_WORLD, "Writing binary...\n");

        // Ensure the vector is assembled before writing
        ierr = VecAssemblyBegin(x); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(x); CHKERRQ(ierr);

        // Write the solution vector `x` to a binary file
        ierr = WriteBinaryFile(x, "temperature.bin"); CHKERRQ(ierr);
    }
    ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
    ierr = MatDestroy(&A); CHKERRQ(ierr);
    ierr = VecDestroy(&x); CHKERRQ(ierr);
    ierr = VecDestroy(&b); CHKERRQ(ierr);
    ierr = IGADestroy(&iga); CHKERRQ(ierr);

    ierr = PetscFree(user.ice); CHKERRQ(ierr);

    ierr = PetscFinalize(); CHKERRQ(ierr);

    return 0;
}
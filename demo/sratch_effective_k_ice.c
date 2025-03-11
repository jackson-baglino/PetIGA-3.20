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

  // Material properties
  PetscReal thcond_ice;  // Thermal conductivity of ice
  PetscReal thcond_air;  // Thermal conductivity of air
  PetscReal cp_ice;      // Specific heat capacity of ice
  PetscReal cp_air;      // Specific heat capacity of air
  PetscReal rho_ice;     // Density of ice
  PetscReal rho_air;     // Density of air

  // Initial conditions
  PetscReal temp0;       // Initial temperature value
  PetscReal grad_temp0[3]; // Initial temperature gradient (x, y, z components)
  PetscReal Qsource;    // Source term for heat equation (could be a constant or function of position)
  PetscReal dQ_dT;      // Derivative of source term with respect to temperature
  PetscReal *ice;        // Ice phase variable (same as phase field definitoin)

  // Domain size and mesh resolution
  PetscInt dim;         // Dimension of the problem (2D or 3D)
  PetscReal Lx, Ly, Lz;  // Physical domain dimensions (length in x, y, z directions)
  PetscInt Nx, Ny, Nz;  // Number of elements in each direction

  // Normal vectors (possibly for boundary conditions or post-processing)
  PetscReal norm0_0; // Normal component in x-direction
  PetscReal norm0_1; // Normal component in y-direction
  PetscReal norm0_2; // Normal component in z-direction

} AppCtx;


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


PetscErrorCode ComputeEffectiveConductivity(IGA iga, Vec U, PetscReal *keff, AppCtx *user) {
    PetscErrorCode ierr;
    IGAElement elem;
    IGAPoint pnt;
    Vec localU;
    const PetscScalar *arrayU;
    PetscReal total_flux = 0.0, avg_gradT = 0.0;
    PetscInt num_elements = 0;
    PetscInt dim;

    PetscFunctionBeginUser;

    ierr = IGAGetDim(iga, &dim);CHKERRQ(ierr);

    // Get local array of U
    ierr = IGACreateLocalVec(iga, &localU);CHKERRQ(ierr);
    ierr = IGAGetLocalVecArray(iga, U, &localU, &arrayU);CHKERRQ(ierr);

    // Iterate over elements
    ierr = IGABeginElement(iga, &elem);CHKERRQ(ierr);
    while (IGANextElement(iga, elem)) {
        ierr = IGAElementBeginPoint(elem, &pnt);CHKERRQ(ierr);
        while (IGAElementNextPoint(elem, pnt)) {
            PetscScalar grad_T[dim];  // Temperature gradient
            PetscReal k_local;        // Thermal conductivity

            // AppCtx *user = (AppCtx*) pnt->user;

            ierr = IGAPointFormGrad(pnt, arrayU, &grad_T[0]);CHKERRQ(ierr);
            PetscScalar ice = user->ice[(PetscInt)pnt->count];
            ThermalCond(user, ice, &k_local, NULL);

            PetscReal q_flux = 0.0;
            for (PetscInt d = 0; d < dim; d++) {
                q_flux -= k_local * grad_T[d];
            }

            total_flux += q_flux;
            avg_gradT += PetscAbsReal(grad_T[0]);
            num_elements++;
        }
        ierr = IGAElementEndPoint(elem, &pnt);CHKERRQ(ierr);
    }
    ierr = IGAEndElement(iga, &elem);CHKERRQ(ierr);

    // Normalize gradient
    if (num_elements > 0) avg_gradT /= num_elements;

    // Compute effective thermal conductivity
    *keff = (avg_gradT > 1e-8) ? (-total_flux / avg_gradT) : 0.0;

    ierr = IGARestoreLocalVecArray(iga, U, &localU, &arrayU);CHKERRQ(ierr);
    ierr = VecDestroy(&localU);CHKERRQ(ierr);

    PetscPrintf(PETSC_COMM_WORLD, "Computed Effective Thermal Conductivity: %.6e W/mK\n", *keff);

    PetscFunctionReturn(0);
}

/*
 * InitializeSolidIceDomain - Sets up the initial ice phase field.
 * Uses a hyperbolic tangent function to smoothly transition ice phase.
 *
 * Inputs:
 *   user - User-defined context with simulation parameters.
 */
PetscErrorCode InitializeSolidIceDomain(IGA iga, AppCtx *user) {
    PetscErrorCode ierr;
    IGAElement elem;
    IGAPoint pnt;

    PetscFunctionBegin;

    // Ensure ice array is allocated
    if (!user->ice) {
        PetscPrintf(PETSC_COMM_WORLD, "ERROR: user->ice is not allocated!\n");
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "user->ice not allocated.");
    }

    PetscPrintf(PETSC_COMM_WORLD, "Initializing ice field...\n");

    // Iterate over elements
    ierr = IGABeginElement(iga, &elem); CHKERRQ(ierr);
    while (IGANextElement(iga, elem)) {
        ierr = IGAElementBeginPoint(elem, &pnt); CHKERRQ(ierr);
        while (IGAElementNextPoint(elem, pnt)) {
            PetscInt idx = pnt->index;  // Global index in the ice array
            user->ice[idx] = 1.0;  // Set ice to 1.0 everywhere
        }
        ierr = IGAElementEndPoint(elem, &pnt); CHKERRQ(ierr);
    }
    ierr = IGAEndElement(iga, &elem); CHKERRQ(ierr);

    PetscPrintf(PETSC_COMM_WORLD, "Ice initialization complete.\n");

    PetscFunctionReturn(0);
}


PetscErrorCode Residual(IGAPoint pnt, const PetscScalar *U, PetscScalar *Re, void *ctx)
{
    // Retrieve user-defined context (contains model parameters)
    AppCtx *user = (AppCtx*) ctx;

    // Get problem dimension from the user context
    PetscInt l, dim = user->dim;

    // Check if we are evaluating at a boundary
    if (pnt->atboundary) {
        // Here, we could impose boundary conditions if needed
        return 0;
    } 

    // Declare scalar values for temperature and its gradient
    // PetscScalar tem, grad_tem[dim];
    PetscScalar tem, *grad_tem;
    PetscMalloc1(dim, &grad_tem);

    // Extract solution values at this integration point
    IGAPointFormValue(pnt, U, &tem);   // Extract temperature T
    IGAPointFormGrad(pnt, U, &grad_tem[0]); // Extract temperature gradient ∇T

    // Declare thermophysical properties
    PetscReal thcond, Qsource;

    // Compute thermal conductivity at this point
    PetscScalar ice = user->ice[pnt->index]; // Get ice phase field value
    ThermalCond(user, ice, &thcond, NULL); // k(T)

    // Get the heat source term
    Qsource = user->Qsource; // This could be spatially varying if needed

    // Get basis function values for this element
    const PetscReal *N0, (*N1)[dim]; 
    IGAPointGetShapeFuns(pnt, 0, (const PetscReal**)&N0); // Basis function values (shape functions)
    IGAPointGetShapeFuns(pnt, 1, (const PetscReal**)&N1); // Gradients of shape functions

    // Pointer to residual array (stores contributions from this point)
    PetscScalar *R = (PetscScalar *) Re;
    PetscInt a, nen = pnt->nen; // Number of shape functions (nodes in this element)

    // Loop over all shape functions
    for (a = 0; a < nen; a++) {
        PetscReal R_tem = 0.0;

        // Diffusion term: ∇ · (k ∇T)
        for (l = 0; l < dim; l++) {
            R_tem += thcond * (N1[a][l] * grad_tem[l]);
        }

        // Add the heat source term: -Q
        R_tem -= Qsource * N0[a];

        // Store contribution to residual array
        R[a] = R_tem;
    }

    // Free memory for gradient array
    PetscFree(grad_tem);

    return 0;
}


PetscErrorCode Jacobian(IGAPoint pnt, const PetscScalar *U, PetscScalar *Je, void *ctx)
{
    PetscErrorCode ierr;

    if (!pnt) {
        PetscPrintf(PETSC_COMM_WORLD, "Error: IGAPoint is NULL in Jacobian function!\n");
        return PETSC_ERR_ARG_NULL;
    }
    if (!U) {
        PetscPrintf(PETSC_COMM_WORLD, "Error: Input U is NULL in Jacobian function!\n");
        return PETSC_ERR_ARG_NULL;
    }
    if (!Je) {
        PetscPrintf(PETSC_COMM_WORLD, "Error: Jacobian matrix Je is NULL in Jacobian function!\n");
        return PETSC_ERR_ARG_NULL;
    }
    if (!ctx) {
        PetscPrintf(PETSC_COMM_WORLD, "Error: Context ctx is NULL in Jacobian function!\n");
        return PETSC_ERR_ARG_NULL;
    }

    // PetscPrintf(PETSC_COMM_WORLD, "Jacobian function: all pointers are valid, proceeding with computation.\n");
    
    AppCtx *user = (AppCtx*) ctx;
    if (!user) {
        PetscPrintf(PETSC_COMM_WORLD, "Error: ctx cast to user struct is NULL!\n");
        return PETSC_ERR_ARG_NULL;
    }

    PetscInt dim = user->dim;
    PetscInt nen = pnt->nen; // Number of shape functions (nodes in this element)
    PetscScalar (*J)[nen] = (PetscScalar (*)[nen]) Je;

    // PetscPrintf(PETSC_COMM_WORLD, "Jacobian Function: dim = %d\n", dim);

    // Check if we are evaluating at a boundary
    if (pnt->atboundary) {
        return 0; // No contribution at boundaries (modify as needed)
    }

    // Declare scalar values for temperature and its gradient
    // PetscScalar tem, grad_tem[dim];
    PetscScalar tem, *grad_tem;
    PetscMalloc1(dim, &grad_tem);

    // Extract solution values at this integration point
    if (!U) {
        PetscPrintf(PETSC_COMM_WORLD, "Error: Input U is NULL in Jacobian function!\n");
        return PETSC_ERR_ARG_NULL;
    }

    IGAPointFormValue(pnt, U, &tem);   // Extract temperature T

    PetscPrintf(PETSC_COMM_WORLD, "Before calling IGAPointFormGrad\n");
    ierr = IGAPointFormGrad(pnt, U, grad_tem);CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD, "Extracted gradient: [%g, %g]\n", grad_tem[0], (dim > 1) ? grad_tem[1] : 0);

    // Declare thermophysical properties and their derivatives
    PetscReal thcond, dthcond_dT, dQ_dT;

    // Compute thermal conductivity and its derivative
    PetscScalar ice = user->ice[(PetscInt)pnt->count];
    ThermalCond(user, ice, &thcond, &dthcond_dT); // k(T) and dk/dT

    // Compute derivative of heat source term Q with respect to temperature
    dQ_dT = user->dQ_dT; // Assuming a constant dQ/dT; modify if needed

    // Get basis function values for this element
    const PetscReal *N0, (*N1)[dim];
    IGAPointGetShapeFuns(pnt, 0, &N0);
    IGAPointGetShapeFuns(pnt, 1, (const PetscReal**)&N1);
    // PetscPrintf(PETSC_COMM_WORLD, "Shape functions!\n");

    // PetscPrintf(PETSC_COMM_WORLD, "Define pointer to J mat.\n");

    // Loop over all shape functions to construct Jacobian entries
    for (PetscInt a = 0; a < nen; a++) {
        for (PetscInt b = 0; b < nen; b++) {

            PetscReal J_tem = 0.0;

            // Diffusion term: ∂(∇ · (k ∇T)) / ∂T
            for (PetscInt l = 0; l < dim; l++) {
                // PetscPrintf(PETSC_COMM_WORLD, "Inside J loop: a=%d, b=%d, l=%d\n", a, b, l);

                // Ensure all pointers are valid before dereferencing
                if (!N0 || !N1[a] || !N1[b]) {
                    PetscPrintf(PETSC_COMM_WORLD, "Error: Invalid shape function pointers (a=%d, b=%d, l=%d)\n", a, b, l);
                    return PETSC_ERR_ARG_NULL;
                }
                if (l >= dim) {
                    PetscPrintf(PETSC_COMM_WORLD, "Error: Index l=%d is out of bounds (dim=%d)\n", l, dim);
                    return PETSC_ERR_ARG_OUTOFRANGE;
                }

                // Check values before computation
                PetscPrintf(PETSC_COMM_WORLD, "thcond=%g, dthcond_dT=%g\n", thcond, dthcond_dT);
                // PetscPrintf(PETSC_COMM_WORLD, "N1[a][l]=%g, N1[b][l]=%g\n", N1[a][l], N1[b][l]);

                // Computation of Jacobian
                J_tem += thcond * (N1[a][l] * N1[b][l]); // Standard term

                if (grad_tem[0]) { // Ensure gradient is valid
                    // PetscPrintf(PETSC_COMM_WORLD, "grad_tem[l]=%g\n", grad_tem[l]);
                    J_tem += dthcond_dT * N0[b] * (N1[a][l] * grad_tem[l]); // dk/dT term
                } else {
                    PetscPrintf(PETSC_COMM_WORLD, "Error: grad_tem is NULL (a=%d, b=%d, l=%d)\n", a, b, l);
                    return PETSC_ERR_ARG_NULL;
                }
            }

            // Source term derivative: ∂(-Q) / ∂T
            J_tem -= dQ_dT * N0[a] * N0[b];

            // Store contribution to Jacobian array
            J[a][b] = J_tem;
        }
    }

    // Free memory for gradient array
    PetscFree(grad_tem);

    return 0;
}


/*
 * WriteVTKFile - Outputs the computed temperature field to a .vtk file for visualization in ParaView.
 *
 * Inputs:
 *   iga       - The IGA object defining the simulation domain.
 *   U         - The solution vector storing the temperature field.
 *   filename  - Name of the output VTK file.
 */
PetscErrorCode WriteVTKFile(IGA iga, Vec U, const char *filename) {
    PetscErrorCode ierr;
    PetscViewer viewer;

    PetscFunctionBegin;

    // Create a VTK file viewer
    ierr = PetscViewerVTKOpen(PETSC_COMM_WORLD, filename, FILE_MODE_WRITE, &viewer);CHKERRQ(ierr);

    // Write the solution vector to the VTK file
    ierr = IGAWriteVec(iga, U, filename);CHKERRQ(ierr);

    // Destroy the viewer after writing
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


/*
 * WriteBinaryFile - Saves the computed temperature field in a binary file.
 *
 * Inputs:
 *   U         - The solution vector storing the temperature field.
 *   filename  - Name of the output binary file.
 */
PetscErrorCode WriteBinaryFile(Vec U, const char *filename) {
    PetscErrorCode ierr;
    PetscViewer viewer;

    PetscFunctionBegin;

    // Create a binary file viewer
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename, FILE_MODE_WRITE, &viewer);CHKERRQ(ierr);

    // Write the solution vector to the binary file
    ierr = VecView(U, viewer);CHKERRQ(ierr);

    // Destroy the viewer after writing
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode IGAComputeFunctionSNES(SNES snes, Vec U, Vec F, void *ctx) {
    AppCtx *user = (AppCtx *)ctx;
    return IGAComputeFunction(user->iga, U, F);
}

PetscErrorCode IGAComputeJacobianSNES(SNES snes, Vec U, Mat J, Mat P, void *ctx) {
    AppCtx *user = (AppCtx *)ctx;

    PetscCheck(user->iga, PETSC_COMM_WORLD, PETSC_ERR_ARG_NULL, "IGA is NULL before calling Jacobian.");
    
    return IGAComputeJacobian(user->iga, U, J);
}

/*
 * main - Main function to set up and solve the thermal diffusion problem.
 */
int main(int argc, char *argv[]) {
    PetscErrorCode ierr;
    ierr = PetscInitialize(&argc, &argv, NULL, NULL); CHKERRQ(ierr);

    /* Define user context */
    AppCtx user;
    PetscPrintf(PETSC_COMM_WORLD, "Initializing thermal diffusion model...\n");

    /* ------------------ Read Environment Variables ------------------ */
    PetscPrintf(PETSC_COMM_WORLD, "Reading simulation parameters...\n");

    const char *Nx_str = getenv("Nx"), *Ny_str = getenv("Ny"), *Nz_str = getenv("Nz");
    const char *Lx_str = getenv("Lx"), *Ly_str = getenv("Ly"), *Lz_str = getenv("Lz");
    const char *temp_str = getenv("temp"), *grad_temp0X_str = getenv("grad_temp0X");
    const char *grad_temp0Y_str = getenv("grad_temp0Y"), *grad_temp0Z_str = getenv("grad_temp0Z");
    const char *dim_str = getenv("dim");

    if (!Nx_str || !Ny_str || !Nz_str || !Lx_str || !Ly_str || !Lz_str || 
        !temp_str || !grad_temp0X_str || !grad_temp0Y_str || !grad_temp0Z_str || !dim_str) {
        PetscPrintf(PETSC_COMM_WORLD, "Error: Missing environment variables.\n");
        PetscFinalize();
        return EXIT_FAILURE;
    }

    /* Convert environment variables */
    char *endptr;
    user.Nx = (PetscInt)strtol(Nx_str, &endptr, 10);
    user.Ny = (PetscInt)strtol(Ny_str, &endptr, 10);
    user.Nz = (PetscInt)strtol(Nz_str, &endptr, 10);
    user.Lx = strtod(Lx_str, &endptr);
    user.Ly = strtod(Ly_str, &endptr);
    user.Lz = strtod(Lz_str, &endptr);
    user.temp0 = strtod(temp_str, &endptr);
    user.grad_temp0[0] = strtod(grad_temp0X_str, &endptr);
    user.grad_temp0[1] = strtod(grad_temp0Y_str, &endptr);
    user.grad_temp0[2] = strtod(grad_temp0Z_str, &endptr);
    user.dim = (PetscInt)strtol(dim_str, &endptr, 10);

    /* ------------------ Initialize IGA ------------------ */
    IGA iga;
    ierr = IGACreate(PETSC_COMM_WORLD, &iga); CHKERRQ(ierr);
    ierr = IGASetDim(iga, user.dim); CHKERRQ(ierr);
    ierr = IGASetDof(iga, 1); CHKERRQ(ierr);
    ierr = IGASetFieldName(iga, 0, "temperature"); CHKERRQ(ierr);

    /* Set up mesh */
    IGAAxis axisX, axisY, axisZ;
    ierr = IGAGetAxis(iga, 0, &axisX); CHKERRQ(ierr);
    ierr = IGAAxisSetDegree(axisX, 2); CHKERRQ(ierr);
    ierr = IGAAxisInitUniform(axisX, user.Nx, 0.0, user.Lx, 1); CHKERRQ(ierr);
    ierr = IGAGetAxis(iga, 1, &axisY); CHKERRQ(ierr);
    ierr = IGAAxisSetDegree(axisY, 2); CHKERRQ(ierr);
    ierr = IGAAxisInitUniform(axisY, user.Ny, 0.0, user.Ly, 1); CHKERRQ(ierr);
    if (user.dim == 3) {
        ierr = IGAGetAxis(iga, 2, &axisZ); CHKERRQ(ierr);
        ierr = IGAAxisSetDegree(axisZ, 2); CHKERRQ(ierr);
        ierr = IGAAxisInitUniform(axisZ, user.Nz, 0.0, user.Lz, 1); CHKERRQ(ierr);
    }

    ierr = IGASetFromOptions(iga); CHKERRQ(ierr);
    ierr = IGASetUp(iga); CHKERRQ(ierr);
    user.iga = iga;

    /* ------------------ Initialize Ice Field ------------------ */
    PetscInt p;
    ierr = IGAGetOrder(iga, &p); CHKERRQ(ierr);
    
    // Compute the number of basis function nodes per element
    PetscInt nmb = iga->elem_width[0] * iga->elem_width[1] * SQ(p+1);
    if (user.dim == 3) {
        nmb = iga->elem_width[0] * iga->elem_width[1] * iga->elem_width[2] * CU(p+1);
    }
    
    // Allocate memory for ice phase field
    ierr = PetscMalloc(sizeof(PetscReal) * nmb, &user.ice); CHKERRQ(ierr);
    
    // Initialize ice field to 1.0 everywhere
    ierr = PetscMemzero(user.ice, sizeof(PetscReal) * nmb); CHKERRQ(ierr);
    for (PetscInt i = 0; i < nmb; i++) {
        user.ice[i] = 1.0;
    }

    ierr = InitializeSolidIceDomain(iga, &user); CHKERRQ(ierr);

    /* ------------------ Define Boundary Conditions ------------------ */
    ierr = IGASetBoundaryValue(iga, 0, 0, 0, user.temp0); CHKERRQ(ierr);  
    ierr = IGASetBoundaryValue(iga, 0, 1, 0, user.temp0 - user.grad_temp0[0] * user.Lx); CHKERRQ(ierr);  

    /* ------------------ Set Up SNES Solver ------------------ */
    SNES snes;
    ierr = SNESCreate(PETSC_COMM_WORLD, &snes); CHKERRQ(ierr);

    PetscMPIInt size;
    MPI_Comm_size(PETSC_COMM_WORLD, &size);
    PetscPrintf(PETSC_COMM_WORLD, "Running solver on %d MPI processes.\n", size);

    /* Create solution vector */
    Vec U;
    ierr = IGACreateVec(iga, &U); CHKERRQ(ierr);

    /* Set up SNES with residual and Jacobian functions */
    ierr = IGASetFormFunction(iga, Residual, &user); CHKERRQ(ierr);
    ierr = IGASetFormJacobian(iga, Jacobian, &user); CHKERRQ(ierr);

    /* Set SNES function and Jacobian correctly */
    ierr = SNESSetFunction(snes, NULL, IGAComputeFunctionSNES, &user); CHKERRQ(ierr);
    ierr = SNESSetJacobian(snes, NULL, NULL, IGAComputeJacobianSNES, &user); CHKERRQ(ierr);

    /* Solve nonlinear system */
    ierr = SNESSolve(snes, NULL, U); CHKERRQ(ierr);

    /* ------------------ Compute Effective Thermal Conductivity ------------------ */
    PetscReal keff;
    ierr = ComputeEffectiveConductivity(iga, U, &keff, &user); CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD, "Computed Effective Thermal Conductivity: %.6e W/mK\n", keff);

    /* ------------------ Cleanup ------------------ */
    ierr = VecDestroy(&U); CHKERRQ(ierr);
    ierr = SNESDestroy(&snes); CHKERRQ(ierr);
    ierr = IGADestroy(&iga); CHKERRQ(ierr);
    ierr = PetscFree(user.ice); CHKERRQ(ierr);
    ierr = PetscFinalize(); CHKERRQ(ierr);

    return 0;
}
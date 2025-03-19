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
    PetscReal T_bottom; // Temperature at y = 0
    PetscReal T_top;    // Temperature at y = Ly

    // **Output options**
    PetscBool outputBinary; // Flag for binary output

} AppCtx;


/* 
 * GetEnvironment - Retrieves environmental variables 
 * 
 * Inputs:
 * 
 */
PetscErrorCode GetEnvironment(AppCtx *user) {
    PetscFunctionBegin;

    PetscPrintf(PETSC_COMM_WORLD, "Reading simulation parameters...\n");

    const char *Nx_str = getenv("Nx"), *Ny_str = getenv("Ny"), *Nz_str = getenv("Nz");
    const char *Lx_str = getenv("Lx"), *Ly_str = getenv("Ly"), *Lz_str = getenv("Lz");
    const char *temp_str = getenv("temp"), *grad_temp0X_str = getenv("grad_temp0X");
    const char *grad_temp0Y_str = getenv("grad_temp0Y"), *grad_temp0Z_str = getenv("grad_temp0Z");
    const char *dim_str = getenv("dim"); const char *outputBinary_str = getenv("OUTPUT_BINARY");

    if (!Nx_str || !Ny_str || !Nz_str || !Lx_str || !Ly_str || !Lz_str || 
        !temp_str || !grad_temp0X_str || !grad_temp0Y_str || !grad_temp0Z_str ||
         !dim_str || !outputBinary_str) {
        PetscPrintf(PETSC_COMM_WORLD, "Error: Missing environment variables.\n");
        PetscFinalize();
        return EXIT_FAILURE;
    }

    /* Convert environment variables */
    char *endptr;
    user->dim = (PetscInt)strtol(dim_str, &endptr, 10);
    user->Nx = (PetscInt)strtol(Nx_str, &endptr, 10);
    user->Ny = (PetscInt)strtol(Ny_str, &endptr, 10);
    if (user->dim == 3) {
        user->Nz = (PetscInt)strtol(Nz_str, &endptr, 10);
    } else {
        user->Nz = 1;
    }
    user->Lx = strtod(Lx_str, &endptr);
    user->Ly = strtod(Ly_str, &endptr);
    user->Lz = strtod(Lz_str, &endptr);
    user->temp0 = strtod(temp_str, &endptr);
    user->grad_temp0[0] = strtod(grad_temp0X_str, &endptr);
    user->grad_temp0[1] = strtod(grad_temp0Y_str, &endptr);
    user->grad_temp0[2] = strtod(grad_temp0Z_str, &endptr);
    user->outputBinary = (PetscBool)strtol(outputBinary_str, &endptr, 10);

    PetscFunctionReturn(0);
}

// Initialize the ice field how Adrian did in the other codes
PetscErrorCode InitializeIceField(AppCtx *user) {
    PetscFunctionBegin;
    
    PetscErrorCode ierr;
    IGAElement    element;
    IGAPoint      point;
    PetscReal     ice, dist;
    PetscInt      idx = 0;
    PetscInt      num_points = user->Nx * user->Ny;  // Store only at element nodes!

    ierr = IGABeginElement(user->iga, &element); CHKERRQ(ierr);
    while (IGANextElement(user->iga, element)) {
        ierr = IGAElementBeginPoint(element, &point); CHKERRQ(ierr);
        while (IGAElementNextPoint(element, point)) {
            ice = 0.0;

            // Compute distance from midline (only using y-coordinate)
            dist = user->Ly / 2.0 - point->mapX[0][1];

            // Compute hyperbolic tangent transition
            // ice = 0.5 - 0.5 * tanh(2.0 / user->eps * dist);
            // ice = PetscMax(0.0, PetscMin(1.0, ice));  // Clamp between 0 and 1
            ice = 1.0;
            // Store ice value at node
            user->ice[idx] = ice;

            PetscPrintf(PETSC_COMM_WORLD, "%0.3f ", ice);
            idx++;

            if (idx % (2*user->Nx) == 0) {
                PetscPrintf(PETSC_COMM_WORLD, "\n");
            }
        }
        ierr = IGAElementEndPoint(element, &point); CHKERRQ(ierr);
    }
    ierr = IGAEndElement(user->iga, &element); CHKERRQ(ierr);

    PetscPrintf(PETSC_COMM_WORLD, "Ice field initialized at %d nodes.\n", idx);

    // ierr = IGABeginElement(user->iga, &element); CHKERRQ(ierr);
    // while(IGANextElement(user->iga, element)) {
    //     ierr = IGAElementBeginPoint(element, &point); CHKERRQ(ierr);
    //     while(IGAElementNextPoint(element, point)) {
    //         ice = 0.0;
    //         // Want to perscribe a hyperbolic tangent profile for the boundary
    //         // between ice and air
    //         dist = user->Ly/2.0 - point->mapX[0][1];
    //         ice += 0.5 - 0.5 * tanh(2.0/user->eps * dist);

    //         if(ice > 1.0) ice = 1.0;
    //         if(ice < 0.0) ice = 0.0;

    //         // idx = point->index + point->count * point->parent->index;
    //         user->ice[idx] = ice;
    //         // PetscPrintf(PETSC_COMM_WORLD, "[%d] Ice field value: %f\n", idx, ice);
    //         PetscPrintf(PETSC_COMM_WORLD, "%f ", dist);
    //         idx++;

    //         if (idx % user->Nx == 0) {
    //             PetscPrintf(PETSC_COMM_WORLD, "\n");
    //         }
    //     }
    //     ierr = IGAElementEndPoint(element, &point); CHKERRQ(ierr);
    // }
    // ierr = IGAEndElement(user->iga, &element); CHKERRQ(ierr);

    PetscPrintf(PETSC_COMM_WORLD, "Ice field initialized after %d iterations.\n", idx);

    // Write Ice Field to Binary File
    PetscViewer viewer;
    Vec iceVec;
    PetscInt *indices;

    // Allocate index array
    ierr = PetscMalloc1(idx, &indices); CHKERRQ(ierr);
    for (PetscInt i = 0; i < idx; i++) {
        indices[i] = i;
    }

    // Create PETSc Vector
    ierr = VecCreate(PETSC_COMM_WORLD, &iceVec); CHKERRQ(ierr);
    ierr = VecSetSizes(iceVec, PETSC_DECIDE, idx); CHKERRQ(ierr);
    ierr = VecSetFromOptions(iceVec); CHKERRQ(ierr);

    // Copy ice field into Vec using explicit indices
    ierr = VecSetValues(iceVec, idx, indices, user->ice, INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(iceVec); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(iceVec); CHKERRQ(ierr);

    // Free indices array
    ierr = PetscFree(indices); CHKERRQ(ierr);

    // Save to binary file
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, "ice_field.bin", FILE_MODE_WRITE, &viewer); CHKERRQ(ierr);
    ierr = VecView(iceVec, viewer); CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
    ierr = VecDestroy(&iceVec); CHKERRQ(ierr);

    PetscPrintf(PETSC_COMM_WORLD, "Ice field successfully written to ice_field.bin\n");

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


PetscErrorCode Residual(IGAPoint pnt, const PetscScalar *U, PetscScalar *Re, void *ctx) 
  {
      PetscFunctionBegin;
      AppCtx *user = (AppCtx *)ctx;
  
      PetscInt dim = user->dim;
      PetscInt nen = pnt->nen;  
      PetscScalar *R = (PetscScalar *)Re; 
  
    //   PetscScalar T;
      PetscScalar grad_T[dim];
      PetscScalar ice;

      // Compute the global index of the current Gauss point
      PetscInt indGP = pnt->index + pnt->count * pnt->parent->index;

      ice = user->ice[indGP];

      PetscReal thcond;
      const PetscReal *N0;
      const PetscReal (*N1)[dim];
  
      PetscInt a, l;
  
      // Get temperature and its gradient
    //   IGAPointFormValue(pnt, U, &T); // Don't need!
      IGAPointFormGrad(pnt, U, grad_T);
      // Right now U 
  
      // Get material properties
      ThermalCond(user, ice, &thcond, NULL);
    //   PetscPrintf(PETSC_COMM_WORLD, "[%d] thcond = %f\n", pnt->index, thcond);
  
      // Get shape function values
      IGAPointGetShapeFuns(pnt, 0, &N0);
      IGAPointGetShapeFuns(pnt, 1, (const PetscReal**)&N1);

    //   PetscPrintf(PETSC_COMM_WORLD, "[RESIDUAL] nen = %d\n", nen);
  
      // Standard residual contribution for diffusion
      for (a = 0; a < nen; a++) {
        PetscScalar R_T = 0.0;
        for (l = 0; l < dim; l++) {
            R_T += thcond * (N1[a][l] * grad_T[l]);
        }
        R[a] += N0[a] * R_T; // Instead of direct assignment
    }
  
      PetscFunctionReturn(0);
}


PetscErrorCode Jacobian(IGAPoint pnt, const PetscScalar *U, PetscScalar *Je, void *ctx) 
{
    // PetscFunctionBegin;
    AppCtx *user = (AppCtx *)ctx;

    PetscInt dim = user->dim;
    PetscInt nen = pnt->nen;
    PetscScalar *J = Je;

    // PetscScalar T; // Don't need
    PetscScalar grad_T[dim];

    PetscScalar ice;
    // Compute the global index of the current Gauss point
    PetscInt indGP = pnt->index + pnt->count * pnt->parent->index;

    ice = user->ice[indGP];
    
    PetscReal thcond;
    PetscReal dthcond_dT;

    const PetscReal *N0;
    const PetscReal (*N1)[dim];

    PetscInt a, b, l;

    // Get temperature and gradient
    IGAPointFormGrad(pnt, U, grad_T);

    // Get material properties
    ThermalCond(user, ice, &thcond, &dthcond_dT);

    // Get shape function values
    IGAPointGetShapeFuns(pnt, 0, &N0);
    IGAPointGetShapeFuns(pnt, 1, (const PetscReal**)&N1);

    // **Modify the Jacobian for diffusion terms**
    for (a = 0; a < nen; a++) {
        for (b = 0; b < nen; b++) {
            PetscScalar J_T = 0.0;
            for (l = 0; l < dim; l++) {
                J_T += thcond * (N1[a][l] * N1[b][l]); // Need to include the  
                // J_T += dthcond_dT * N0[b] * (N1[a][l] * grad_T[l]); // Thermal conductivity does not depend on temperature dK_dT = 0
            }
            J[a * nen + b] += J_T; // Instead of direct assignment
        }
    }

    return 0;
}

/*
 * StiffnessMatirx - Computes the stiffness matrix for the heat conduction problem.
    * Inputs:
    *   iga  - IGA object defining the domain.
    *   user - Pointer to the user-defined structure containing boundary condition settings. 
 */
PetscErrorCode AssembleStiffnessMatrix(IGAPoint pnt,PetscScalar *K,PetscScalar *F,void *ctx) {
    PetscErrorCode ierr;

    PetscFunctionBegin;
    AppCtx *user = (AppCtx *)ctx;

    PetscPrintf(PETSC_COMM_WORLD, "Assembling stiffness matrix...\n\n");

    // PetscInt dim = user->dim;
    PetscPrintf(PETSC_COMM_WORLD, "dim = %d\n", user->dim);
    PetscInt nen = pnt->nen;
    PetscPrintf(PETSC_COMM_WORLD, "nen = %d\n", nen);
    PetscInt a, b;
    PetscInt dim = user->dim;



    const PetscScalar *N0, (*N1)[dim];


    // Get Thermal conductivity
    PetscReal ice, thcond;
    PetscInt indGP = pnt->index + pnt->count * pnt->parent->index;
    PetscPrintf(PETSC_COMM_WORLD, "indGP = %d\n", indGP);
    ice = user->ice[indGP];
    ThermalCond(user, ice, &thcond, NULL);

    ierr = IGAPointGetShapeFuns(pnt, 0, &N0); CHKERRQ(ierr);
    ierr = IGAPointGetShapeFuns(pnt, 1, (const PetscReal**)&N1); CHKERRQ(ierr);

    for (a = 0; a < nen; a++) {
        PetscReal Na   = N0[a];
        PetscReal Na_x = N1[a][0];
        PetscReal Na_y = N1[a][1];
        for (b = 0; b < nen; b++) {
            PetscReal Nb_x = N1[b][0];
            PetscReal Nb_y = N1[b][1];

            // Stiffness matrix
            K[a*nen+b] = thcond * (Na_x*Nb_x + Na_y*Nb_y);
            PetscPrintf(PETSC_COMM_WORLD, "[a=%d, b=%d] thcond = %f\n", a, b, thcond);
        }
        // Load vector
        F[a] = Na * 0.0; 
    }
    
    return 0;
}


/*
 * ApplyBoundaryConditions - Sets Dirichlet and Neumann boundary conditions 
 *                           based on user-defined flags.
 *
 * Inputs:
 *   iga  - IGA object defining the domain.
 *   user - Pointer to the user-defined structure containing boundary condition settings.
 */
PetscErrorCode ApplyBoundaryConditions(IGA iga, AppCtx *user) {
    PetscErrorCode ierr;
    PetscFunctionBegin;

    // **Fixed Temperature at y = 0 (Bottom)**
    ierr = IGASetBoundaryValue(iga, 1, 0, 0, user->T_bottom); CHKERRQ(ierr);

    // **Fixed Temperature at y = Ly (Top)**
    ierr = IGASetBoundaryValue(iga, 1, 1, 0, user->T_top); CHKERRQ(ierr);

    // **Zero-Flux BCs on x = 0 and x = Lx**
    ierr = IGASetBoundaryForm(iga, 0, 0, PETSC_TRUE); CHKERRQ(ierr);
    ierr = IGASetBoundaryForm(iga, 0, 1, PETSC_TRUE); CHKERRQ(ierr);

    // **Zero-Flux BCs on z = 0 and z = Lz (Only in 3D)**
    if (user->dim == 3) {
        ierr = IGASetBoundaryForm(iga, 2, 0, PETSC_TRUE); CHKERRQ(ierr);
        ierr = IGASetBoundaryForm(iga, 2, 1, PETSC_TRUE); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}
  

PetscErrorCode IGAComputeFunctionSNES(SNES snes, Vec U, Vec F, void *ctx) 
{
    AppCtx *user = (AppCtx *)ctx;
    return IGAComputeFunction(user->iga, U, F);
}


PetscErrorCode IGAComputeJacobianSNES(SNES snes, Vec U, Mat J, Mat P, void *ctx) 
{
    AppCtx *user = (AppCtx *)ctx;
    return IGAComputeJacobian(user->iga, U, J);
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


/*
 * WriteBinaryFile - Saves the computed temperature field in a binary file.
 *
 * Inputs:
 *   U         - The solution vector storing the temperature field.
 *   filename  - Name of the output binary file.
 */
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

/* */
int main (int argc, char *argv[]) {
    // Vec                 U;
    AppCtx              user;
    PetscErrorCode      ierr;
    // Vec                 K;                      // Thermal conductivity field (debugging)
    // PetscReal           keff;

    /* ------------------ Initialize PETSc ------------------ */
    ierr = PetscInitialize(&argc, &argv, NULL, NULL); CHKERRQ(ierr);

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

    user.eps                    = 2.0*user.Ly/user.Ny; // Interface width parameter for phase field method

    // Basis functions order and continuity
    PetscInt p = 1, C = 0;
    user.p                      = p;        // Polynomial order for basis functions      
    user.C                      = C;        // Continuity of basis functions

    // Initial conditions
    // user.Qsource                = 0.0;      // Source term for heat equation (could be a constant or function of position)
    // user.dQ_dT                  = 0.0;      // Derivative of source term with respect to temperature
    // user.q_flux                 = 100.0;      // Prescribed heat flux at the boundary

    // Set boundary conditions for a simple test case
    user.T_bottom = 265.15;  // -8°C at y = 0  // NEED TO FIX INITAL GUESS!!!!å
    user.T_top    = 270.15;  // -3°C at y = Ly

    /* Define user context */
    PetscPrintf(PETSC_COMM_WORLD, "Initializing thermal diffusion model...\n");

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

    user.iga = iga;

    /* ------------------ Initialize Field Variables ------------------ */
    PetscInt nmb = iga->elem_width[0] * iga->elem_width[1] * SQ(p + 1);
    ierr = PetscMalloc1(nmb, &user.ice); CHKERRQ(ierr);
    ierr = PetscMemzero(user.ice, nmb); CHKERRQ(ierr);
    ierr = InitializeIceField(&user); CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD, "Ice field initialized.\n");

    // Display Ice Field
    // PetscPrintf(PETSC_COMM_WORLD, "Ice Field:\n");
    // for (int i = 0; i < 4.0 * user.Nx * user.Ny; i++) {
    //     PetscPrintf(PETSC_COMM_WORLD, "%f ", user.ice[i]);
    //     if ((i+1) % user.Nx == 0) {
    //         PetscPrintf(PETSC_COMM_WORLD, "\n");
    //     }
    // }
    // PetscPrintf(PETSC_COMM_WORLD, "nmb = %d.\n", nmb);


    /* ------------------ Define Boundary Conditions ------------------ */
    // Apply the boundary conditions
    ierr = ApplyBoundaryConditions(iga, &user); CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD, "Boundary conditions applied.\n");

    /* ------------------ Set Up KSP Solver ------------------ */
    // Creat KSP solver
    Mat A;
    Vec x, b;

    ierr = IGACreateMat(iga, &A); CHKERRQ(ierr);
    ierr = IGACreateVec(iga, &x); CHKERRQ(ierr);
    ierr = IGACreateVec(iga, &b); CHKERRQ(ierr);
    ierr = IGASetFormSystem(iga, AssembleStiffnessMatrix, &user); CHKERRQ(ierr);
    ierr = IGAComputeSystem(iga, A, b); CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD, "System assembled.\n");

    KSP ksp;
    ierr = IGACreateKSP(iga, &ksp); CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp, A, A); CHKERRQ(ierr);
    ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);
    ierr = KSPSolve(ksp, b, x); CHKERRQ(ierr);

    ierr = VecView(x,PETSC_VIEWER_DRAW_WORLD);CHKERRQ(ierr);

    ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
    ierr = MatDestroy(&A); CHKERRQ(ierr);
    ierr = VecDestroy(&x); CHKERRQ(ierr);
    ierr = VecDestroy(&b); CHKERRQ(ierr);
    ierr = IGADestroy(&iga); CHKERRQ(ierr);

    ierr = PetscFinalize(); CHKERRQ(ierr);

    return 0;
}





    // // Create SNES solver
    // ierr = SNESCreate(PETSC_COMM_WORLD, &snes); CHKERRQ(ierr);
    // ierr = SNESSetType(snes, SNESNEWTONLS); CHKERRQ(ierr);

    // // Create solution vector with initial guess
    // ierr = IGACreateVec(iga, &U); CHKERRQ(ierr);

    // // Create the thermal conductivity field for debugging
    // ierr = IGACreateVec(iga, &K); CHKERRQ(ierr);

    // PetscScalar T_init;
    // PetscInt i, size;
    // ierr = VecGetSize(U, &size); CHKERRQ(ierr);

    // // This is not doing what I think it is doing! Size is related to the number of elements in the mesh per core.
    // // Look for initial condition funtion
    // for (i = 0; i < size; i++) {
    //     T_init = user.T_bottom + (user.T_top - user.T_bottom) * ((PetscReal)i / (size - 1));
    //     ierr = VecSetValue(U, i, T_init, INSERT_VALUES); CHKERRQ(ierr);
    // }
    // ierr = VecAssemblyBegin(U); CHKERRQ(ierr);
    // ierr = VecAssemblyEnd(U); CHKERRQ(ierr);

    // // Attach Residual and Jacobian functions to IGA
    // ierr = IGASetFormFunction(iga, Residual, &user); CHKERRQ(ierr);
    // ierr = IGASetFormJacobian(iga, Jacobian, &user); CHKERRQ(ierr);

    // // Link Residual and Jacobian functions to SNES
    // ierr = SNESSetFunction(snes, NULL, IGAComputeFunctionSNES, &user); CHKERRQ(ierr);
    // ierr = SNESSetJacobian(snes, NULL, NULL, IGAComputeJacobianSNES, &user); CHKERRQ(ierr);

    // // Set SNES solver from command line options
    // ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);

    // /* ------------------ Solve System ------------------ */
    // PetscPrintf(PETSC_COMM_WORLD, "Solving nonlinear system...\n");
    // ierr = SNESSolve(snes, NULL, U); CHKERRQ(ierr);

    /* ------------------ Write Output ------------------ */
    // PetscPrintf(PETSC_COMM_WORLD, "Writing VTK output...\n");
    // ierr = WriteVTKFile(iga, U, "temperature.vtk"); CHKERRQ(ierr);
    
    // if (user.outputBinary) {
    //     PetscPrintf(PETSC_COMM_WORLD, "Writing binary...\n");
    //     ierr = WriteBinaryFile(U, "temperature.bin"); CHKERRQ(ierr);
    // }

    // PetscPrintf(PETSC_COMM_WORLD, "Computing and storing thermal conductivity field...\n");
    // ComputeAndStoreThermalConductivity(&user, K);

    // if (user.outputBinary) {
    //     PetscPrintf(PETSC_COMM_WORLD, "Writing thermal conductivity to binary...\n");
    //     WriteBinaryFile(K, "thermal_conductivity.bin");
    // }


    /* ------------------ Compute Effective Conductivity ------------------ */
    // ierr = ComputeEffectiveConductivity(iga, U, &keff, &user); CHKERRQ(ierr);
    // PetscPrintf(PETSC_COMM_WORLD, "Computed Effective Thermal Conductivity: %.6e W/mK\n", keff);

    /* ------------------ Cleanup ------------------ */
    // ierr = VecDestroy(&U); CHKERRQ(ierr);
    // ierr = SNESDestroy(&snes); CHKERRQ(ierr);
    // ierr = IGADestroy(&iga); CHKERRQ(ierr);
    // ierr = PetscFree(user.ice); CHKERRQ(ierr);
    // ierr = PetscFinalize(); CHKERRQ(ierr);
// }
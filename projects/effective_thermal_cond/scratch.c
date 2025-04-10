#include "user_context.h"

static PetscErrorCode ComputeCircleIceField(AppCtx *user,Vec vec_phase)
{
    PetscErrorCode ierr;
    IGAElement element;
    IGAPoint point;
    PetscReal dist;
    PetscInt idx = 0;
    PetscInt total = user->Nx * user->Ny * (user->dim == 3 ? user->Nz : 1);
    Vec               localU;
    const PetscScalar *arrayU;
    PetscScalar       *U;
    /* Create a PETSc Vec to hold the ice field at integration (Gauss) points */
    Vec ice;
    ierr = VecCreate(PETSC_COMM_WORLD, &ice); CHKERRQ(ierr);
    ierr = VecSetSizes(ice, PETSC_DECIDE, total); CHKERRQ(ierr);
    ierr = VecSetFromOptions(ice); CHKERRQ(ierr);
    /* Get array pointer for direct access */
    PetscScalar *iceArray;
    // ierr = VecGetArray(ice, &iceArray); CHKERRQ(ierr);
    ierr = IGAGetLocalVecArray(iga_1dof,vec_phase,&localU,&arrayU);CHKERRQ(ierr);
    ierr = IGABeginElement(user->iga, &element); CHKERRQ(ierr);
    while (IGANextElement(user->iga, element)) {
      ierr = IGAElementGetValues(element,arrayU,&U);CHKERRQ(ierr);
        ierr = IGAElementBeginPoint(element, &point); CHKERRQ(ierr);
        while (IGAElementNextPoint(element, point)) {
            PetscScalar sol[3];
            IGAPointFormValue(point,U,&sol[0]);
            PetscScalar phi = sol[0];
            PetscScalar T = sol[1];
            PetscScalar rho = sol[2];
            iceArray[idx] = phi;
            // PetscReal radius = PetscMin(user->Lx / 6.0, user->Ly / 6.0);
            // /* Compute distance from the circle center */
            // dist = PetscSqrtReal(SQ(point->mapX[0][0] - user->Lx / 2.0) +
            //                      SQ(point->mapX[0][1] - user->Ly / 2.0)) - radius;
            // /* Compute ice phase using a hyperbolic tangent transition */
            // iceArray[idx] = 0.5 - 0.5 * PetscTanhReal(0.5 / user->eps * dist);
            // iceArray[idx] = PetscMax(0.0, PetscMin(1.0, iceArray[idx])); /* Clamp [0,1] */
            /* For debugging */
            // if (idx > total) {
            //     PetscPrintf(PETSC_COMM_WORLD, "There are %d points in ice Vec. Originally set to %d.\n", idx, total);
            // }
            idx++;
        }
        ierr = IGAElementEndPoint(element, &point); CHKERRQ(ierr);
    }
    ierr = IGAEndElement(user->iga, &element); CHKERRQ(ierr);
    ierr = IGARestoreLocalVecArray(iga_1dof,vec_phase,&localU,&arrayU);CHKERRQ(ierr);
    // ierr = VecRestoreArray(ice, &iceArray); CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD, "Circle mode: ice field computed with %d points.\n", idx);
    PetscPrintf(PETSC_COMM_WORLD, "Ice field is %g times larger than expected.\n", (PetscReal)idx / (PetscReal)total);
    /* Assign the computed Vec to user->ice */
    user->ice = ice;
    return 0;
}

/* -------------------------------------------------------------------------- */
IGA iga_1dof;
ierr = IGACreate(PETSC_COMM_WORLD,&iga_1dof);CHKERRQ(ierr);
ierr = IGASetDim(iga_1dof,dim);CHKERRQ(ierr);
ierr = IGASetDof(iga_1dof,1);CHKERRQ(ierr);
ierr = IGASetFieldName(iga_1dof,0,"phaseice"); CHKERRQ(ierr);
IGAAxis axis0, axis1, axis2;
ierr = IGAGetAxis(iga_1dof,0,&axis0);CHKERRQ(ierr);
if(user.periodic==1) {ierr = IGAAxisSetPeriodic(axis0,PETSC_TRUE);CHKERRQ(ierr);}
ierr = IGAAxisSetDegree(axis0,p);CHKERRQ(ierr);
ierr = IGAAxisInitUniform(axis0,Nx,0.0,Lx,C);CHKERRQ(ierr);
ierr = IGAGetAxis(iga_1dof,1,&axis1);CHKERRQ(ierr);
if(user.periodic==1) {ierr = IGAAxisSetPeriodic(axis1,PETSC_TRUE);CHKERRQ(ierr);}
ierr = IGAAxisSetDegree(axis1,p);CHKERRQ(ierr);
ierr = IGAAxisInitUniform(axis1,Ny,0.0,Ly,C);CHKERRQ(ierr);
if(dim==3){
    ierr = IGAGetAxis(iga_1dof,2,&axis2);CHKERRQ(ierr);
    if(user.periodic==1) {ierr = IGAAxisSetPeriodic(axis2,PETSC_TRUE);CHKERRQ(ierr);}
    ierr = IGAAxisSetDegree(axis2,p);CHKERRQ(ierr);
    ierr = IGAAxisInitUniform(axis2,Nz,0.0,Lz,C);CHKERRQ(ierr);
}
ierr = IGASetFromOptions(iga_1dof);CHKERRQ(ierr);
ierr = IGASetUp(iga_1dof);CHKERRQ(ierr);
user.iga_1dof = iga_1dof;
Vec vec_phase;
// ierr = IGACreateVec(iga_1dof,&vec_phase);CHKERRQ(ierr);
ierr = IGAReadVec(iga_1dof,&vec_phase,"path/name.dat");
AssignGaussPoints(iga_1dof,vec_phase,user);
IGADeleteVec(vec_phase);
IGAdelete(iga_1dof);

/* -------------------------------------------------------------------------- */
// Inside the monitor function
Vec localU;
        const PetscScalar *arrayU;
        IGAElement element;
        IGAPoint point;
        PetscScalar *UU;
        PetscInt indd=0;
        ierr = IGAGetLocalVecArray(user->iga,U,&localU,&arrayU);CHKERRQ(ierr);
        ierr = IGABeginElement(user->iga,&element);CHKERRQ(ierr);
        while (IGANextElement(user->iga,element)) {
            ierr = IGAElementGetValues(element,arrayU,&UU);CHKERRQ(ierr);
            ierr = IGAElementBeginPoint(element,&point);CHKERRQ(ierr);
            while (IGAElementNextPoint(element,point)) {
                PetscScalar solS[5];
                ierr = IGAPointFormValue(point,UU,&solS[0]);CHKERRQ(ierr);
                if(user->flag_rad_hcap==1){
                    PetscScalar por0;
                    Porosity0(user,point->mapX[0][dim-1],&por0);
                    if(solS[0]<0.01) solS[0] = 0.01;
                    user->h_c[indd] = user->h_cap*pow(por0/solS[0],0.4);
                }
                if(user->beta_sol<=1.0e-4){
                    PetscScalar Tint, beta_sol = 800.0;
                    if(step>0) beta_sol = user->beta_s[indd];
                    InterfaceTemp(user,beta_sol,solS[3],solS[4],&Tint,NULL,NULL);
                    Tint=fabs(Tint);
                    if (Tint<1.0e-5) Tint=1.0e-5;
                    user->beta_s[indd] = user->cp_wat/user->latheat/1.51e-3/pow(Tint,0.67);
                }
                indd ++;
            }
            ierr = IGAElementEndPoint(element,&point);CHKERRQ(ierr);
        }
        ierr = IGAEndElement(user->iga,&element);CHKERRQ(ierr);
        ierr = IGARestoreLocalVecArray(user->iga,U,&localU,&arrayU);CHKERRQ(ierr);
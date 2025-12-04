#include "initial_conditions.h"
#include "material_properties.h"

PetscErrorCode FormInitialSoil2D(IGA igaS,Vec S,AppCtx *user)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;
    DM da;
    ierr = IGACreateNodeDM(igaS,1,&da);CHKERRQ(ierr);
    FieldS **u;
    ierr = DMDAVecGetArray(da,S,&u);CHKERRQ(ierr);
    DMDALocalInfo info;
    ierr = DMDAGetLocalInfo(da,&info);CHKERRQ(ierr);
    PetscReal dist=0.0,value;
    PetscInt i,j,kk, l=-1;
    if(user->periodic==1) l=user->p-1;
    for(i=info.xs;i<info.xs+info.xm;i++){
        for(j=info.ys;j<info.ys+info.ym;j++){
            PetscReal x = user->Lx*(PetscReal)i / ( (PetscReal)(info.mx+l) );
            PetscReal y = user->Ly*(PetscReal)j / ( (PetscReal)(info.my+l) );
            value=0.0;
            for(kk=0;kk<user->n_actsed;kk++){
                dist = sqrt(SQ(x-user->centsed[0][kk])+SQ(y-user->centsed[1][kk]));
                value += 0.5-0.5*tanh(0.5/user->eps*(dist-user->radiussed[kk]));
            }
            if(value>1.0) value=1.0;

            u[j][i].soil = value;
        }
    }
    ierr = DMDAVecRestoreArray(da,S,&u);CHKERRQ(ierr);
    ierr = DMDestroy(&da);;CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

PetscErrorCode FormInitialSoil3D(IGA igaS,Vec S,AppCtx *user)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;
    DM da;
    ierr = IGACreateNodeDM(igaS,1,&da);CHKERRQ(ierr);
    FieldS ***u;
    ierr = DMDAVecGetArray(da,S,&u);CHKERRQ(ierr);
    DMDALocalInfo info;
    ierr = DMDAGetLocalInfo(da,&info);CHKERRQ(ierr);
    PetscReal dist=0.0,value;
    PetscInt i,j,k, kk, l=-1;
    if(user->periodic==1) l=user->p-1;
    for(i=info.xs;i<info.xs+info.xm;i++){
        for(j=info.ys;j<info.ys+info.ym;j++){
            for(k=info.zs;k<info.zs+info.zm;k++){
                PetscReal x = user->Lx*(PetscReal)i / ( (PetscReal)(info.mx+l) );
                PetscReal y = user->Ly*(PetscReal)j / ( (PetscReal)(info.my+l) );
                PetscReal z = user->Lz*(PetscReal)k / ( (PetscReal)(info.mz+l) );
                value=0.0;
                for(kk=0;kk<user->n_actsed;kk++){
                    dist = sqrt(SQ(x-user->centsed[0][kk])+SQ(y-user->centsed[1][kk])+SQ(z-user->centsed[2][kk]));
                    value += 0.5-0.5*tanh(0.5/user->eps*(dist-user->radiussed[kk]));
                }
                if(value>1.0) value=1.0;
                
                u[k][j][i].soil = value;
            }
        }
    }
    ierr = DMDAVecRestoreArray(da,S,&u);CHKERRQ(ierr);
    ierr = DMDestroy(&da);;CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

PetscErrorCode FormInitialCondition2D(IGA iga, PetscReal t, Vec U,AppCtx *user, 
                                    const char datafile[],const char dataPF[])
{
  PetscErrorCode ierr;
  PetscFunctionBegin;

  if (datafile[0] != 0) { /* initial condition from datafile */
    MPI_Comm comm;
    PetscViewer viewer;
    ierr = PetscObjectGetComm((PetscObject)U,&comm);CHKERRQ(ierr);
    ierr = PetscViewerBinaryOpen(comm,datafile,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
    ierr = VecLoad(U,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);

  } else if (dataPF[0] != 0){
    IGA igaPF;
    ierr = IGACreate(PETSC_COMM_WORLD,&igaPF);CHKERRQ(ierr);
    ierr = IGASetDim(igaPF,2);CHKERRQ(ierr);
    ierr = IGASetDof(igaPF,1);CHKERRQ(ierr);

    IGAAxis axisPF0,axisPF1;
    ierr = IGAGetAxis(igaPF,0,&axisPF0);CHKERRQ(ierr);
    if(user->periodic==1) {ierr = IGAAxisSetPeriodic(axisPF0,PETSC_TRUE);CHKERRQ(ierr);}
    ierr = IGAAxisSetDegree(axisPF0,user->p);CHKERRQ(ierr);
    ierr = IGAAxisInitUniform(axisPF0,user->Nx,0.0,user->Lx,user->C);CHKERRQ(ierr);
    ierr = IGAGetAxis(igaPF,1,&axisPF1);CHKERRQ(ierr);
    if(user->periodic==1) {ierr = IGAAxisSetPeriodic(axisPF1,PETSC_TRUE);CHKERRQ(ierr);}
    ierr = IGAAxisSetDegree(axisPF1,user->p);CHKERRQ(ierr);
    ierr = IGAAxisInitUniform(axisPF1,user->Ny,0.0,user->Ly,user->C);CHKERRQ(ierr);
    ierr = IGASetFromOptions(igaPF);CHKERRQ(ierr);
    ierr = IGASetUp(igaPF);CHKERRQ(ierr);

    Vec PF;
    ierr = IGACreateVec(igaPF,&PF);CHKERRQ(ierr);
    
    MPI_Comm comm;
    PetscViewer viewer;
    ierr = PetscObjectGetComm((PetscObject)PF,&comm);CHKERRQ(ierr);
    ierr = PetscViewerBinaryOpen(comm,dataPF,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
    ierr = VecLoad(PF,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);
    ierr = VecStrideScatter(PF,0,U,INSERT_VALUES);
    ierr = VecDestroy(&PF);CHKERRQ(ierr);
    ierr = IGADestroy(&igaPF);CHKERRQ(ierr);

    DM da;
    ierr = IGACreateNodeDM(iga,3,&da);CHKERRQ(ierr);
    Field **u;
    ierr = DMDAVecGetArray(da,U,&u);CHKERRQ(ierr);
    DMDALocalInfo info;
    ierr = DMDAGetLocalInfo(da,&info);CHKERRQ(ierr);
    PetscInt i,j, k=-1;
    if(user->periodic==1) k=user->p -1;
    for(i=info.xs;i<info.xs+info.xm;i++){
      for(j=info.ys;j<info.ys+info.ym;j++){
        PetscReal x = user->Lx*(PetscReal)i / ( (PetscReal)(info.mx+k) );
        PetscReal y = user->Ly*(PetscReal)j / ( (PetscReal)(info.my+k) );

        u[j][i].tem = user->temp0 + user->grad_temp0[0]*(x-0.5*user->Lx) + user->grad_temp0[1]*(y-0.5*user->Ly);
        PetscScalar rho_vs, temp=u[j][i].tem;
        RhoVS_I(user,temp,&rho_vs,NULL);
        u[j][i].rhov = user->hum0*rho_vs;
      }
    }
    ierr = DMDAVecRestoreArray(da,U,&u);CHKERRQ(ierr); 
    ierr = DMDestroy(&da);;CHKERRQ(ierr); 

  } else {
    DM da;
    ierr = IGACreateNodeDM(iga,3,&da);CHKERRQ(ierr);
    Field **u;
    ierr = DMDAVecGetArray(da,U,&u);CHKERRQ(ierr);
    DMDALocalInfo info;
    ierr = DMDAGetLocalInfo(da,&info);CHKERRQ(ierr);

    PetscInt i,j,k=-1;
    if(user->periodic==1) k=user->p -1;
    for(i=info.xs;i<info.xs+info.xm;i++){
      for(j=info.ys;j<info.ys+info.ym;j++){
        PetscReal x = user->Lx*(PetscReal)i / ( (PetscReal)(info.mx+k) );
        PetscReal y = user->Ly*(PetscReal)j / ( (PetscReal)(info.my+k) );

        PetscReal dist,ice=0.0;
        PetscInt aa;
        for(aa=0;aa<user->n_act;aa++){
          dist=sqrt(SQ(x-user->cent[0][aa])+SQ(y-user->cent[1][aa]));
          ice += 0.5-0.5*tanh(0.5/user->eps*(dist-user->radius[aa]));
        }
        if(ice>1.0) ice=1.0;

        u[j][i].ice = ice;    
        u[j][i].tem = user->temp0 + user->grad_temp0[0]*(x-0.5*user->Lx) + user->grad_temp0[1]*(y-0.5*user->Ly);
        PetscScalar rho_vs, temp=u[j][i].tem;
        RhoVS_I(user,temp,&rho_vs,NULL);
        u[j][i].rhov = user->hum0*rho_vs;
      }
    }
    ierr = DMDAVecRestoreArray(da,U,&u);CHKERRQ(ierr); 
    ierr = DMDestroy(&da);;CHKERRQ(ierr); 
  }
  PetscFunctionReturn(0); 
}

PetscErrorCode FormInitialCondition3D(IGA iga, PetscReal t, Vec U,AppCtx *user, 
                                    const char datafile[],const char dataPF[])
{
  PetscErrorCode ierr;
  PetscFunctionBegin;

  if (datafile[0] != 0) { /* initial condition from datafile */
    MPI_Comm comm;
    PetscViewer viewer;
    ierr = PetscObjectGetComm((PetscObject)U,&comm);CHKERRQ(ierr);
    ierr = PetscViewerBinaryOpen(comm,datafile,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
    ierr = VecLoad(U,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);

  } else if (dataPF[0] != 0){
    IGA igaPF;
    ierr = IGACreate(PETSC_COMM_WORLD,&igaPF);CHKERRQ(ierr);
    ierr = IGASetDim(igaPF,3);CHKERRQ(ierr);
    ierr = IGASetDof(igaPF,1);CHKERRQ(ierr);
    IGAAxis axisPF0, axisPF1, axisPF2;
    ierr = IGAGetAxis(igaPF,0,&axisPF0);CHKERRQ(ierr);
    if(user->periodic==1) {ierr = IGAAxisSetPeriodic(axisPF0,PETSC_TRUE);CHKERRQ(ierr);}
    ierr = IGAAxisSetDegree(axisPF0,user->p);CHKERRQ(ierr);
    ierr = IGAAxisInitUniform(axisPF0,user->Nx,0.0,user->Lx,user->C);CHKERRQ(ierr);
    ierr = IGAGetAxis(igaPF,1,&axisPF1);CHKERRQ(ierr);
    if(user->periodic==1) {ierr = IGAAxisSetPeriodic(axisPF1,PETSC_TRUE);CHKERRQ(ierr);}
    ierr = IGAAxisSetDegree(axisPF1,user->p);CHKERRQ(ierr);
    ierr = IGAAxisInitUniform(axisPF1,user->Ny,0.0,user->Ly,user->C);CHKERRQ(ierr);
    ierr = IGAGetAxis(igaPF,2,&axisPF2);CHKERRQ(ierr);
    if(user->periodic==1) {ierr = IGAAxisSetPeriodic(axisPF2,PETSC_TRUE);CHKERRQ(ierr);}
    ierr = IGAAxisSetDegree(axisPF2,user->p);CHKERRQ(ierr);
    ierr = IGAAxisInitUniform(axisPF2,user->Nz,0.0,user->Lz,user->C);CHKERRQ(ierr);
    ierr = IGASetFromOptions(igaPF);CHKERRQ(ierr);
    ierr = IGASetUp(igaPF);CHKERRQ(ierr);
    Vec PF;
    ierr = IGACreateVec(igaPF,&PF);CHKERRQ(ierr);
    MPI_Comm comm;
    PetscViewer viewer;
    ierr = PetscObjectGetComm((PetscObject)PF,&comm);CHKERRQ(ierr);
    ierr = PetscViewerBinaryOpen(comm,dataPF,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
    ierr = VecLoad(PF,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);
    ierr = VecStrideScatter(PF,0,U,INSERT_VALUES);
    ierr = VecDestroy(&PF);CHKERRQ(ierr);
    ierr = IGADestroy(&igaPF);CHKERRQ(ierr);

    DM da;
    ierr = IGACreateNodeDM(iga,3,&da);CHKERRQ(ierr);
    Field ***u;
    ierr = DMDAVecGetArray(da,U,&u);CHKERRQ(ierr);
    DMDALocalInfo info;
    ierr = DMDAGetLocalInfo(da,&info);CHKERRQ(ierr);
    PetscInt i,j,k, l=-1;
    if(user->periodic==1) l=user->p -1;
    for(i=info.xs;i<info.xs+info.xm;i++){
      for(j=info.ys;j<info.ys+info.ym;j++){
        for(k=info.zs;k<info.zs+info.zm;k++){
          PetscReal x = user->Lx*(PetscReal)i / ( (PetscReal)(info.mx+l) );
          PetscReal y = user->Ly*(PetscReal)j / ( (PetscReal)(info.my+l) );
          PetscReal z = user->Lz*(PetscReal)k / ( (PetscReal)(info.mz+l) );

          u[k][j][i].tem = user->temp0 + user->grad_temp0[0]*(x-0.5*user->Lx) + user->grad_temp0[1]*(y-0.5*user->Ly) + user->grad_temp0[2]*(z-0.5*user->Lz);
          PetscScalar rho_vs, temp=u[k][j][i].tem;
          RhoVS_I(user,temp,&rho_vs,NULL);
          u[k][j][i].rhov = user->hum0*rho_vs;
        }
      }
    }
    ierr = DMDAVecRestoreArray(da,U,&u);CHKERRQ(ierr); 
    ierr = DMDestroy(&da);;CHKERRQ(ierr); 

  } else {
    DM da;
    ierr = IGACreateNodeDM(iga,3,&da);CHKERRQ(ierr);
    Field ***u;
    ierr = DMDAVecGetArray(da,U,&u);CHKERRQ(ierr);
    DMDALocalInfo info;
    ierr = DMDAGetLocalInfo(da,&info);CHKERRQ(ierr);

    PetscInt i,j,k, l=-1;
    if(user->periodic==1) k=user->p -1;
    for(i=info.xs;i<info.xs+info.xm;i++){
      for(j=info.ys;j<info.ys+info.ym;j++){
        for(k=info.zs;k<info.zs+info.zm;k++){
          PetscReal x = user->Lx*(PetscReal)i / ( (PetscReal)(info.mx+l) );
          PetscReal y = user->Ly*(PetscReal)j / ( (PetscReal)(info.my+l) );
          PetscReal z = user->Lz*(PetscReal)k / ( (PetscReal)(info.mz+l) );

          PetscReal dist,ice=0.0;
          PetscInt aa;
          for(aa=0;aa<user->n_act;aa++){
            dist=sqrt(SQ(x-user->cent[0][aa])+SQ(y-user->cent[1][aa])+SQ(z-user->cent[2][aa]));
            ice += 0.5-0.5*tanh(0.5/user->eps*(dist-user->radius[aa]));
          }
          if(ice>1.0) ice=1.0;

          u[k][j][i].ice = ice;    
          u[k][j][i].tem = user->temp0 + user->grad_temp0[0]*(x-0.5*user->Lx) + user->grad_temp0[1]*(y-0.5*user->Ly) + user->grad_temp0[2]*(z-0.5*user->Lz);
          PetscScalar rho_vs, temp=u[k][j][i].tem;
          RhoVS_I(user,temp,&rho_vs,NULL);
          u[k][j][i].rhov = user->hum0*rho_vs;
        }
      }
    }
    ierr = DMDAVecRestoreArray(da,U,&u);CHKERRQ(ierr); 
    ierr = DMDestroy(&da);;CHKERRQ(ierr); 
  }
  PetscFunctionReturn(0); 
}

/* Create a 2D System that has an ice grain-pack sitting above a solid block of 
   ice with air inclusions inside of it.
   Note: If we call this function, we assume that the ice grain-pack is already 
   intialized only in the top half of the domain. Here we will assign the 
   phase-field values for those grains and for the top of the domain. The last 
   thing we do is remove ice 'grains' for the inclusions.
   The best way to initialize ice grains in the top half will be to read it form
   a file. We should do the same for the inclusions. The data for the ice grains
   and the data for the inclusions should be in the same file. The defining 
   difference will be that the inclusions will have coordinates that are below 
   Ly/2 and the ice grains will have coordinates that are above Ly/2.
*/
PetscErrorCode FormLayeredInitialCondition2D(IGA iga, PetscReal t, Vec U, 
                                            AppCtx *user, const char datafile[],
                                            const char dataPF[])
{
  PetscErrorCode ierr;
  PetscFunctionBegin;


  if (datafile[0] != 0) { /* initial condition from datafile */
    MPI_Comm comm;
    PetscViewer viewer;
    ierr = PetscObjectGetComm((PetscObject)U,&comm);CHKERRQ(ierr);
    ierr = PetscViewerBinaryOpen(comm,datafile,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
    ierr = VecLoad(U,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);

  } else if (dataPF[0] != 0){
    IGA igaPF;
    ierr = IGACreate(PETSC_COMM_WORLD,&igaPF);CHKERRQ(ierr);
    ierr = IGASetDim(igaPF,2);CHKERRQ(ierr);
    ierr = IGASetDof(igaPF,1);CHKERRQ(ierr);
    IGAAxis axisPF0,axisPF1;
    ierr = IGAGetAxis(igaPF,0,&axisPF0);CHKERRQ(ierr);
    if(user->periodic==1) {ierr = IGAAxisSetPeriodic(axisPF0,PETSC_TRUE);CHKERRQ(ierr);}
    ierr = IGAAxisSetDegree(axisPF0,user->p);CHKERRQ(ierr);
    ierr = IGAAxisInitUniform(axisPF0,user->Nx,0.0,user->Lx,user->C);CHKERRQ(ierr);
    ierr = IGAGetAxis(igaPF,1,&axisPF1);CHKERRQ(ierr);
    if(user->periodic==1) {ierr = IGAAxisSetPeriodic(axisPF1,PETSC_TRUE);CHKERRQ(ierr);}
    ierr = IGAAxisSetDegree(axisPF1,user->p);CHKERRQ(ierr);
    ierr = IGAAxisInitUniform(axisPF1,user->Ny,0.0,user->Ly,user->C);CHKERRQ(ierr);
    ierr = IGASetFromOptions(igaPF);CHKERRQ(ierr);
    ierr = IGASetUp(igaPF);CHKERRQ(ierr);
    Vec PF;
    ierr = IGACreateVec(igaPF,&PF);CHKERRQ(ierr);
    MPI_Comm comm;
    PetscViewer viewer;
    ierr = PetscObjectGetComm((PetscObject)PF,&comm);CHKERRQ(ierr);
    ierr = PetscViewerBinaryOpen(comm,dataPF,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
    ierr = VecLoad(PF,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);
    ierr = VecStrideScatter(PF,0,U,INSERT_VALUES);
    ierr = VecDestroy(&PF);CHKERRQ(ierr);
    ierr = IGADestroy(&igaPF);CHKERRQ(ierr);

    DM da;
    ierr = IGACreateNodeDM(iga,3,&da);CHKERRQ(ierr);
    Field **u;
    ierr = DMDAVecGetArray(da,U,&u);CHKERRQ(ierr);
    DMDALocalInfo info;
    ierr = DMDAGetLocalInfo(da,&info);CHKERRQ(ierr);
    PetscInt i,j, k=-1;
    if(user->periodic==1) k=user->p -1;
    for(i=info.xs;i<info.xs+info.xm;i++){
      for(j=info.ys;j<info.ys+info.ym;j++){
        PetscReal x = user->Lx*(PetscReal)i / ( (PetscReal)(info.mx+k) );
        PetscReal y = user->Ly*(PetscReal)j / ( (PetscReal)(info.my+k) );

        u[j][i].tem = user->temp0 + user->grad_temp0[0]*(x-0.5*user->Lx) + user->grad_temp0[1]*(y-0.5*user->Ly);
        PetscScalar rho_vs, temp=u[j][i].tem;
        RhoVS_I(user,temp,&rho_vs,NULL);
        u[j][i].rhov = user->hum0*rho_vs;
      }
    }
    ierr = DMDAVecRestoreArray(da,U,&u);CHKERRQ(ierr); 
    ierr = DMDestroy(&da);;CHKERRQ(ierr); 



  } else {
    /* Permafrost Intializatoin*/
    DM da;
    ierr = IGACreateNodeDM(iga,3,&da);CHKERRQ(ierr);
    Field **u;
    ierr = DMDAVecGetArray(da,U,&u);CHKERRQ(ierr);
    DMDALocalInfo info;
    ierr = DMDAGetLocalInfo(da,&info);CHKERRQ(ierr);

    PetscInt        i, j, k = -1, aa;
    PetscReal       dist, ice = 0.0;

    ierr = IGACreateNodeDM(iga, 3, &da); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da, U, &u); CHKERRQ(ierr);
    ierr = DMDAGetLocalInfo(da, &info); CHKERRQ(ierr);

    if (user->periodic == 1) k = user->p - 1;

    for (i = info.xs; i < info.xs + info.xm; i++) {
      for (j = info.ys; j < info.ys + info.ym; j++) {
        PetscReal x = user->Lx * (PetscReal)i / ((PetscReal)(info.mx + k));
        PetscReal y = user->Ly * (PetscReal)j / ((PetscReal)(info.my + k));

        // Initialize temperature and density fields
        u[j][i].tem = user->temp0 + user->grad_temp0[0] * (x - 0.5 * user->Lx) + user->grad_temp0[1] * (y - 0.5 * user->Ly);
        PetscScalar rho_vs, temp = u[j][i].tem;
        RhoVS_I(user, temp, &rho_vs, NULL);
        u[j][i].rhov = user->hum0 * rho_vs;

        // Initialize phase-field variable for ice
        ice = 0.0;
        for(aa = 0; aa < user->n_act; aa++) {
          dist = sqrt(SQ(x - user->cent[0][aa]) + SQ(y - user->cent[1][aa]));
          ice += 0.5 - 0.5 * tanh(0.5 / user->eps * (dist - user->radius[aa]));
        }
        
        ice = PetscMin(PetscMax(ice, 0.0), 1.0); // Clamp ice value between 0 and 1

        u[j][i].ice = ice;
      }
    }

    ierr = DMDAVecRestoreArray(da,U,&u);CHKERRQ(ierr); 
    ierr = DMDestroy(&da);;CHKERRQ(ierr); 
  } 

  /* The top half of the domain is initialized with ice grains. The bottom half of */
  // } else {
  //   DM da;
  //   ierr = IGACreateNodeDM(iga,3,&da);CHKERRQ(ierr);
  //   Field **u;
  //   ierr = DMDAVecGetArray(da,U,&u);CHKERRQ(ierr);
  //   DMDALocalInfo info;
  //   ierr = DMDAGetLocalInfo(da,&info);CHKERRQ(ierr);

  //   /* Initialize the bottom half of the domain with solid ice. */
  //   PetscInt i,j,k=-1;
  //   if(user->periodic==1) k=user->p -1;
  //   // Loop over the domain
  //   for(i=info.xs;i<info.xs+info.xm;i++){
  //     for(j=info.ys;j<info.ys+info.ym;j++){
  //       // Define the coordinates
  //       PetscReal x = user->Lx*(PetscReal)i / ( (PetscReal)(info.mx+k) );
  //       PetscReal y = user->Ly*(PetscReal)j / ( (PetscReal)(info.my+k) );

  //       // Initialize the ice phase-field variable for the top (air) and bottom 
  //       // (ice) layers.
  //       PetscReal dist, ice=0.0;
  //       dist = y - (user->Ly/2.0);
  //       ice = 0.5-0.5*tanh(0.5/user->eps*dist);

  //       // Remove the air inclusions
  //       PetscInt aa;
  //       for(aa=0;aa<user->n_act;aa++){
  //         dist=sqrt(SQ(x-user->cent[0][aa])+SQ(y-user->cent[1][aa]));

  //         if((user->cent[1][aa] - 0.9*user->radius[aa]) < user->Ly/2.0){
  //           // Remove the air inclusion
  //           // ice -= 0.5-0.5*tanh(0.5/user->eps*(dist-user->radius[aa]));
  //           // PetscPrintf(PETSC_COMM_WORLD, "Intialized air inclusion %d at (%.2e, %.2e)\n", aa, user->cent[0][aa], user->cent[1][aa]);
  //         } else {
  //           ice += 0.5-0.5*tanh(0.5/user->eps*(dist-user->radius[aa]));
  //           // PetscPrintf(PETSC_COMM_WORLD, "Intialized ice grain %d at (%.2e, %.2e)\n", aa, user->cent[0][aa], user->cent[1][aa]);
  //         }
  //       }
  //       if(ice>1.0) ice=1.0;
  //       if(ice<0.0) ice=0.0;

  //       u[j][i].ice = ice;    
  //       u[j][i].tem = user->temp0 + user->grad_temp0[0]*(x-0.5*user->Lx) + user->grad_temp0[1]*(y-0.5*user->Ly);
  //       PetscScalar rho_vs, temp=u[j][i].tem;
  //       RhoVS_I(user,temp,&rho_vs,NULL);
  //       u[j][i].rhov = user->hum0*rho_vs;
  //     }
  //   }

  //   ierr = DMDAVecRestoreArray(da,U,&u);CHKERRQ(ierr); 
  //   ierr = DMDestroy(&da);;CHKERRQ(ierr); 
  // }
  PetscFunctionReturn(0); 
}















/* 
   LoadInputSolutionVec:
   - Loads Sthavishtha's solution Vec from file
*/
PetscErrorCode LoadInputSolutionVec(const char *filename, Vec *U_in_seq)
{
  PetscErrorCode ierr;
  PetscViewer    viewer;
  PetscFunctionBegin;

  if (!filename || filename[0] == '\0') {
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG,
            "LoadInputSolutionVec: filename is empty");
  }

  /* This is rank-0 only, so use PETSC_COMM_SELF */
  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF, filename,
                               FILE_MODE_READ, &viewer);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_SELF, U_in_seq);CHKERRQ(ierr);
  ierr = VecLoad(*U_in_seq, viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/* 
   InitializeFromInputSolution:
   - Reads Sthavishtha's solution Vec (ice, sediment, temp)
   - Sets your U (ice, temp, rhov)
   - Fills user->Phi_sed with sediment
*/
PetscErrorCode InitializeFromInputSolution(IGA iga, Vec U, Vec S, AppCtx *user)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;

    const char *solutionFile = user->initial_cond;
    if (!solutionFile || solutionFile[0] == '\0') {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
                "InitializeFromInputSolution: user->initial_cond is empty");
    }

    /* 1) Load Sthavishtha's solution Vec in parallel (same IGA, same layout) */
    MPI_Comm   comm;
    Vec        U_in;
    PetscInt   Nloc_in, Nloc_out, Nloc_s, dof_out;

    ierr = PetscObjectGetComm((PetscObject)U, &comm);CHKERRQ(ierr);

    ierr = IGACreateVec(iga, &U_in);CHKERRQ(ierr);

    PetscViewer viewer;
    ierr = PetscViewerBinaryOpen(comm, solutionFile, FILE_MODE_READ, &viewer);CHKERRQ(ierr);
    ierr = VecLoad(U_in, viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

    /* 2) Sanity checks on sizes and dofs */
    ierr = VecGetLocalSize(U_in, &Nloc_in);CHKERRQ(ierr);
    ierr = VecGetLocalSize(U,    &Nloc_out);CHKERRQ(ierr);
    ierr = IGAGetDof(iga,        &dof_out);CHKERRQ(ierr); /* should be 3 (ice, tem, rhov) */

    if (Nloc_in != Nloc_out) {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_SIZ,
                "InitializeFromInputSolution: local sizes mismatch: in = %d, out = %d",
                (int)Nloc_in, (int)Nloc_out);
    }
    if (Nloc_out % dof_out != 0) {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_SIZ,
                "InitializeFromInputSolution: Nloc_out = %d not divisible by dof_out = %d",
                (int)Nloc_out, (int)dof_out);
    }

    PetscInt nNodeLocal = Nloc_out / dof_out;

    ierr = VecGetLocalSize(S, &Nloc_s);CHKERRQ(ierr);
    if (Nloc_s != nNodeLocal) {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_SIZ,
                "InitializeFromInputSolution: S local size %d != nNodeLocal %d",
                (int)Nloc_s, (int)nNodeLocal);
    }

    /* 3) Access arrays: input solution, output solution, sediment Vec */
    const PetscScalar *ain  = NULL;
    PetscScalar       *aout = NULL;
    PetscScalar       *s_arr = NULL;

    ierr = VecGetArrayRead(U_in, &ain);CHKERRQ(ierr);
    ierr = VecGetArray(U,       &aout);CHKERRQ(ierr);
    ierr = VecGetArray(S,       &s_arr);CHKERRQ(ierr);

    /* 4) Precompute uniform rhov0 used as initial guess */
    PetscReal rho_vs, d_rhovs;
    RhoVS_I(user, user->temp0, &rho_vs, &d_rhovs);
    PetscReal rhov0 = user->hum0 * rho_vs;

    /* Input component mapping (Sthavishtha): [0]=ice, [1]=sediment, [2]=temperature */
    const PetscInt ice_in_comp = 0;
    const PetscInt sed_in_comp = 1;
    /* const PetscInt temp_in_comp = 2;  <-- currently ignored */

    /* Output mapping: [0]=ice, [1]=temp, [2]=rhov */
    const PetscInt ice_out_comp  = 0;
    const PetscInt rhov_out_comp = 2;

    /* 5) Loop over local nodes and map phases */
    for (PetscInt inode = 0; inode < nNodeLocal; inode++) {
        PetscInt base = inode * dof_out; /* local base index in both ain and aout */

        PetscReal ice_in = PetscRealPart(ain[base + ice_in_comp]);
        PetscReal sed_in = PetscRealPart(ain[base + sed_in_comp]);

        /* Clamp sediment to [0,1] */
        PetscReal sed_clamped = PetscMin(PetscMax(sed_in, 0.0), 1.0);

        /* Copy ice phase */
        aout[base + ice_out_comp]  = (PetscScalar)ice_in;

        /* Set vapor density to uniform rhov0 for now;
           we will overwrite based on T(x) below. */
        aout[base + rhov_out_comp] = (PetscScalar)rhov0;

        /* Store sediment in S and in user->Phi_sed */
        s_arr[inode]        = (PetscScalar)sed_clamped;
        user->Phi_sed[inode] = sed_clamped;
    }

    ierr = VecRestoreArrayRead(U_in, &ain);CHKERRQ(ierr);
    ierr = VecRestoreArray(U,       &aout);CHKERRQ(ierr);
    ierr = VecRestoreArray(S,       &s_arr);CHKERRQ(ierr);

    ierr = VecDestroy(&U_in);CHKERRQ(ierr);

    user->n_actsed = user->NCsed;

    /* 6) Overwrite temperature and rhov based on temp0 + grad*(x-0.5L) */
    DM            da;
    DMDALocalInfo info;

    if (user->dim == 2) {
        Field **u;
        PetscInt k = -1;

        ierr = IGACreateNodeDM(iga, 3, &da);CHKERRQ(ierr);
        ierr = DMDAGetLocalInfo(da, &info);CHKERRQ(ierr);
        ierr = DMDAVecGetArray(da, U, &u);CHKERRQ(ierr);

        if (user->periodic == 1) k = user->p - 1;

        for (PetscInt i = info.xs; i < info.xs + info.xm; i++) {
            for (PetscInt j = info.ys; j < info.ys + info.ym; j++) {
                PetscReal x = user->Lx * (PetscReal)i / (PetscReal)(info.mx + k);
                PetscReal y = user->Ly * (PetscReal)j / (PetscReal)(info.my + k);

                u[j][i].tem =
                    user->temp0
                    + user->grad_temp0[0] * (x - 0.5 * user->Lx)
                    + user->grad_temp0[1] * (y - 0.5 * user->Ly);

                PetscScalar rho_vs_loc, temp_loc = u[j][i].tem;
                RhoVS_I(user, temp_loc, &rho_vs_loc, NULL);
                u[j][i].rhov = user->hum0 * rho_vs_loc;
            }
        }

        ierr = DMDAVecRestoreArray(da, U, &u);CHKERRQ(ierr);
        ierr = DMDestroy(&da);CHKERRQ(ierr);

    } else if (user->dim == 3) {
        Field ***u;
        PetscInt l = -1;

        ierr = IGACreateNodeDM(iga, 3, &da);CHKERRQ(ierr);
        ierr = DMDAGetLocalInfo(da, &info);CHKERRQ(ierr);
        ierr = DMDAVecGetArray(da, U, &u);CHKERRQ(ierr);

        if (user->periodic == 1) l = user->p - 1;

        for (PetscInt i = info.xs; i < info.xs + info.xm; i++) {
            for (PetscInt j = info.ys; j < info.ys + info.ym; j++) {
                for (PetscInt k = info.zs; k < info.zs + info.zm; k++) {
                    PetscReal x = user->Lx * (PetscReal)i / (PetscReal)(info.mx + l);
                    PetscReal y = user->Ly * (PetscReal)j / (PetscReal)(info.my + l);
                    PetscReal z = user->Lz * (PetscReal)k / (PetscReal)(info.mz + l);

                    u[k][j][i].tem =
                        user->temp0
                        + user->grad_temp0[0] * (x - 0.5 * user->Lx)
                        + user->grad_temp0[1] * (y - 0.5 * user->Ly)
                        + user->grad_temp0[2] * (z - 0.5 * user->Lz);

                    PetscScalar rho_vs_loc, temp_loc = u[k][j][i].tem;
                    RhoVS_I(user, temp_loc, &rho_vs_loc, NULL);
                    u[k][j][i].rhov = user->hum0 * rho_vs_loc;
                }
            }
        }

        ierr = DMDAVecRestoreArray(da, U, &u);CHKERRQ(ierr);
        ierr = DMDestroy(&da);CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}
#include <petsc/private/tsimpl.h>
#include <petsc/private/snesimpl.h>
#include "petiga.h"
#define SQ(x) ((x)*(x))
#define CU(x) ((x)*(x)*(x))

typedef struct {
  IGA       iga;
  // problem parameters
  PetscReal eps,nucleat;
  PetscReal mob_sol,mob_sub,mob_eva,mav,Etai,Etaw,Etaa,alph_sol,alph_sub,alph_eva,Lambd;
  PetscReal thcond_ice,thcond_wat,thcond_air,cp_ice,cp_wat,cp_air,rho_ice,rho_wat,rho_air,dif_vap,lat_sol,lat_sub;
  PetscReal phi_L,air_lim,xi_v,xi_T;
  PetscReal T_melt,temp0,grad_temp0[2],tem_nucl,temp_m_ampl,temp_m_fre,costhet;
  PetscReal Lx,Ly,Lz,Nx,Ny,Nz;
  PetscReal norm0_0,norm0_1,norm0_2,norm0_3;
  PetscInt  flag_it0, flag_tIC, outp, nsteps_IC,flag_xiT,flag_contang;
  PetscInt  xiT_count;
  PetscInt  p,C, periodic,BC_Tfix,BC_Tvar;
  PetscReal t_out,t_interv,t_IC;
  PetscInt  *FlagMob;
  PetscInt  NCice, n_act, seed;
  PetscReal cent[2][200], radius[200], overl, RCice, RCice_dev;
} AppCtx;

PetscErrorCode SNESDOFConvergence(SNES snes,PetscInt it_number,PetscReal xnorm,PetscReal gnorm,PetscReal fnorm,SNESConvergedReason *reason,void *cctx)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;
  AppCtx *user = (AppCtx *)cctx;

  Vec Res,Sol,Sol_upd;
  PetscScalar n2dof0,n2dof1,n2dof2,n2dof3;
  PetscScalar solv,solupdv;

  ierr = SNESGetFunction(snes,&Res,0,0);CHKERRQ(ierr);
  ierr = VecStrideNorm(Res,0,NORM_2,&n2dof0);CHKERRQ(ierr);
  ierr = VecStrideNorm(Res,1,NORM_2,&n2dof1);CHKERRQ(ierr);
  ierr = VecStrideNorm(Res,2,NORM_2,&n2dof2);CHKERRQ(ierr);
  ierr = VecStrideNorm(Res,3,NORM_2,&n2dof3);CHKERRQ(ierr);
  if(user->flag_tIC == 1){
    ierr = SNESGetSolution(snes,&Sol);CHKERRQ(ierr);
    ierr = VecStrideNorm(Sol,3,NORM_2,&solv);CHKERRQ(ierr);
    ierr = SNESGetSolutionUpdate(snes,&Sol_upd);CHKERRQ(ierr);
    ierr = VecStrideNorm(Sol_upd,3,NORM_2,&solupdv);CHKERRQ(ierr);
  }

  if(it_number==0) {
    user->norm0_0 = n2dof0;
    user->norm0_1 = n2dof1;
    user->norm0_2 = n2dof2;
    user->norm0_3 = n2dof3;
    if(user->flag_tIC == 1) solupdv = solv;  
  }

  PetscPrintf(PETSC_COMM_WORLD,"    IT_NUMBER: %d ", it_number);
  PetscPrintf(PETSC_COMM_WORLD,"    fnorm: %.4e \n", fnorm);
  PetscPrintf(PETSC_COMM_WORLD,"    n0: %.2e r %.1e ", n2dof0, n2dof0/user->norm0_0);
  PetscPrintf(PETSC_COMM_WORLD,"  n1: %.2e r %.1e ", n2dof1, n2dof1/user->norm0_1);
  PetscPrintf(PETSC_COMM_WORLD,"  n2: %.2e r %.1e ", n2dof2, n2dof2/user->norm0_2);
  if(user->flag_tIC == 1) PetscPrintf(PETSC_COMM_WORLD,"  x3: %.2e s %.1e \n", solv, solupdv/solv);
  else PetscPrintf(PETSC_COMM_WORLD,"  n3: %.2e r %.1e \n", n2dof3, n2dof3/user->norm0_3);

  PetscScalar atol, rtol, stol;
  PetscInt maxit, maxf;
  ierr = SNESGetTolerances(snes,&atol,&rtol,&stol,&maxit,&maxf);
  if(snes->prev_dt_red ==1) rtol *= 10.0;

  if(user->flag_it0 == 1){
    atol = 1.0e-12;
    if ( (n2dof0 <= rtol*user->norm0_0 || n2dof0 < atol) 
      && (n2dof1 <= rtol*user->norm0_1 || n2dof1 < atol)  
      && (n2dof2 <= rtol*user->norm0_2 || n2dof2 < atol) 
      && (n2dof3 <= rtol*user->norm0_3 || n2dof3 < atol) ) {

    *reason = SNES_CONVERGED_FNORM_RELATIVE;
    }    
  } else {
    atol = 1.0e-25;
    if ( (n2dof0 <= rtol*user->norm0_0 || n2dof0 < atol) 
      && (n2dof1 <= rtol*user->norm0_1 || n2dof1 < atol)  
      && (n2dof2 <= rtol*user->norm0_2 || n2dof2 < atol) 
      && (n2dof3 <= rtol*user->norm0_3 || n2dof3 < atol) ) {

    *reason = SNES_CONVERGED_FNORM_RELATIVE;
    }     
  }

  PetscFunctionReturn(0);
}


void ThermalCond(AppCtx *user, PetscScalar ice, PetscScalar air, PetscScalar *cond, PetscScalar *dcond_ice, PetscScalar *dcond_air)
{
  PetscReal dice=1.0, dair=1.0, dwat=1.0;
  PetscReal wat = 1.0-air-ice;
  if(wat<0.0) {wat=0.0;dwat=0.0;}
  if(ice<0.0) {ice=0.0;dice=0.0;}
  if(air<0.0) {air=0.0;dair=0.0;}
  PetscReal cond_ice = user->thcond_ice;
  PetscReal cond_wat = user->thcond_wat;
  PetscReal cond_air = user->thcond_air;
  if(cond)      (*cond)  = ice*cond_ice + wat*cond_wat + air*cond_air;
  if(dcond_ice)    (*dcond_ice) = cond_ice*dice-cond_wat*dwat;
  if(dcond_air)    (*dcond_air) = cond_air*dair-cond_wat*dwat;

  return;
}

void HeatCap(AppCtx *user, PetscScalar ice, PetscScalar air, PetscScalar *cp, PetscScalar *dcp_ice, PetscScalar *dcp_air)
{
  PetscReal dice=1.0, dair=1.0, dwat=1.0;
  PetscReal wat = 1.0-air-ice;
  if(wat<0.0) {wat=0.0;dwat=0.0;}
  if(ice<0.0) {ice=0.0;dice=0.0;}
  if(air<0.0) {air=0.0;dair=0.0;}
  PetscReal cp_ice = user->cp_ice;
  PetscReal cp_wat = user->cp_wat;
  PetscReal cp_air = user->cp_air;
  if(cp)     (*cp)  = ice*cp_ice + wat*cp_wat + air*cp_air;
  if(dcp_ice)    (*dcp_ice) = cp_ice*dice-cp_wat*dwat;
  if(dcp_air)    (*dcp_air) = cp_air*dair-cp_wat*dwat;

  return;
}

void Density(AppCtx *user, PetscScalar ice, PetscScalar air, PetscScalar *rho, PetscScalar *drho_ice, PetscScalar *drho_air)
{
  PetscReal dice=1.0, dair=1.0, dwat=1.0;
  PetscReal wat = 1.0-air-ice;
  if(wat<0.0) {wat=0.0;dwat=0.0;}
  if(ice<0.0) {ice=0.0;dice=0.0;}
  if(air<0.0) {air=0.0;dair=0.0;}
  PetscReal rho_ice = user->rho_ice;
  PetscReal rho_wat = user->rho_wat;
  PetscReal rho_air = user->rho_air;
  if(rho)     (*rho)  = ice*rho_ice + wat*rho_wat + air*rho_air;
  if(drho_ice)    (*drho_ice) = rho_ice*dice-rho_wat*dwat;
  if(drho_air)    (*drho_air) = rho_air*dair-rho_wat*dwat;
  
  return;
}

void VaporDiffus(AppCtx *user, PetscScalar tem, PetscScalar *difvap, PetscScalar *d_difvap)
{

  PetscReal dif_vap = user->dif_vap;
  PetscReal Kratio = (tem+273.15)/273.15;
  PetscReal aa = 1.81;
  if(difvap)     (*difvap)  = dif_vap*pow(Kratio,aa);
  if(d_difvap)    (*d_difvap) = dif_vap*aa*pow(Kratio,aa-1.0)/273.15;
  
  return;
}

void RhoVS_I(AppCtx *user, PetscScalar tem, PetscScalar *rho_vs, PetscScalar *d_rhovs)
{

  PetscReal rho_air = user->rho_air;
  PetscReal K0,K1,K2,K3,K4,K5;
  K0 = -0.5865*1.0e4;   K1 = 0.2224*1.0e2;    K2 = 0.1375*1.0e-1;
  K3 = -0.3403*1.0e-4;  K4 = 0.2697*1.0e-7;   K5 = 0.6918;
  PetscReal Patm = 101325.0;
  PetscReal bb = 0.62;
  PetscReal temK = tem+273.15;
  PetscReal Pvs = exp(K0*pow(temK,-1.0)+K1+K2*pow(temK,1.0)+K3*pow(temK,2.0)+K4*pow(temK,3.0)+K5*log(temK));
  PetscReal Pvs_T = Pvs*(-K0*pow(temK,-2.0)+K2+2.0*K3*pow(temK,1.0)+3.0*K4*pow(temK,2.0)+K5/temK);

  if(rho_vs)     (*rho_vs)  = rho_air*bb*Pvs/(Patm-Pvs);
  if(d_rhovs)  (*d_rhovs) = rho_air*bb*(Pvs_T*(Patm-Pvs)+Pvs*Pvs_T)/(Patm-Pvs)/(Patm-Pvs);
  
  return;
}

void RhoVS_W(AppCtx *user, PetscScalar tem, PetscScalar *rho_vs, PetscScalar *d_rhovs)
{

  PetscReal rho_air = user->rho_air;
  PetscReal K0,K1,K2,K3,K4,K5,K6,K7;
  K0 = -0.2991*1.0e4;   K1 = -0.6017*1.0e4;    K2 = 0.18876*1.0e2;      K3 = -0.28355*1.0e-1;
  K4 = 0.1784*1.0e-4;   K5 = -0.8415*1.0e-9;   K6 = 0.444125*1.0e-12;   K7 = 0.28585*1.0e1;
  PetscReal Patm = 101325.0;
  PetscReal bb = 0.62;
  PetscReal temK = tem+273.15;
  PetscReal Pvs = exp(K0*pow(temK,-2.0)+K1*pow(temK,-1.0)+K2+K3*pow(temK,1.0)+K4*pow(temK,2.0)+K5*pow(temK,3.0)+K6*pow(temK,4.0)+K7*log(temK));
  PetscReal Pvs_T = Pvs*(-2.0*K0*pow(temK,-3.0)-K1*pow(temK,-2.0)+K3+2.0*K4*pow(temK,1.0)+3.0*K5*pow(temK,2.0)+4.0*K6*pow(temK,3.0)+K7/temK);

  if(rho_vs)     (*rho_vs)  = rho_air*bb*Pvs/(Patm-Pvs);
  if(d_rhovs)  (*d_rhovs) = rho_air*bb*(Pvs_T*(Patm-Pvs)+Pvs*Pvs_T)/(Patm-Pvs)/(Patm-Pvs);
  
  return;
}


void Density_air_PhCh(AppCtx *user, PetscScalar ice, PetscScalar air, PetscScalar *rhoairPC, PetscScalar *rhoairPC_ice, PetscScalar *rhoairPC_air)
{

  PetscReal phi_L = 1.0-user->phi_L;
  PetscReal rho_ice = user->rho_ice;
  PetscReal rho_wat = user->rho_wat;
  PetscReal wat = 1.0-air-ice;
  PetscReal x=0.0, dice=0.0, dair=0.0;

  if (air>=0.0 && ice>=0.0 && wat>=0.0 && air<=phi_L){ //1
    x = (1.0-air-2.0*ice)/(1.0-air);    dair = -2.0*ice/(1.0-air)/(1.0-air);    dice = -2.0/(1.0-air);
    if(rhoairPC)         (*rhoairPC) = 0.5*(rho_ice+rho_wat) + 0.5*(rho_wat-rho_ice)*0.5*x*(3.0-x*x);
    if(rhoairPC_ice)    (*rhoairPC_ice) = 0.5*(rho_wat-rho_ice)*0.5*3.0*(1.0-x*x)*dice;
    if(rhoairPC_air)    (*rhoairPC_air) = 0.5*(rho_wat-rho_ice)*0.5*3.0*(1.0-x*x)*dair;
  } else if (ice>=0.0 && wat>=0.0 && air>phi_L) { //2
    x = (1.0-air-2.0*ice)/(1.0-phi_L);    dair = -1.0/(1.0-phi_L);    dice = -2.0/(1.0-phi_L);
    if(rhoairPC)         (*rhoairPC) =  0.5*(rho_ice+rho_wat) + 0.5*(rho_wat-rho_ice)*0.5*x*(3.0-x*x);
    if(rhoairPC_ice)    (*rhoairPC_ice) = 0.5*(rho_wat-rho_ice)*0.5*3.0*(1.0-x*x)*dice;
    if(rhoairPC_air)    (*rhoairPC_air) = 0.5*(rho_wat-rho_ice)*0.5*3.0*(1.0-x*x)*dair;    
  } else if (air<0.0 && wat-ice <=1.0 && wat-ice >=-1) { //3
    x = (1.0-air-2.0*ice);    dair = -1.0;    dice = -2.0;
    if(rhoairPC)         (*rhoairPC) = 0.5*(rho_ice+rho_wat) + 0.5*(rho_wat-rho_ice)*0.5*x*(3.0-x*x);
    if(rhoairPC_ice)    (*rhoairPC_ice) = 0.5*(rho_wat-rho_ice)*0.5*3.0*(1.0-x*x)*dice;
    if(rhoairPC_air)    (*rhoairPC_air) = 0.5*(rho_wat-rho_ice)*0.5*3.0*(1.0-x*x)*dair;
  } else if (wat<0.0 && (air-ice)>=(2.0*phi_L-1.0) && (air-ice)<=1.0 ){ //6
    x = (air-ice-1.0)/(2.0*phi_L-2.0);    dice = -1.0/(2.0*phi_L-2.0);      dair = 1.0/(2.0*phi_L-2.0);
    if(rhoairPC)         (*rhoairPC) = 0.5*(rho_ice+rho_wat) + 0.5*(rho_wat-rho_ice)*0.5*x*(3.0-x*x);
    if(rhoairPC_ice)    (*rhoairPC_ice) = 0.5*(rho_wat-rho_ice)*0.5*3.0*(1.0-x*x)*dice;
    if(rhoairPC_air)    (*rhoairPC_air) = 0.5*(rho_wat-rho_ice)*0.5*3.0*(1.0-x*x)*dair;
  } else if (ice<0.0 && (air-wat)<=1.0 && (air-wat)>=(2.0*phi_L-1.0) ){ //7
    x = (2.0*air+ice-2.0)/(2.0*phi_L-2.0);    dice = 1.0/(2.0*phi_L-2.0);      dair = 2.0/(2.0*phi_L-2.0);
    if(rhoairPC)         (*rhoairPC) = 0.5*(rho_ice+rho_wat) - 0.5*(rho_wat-rho_ice)*0.5*x*(3.0-x*x);
    if(rhoairPC_ice)    (*rhoairPC_ice) = -0.5*(rho_wat-rho_ice)*0.5*3.0*(1.0-x*x)*dice;
    if(rhoairPC_air)    (*rhoairPC_air) = -0.5*(rho_wat-rho_ice)*0.5*3.0*(1.0-x*x)*dair;
  } else if ( (air-ice)>1.0 && (air-wat)>1.0 ){ //8
    if(rhoairPC)         (*rhoairPC) = 0.5*(rho_ice+rho_wat);
    if(rhoairPC_ice)    (*rhoairPC_ice) = 0.0;
    if(rhoairPC_air)    (*rhoairPC_air) = 0.0; 
  } else if ( ( wat<0.0 && (air-ice)<(2.0*phi_L-1.0) ) || ( (ice-air)>1.0 && (wat-ice) < -1.0 ) ) { //4
    if(rhoairPC)         (*rhoairPC) = rho_ice;
    if(rhoairPC_ice)    (*rhoairPC_ice) = 0.0;
    if(rhoairPC_air)    (*rhoairPC_air) = 0.0;     
  } else { //5
    if(rhoairPC)         (*rhoairPC) = rho_wat;
    if(rhoairPC_ice)    (*rhoairPC_ice) = 0.0;
    if(rhoairPC_air)    (*rhoairPC_air) = 0.0;
  }

  return;
}


void Fice(AppCtx *user, PetscScalar ice, PetscScalar air, PetscScalar *fice, PetscScalar *dfice_ice, PetscScalar *dfice_air)
{

  PetscReal Lambd = user->Lambd;
  PetscReal etai  = user->Etai;
  if(fice)     (*fice)  = etai*ice*(1.0-ice)*(1.0-2.0*ice) + 2.0*Lambd*ice*(1.0-air-ice)*(1.0-air-ice)*air*air;
  if(dfice_ice)    (*dfice_ice) = etai*(1.0-6.0*ice+6.0*ice*ice) + 2.0*Lambd*air*air*((1.0-air-ice)*(1.0-air-ice)-ice*2.0*(1.0-air-ice));
  if(dfice_air)    (*dfice_air) = 2.0*Lambd*ice*(2.0*air*(1.0-air-ice)*(1.0-air-ice)-air*air*2.0*(1.0-air-ice));
  
  return;
}

void Fwat(AppCtx *user, PetscScalar ice, PetscScalar air, PetscScalar *fwat, PetscScalar *dfwat_ice, PetscScalar *dfwat_air)
{

  PetscReal Lambd = user->Lambd;
  PetscReal etaw  = user->Etaw;
  if(fwat)     (*fwat)  = etaw*(1.0-air-ice)*(ice+air)*(2.0*ice+2.0*air-1.0) + 2.0*Lambd*ice*ice*(1.0-air-ice)*air*air;
  if(dfwat_ice)    {
    (*dfwat_ice)  = etaw*(2.0*(1.0-air-ice)*(ice+air) + (1.0-air-ice)*(2.0*ice+2.0*air-1.0) - (ice+air)*(2.0*ice+2.0*air-1.0));
    (*dfwat_ice) += 2.0*Lambd*air*air*(2.0*ice*(1.0-air-ice) - ice*ice);
  }
  if(dfwat_air)   {
    (*dfwat_air)  = etaw*(2.0*(1.0-air-ice)*(ice+air) + (1.0-air-ice)*(2.0*ice+2.0*air-1.0) - (ice+air)*(2.0*ice+2.0*air-1.0));
    (*dfwat_air) += 2.0*Lambd*ice*ice*(2.0*air*(1.0-air-ice) - air*air);
  }
  
  return;
}

void Fair(AppCtx *user, PetscScalar ice, PetscScalar air, PetscScalar *fair, PetscScalar *dfair_ice, PetscScalar *dfair_air)
{

  PetscReal Lambd = user->Lambd;
  PetscReal etaa  = user->Etaa;
  if(fair)     (*fair)  = etaa*air*(1.0-air)*(1.0-2.0*air) + 2.0*Lambd*ice*ice*(1.0-air-ice)*(1.0-air-ice)*air;
  if(dfair_ice)    (*dfair_ice) = 2.0*Lambd*air*(2.0*ice*(1.0-air-ice)*(1.0-air-ice) - ice*ice*2.0*(1.0-air-ice));
  if(dfair_air)    (*dfair_air) = etaa*(1.0-6.0*air+6.0*air*air) +2.0*Lambd*ice*ice*((1.0-air-ice)*(1.0-air-ice) - air*2.0*(1.0-air-ice));
  
  return;
}

void Nucl_funct(AppCtx *user, PetscScalar tem, PetscScalar *nucI, PetscScalar *nucW, PetscScalar *dnucI, PetscScalar *dnucW)
{

  PetscReal tem_nucl = user->tem_nucl;
  if(nucI)     (*nucI)  = 0.5-0.5*tanh(20.0*(tem-tem_nucl+0.1));
  if(nucW)     (*nucW)  = 0.5+0.5*tanh(20.0*(tem-0.1));
  if(dnucI)    (*dnucI) = -0.5*(1.0-tanh(20.0*(tem-tem_nucl+0.1))*tanh(20.0*(tem-tem_nucl+0.1)))*20.0;
  if(dnucW)    (*dnucW) = 0.5*(1.0-tanh(20.0*(tem-0.1))*tanh(20.0*(tem-0.1)))*20.0;
  
  return;
}


PetscErrorCode Residual(IGAPoint pnt,
                        PetscReal shift,const PetscScalar *V,
                        PetscReal t,const PetscScalar *U,
                        PetscScalar *Re,void *ctx)
{
  AppCtx *user = (AppCtx*) ctx;

  PetscReal eps = user->eps;
  PetscReal Etai = user->Etai;
  PetscReal Etaw = user->Etaw;
  PetscReal Etaa = user->Etaa;
  PetscReal ETA = Etaa*Etai + Etaa*Etaw + Etaw*Etai; 
  PetscReal alph_sol = user->alph_sol;
  PetscReal alph_sub = user->alph_sub;
  PetscReal alph_eva = user->alph_eva;
  PetscReal rho_ice = user->rho_ice;
  PetscReal rho_wat = user->rho_wat;
  PetscReal lat_sol = user->lat_sol;
  PetscReal lat_sub = user->lat_sub;
  PetscReal cp_wat = user->cp_wat;
  PetscReal T_melt = user->T_melt;
  PetscReal air_lim = user->air_lim;
  PetscReal nucleat = user->nucleat/eps;
  PetscReal xi_v = user->xi_v;
  PetscReal xi_T = user->xi_T;
  if(user->flag_xiT==0) xi_T=1.0;
  PetscInt pntind = pnt->parent->index*SQ(user->p+1)+pnt->index;
  PetscInt flag_mob = user->FlagMob[pntind];
  PetscReal costhet = user->costhet;

  if(pnt->atboundary){

    PetscScalar grad_sol[4][2];
    IGAPointFormGrad (pnt,U,&grad_sol[0][0]);

    PetscScalar grad_ice[2],modgradice;
    grad_ice[0]  = grad_sol[0][0];
    grad_ice[1]  = grad_sol[0][1];
    modgradice = sqrt(SQ(grad_ice[0])+SQ(grad_ice[1]));

    PetscScalar grad_air[2],modgradair;
    grad_air[0]  = grad_sol[1][0];
    grad_air[1]  = grad_sol[1][1];
    modgradair = sqrt(SQ(grad_air[0])+SQ(grad_air[1]));

    PetscReal mob;
    if (flag_mob==1) mob=user->mav;
    else if (flag_mob==2) mob=user->mob_eva;
    else if (flag_mob==3) mob=user->mob_sol;
    else if (flag_mob==4) mob=user->mob_sub;
    else {
      mob=user->mav;
      PetscPrintf(PETSC_COMM_SELF,"ERROR, NO MOBILITY DEFINED \n");
    }

    const PetscReal *N0; 
    IGAPointGetShapeFuns(pnt,0,(const PetscReal**)&N0);
    
    PetscScalar (*R)[4] = (PetscScalar (*)[4])Re;
    PetscInt a,nen=pnt->nen;
    for(a=0; a<nen; a++) {
      R[a][0] = -N0[a]*3.0*eps*mob*modgradice*costhet;
      R[a][1] = N0[a]*3.0*eps*mob*modgradair*costhet;
      R[a][2] = 0.0;
      R[a][3] = 0.0;
    }
      
  } else  {
    
    PetscScalar sol_t[4],sol[4];
    PetscScalar grad_sol[4][2];
    IGAPointFormValue(pnt,V,&sol_t[0]);
    IGAPointFormValue(pnt,U,&sol[0]);
    IGAPointFormGrad (pnt,U,&grad_sol[0][0]);

    PetscScalar ice, ice_t, grad_ice[2];
    ice          = sol[0]; 
    ice_t        = sol_t[0]; 
    grad_ice[0]  = grad_sol[0][0];
    grad_ice[1]  = grad_sol[0][1];

    PetscScalar air, air_t, grad_air[2];
    air          = sol[1]; 
    air_t        = sol_t[1]; 
    grad_air[0]  = grad_sol[1][0];
    grad_air[1]  = grad_sol[1][1];

    PetscScalar tem, tem_t, grad_tem[2];
    tem          = sol[2];
    tem_t        = sol_t[2];
    grad_tem[0]  = grad_sol[2][0];
    grad_tem[1]  = grad_sol[2][1]; 

    PetscScalar rhov, rhov_t, grad_rhov[2];
    rhov           = sol[3];
    rhov_t         = sol_t[3];
    grad_rhov[0]   = grad_sol[3][0];
    grad_rhov[1]   = grad_sol[3][1];

    PetscReal thcond,cp,rho,difvap,rhoI_vs,rhoW_vs,mob,rhoSE,fice,fwat,fair,nucI,nucW;
    ThermalCond(user,ice,air,&thcond,NULL,NULL);
    HeatCap(user,ice,air,&cp,NULL,NULL);
    Density(user,ice,air,&rho,NULL,NULL);
    VaporDiffus(user,tem,&difvap,NULL);
    RhoVS_I(user,tem,&rhoI_vs,NULL);
    RhoVS_W(user,tem,&rhoW_vs,NULL);
    Density_air_PhCh(user,ice,air,&rhoSE,NULL,NULL);
    Fice(user,ice,air,&fice,NULL,NULL);
    Fwat(user,ice,air,&fwat,NULL,NULL);
    Fair(user,ice,air,&fair,NULL,NULL);
    Nucl_funct(user,tem,&nucI,&nucW,NULL,NULL);

    if (flag_mob==1)  mob=user->mav; 
    else if (flag_mob==2) mob=user->mob_eva; 
    else if (flag_mob==3) mob=user->mob_sol; 
    else if (flag_mob==4) mob=user->mob_sub; 
    else {
      mob=user->mav; 
      PetscPrintf(PETSC_COMM_SELF,"ERROR, NO MOBILITY DEFINED \n");
    }

    const PetscReal *N0,(*N1)[2]; 
    IGAPointGetShapeFuns(pnt,0,(const PetscReal**)&N0);
    IGAPointGetShapeFuns(pnt,1,(const PetscReal**)&N1);
    
    PetscScalar (*R)[4] = (PetscScalar (*)[4])Re;
    PetscInt a,nen=pnt->nen;
    for(a=0; a<nen; a++) {
      
      PetscReal R_ice,R_air,R_tem,R_vap;

      if(user->flag_tIC==1){

        R_ice  = N0[a]*ice_t;
        R_air  = N0[a]*air_t;

        R_tem  =  rho*cp*N0[a]*tem_t;
        R_tem +=  thcond*(N1[a][0]*grad_tem[0] + N1[a][1]*grad_tem[1]);

        R_vap  =  N0[a]*rhov;
        R_vap -=  N0[a]*rhoI_vs;

      } else {

        R_ice  = N0[a]*ice_t; 
        R_ice += 3.0*mob*eps*(N1[a][0]*grad_ice[0] + N1[a][1]*grad_ice[1]);
        R_ice += N0[a]*mob*3.0/eps/ETA*((Etaw+Etaa)*fice - Etaa*fwat - Etaw*fair);
        R_ice += N0[a]*alph_sol*ice*ice*(1.0-ice-air)*(1.0-ice-air)*(tem-T_melt)/lat_sol*cp_wat;
        R_ice -= N0[a]*alph_sub*ice*ice*air*air*(rhov-rhoI_vs)/rho_ice;
        R_ice -= mob*nucleat*N0[a]*air*air*(1.0-air-ice)*(1.0-air-ice)*nucI;
        R_ice += mob*nucleat*N0[a]*air*air*ice*ice*nucW;

        R_air  = N0[a]*air_t;
        R_air += 3.0*mob*eps*(N1[a][0]*grad_air[0] + N1[a][1]*grad_air[1]);
        R_air += N0[a]*mob*3.0/eps/ETA*((Etaw+Etai)*fair - Etaw*fice - Etai*fwat);
        R_air += N0[a]*alph_sub*ice*ice*air*air*(rhov-rhoI_vs)/rho_ice;
        R_air += N0[a]*alph_eva*(1.0-ice-air)*(1.0-ice-air)*air*air*(rhov-rhoW_vs)/rho_wat;

        R_tem  = rho*cp*N0[a]*tem_t;
        R_tem += xi_T*thcond*(N1[a][0]*grad_tem[0] + N1[a][1]*grad_tem[1]);
        R_tem -= xi_T*rho*lat_sol*N0[a]*(air_t+ice_t);
        R_tem += xi_T*rho*lat_sub*N0[a]*air_t;
	//R_tem -= N0[a]*1.0e8;

        R_vap  = N0[a]*rhov*air_t;
        if(air>air_lim){
          R_vap += N0[a]*air*rhov_t;
          R_vap += xi_v*difvap*air*(N1[a][0]*grad_rhov[0] + N1[a][1]*grad_rhov[1]);
        } else {
          R_vap += N0[a]*air_lim*rhov_t;
          R_vap += xi_v*difvap*air_lim*(N1[a][0]*grad_rhov[0] + N1[a][1]*grad_rhov[1]);
        }
        R_vap -=  xi_v*N0[a]*rhoSE*air_t;        

      }

      R[a][0] = R_ice;
      R[a][1] = R_air;
      R[a][2] = R_tem;
      R[a][3] = R_vap;
    }
    
  }

  return 0;
}

PetscErrorCode Jacobian(IGAPoint pnt,
                        PetscReal shift,const PetscScalar *V,
                        PetscReal t,const PetscScalar *U,
                        PetscScalar *Je,void *ctx)
{
  AppCtx *user = (AppCtx*) ctx;

  PetscReal eps = user->eps;
  PetscReal Etai = user->Etai;
  PetscReal Etaw = user->Etaw;
  PetscReal Etaa = user->Etaa;
  PetscReal ETA = Etaa*Etai + Etaa*Etaw + Etaw*Etai; 
  PetscReal alph_sol = user->alph_sol;
  PetscReal alph_sub = user->alph_sub;
  PetscReal alph_eva = user->alph_eva;
  PetscReal rho_ice = user->rho_ice;
  PetscReal rho_wat = user->rho_wat;
  PetscReal lat_sol = user->lat_sol;
  PetscReal lat_sub = user->lat_sub;
  PetscReal cp_wat = user->cp_wat;
  PetscReal T_melt = user->T_melt;
  PetscReal air_lim = user->air_lim;
  PetscReal nucleat = user->nucleat/eps;
  PetscReal xi_v = user->xi_v;
  PetscReal xi_T = user->xi_T;
  if(user->flag_xiT==0) xi_T=1.0;
  PetscInt pntind = pnt->parent->index*SQ(user->p+1)+pnt->index;
  PetscInt flag_mob = user->FlagMob[pntind];
  PetscReal costhet = user->costhet;

 if(pnt->atboundary){

    PetscScalar grad_sol[4][2];
    IGAPointFormGrad (pnt,U,&grad_sol[0][0]);

    PetscScalar grad_ice[2], modgradice;
    grad_ice[0]  = grad_sol[0][0];
    grad_ice[1]  = grad_sol[0][1];
    modgradice = sqrt(SQ(grad_ice[0])+SQ(grad_ice[1]));
    if(modgradice<1.0e-5) modgradice=1.0e-5;

    PetscScalar grad_air[2], modgradair;
    grad_air[0]  = grad_sol[1][0];
    grad_air[1]  = grad_sol[1][1];
    modgradair = sqrt(SQ(grad_air[0])+SQ(grad_air[1]));
    if(modgradair<1.0e-5) modgradair=1.0e-5;

    PetscReal mob;
    if (flag_mob==1) mob=user->mav;
    else if (flag_mob==2) mob=user->mob_eva;
    else if (flag_mob==3) mob=user->mob_sol;
    else if (flag_mob==4) mob=user->mob_sub;
    else {
      mob=user->mav;
      PetscPrintf(PETSC_COMM_SELF,"ERROR, NO MOBILITY DEFINED \n");
    }

    const PetscReal *N0,(*N1)[2]; 
    IGAPointGetShapeFuns(pnt,0,(const PetscReal**)&N0);
    IGAPointGetShapeFuns(pnt,1,(const PetscReal**)&N1);

    PetscInt a,b,nen=pnt->nen;
    PetscScalar (*J)[4][nen][4] = (PetscScalar (*)[4][nen][4])Je;
    for(a=0; a<nen; a++) {
      for(b=0; b<nen; b++) {

        J[a][0][b][0] -= N0[a]*3.0*eps*mob*(grad_ice[0]*N1[b][0]+grad_ice[1]*N1[b][1])/modgradice*costhet;
        J[a][1][b][1] += N0[a]*3.0*eps*mob*(grad_air[0]*N1[b][0]+grad_air[1]*N1[b][1])/modgradair*costhet;

      }
    }
    
  } else {

    PetscScalar sol_t[4],sol[4];
    PetscScalar grad_sol[4][2];
    IGAPointFormValue(pnt,V,&sol_t[0]);
    IGAPointFormValue(pnt,U,&sol[0]);
    IGAPointFormGrad (pnt,U,&grad_sol[0][0]);

    PetscScalar ice, ice_t;
    ice          = sol[0]; 
    ice_t        = sol_t[0]; 

    PetscScalar air, air_t; 
    air          = sol[1]; 
    air_t        = sol_t[1]; 

    PetscScalar tem, tem_t, grad_tem[2];
    tem          = sol[2];
    tem_t        = sol_t[2];
    grad_tem[0]  = grad_sol[2][0];
    grad_tem[1]  = grad_sol[2][1]; 

    PetscScalar rhov, rhov_t, grad_rhov[2];
    rhov           = sol[3];
    rhov_t         = sol_t[3];
    grad_rhov[0]   = grad_sol[3][0];
    grad_rhov[1]   = grad_sol[3][1];

    PetscReal thcond,dthcond_ice,dthcond_air;
    ThermalCond(user,ice,air,&thcond,&dthcond_ice,&dthcond_air);
    PetscReal cp,dcp_ice,dcp_air;
    HeatCap(user,ice,air,&cp,&dcp_ice,&dcp_air);
    PetscReal rho,drho_ice,drho_air;
    Density(user,ice,air,&rho,&drho_ice,&drho_air);
    PetscReal difvap,d_difvap;
    VaporDiffus(user,tem,&difvap,&d_difvap);
    PetscReal rhoI_vs,drhoI_vs;
    RhoVS_I(user,tem,&rhoI_vs,&drhoI_vs);
    PetscReal rhoW_vs,drhoW_vs;
    RhoVS_W(user,tem,&rhoW_vs,&drhoW_vs);
    PetscReal rhoSE,rhoSE_air,rhoSE_ice;
    Density_air_PhCh(user,ice,air,&rhoSE,&rhoSE_ice,&rhoSE_air);
    PetscReal fice,fice_ice,fice_air;
    Fice(user,ice,air,&fice,&fice_ice,&fice_air);
    PetscReal fwat,fwat_ice,fwat_air;
    Fwat(user,ice,air,&fwat,&fwat_ice,&fwat_air);
    PetscReal fair,fair_ice,fair_air;
    Fair(user,ice,air,&fair,&fair_ice,&fair_air);
    PetscReal nucI,nucW,dnucI,dnucW;
    Nucl_funct(user,tem,&nucI,&nucW,&dnucI,&dnucW);

    PetscReal mob;
    if (flag_mob==1)  mob=user->mav; 
    else if (flag_mob==2) mob=user->mob_eva; 
    else if (flag_mob==3) mob=user->mob_sol; 
    else if (flag_mob==4) mob=user->mob_sub; 
    else {
      mob=user->mav; 
      PetscPrintf(PETSC_COMM_SELF,"ERROR, NO MOBILITY DEFINED \n");
    }

    const PetscReal *N0,(*N1)[2]; 
    IGAPointGetShapeFuns(pnt,0,(const PetscReal**)&N0);
    IGAPointGetShapeFuns(pnt,1,(const PetscReal**)&N1);
 
    PetscInt a,b,nen=pnt->nen;
    PetscScalar (*J)[4][nen][4] = (PetscScalar (*)[4][nen][4])Je;
    for(a=0; a<nen; a++) {
      for(b=0; b<nen; b++) {

        if(user->flag_tIC==1){

          J[a][0][b][0] += shift*N0[a]*N0[b];
          J[a][1][b][1] += shift*N0[a]*N0[b];

          J[a][2][b][2] += shift*rho*cp*N0[a]*N0[b];
          J[a][2][b][0] += drho_ice*N0[b]*cp*N0[a]*tem_t;
          J[a][2][b][1] += drho_air*N0[b]*cp*N0[a]*tem_t;
          J[a][2][b][0] += rho*dcp_ice*N0[b]*N0[a]*tem_t;
          J[a][2][b][1] += rho*dcp_air*N0[b]*N0[a]*tem_t;
          J[a][3][b][0] += dthcond_ice*N0[b]*(N1[a][0]*grad_tem[0]+N1[a][1]*grad_tem[1]);
          J[a][2][b][1] += dthcond_air*N0[b]*(N1[a][0]*grad_tem[0]+N1[a][1]*grad_tem[1]);
          J[a][2][b][2] += thcond*(N1[a][0]*N1[b][0]+N1[a][1]*N1[b][1]);

          J[a][3][b][3] += N0[a]*N0[b];
          J[a][3][b][2] -= N0[a]*drhoI_vs*N0[b];

        } else {
        
        //ice
          J[a][0][b][0] += shift*N0[a]*N0[b];
          J[a][0][b][0] += 3.0*mob*eps*(N1[a][0]*N1[b][0] + N1[a][1]*N1[b][1]);

          J[a][0][b][0] += N0[a]*mob*3.0/eps/ETA*((Etaw+Etaa)*fice_ice - Etaa*fwat_ice - Etaw*fair_ice)*N0[b];
          J[a][0][b][1] += N0[a]*mob*3.0/eps/ETA*((Etaw+Etaa)*fice_air - Etaa*fwat_air - Etaw*fair_air)*N0[b];
          J[a][0][b][0] += N0[a]*alph_sol*2.0*ice*N0[b]*(1.0-ice-air)*(1.0-ice-air)*(tem-T_melt)/lat_sol*cp_wat;
          J[a][0][b][0] -= N0[a]*alph_sol*ice*ice*2.0*(1.0-ice-air)*N0[b]*(tem-T_melt)/lat_sol*cp_wat;
          J[a][0][b][1] -= N0[a]*alph_sol*ice*ice*2.0*(1.0-ice-air)*N0[b]*(tem-T_melt)/lat_sol*cp_wat;
          J[a][0][b][2] += N0[a]*alph_sol*ice*ice*(1.0-ice-air)*(1.0-ice-air)*N0[b]/lat_sol*cp_wat;
          J[a][0][b][0] -= N0[a]*alph_sub*2.0*ice*N0[b]*air*air*(rhov-rhoI_vs)/rho_ice;
          J[a][0][b][1] -= N0[a]*alph_sub*ice*ice*2.0*air*N0[b]*(rhov-rhoI_vs)/rho_ice;
          J[a][0][b][2] += N0[a]*alph_sub*ice*ice*air*air*drhoI_vs*N0[b]/rho_ice;
          J[a][0][b][3] -= N0[a]*alph_sub*ice*ice*air*air*N0[b]/rho_ice;
          J[a][0][b][0] += mob*nucleat*N0[a]*air*air*2.0*(1.0-air-ice)*N0[b]*nucI;
          J[a][0][b][1] -= mob*nucleat*N0[a]*2.0*air*N0[b]*(1.0-air-ice)*(1.0-air-ice)*nucI;
          J[a][0][b][1] += mob*nucleat*N0[a]*air*air*2.0*(1.0-air-ice)*N0[b]*nucI;
          J[a][0][b][2] -= mob*nucleat*N0[a]*air*air*(1.0-air-ice)*(1.0-air-ice)*dnucI*N0[b];
          J[a][0][b][0] += mob*nucleat*N0[a]*air*air*2.0*ice*N0[b]*nucW;
          J[a][0][b][1] += mob*nucleat*N0[a]*2.0*air*N0[b]*ice*ice*nucW;
          J[a][0][b][2] += mob*nucleat*N0[a]*air*air*ice*ice*dnucW*N0[b];

        //air
          J[a][1][b][1] += shift*N0[a]*N0[b];
          J[a][1][b][1] += 3.0*mob*eps*(N1[a][0]*N1[b][0] + N1[a][1]*N1[b][1]);

          J[a][1][b][0] += N0[a]*mob*3.0/eps/ETA*((Etaw+Etai)*fair_ice - Etaw*fice_ice - Etai*fwat_ice)*N0[b];
          J[a][1][b][1] += N0[a]*mob*3.0/eps/ETA*((Etaw+Etai)*fair_air - Etaw*fice_air - Etai*fwat_air)*N0[b];
          J[a][1][b][0] += N0[a]*alph_sub*2.0*ice*N0[b]*air*air*(rhov-rhoI_vs)/rho_ice;
          J[a][1][b][1] += N0[a]*alph_sub*ice*ice*2.0*air*N0[b]*(rhov-rhoI_vs)/rho_ice;
          J[a][1][b][2] -= N0[a]*alph_sub*ice*ice*air*air*drhoI_vs*N0[b]/rho_ice;
          J[a][1][b][3] += N0[a]*alph_sub*ice*ice*air*air*N0[b]/rho_ice;
          J[a][1][b][0] -= N0[a]*alph_eva*2.0*(1.0-ice-air)*N0[b]*air*air*(rhov-rhoW_vs)/rho_wat;
          J[a][1][b][1] -= N0[a]*alph_eva*2.0*(1.0-ice-air)*N0[b]*air*air*(rhov-rhoW_vs)/rho_wat;
          J[a][1][b][1] += N0[a]*alph_eva*(1.0-ice-air)*(1.0-ice-air)*2.0*air*N0[b]*(rhov-rhoW_vs)/rho_wat;
          J[a][1][b][2] -= N0[a]*alph_eva*(1.0-ice-air)*(1.0-ice-air)*air*air*drhoW_vs*N0[b]/rho_wat;
          J[a][1][b][3] += N0[a]*alph_eva*(1.0-ice-air)*(1.0-ice-air)*air*air*N0[b]/rho_wat;

        //temperature
          J[a][2][b][2] += shift*rho*cp*N0[a]*N0[b];
          J[a][2][b][0] += drho_ice*N0[b]*cp*N0[a]*tem_t;
          J[a][2][b][1] += drho_air*N0[b]*cp*N0[a]*tem_t;
          J[a][2][b][0] += rho*dcp_ice*N0[b]*N0[a]*tem_t;
          J[a][2][b][1] += rho*dcp_air*N0[b]*N0[a]*tem_t;
          J[a][2][b][0] += xi_T*dthcond_ice*N0[b]*(N1[a][0]*grad_tem[0] + N1[a][1]*grad_tem[1]);
          J[a][2][b][1] += xi_T*dthcond_air*N0[b]*(N1[a][0]*grad_tem[0] + N1[a][1]*grad_tem[1]);
          J[a][2][b][2] += xi_T*thcond*(N1[a][0]*N1[b][0] + N1[a][1]*N1[b][1]);
          J[a][2][b][0] -= xi_T*drho_ice*N0[b]*lat_sol*N0[a]*(air_t+ice_t);
          J[a][2][b][1] -= xi_T*drho_air*N0[b]*lat_sol*N0[a]*(air_t+ice_t);
          J[a][2][b][0] -= xi_T*rho*lat_sol*N0[a]*shift*N0[b];
          J[a][2][b][1] -= xi_T*rho*lat_sol*N0[a]*shift*N0[b];
          J[a][2][b][0] += xi_T*drho_ice*N0[b]*lat_sub*N0[a]*air_t;
          J[a][2][b][1] += xi_T*drho_air*N0[b]*lat_sub*N0[a]*air_t;
          J[a][2][b][1] += xi_T*rho*lat_sub*N0[a]*shift*N0[b];

        //vapor density
          J[a][3][b][1] += N0[a]*rhov*shift*N0[b];
          J[a][3][b][3] += N0[a]*N0[b]*air_t;
          if(air>air_lim){
            J[a][3][b][1] += N0[a]*N0[b]*rhov_t;
            J[a][3][b][3] += N0[a]*air*shift*N0[b];
            J[a][3][b][1] += xi_v*difvap*N0[b]*(N1[a][0]*grad_rhov[0] + N1[a][1]*grad_rhov[1]);
            J[a][3][b][2] += xi_v*d_difvap*N0[b]*air*(N1[a][0]*grad_rhov[0] + N1[a][1]*grad_rhov[1]);
            J[a][3][b][3] += xi_v*difvap*air*(N1[a][0]*N1[b][0] + N1[a][1]*N1[b][1]);        
          } else {
            J[a][3][b][3] += N0[a]*air_lim*shift*N0[b];
            J[a][3][b][2] += xi_v*d_difvap*N0[b]*air_lim*(N1[a][0]*grad_rhov[0] + N1[a][1]*grad_rhov[1]);
            J[a][3][b][3] += xi_v*difvap*air_lim*(N1[a][0]*N1[b][0] + N1[a][1]*N1[b][1]);
          }
          J[a][3][b][0] -= xi_v*N0[a]*rhoSE_ice*N0[b]*air_t;
          J[a][3][b][1] -= xi_v*N0[a]*rhoSE_air*N0[b]*air_t;
          J[a][3][b][1] -= xi_v*N0[a]*rhoSE*shift*N0[b];
        }

      }
    }

  //return 0;
  }
  return 0;
}


PetscErrorCode Integration(IGAPoint pnt,const PetscScalar *U,PetscInt n,PetscScalar *S,void *ctx)
{
  PetscFunctionBegin;
  AppCtx *user = (AppCtx *)ctx;

  PetscScalar sol[4],grad_sol[4][2];
  IGAPointFormValue(pnt,U,&sol[0]);
  IGAPointFormGrad(pnt,U,&grad_sol[0][0]);

  PetscReal grad_wat[2];
  grad_wat[0] = -grad_sol[0][0]-grad_sol[1][0];
  grad_wat[1] = -grad_sol[0][1]-grad_sol[1][1];
  PetscReal ice     = sol[0]; 
  PetscReal air     = sol[1];
  PetscReal wat     = 1.0-ice-air;
  PetscReal sol_interf = SQ(wat)*SQ(ice);
  PetscReal ice_SSA = user->eps*(SQ(grad_sol[0][0])+SQ(grad_sol[0][1]));
  PetscReal W_SSA = user->eps*fabs(grad_sol[0][0]*grad_wat[0]+grad_sol[0][1]*grad_wat[1]);

  S[0]  = ice;
  S[1]  = wat;
  S[2]  = air;
  S[3]  = sol_interf;
  S[4]  = ice_SSA;
  S[5]  = W_SSA;

  PetscFunctionReturn(0);
}

PetscErrorCode Monitor(TS ts,PetscInt step,PetscReal t,Vec U,void *mctx)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  AppCtx *user = (AppCtx *)mctx;

//------Motility region
  PetscInt ii,act;
  if(step % 1 == 0) {
    Vec localU;
    const PetscScalar *arrayU;
    IGAElement element;
    IGAPoint point;
    PetscScalar *UU;
    act=0; ii=0;
    PetscReal air,ice,wat,lim=0.01;

    ierr = IGAGetLocalVecArray(user->iga,U,&localU,&arrayU);CHKERRQ(ierr);
    ierr = IGABeginElement(user->iga,&element);CHKERRQ(ierr);
    while (IGANextElement(user->iga,element)) {
      act=0;
      ierr = IGAElementGetValues(element,arrayU,&UU);CHKERRQ(ierr);
      ierr = IGAElementBeginPoint(element,&point);CHKERRQ(ierr);
      while (IGAElementNextPoint(element,point)) {
        PetscScalar sol[4];
        ierr = IGAPointFormValue(point,UU,&sol[0]);CHKERRQ(ierr);
        ice = sol[0];
        air = sol[1];
        wat = 1.0-air-ice;
        if(air>lim && ice>lim && wat>lim) {act=1; }
        else if(ice<=lim && wat>=ice && air>ice) {act=2; }
        else if(air<=lim && ice>=air && wat>=air) {act=3; }
        else if(wat<=lim && air>wat && ice>wat) {act=4; }
        else {act=1; }
        user->FlagMob[ii] = act;

        ii++;
      }
      ierr = IGAElementEndPoint(element,&point);CHKERRQ(ierr);
    }
    ierr = IGAEndElement(user->iga,&element);CHKERRQ(ierr);
    ierr = IGARestoreLocalVecArray(user->iga,U,&localU,&arrayU);CHKERRQ(ierr);
  }

//--------Integration on the domain
  PetscScalar stats[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
  ierr = IGAComputeScalar(user->iga,U,6,&stats[0],Integration,mctx);CHKERRQ(ierr);
  PetscReal tot_ice     = PetscRealPart(stats[0]);
  PetscReal tot_wat     = PetscRealPart(stats[1]);
  PetscReal tot_air     = PetscRealPart(stats[2]);
  PetscReal sol_interf  = PetscRealPart(stats[3]);
  PetscReal ice_SSA     = PetscRealPart(stats[4]);
  PetscReal W_SSA       = PetscRealPart(stats[5]); 

  if(step==0) PetscPrintf(PETSC_COMM_WORLD,"\nPOROSITY %.3f\n\n",1.0 - tot_ice/user->Lx/user->Ly);

  if(step==0) {user->xiT_count=10;}

//-------- temperature time-scaling
  PetscReal  solinterf_value;       // for 2D
  PetscInt temp_xiT=user->flag_xiT; 

  if(user->xiT_count>2){
    solinterf_value=3.0*user->eps*SQ(0.5)*SQ(0.5)*0.01*sqrt(user->Lx*user->Ly);        //1.87e-13  // (1% of domain length) times (1D-value)
    if(sol_interf<0.001*solinterf_value) user->flag_xiT = 1; //1.87e-16
    else user->flag_xiT=0;
    user->xiT_count = 0;
  }
  user->xiT_count ++;
  if(temp_xiT != user->flag_xiT) PetscPrintf(PETSC_COMM_WORLD,"Temp time scale changes from (%d) to (%d) -- (0):xi_T=1.0   (1):xi_T=0.01  \n",temp_xiT,user->flag_xiT);


//-----------
  PetscReal dt;
  TSGetTimeStep(ts,&dt);
  if(step==1) user->flag_it0 = 0;

//--------- Boundary Conditions
  if(user->BC_Tvar==1){
    //if(t<0.5288) user->iga->Temp_metam = user->temp0;
    //else user->iga->Temp_metam = -0.01;
    //user->iga->Temp_metam = 0.1*sin(2.0*3.14159*(t+0.6666*dt)/2.0);
    user->iga->Temp_metam = user->temp0;
  }

//--------- Initial condition
  if(user->flag_tIC==1) if(step==user->nsteps_IC) {
    user->flag_tIC = 0; user->t_IC = t; 
    PetscPrintf(PETSC_COMM_WORLD,"INITIAL_CONDITION!!! \n");
  }
  
  if(step%10==0) PetscPrintf(PETSC_COMM_WORLD,"\nTIME(s)       TIME_STEP     TOT_ICE      TOT_WAT       TOT_AIR      SOL_INTERF ICE_SSA    W_SSA \n");
              PetscPrintf(PETSC_COMM_WORLD,"\n(%.0f) %.3e    %.3e   %.3e   %.3e   %.3e   %.3e   %.3e   %.3e \n\n",
                t,t,dt,tot_ice,tot_wat,tot_air,sol_interf,ice_SSA,W_SSA);


  if(step % 2 == 0) {
    const char *env = "folder"; char *dir; dir = getenv(env);
    char  filedata[256];
    sprintf(filedata,"%s/Data.dat",dir);
    PetscViewer       view;
    PetscViewerCreate(PETSC_COMM_WORLD,&view);
    PetscViewerSetType(view,PETSCVIEWERASCII);
    if (step==0) PetscViewerFileSetMode(view,FILE_MODE_WRITE); else PetscViewerFileSetMode(view,FILE_MODE_APPEND);
    PetscViewerFileSetName(view,filedata);
    PetscViewerASCIIPrintf(view," %d %e %e %e %e %e %e \n",step,t,dt,tot_ice,tot_wat,ice_SSA,W_SSA);
    PetscViewerDestroy(&view);
  }


  PetscFunctionReturn(0);
}

PetscErrorCode OutputMonitor(TS ts, PetscInt step, PetscReal t, Vec U,void *mctx)
{
  PetscFunctionBegin;
  PetscErrorCode ierr;
  AppCtx *user = (AppCtx *)mctx;

  if(step==0) {
   
    const char *env = "folder"; char *dir; dir = getenv(env);
    char  fileiga[256],  filename[256];
    PetscPrintf(PETSC_COMM_WORLD,"folder %s \n",dir);
    sprintf(fileiga,"%s/igasol.dat",dir);
    sprintf(filename,"%s/sol%d.dat",dir,step); 
    ierr = IGAWrite(user->iga,fileiga);CHKERRQ(ierr);
    ierr = IGAWriteVec(user->iga,U,filename);CHKERRQ(ierr);
  }

  PetscInt print=0;
  if(user->outp > 0) {
    if(step % user->outp == 0) print=1;
  } else {
    if (t>= user->t_out) print=1;
  }

  if(print == 1) {
    PetscPrintf(PETSC_COMM_WORLD,"OUTPUT print!\n");

    user->t_out += user->t_interv;

    const char *env = "folder"; char *dir; dir = getenv(env);
    char   filename[256];
    sprintf(filename,"%s/sol%d.dat",dir,step);
    ierr = IGAWriteVec(user->iga,U,filename);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode InitialIceGrains(IGA iga,AppCtx *user)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;

  PetscPrintf(PETSC_COMM_WORLD,"--------------------- ICE GRAINS --------------------------\n");

  if(user->NCice==0) {
    user->n_act = 0;
    PetscPrintf(PETSC_COMM_WORLD,"No ice grains\n\n");
    PetscFunctionReturn(0);
  }

  PetscReal rad = user->RCice;
  PetscReal rad_dev = user->RCice_dev;
  PetscInt  numb_clust = user->NCice, ii,jj,tot=10000;
  PetscInt  l, dim=2, n_act=0,flag,seed=user->seed;

  PetscReal centX[dim][numb_clust], radius[numb_clust];
  PetscRandom randcX,randcY,randcR;
  ierr = PetscRandomCreate(PETSC_COMM_WORLD,&randcX);CHKERRQ(ierr);
  ierr = PetscRandomCreate(PETSC_COMM_WORLD,&randcY);CHKERRQ(ierr);
  ierr = PetscRandomCreate(PETSC_COMM_WORLD,&randcR);CHKERRQ(ierr);
  ierr = PetscRandomSetInterval(randcX,0.0,user->Lx);CHKERRQ(ierr);
  ierr = PetscRandomSetInterval(randcY,0.0,user->Ly);CHKERRQ(ierr);
  ierr = PetscRandomSetInterval(randcR,sqrt(rad*(1.0-rad_dev)),sqrt(rad*(1.0+rad_dev)));CHKERRQ(ierr);
  ierr = PetscRandomSetSeed(randcX,seed+24+9*iga->elem_start[0]+11*iga->elem_start[1]);CHKERRQ(ierr);
  ierr = PetscRandomSetSeed(randcY,seed+numb_clust*35+5*iga->elem_start[1]+3*iga->elem_start[0]);CHKERRQ(ierr);
  ierr = PetscRandomSetSeed(randcR,seed*numb_clust+6*iga->proc_ranks[1]+5*iga->elem_start[0]+9);CHKERRQ(ierr);
  ierr = PetscRandomSeed(randcX);CHKERRQ(ierr);
  ierr = PetscRandomSeed(randcY);CHKERRQ(ierr);
  ierr = PetscRandomSeed(randcR);CHKERRQ(ierr);
  ierr = PetscRandomSetFromOptions(randcX);CHKERRQ(ierr);
  ierr = PetscRandomSetFromOptions(randcY);CHKERRQ(ierr);
  ierr = PetscRandomSetFromOptions(randcR);CHKERRQ(ierr);

  PetscReal xc[dim], rc=0.0, dist=0.0;
  xc[0] = xc[1] =  0.0;

  for(ii=0;ii<tot*numb_clust;ii++)
  {
    ierr=PetscRandomGetValue(randcX,&xc[0]);CHKERRQ(ierr);
    ierr=PetscRandomGetValue(randcY,&xc[1]);CHKERRQ(ierr);
    ierr=PetscRandomGetValue(randcR,&rc);CHKERRQ(ierr);
    rc = SQ(rc);
    flag=1;
    
    if(flag==1){
      for(jj=0;jj<n_act;jj++){
        dist = 0.0;
        for(l=0;l<dim;l++) dist += SQ(xc[l]-centX[l][jj]);
        dist = sqrt(dist);
        if(dist< user->overl*(rc+radius[jj]) ) flag = 0;
        if(dist> 1.02*(rc+radius[jj]) && dist< 1.15*(rc+radius[jj]) ) flag = 0;
      }
    }
    if(flag==1){
      if(dim==3) PetscPrintf(PETSC_COMM_WORLD," new ice grain %d!!  x %.2e  y %.2e  z %.2e  r %.2e \n",n_act,xc[0],xc[1],xc[2],rc);
      else PetscPrintf(PETSC_COMM_WORLD," new ice grain %d!!  x %.2e  y %.2e  r %.2e \n",n_act,xc[0],xc[1],rc);
      for(l=0;l<dim;l++) centX[l][n_act] = xc[l];
      radius[n_act] = rc;
      n_act++;
    }
    if(n_act==numb_clust) {
      PetscPrintf(PETSC_COMM_WORLD," %d ice grains in %d iterations \n\n", n_act,ii+1);
      ii=tot*numb_clust;
    }
  }
  if(n_act != numb_clust) PetscPrintf(PETSC_COMM_WORLD," %d ice grains in maximum number of iterations allowed (%d) \n\n", n_act, ii);

  ierr = PetscRandomDestroy(&randcX);CHKERRQ(ierr);
  ierr = PetscRandomDestroy(&randcY);CHKERRQ(ierr);
  ierr = PetscRandomDestroy(&randcR);CHKERRQ(ierr);

  for(l=0;l<dim;l++) {ierr = MPI_Bcast(centX[l],numb_clust,MPI_DOUBLE,0,PETSC_COMM_WORLD);CHKERRQ(ierr);}
  ierr = MPI_Bcast(radius,numb_clust,MPI_DOUBLE,0,PETSC_COMM_WORLD);CHKERRQ(ierr);
  ierr = MPI_Bcast(&n_act,1,MPI_INT,0,PETSC_COMM_WORLD);CHKERRQ(ierr);

  user->n_act = n_act;
  for(jj=0;jj<n_act;jj++){
    for(l=0;l<dim;l++) user->cent[l][jj] = centX[l][jj];
    user->radius[jj] = radius[jj];
  }

  PetscFunctionReturn(0); 
}


typedef struct {
  PetscScalar ice,air,tem,rhov;
} Field;

PetscErrorCode FormInitialCondition(IGA iga,PetscReal t,Vec U,AppCtx *user,const char datafile[],const char dataPF[])
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

  } else if (dataPF[0]){
    IGA igaPF;
    ierr = IGACreate(PETSC_COMM_WORLD,&igaPF);CHKERRQ(ierr);
    ierr = IGASetDim(igaPF,2);CHKERRQ(ierr);
    ierr = IGASetDof(igaPF,1);CHKERRQ(ierr);
    IGAAxis axisPF0;
    ierr = IGAGetAxis(igaPF,0,&axisPF0);CHKERRQ(ierr);
    ierr = IGAAxisSetDegree(axisPF0,user->p);CHKERRQ(ierr);
    ierr = IGAAxisInitUniform(axisPF0,user->Nx,0.0,user->Lx,user->C);CHKERRQ(ierr);
    IGAAxis axisPF1;
    ierr = IGAGetAxis(igaPF,1,&axisPF1);CHKERRQ(ierr);
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
    ierr = IGACreateNodeDM(iga,4,&da);CHKERRQ(ierr);
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

        PetscReal ice0 = u[j][i].ice;
        u[j][i].air = 1.0-ice0;
        u[j][i].tem = user->temp0 + user->grad_temp0[0]*(x-0.5*user->Lx) + user->grad_temp0[1]*(y-0.5*user->Ly);
        PetscScalar rho_vs, temp=u[j][i].tem;
        RhoVS_I(user,temp,&rho_vs,NULL);
        u[j][i].rhov = rho_vs;
      }
    }
    ierr = DMDAVecRestoreArray(da,U,&u);CHKERRQ(ierr); 
    ierr = DMDestroy(&da);;CHKERRQ(ierr);

  } else {

    DM da;
    ierr = IGACreateNodeDM(iga,4,&da);CHKERRQ(ierr);
    Field **u;
    ierr = DMDAVecGetArray(da,U,&u);CHKERRQ(ierr);
    DMDALocalInfo info;
    ierr = DMDAGetLocalInfo(da,&info);CHKERRQ(ierr);

    PetscInt i,j,m,k=-1;
    if(user->periodic==1) k=user->p -1;
    for(i=info.xs;i<info.xs+info.xm;i++){
      for(j=info.ys;j<info.ys+info.ym;j++){
        PetscReal x = user->Lx*(PetscReal)i / ( (PetscReal)(info.mx+k) );
        PetscReal y = user->Ly*(PetscReal)j / ( (PetscReal)(info.my+k) );

        //PetscReal bot_hor = 0.5 - 0.5*tanh(0.5/user->eps*(y-0.5*user->Ly));
        //PetscReal top_hor = 0.5 + 0.5*tanh(0.5/user->eps*(y-0.5*user->Ly));
        //PetscReal lef_ver = 0.5 - 0.5*tanh(0.5/user->eps*(x-0.5*user->Lx));
        //PetscReal rig_ver = 0.5 - 0.5*tanh(0.5/user->eps*(x-0.5*user->Lx));

 /*     	PetscReal xc[4],yc[4],Rc[4],arg, ice=0.0, wat=0.0;
      	xc[0]=6.5e-5 ; yc[0]=7.0e-5 ; Rc[0]=4.2e-5 ;
      	xc[1]=1.35e-4 ; yc[1]=6.5e-5 ; Rc[1]=2.7e-5 ;
      	xc[2]=1.2e-4 ; yc[2]=1.23e-4 ; Rc[2]=3.2e-5 ;
      	xc[3]=5.0e-5 ; yc[3]=1.34e-4 ; Rc[3]=2.3e-5 ;
      	for(m=0;m<4;m++){
      	  arg = sqrt(SQ(x-xc[m])+SQ(y-yc[m]))-(Rc[m]- 0.3e-5 );
      	  ice += 0.5-0.5*tanh(0.5/user->eps*arg);
	  arg = sqrt(SQ(x-xc[m])+SQ(y-yc[m]))-Rc[m];
	  wat += 0.5-0.5*tanh(0.5/user->eps*arg);
      	}
*/
	PetscReal dist, small=0.0, big=0.0, alp_i=0.0, alp_e=0.16;
	for(m=0;m<user->n_act;m++){
	   dist=sqrt(SQ(x-user->cent[0][m])+SQ(y-user->cent[1][m]));
           small += 0.5-0.5*tanh(0.5/user->eps*(dist-(user->radius[m] - alp_i*user->RCice)));
	   big += 0.5-0.5*tanh(0.5/user->eps*(dist-(user->radius[m] + alp_e*user->RCice)));
	}

      	if(small>1.0) small=1.0;
      	if(small<0.0) small=0.0;
	if(big>1.0) big=1.0;
	if(big<0.0) big=0.0;

        u[j][i].ice = small;     //bot_hor;
        u[j][i].air = 1.0-big; //top_hor*lef_ver;
        u[j][i].tem = user->temp0 + user->grad_temp0[0]*(x-0.5*user->Lx) + user->grad_temp0[1]*(y-0.5*user->Ly);
        PetscScalar rho_vs, temp=u[j][i].tem;
        RhoVS_I(user,temp,&rho_vs,NULL);
        u[j][i].rhov = rho_vs;
      }
    }
    ierr = DMDAVecRestoreArray(da,U,&u);CHKERRQ(ierr); 
    ierr = DMDestroy(&da);;CHKERRQ(ierr); 
  }
  PetscFunctionReturn(0); 
}

int main(int argc, char *argv[]) {

  // Petsc Initialization rite of passage 
  PetscErrorCode  ierr;
  ierr = PetscInitialize(&argc,&argv,0,0);CHKERRQ(ierr);

  PetscLogDouble itim;
  ierr = PetscTime(&itim); CHKERRQ(ierr);

  AppCtx user;
 
  PetscReal d0_sol = 4.0e-9, d0_sub0 = 1.0e-7, d0_eva0 = 6.0e-8; //d0   sol 4e-10  sub 9.5e-10  eva 6e-10
    //note: if minimum ice radius > 0.05mm, d0_sub d0_eva can be increased 2 orders of magnitude (*100) to speed numerics
  PetscReal beta_sol = 125.0, beta_sub0 = 1.4e5, beta_eva0 = 1.53e4;   //bet  sol 12     sub 1.4e5    eva 1.53e4
  PetscReal a1=5.0, a2=0.1581, rho_rhovs = 2.0e5;   // at 0C,  rho_rhovs=5e5 at -10C
  PetscReal gamma_iw = 0.033, gamma_iv = 0.109, gamma_wv = 0.076;

  PetscReal d0_sub, d0_eva, beta_sub, beta_eva;
  d0_sub = d0_sub0/rho_rhovs; d0_eva = d0_eva0/rho_rhovs; beta_sub = beta_sub0/rho_rhovs; beta_eva = beta_eva0/rho_rhovs;

  // PetscInt angle = 0; // 0:real,  1:120deg,
  // PetscInt mmm   = 0; // 1 if constant mobility
  // PetscInt aaa   = 0; // 1 if no alpha (sol,sub,eva)

// Import all environmnetal variables
  PetscPrintf(PETSC_COMM_WORLD, "Unpacking environment variables...\n");

  const char *angle_str       = getenv("angl");
  const char *mmm_str         = getenv("mm");
  const char *aaa_str         = getenv("aa");

  const char *Nx_str          = getenv("Nx");
  const char *Ny_str          = getenv("Ny");
  const char *Nz_str          = getenv("Nz");

  const char *Lx_str          = getenv("Lx");
  const char *Ly_str          = getenv("Ly");
  const char *Lz_str          = getenv("Lz");

  const char *delt_t_str      = getenv("delt_t");
  const char *t_final_str     = getenv("t_final");
  const char *n_out_str       = getenv("n_out");

  const char *temp_str        = getenv("temp");

  const char *grad_temp0X_str = getenv("grad_temp0X");
  const char *grad_temp0Y_str = getenv("grad_temp0Y");
  const char *grad_temp0Z_str = getenv("grad_temp0Z");

  const char *dim_str         = getenv("dim");
	
	const char *eps_str 				= getenv("eps");
	
  // Verify that all environment variables are set
  if (!angle_str) {
    PetscPrintf(PETSC_COMM_WORLD, "Error: angl environment variable is not set.\n");
    PetscFinalize();
    return EXIT_FAILURE;
  } else if (!mmm_str) {
    PetscPrintf(PETSC_COMM_WORLD, "Error: mm environment variable is not set.\n");
    PetscFinalize();
    return EXIT_FAILURE;
  } else if (!aaa_str) {
    PetscPrintf(PETSC_COMM_WORLD, "Error: aa environment variable is not set.\n");
    PetscFinalize();
    return EXIT_FAILURE;
  } else if (!Nx_str) {
      PetscPrintf(PETSC_COMM_WORLD, "Error: Nx_str environment variable is not set.\n");
      PetscFinalize();
      return EXIT_FAILURE;
  } else if (!Ny_str) {
      PetscPrintf(PETSC_COMM_WORLD, "Error: Ny_str environment variable is not set.\n");
      PetscFinalize();
      return EXIT_FAILURE;
  } else if (!Nz_str) {
      PetscPrintf(PETSC_COMM_WORLD, "Error: Nz_str environment variable is not set.\n");
      PetscFinalize();
      return EXIT_FAILURE;
  } else if (!Lx_str) {
      PetscPrintf(PETSC_COMM_WORLD, "Error: Lx_str environment variable is not set.\n");
      PetscFinalize();
      return EXIT_FAILURE;
  } else if (!Ly_str) {
      PetscPrintf(PETSC_COMM_WORLD, "Error: Ly_str environment variable is not set.\n");
      PetscFinalize();
      return EXIT_FAILURE;
  } else if (!Lz_str) {
      PetscPrintf(PETSC_COMM_WORLD, "Error: Lz_str environment variable is not set.\n");
      PetscFinalize();
      return EXIT_FAILURE;
  } else if (!delt_t_str) {
      PetscPrintf(PETSC_COMM_WORLD, "Error: delt_t_str environment variable is not set.\n");
      PetscFinalize();
      return EXIT_FAILURE;
  } else if (!t_final_str) {
      PetscPrintf(PETSC_COMM_WORLD, "Error: t_final_str environment variable is not set.\n");
      PetscFinalize();
      return EXIT_FAILURE;
  } else if (!temp_str) {
      PetscPrintf(PETSC_COMM_WORLD, "Error: temp_str environment variable is not set.\n");
      PetscFinalize();
      return EXIT_FAILURE;
  } else if (!grad_temp0X_str) {
      PetscPrintf(PETSC_COMM_WORLD, "Error: grad_temp0X_str environment variable is not set.\n");
      PetscFinalize();
      return EXIT_FAILURE;
  } else if (!grad_temp0Y_str) {
      PetscPrintf(PETSC_COMM_WORLD, "Error: grad_temp0Y_str environment variable is not set.\n");
      PetscFinalize();
      return EXIT_FAILURE;
  } else if (!grad_temp0Z_str) {
      PetscPrintf(PETSC_COMM_WORLD, "Error: grad_temp0Z_str environment variable is not set.\n");
      PetscFinalize();
      return EXIT_FAILURE;
  } else if (!dim_str) {
      PetscPrintf(PETSC_COMM_WORLD, "Error: dim_str environment variable is not set.\n");
      PetscFinalize();
      return EXIT_FAILURE;
  } else if (!eps_str) {
      PetscPrintf(PETSC_COMM_WORLD, "Error: eps_str environment variable is not set.\n");
      PetscFinalize();
      return EXIT_FAILURE;
  }

  // Convert environment variables to appropriate types
  char *endptr;
  PetscInt angle       = strtod(angle_str, &endptr);
  PetscInt mmm         = strtod(mmm_str, &endptr);
  PetscInt aaa         = strtod(aaa_str, &endptr);

  PetscInt Nx          = strtod(Nx_str, &endptr);
  PetscInt Ny          = strtod(Ny_str, &endptr);
  PetscInt Nz          = strtod(Nz_str, &endptr);

  PetscReal Lx          = strtod(Lx_str, &endptr);
  PetscReal Ly          = strtod(Ly_str, &endptr);
  PetscReal Lz          = strtod(Lz_str, &endptr);

  PetscReal delt_t      = strtod(delt_t_str, &endptr);
  PetscReal t_final     = strtod(t_final_str, &endptr);
  PetscInt n_out        = strtod(n_out_str, &endptr);

  // PetscReal humidity    = strtod(humidity_str, &endptr);
  PetscReal temp        = strtod(temp_str, &endptr);

  PetscReal grad_temp0X = strtod(grad_temp0X_str, &endptr);
  PetscReal grad_temp0Y = strtod(grad_temp0Y_str, &endptr);
  PetscReal grad_temp0Z = strtod(grad_temp0Z_str, &endptr);

  PetscInt dim          = strtod(dim_str, &endptr);
  
	PetscReal eps         = strtod(eps_str, &endptr);

  // Verify that conversion was successful
  if (*endptr != '\0') {
      PetscPrintf(PETSC_COMM_WORLD, "Error: One or more environment variables contain invalid values.\n");
      PetscFinalize();
      return EXIT_FAILURE;
  }

  // Print out the environment variables
  PetscPrintf(PETSC_COMM_WORLD, "Environment variables successfully set.\n");
  PetscPrintf(PETSC_COMM_WORLD, "angle: %d\n", angle);
  PetscPrintf(PETSC_COMM_WORLD, "mmm: %d\n", mmm);
  PetscPrintf(PETSC_COMM_WORLD, "aaa: %d\n", aaa);
  PetscPrintf(PETSC_COMM_WORLD, "Nx: %d\n", Nx);
  PetscPrintf(PETSC_COMM_WORLD, "Ny: %d\n", Ny);
  PetscPrintf(PETSC_COMM_WORLD, "Nz: %d\n", Nz);
  PetscPrintf(PETSC_COMM_WORLD, "Lx: %f\n", Lx);
  PetscPrintf(PETSC_COMM_WORLD, "Ly: %f\n", Ly);
  PetscPrintf(PETSC_COMM_WORLD, "Lz: %f\n", Lz);
  PetscPrintf(PETSC_COMM_WORLD, "delt_t: %f\n", delt_t);
  PetscPrintf(PETSC_COMM_WORLD, "t_final: %f\n", t_final);
  PetscPrintf(PETSC_COMM_WORLD, "temp: %f\n", user.temp0);
  PetscPrintf(PETSC_COMM_WORLD, "grad_temp0X: %f\n", user.grad_temp0[0]);
  PetscPrintf(PETSC_COMM_WORLD, "grad_temp0Y: %f\n", user.grad_temp0[1]);
  PetscPrintf(PETSC_COMM_WORLD, "dim: %d\n", dim);
  PetscPrintf(PETSC_COMM_WORLD, "eps: %e\n", eps);

  // Assign necessary variables to user struct
  user.Nx = Nx;
  user.Ny = Ny;
  user.Nz = Nz;

  user.Lx = Lx;
  user.Ly = Ly;
  user.Lz = Lz;

  user.temp0 = temp;
  user.grad_temp0[0] = grad_temp0X;
  user.grad_temp0[1] = grad_temp0Y;

  user.eps = eps;


  user.xi_v       = 1.0e-3; //can be used safely all time
  user.xi_T       = 1.0e-2;
  user.flag_xiT   = 0;
  user.flag_it0   = 1; 
  user.flag_tIC   = 0;

  // user.eps        = 2.0e-7;
  user.nucleat    = 2.5;//5.0e6;
  user.Lambd      = 1.0;
  user.air_lim    = 1.0e-6;
  user.phi_L      = 1.0e-8;
  user.nsteps_IC  = 0; //16;

  user.lat_sol    = 3.34e5;
  user.lat_sub    = 2.83e6;
  user.thcond_ice = 2.29;
  user.thcond_wat = 0.554;
  user.thcond_air = 0.02;
  user.cp_ice     = 1.96e3;
  user.cp_wat     = 4.2e3;
  user.cp_air     = 1.044e3;
  user.rho_ice    = 917.0;
  user.rho_wat    = 1000.0;
  user.rho_air    = 1.341;
  user.dif_vap    = 2.178e-5;
  user.T_melt     = 0.0;
  user.tem_nucl   = -5.0;

  PetscReal lambda_sol, lambda_sub, lambda_eva, tau_sol, tau_sub, tau_eva;
  PetscReal diff_sol = 0.5*(user.thcond_wat/user.rho_wat/user.cp_wat + user.thcond_ice/user.rho_ice/user.cp_ice);
  PetscReal diff_sub = 0.5*(user.thcond_air/user.rho_air/user.cp_air + user.thcond_ice/user.rho_ice/user.cp_ice);
  PetscReal diff_eva = 0.5*(user.thcond_wat/user.rho_wat/user.cp_wat + user.thcond_air/user.rho_air/user.cp_air);
  lambda_sol    = a1*user.eps/d0_sol;
  tau_sol       = user.eps*lambda_sol*(beta_sol/a1 + a2*user.eps/diff_sol );
  lambda_sub    = a1*user.eps/d0_sub;
  tau_sub       = user.eps*lambda_sub*(beta_sub/a1 + a2*user.eps/diff_sub + a2*user.eps/user.dif_vap);
  lambda_eva    = a1*user.eps/d0_eva;
  tau_eva       = user.eps*lambda_eva*(beta_eva/a1 + a2*user.eps/diff_eva + a2*user.eps/user.dif_vap);

  user.mob_sol    = user.eps/3.0/tau_sol;
  user.mob_sub    = user.eps/3.0/tau_sub; 
  user.mob_eva    = user.eps/3.0/tau_eva; 
  user.mav = cbrt(user.mob_sub*user.mob_eva*user.mob_sol);
  user.alph_sol   = lambda_sol/tau_sol; 
  user.alph_sub   = lambda_sub/tau_sub;
  user.alph_eva   = lambda_eva/tau_eva; 

  // const char *env1 = "angl"; const char *env2 = "aa"; const char *env3 = "mm";
  // char *angl1, *aa1, *mm1; 
  // angl1 = getenv(env1); aa1 = getenv(env2); mm1 = getenv(env3);
  // angle = atoi(angl1); aaa = atoi(aa1); mmm = atoi(mm1);
  // PetscPrintf(PETSC_COMM_WORLD,"Options: angle:%d alph:%d mobil:%d \n",angle,aaa,mmm);

  if(mmm==1) user.mob_sol=user.mob_sub=user.mob_eva=user.mav;
  if(aaa==1) user.alph_sol = user.alph_sub = user.alph_eva = 0.0;

  if(angle==1) gamma_iv = gamma_iw = gamma_wv; // 120 degrees
  user.Etai       = gamma_iv + gamma_iw - gamma_wv;
  user.Etaw       = gamma_wv + gamma_iw - gamma_iv;
  user.Etaa       = gamma_iv + gamma_wv - gamma_iw;

  PetscPrintf(PETSC_COMM_WORLD,"SOLID: tau %.4e  lambda %.4e  M0 %.4e  alpha %.4e \n",tau_sol,lambda_sol,user.mob_sol,user.alph_sol);
  PetscPrintf(PETSC_COMM_WORLD,"SUBLI: tau %.4e  lambda %.4e  M0 %.4e  alpha %.4e \n",tau_sub,lambda_sub,user.mob_sub,user.alph_sub);
  PetscPrintf(PETSC_COMM_WORLD,"EVAPO: tau %.4e  lambda %.4e  M0 %.4e  alpha %.4e \n",tau_eva,lambda_eva,user.mob_eva,user.alph_eva);

//ice grains
  user.NCice      = 35;//9; //less than 200, otherwise update in user
  user.RCice      = 0.2e-4;//0.18e-4;
  user.RCice_dev  = 0.5;
  user.overl      = 0.9;
  user.seed       = 108;//129;

  //boundary conditions : "periodic" >> "fixed-T" >> "variable-T"
  user.BC_Tfix      = 1;    //fixed T on boundary
  user.BC_Tvar      = 0;    //variable T on boundary
  user.flag_contang = 0;    //wall-wat, wall-ice contact angle
  user.periodic     = 0;    //periodic BC
  if(user.periodic==1) { user.BC_Tvar = user.BC_Tfix = 0; user.flag_contang = 0; } //"periodic" >> "fixed-T" >> "variable-T"
  if(user.BC_Tfix==1) user.BC_Tvar = 0;

  user.temp_m_ampl = 2.0; 	user.temp_m_fre = 1.0; // must be != 0
  user.costhet = 0.0; // wall-wat wall-ice contact angle; activate flag_contan; (sin, not cos)

  //domain and mesh characteristics
  PetscInt  p=1, C=0;
  user.p=p; user.C=C;

  //initial conditions
  // user.temp0      = temp;
  // user.grad_temp0[0] = grad_temp0X;
  // user.grad_temp0[1] = grad_temp0Y;
  // user.grad_temp0[2] = grad_temp0Z;

  //time stepping
  //output
  user.outp = 0;    //--------------------- if !=0 : output save every 'outp'
  user.t_out = 0.0;  user.t_interv = t_final/(n_out-1); //output every t_interv

  PetscInt adap = 1;
  PetscInt NRmin = 2, NRmax = 4;
  PetscReal factor = pow(10.0,1.0/8.0);
  PetscReal dtmin = 0.1*delt_t, dtmax = 0.5*user.t_interv;
  if(dtmax>0.5*user.t_interv) PetscPrintf(PETSC_COMM_WORLD,"OUTPUT DATA ERROR: Reduce maximum time step, or increase t_interval \n\n");
  PetscInt max_rej = 10;
  if(adap==1) PetscPrintf(PETSC_COMM_WORLD,"Adapative time stepping scheme: NR_iter %d-%d  factor %.3f  dt0 %.2e  dt_range %.2e-%.2e  \n\n",NRmin,NRmax,factor,delt_t,dtmin,dtmax);

  PetscBool output=PETSC_TRUE,monitor=PETSC_TRUE;
  char initial[PETSC_MAX_PATH_LEN] = {0};
  char PFgeom[PETSC_MAX_PATH_LEN] = {0};
  PetscOptionsBegin(PETSC_COMM_WORLD, "", "Metamorph Options", "IGA");//CHKERRQ(ierr);
  ierr = PetscOptionsInt("-Nx", "number of elements along x dimension", __FILE__, Nx, &Nx, NULL);CHKERRQ(ierr);
  ierr = PetscOptionsInt("-Ny", "number of elements along y dimension", __FILE__, Ny, &Ny, NULL);CHKERRQ(ierr);
  ierr = PetscOptionsInt("-p", "polynomial order", __FILE__, p, &p, NULL);CHKERRQ(ierr);
  ierr = PetscOptionsInt("-C", "global continuity order", __FILE__, C, &C, NULL);CHKERRQ(ierr);
  ierr = PetscOptionsString("-initial_cond","Load initial solution from file",__FILE__,initial,initial,sizeof(initial),NULL);CHKERRQ(ierr);
  ierr = PetscOptionsString("-initial_PFgeom","Load initial ice geometry from file",__FILE__,PFgeom,PFgeom,sizeof(PFgeom),NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool("-meta_output","Enable output files",__FILE__,output,&output,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool("-meta_monitor","Monitor the solution",__FILE__,monitor,&monitor,NULL);CHKERRQ(ierr);
  PetscOptionsEnd();//CHKERRQ(ierr);

  IGA iga;
  ierr = IGACreate(PETSC_COMM_WORLD,&iga);CHKERRQ(ierr);
  ierr = IGASetDim(iga,dim);CHKERRQ(ierr);
  ierr = IGASetDof(iga,4);CHKERRQ(ierr);
  ierr = IGASetFieldName(iga,0,"phaseice"); CHKERRQ(ierr);
  ierr = IGASetFieldName(iga,1,"phaseair"); CHKERRQ(ierr);
  ierr = IGASetFieldName(iga,2,"temperature"); CHKERRQ(ierr);
  ierr = IGASetFieldName(iga,3,"vap_density"); CHKERRQ(ierr);

  IGAAxis axis0;
  ierr = IGAGetAxis(iga,0,&axis0);CHKERRQ(ierr);
  if(user.periodic==1) {ierr = IGAAxisSetPeriodic(axis0,PETSC_TRUE);CHKERRQ(ierr);}
  ierr = IGAAxisSetDegree(axis0,p);CHKERRQ(ierr);
  ierr = IGAAxisInitUniform(axis0,Nx,0.0,Lx,C);CHKERRQ(ierr);
  IGAAxis axis1;
  ierr = IGAGetAxis(iga,1,&axis1);CHKERRQ(ierr);
  if(user.periodic==1) {ierr = IGAAxisSetPeriodic(axis1,PETSC_TRUE);CHKERRQ(ierr);}
  ierr = IGAAxisSetDegree(axis1,p);CHKERRQ(ierr);
  ierr = IGAAxisInitUniform(axis1,Ny,0.0,Ly,C);CHKERRQ(ierr);

  ierr = IGASetFromOptions(iga);CHKERRQ(ierr);
  ierr = IGASetUp(iga);CHKERRQ(ierr);
  user.iga = iga;

  PetscInt ngp = iga->elem_width[0]*iga->elem_width[1]*SQ(p+1); // #Gauss points local
  ierr = PetscMalloc(sizeof(PetscInt)*(ngp),&user.FlagMob);CHKERRQ(ierr);
  ierr = PetscMemzero(user.FlagMob,sizeof(PetscInt)*(ngp));CHKERRQ(ierr);

  //Residual and Tangent
  ierr = IGASetFormIFunction(iga,Residual,&user);CHKERRQ(ierr);
  ierr = IGASetFormIJacobian(iga,Jacobian,&user);CHKERRQ(ierr);

//Boundary Conditions
  iga->BCmetam = user.BC_Tvar;
  iga->Temp_metam = user.temp0;

  if(user.periodic==1) PetscPrintf(PETSC_COMM_WORLD,"Periodic Boundary Conditions!\n");
  else if(user.BC_Tfix==1) PetscPrintf(PETSC_COMM_WORLD,"Fixed temperature on the boundary!\n");
  else if(user.BC_Tvar==1) PetscPrintf(PETSC_COMM_WORLD,"Time-dependent temperature on the boundary!\n");

  if(user.flag_contang==1){
    ierr = IGASetBoundaryForm(iga,0,0,PETSC_TRUE);CHKERRQ(ierr);
    ierr = IGASetBoundaryForm(iga,0,1,PETSC_TRUE);CHKERRQ(ierr);
    ierr = IGASetBoundaryForm(iga,1,0,PETSC_TRUE);CHKERRQ(ierr);
    ierr = IGASetBoundaryForm(iga,1,1,PETSC_TRUE);CHKERRQ(ierr);
  }
  if(user.BC_Tfix==1){
    PetscReal Ttop,Tbot,Tlef,Trig;
    Tlef = user.temp0 - user.grad_temp0[0]*0.5*Lx;
    Trig = user.temp0 + user.grad_temp0[0]*0.5*Lx;  
    Tbot = user.temp0 - user.grad_temp0[1]*0.5*Ly;
    Ttop = user.temp0 + user.grad_temp0[1]*0.5*Ly;
    ierr = IGASetBoundaryValue(iga,0,0,2,Tlef);CHKERRQ(ierr);
    ierr = IGASetBoundaryValue(iga,0,1,2,Trig);CHKERRQ(ierr);
    ierr = IGASetBoundaryValue(iga,1,0,2,Tbot);CHKERRQ(ierr);
    ierr = IGASetBoundaryValue(iga,1,1,2,Ttop);CHKERRQ(ierr);
  }

  TS ts;
  ierr = IGACreateTS(iga,&ts);CHKERRQ(ierr);
  ierr = TSSetMaxTime(ts,t_final);CHKERRQ(ierr);
  ierr = TSSetExactFinalTime(ts,TS_EXACTFINALTIME_MATCHSTEP);CHKERRQ(ierr);
  ierr = TSSetTimeStep(ts,delt_t);CHKERRQ(ierr);
  ierr = TSSetType(ts,TSALPHA);CHKERRQ(ierr);
  ierr = TSAlphaSetRadius(ts,0.5);CHKERRQ(ierr);
  if (monitor) {ierr = TSMonitorSet(ts,Monitor,&user,NULL);CHKERRQ(ierr);}
  if (output)  {ierr = TSMonitorSet(ts,OutputMonitor,&user,NULL);CHKERRQ(ierr);}
  ierr = TSSetFromOptions(ts);CHKERRQ(ierr);

  ts->adap = adap;
  ts->NRmin = NRmin;
  ts->NRmax = NRmax;
  ts->factor = factor;
  ts->dtmax = dtmax;
  ts->dtmin = dtmin;
  ts->max_reject = max_rej;
  ts->max_snes_failures = -1;

  SNES nonlin;
  ierr = TSGetSNES(ts,&nonlin);CHKERRQ(ierr);
  ierr = SNESSetConvergenceTest(nonlin,SNESDOFConvergence,&user,NULL);CHKERRQ(ierr);

  ierr = InitialIceGrains(iga,&user);CHKERRQ(ierr);

  PetscReal t=0; Vec U;
  ierr = IGACreateVec(iga,&U);CHKERRQ(ierr);
  ierr = FormInitialCondition(iga,t,U,&user,initial,PFgeom);CHKERRQ(ierr);
  ierr = TSSolve(ts,U);CHKERRQ(ierr);

  ierr = VecDestroy(&U);CHKERRQ(ierr);
  ierr = TSDestroy(&ts);CHKERRQ(ierr);
  ierr = IGADestroy(&iga);CHKERRQ(ierr);

  ierr = PetscFree(user.FlagMob);CHKERRQ(ierr);

  PetscLogDouble ltim,tim;
  ierr = PetscTime(&ltim); CHKERRQ(ierr);
  tim = ltim-itim;
  PetscPrintf(PETSC_COMM_WORLD," comp time %e sec  =  %.2f min \n\n",tim,tim/60.0);

  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}


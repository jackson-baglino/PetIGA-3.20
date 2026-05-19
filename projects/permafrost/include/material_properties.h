#ifndef MATERIAL_PROPERTIES_H
#define MATERIAL_PROPERTIES_H

#include "NASA_types.h"

/* Computes effective thermal conductivity and its derivative with respect to ice */
void ThermalCond(AppCtx *user, PetscScalar ice, PetscScalar sed, PetscScalar *cond, PetscScalar *dcond_ice);

/* Computes effective heat capacity and its derivative with respect to ice */
void HeatCap(AppCtx *user, PetscScalar ice, PetscScalar sed, PetscScalar *cp, PetscScalar *dcp_ice);

/* Computes effective density and its derivative with respect to ice */
void Density(AppCtx *user, PetscScalar ice, PetscScalar sed, PetscScalar *rho, PetscScalar *drho_ice);

/* Computes vapor diffusivity and its temperature derivative */
void VaporDiffus(AppCtx *user, PetscScalar tem, PetscScalar *difvap, PetscScalar *d_difvap);

/* Computes the smooth Heaviside function and its derivative*/
void SmoothHeavisidePoly(PetscScalar phi, PetscScalar *g, PetscScalar *dg_dphi);

/* Penalty weight: smooth Heaviside concentrated near phi = 1.
 * Zero for phi <= PENALTY_PHI_LO, unity for phi >= PENALTY_PHI_HI,
 * smooth ramp in between (via SmoothHeavisidePoly of the shifted variable).
 * Used so the vapor interface-equilibrium penalty is active only deep in
 * solid (phi = ice+sed close to 1) and OFF at the diffuse ice-air interface
 * so Gibbs-Thomson curvature dependence can emerge. */
void PenaltyWeight(PetscScalar phi, PetscScalar *g, PetscScalar *dg_dphi);

/* Computes the saturation vapor density and its derivative */
void RhoVS_I(AppCtx *user, PetscScalar tem, PetscScalar *rho_vs, PetscScalar *d_rhovs);

/* Computes the free energy function for the ice phase and its derivative */
void Fice(AppCtx *user, PetscScalar ice, PetscScalar sed, PetscScalar *fice, PetscScalar *dfice_ice);

/* Computes the phase evolution function for the water phase and its derivative */
void Fsed(AppCtx *user, PetscScalar ice, PetscScalar sed, PetscScalar *fsed, PetscScalar *dfsed_ice);

/* Computes the phase evolution function for the air phase and its derivative */
void Fair(AppCtx *user, PetscScalar ice, PetscScalar sed, PetscScalar *fair, PetscScalar *dfair_ice);

/* computes the phase dependent mobility*/
void Mobility(AppCtx *user, PetscScalar ice, PetscScalar sed, PetscScalar *mob);

/* Computes sigma0 using logarithmic interpolation based on temperature */
void Sigma0(PetscScalar temp, PetscScalar *sigm0);

#endif // MATERIAL_PROPERTIES_H
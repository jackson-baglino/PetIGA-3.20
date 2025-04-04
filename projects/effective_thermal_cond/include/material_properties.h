#ifndef MATERIAL_PROPERTIES_H
#define MATERIAL_PROPERTIES_H

#include "user_context.h"

// ========================= Function Declarations =========================
void ThermalCond(AppCtx *user, PetscScalar ice, PetscScalar *cond, PetscScalar *dcond_ice);
void HeatCap(AppCtx *user, PetscScalar ice, PetscScalar *cp, PetscScalar *dcp_ice);
void Density(AppCtx *user, PetscScalar ice, PetscScalar *rho, PetscScalar *drho_ice);

#endif // ENV_CONFIG_H
#include "material_properties.h"

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

/* Function: HeatCap  
   Computes heat capacity and its derivative w.r.t. ice content. */
void HeatCap(AppCtx *user, PetscScalar ice, PetscScalar *cp, PetscScalar *dcp_ice)
{
    PetscReal dice = 1.0, dair = 1.0;
    PetscReal air = 1.0 - ice;
    if (ice < 0.0) { ice = 0.0; dice = 0.0; }
    if (air < 0.0) { air = 0.0; dair = 0.0; }

    PetscReal cp_ice = user->cp_ice;
    PetscReal cp_air = user->cp_air;

    if (cp)       (*cp)      = ice * cp_ice + air * cp_air;
    if (dcp_ice)  (*dcp_ice) = cp_ice * dice - cp_air * dair;

    return;
}

/* Function: Density  
   Computes density and its derivative w.r.t. ice content. */
void Density(AppCtx *user, PetscScalar ice, PetscScalar *rho, PetscScalar *drho_ice)
{
    PetscReal dice = 1.0, dair = 1.0;
    PetscReal air = 1.0 - ice;
    if (ice < 0.0) { ice = 0.0; dice = 0.0; }
    if (air < 0.0) { air = 0.0; dair = 0.0; }

    PetscReal rho_ice = user->rho_ice;
    PetscReal rho_air = user->rho_air;

    if (rho)        (*rho)      = ice * rho_ice + air * rho_air;
    if (drho_ice)   (*drho_ice) = rho_ice * dice - rho_air * dair;

    return;
}
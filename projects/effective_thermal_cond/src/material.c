#include "material.h"

/*-----------------------------------------------------------------------------
  ThermalCond
  Linear mixture rule: k(x) = phi_ice * k_ice + (1 - phi_ice) * k_air
  Clamps ice and air fractions to [0, 1] before computing.
  Pass NULL for outputs you do not need.
-----------------------------------------------------------------------------*/
void ThermalCond(AppCtx *user, PetscScalar ice,
                 PetscScalar *cond, PetscScalar *dcond_ice)
{
  PetscReal dice = 1.0, dair = 1.0;
  PetscReal air  = 1.0 - ice;

  if (ice < 0.0) { ice  = 0.0; dice = 0.0; }
  if (air < 0.0) { air  = 0.0; dair = 0.0; }

  if (cond)      *cond      = ice * user->thcond_ice + air * user->thcond_air;
  if (dcond_ice) *dcond_ice = user->thcond_ice * dice - user->thcond_air * dair;
}

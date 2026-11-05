#ifndef MONITORING_H
#define MONITORING_H

#include "NASA_types.h"

/* Monitor function for TS solver */
PetscErrorCode Monitor(TS ts, PetscInt step, PetscReal t, Vec U, void *mctx);

/* Output monitor function for writing solution data */
PetscErrorCode OutputMonitor(TS ts, PetscInt step, PetscReal t, Vec U, void *mctx);

/* Interface-CFL timestep limiter: clamps the next dt from the last accepted
   step's max pointwise phase-change rate (see monitoring.c) */
PetscErrorCode InterfaceCFLMonitor(TS ts, PetscInt step, PetscReal t, Vec U, void *mctx);

/* Physical (unscaled) vapor mass flux and heat flux out through every face
 * that carries a T or rho_v Dirichlet value; one row per step in bflux.dat. */
PetscErrorCode BoundaryFluxMonitor(TS ts, PetscInt step, PetscReal t, Vec U, void *mctx);

/* TS pre-step callback: consumes deferred bounds-rollback requests set by Monitor() */
PetscErrorCode BoundsRollbackPreStep(TS ts);

#endif // MONITORING_H
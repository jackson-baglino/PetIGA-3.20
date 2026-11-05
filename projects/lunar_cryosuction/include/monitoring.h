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

/* Resolution checker: measured interface speed against the K&P Eq. 45 bound
 * eps < d0/(beta0 v_n), and ice on open faces. Read-only; resolution.dat. */
PetscErrorCode ResolutionMonitor(TS ts, PetscInt step, PetscReal t, Vec U, void *mctx);

/* Start-up static check (mesh spacing against eps/sqrt(2), eps against the
 * smallest IC grain radius) and the end-of-run summary of ResolutionMonitor. */
PetscErrorCode ResolutionStaticCheck(AppCtx *user);
PetscErrorCode ResolutionSummary(AppCtx *user);

/* TS pre-step callback: consumes deferred bounds-rollback requests set by Monitor() */
PetscErrorCode BoundsRollbackPreStep(TS ts);

#endif // MONITORING_H
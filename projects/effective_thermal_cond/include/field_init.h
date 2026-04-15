#ifndef FIELD_INIT_H
#define FIELD_INIT_H

#include "app_ctx.h"

/*
 * AllocateAppCtxFields — allocate user->ice for the local Gauss-point count.
 * Uses the local element widths from the given IGA object.
 */
PetscErrorCode AllocateAppCtxFields(IGA iga, AppCtx *user, PetscScalar **field);

/*
 * EvaluateFieldAtGaussPoints — project a PETSc solution vector onto the local
 * Gauss-point array stored in user->ice.
 */
PetscErrorCode EvaluateFieldAtGaussPoints(AppCtx *user, IGA iga, Vec vec_phase);

/*
 * ReadSolutionVec — read an IGA + solution vector from file and populate user->ice.
 *   iga_file : path to igasol.dat
 *   vec_file : path to sol_NNNNN.dat
 */
PetscErrorCode ReadSolutionVec(const char *iga_file, const char *vec_file,
                               IGA *iga_out, AppCtx *user);

/*
 * ComputeCircleIceField  — populate user->ice with a circular grain centred in
 * the domain (tanh interface, radius = min(Lx,Ly)/16).
 */
PetscErrorCode ComputeCircleIceField(AppCtx *user);

/*
 * ComputeLayeredIceField — populate user->ice with a horizontal layer: ice below
 * y = Ly/2, air above (tanh interface).
 */
PetscErrorCode ComputeLayeredIceField(AppCtx *user);

/*
 * FormInitialCondition — dispatch to the appropriate initialiser based on
 * user->init_mode ("circle", "layered", or file-based).
 */
PetscErrorCode FormInitialCondition(AppCtx *user);

#endif /* FIELD_INIT_H */

#include "snes_convergence.h"
#include <petscstring.h>

PetscErrorCode SNESDOFConvergence(SNES snes, PetscInt it_number, PetscReal xnorm,
  PetscReal gnorm, PetscReal fnorm, SNESConvergedReason *reason, void *cctx)
{
    PetscFunctionBegin;
    PetscErrorCode ierr;
    AppCtx *user = (AppCtx *)cctx;

    Vec Res, Sol, Sol_upd;
    PetscInt  i, dof = user->iga->dof;
    PetscReal n2dof[dof], sol_n2dof[dof], sol_update_n2dof[dof];

    // Retrieve residual vector and compute per-DOF 2-norms
    ierr = SNESGetFunction(snes, &Res, 0, 0); CHKERRQ(ierr);
    for (i = 0; i < dof; i++) {
        ierr = VecStrideNorm(Res, i, NORM_2, &n2dof[i]); CHKERRQ(ierr);
    }

    // Retrieve solution and solution update vectors
    ierr = SNESGetSolution(snes, &Sol); CHKERRQ(ierr);
    ierr = SNESGetSolutionUpdate(snes, &Sol_upd); CHKERRQ(ierr);

    for (i = 0; i < dof; i++) {
        ierr = VecStrideNorm(Sol,     i, NORM_2, &sol_n2dof[i]);     CHKERRQ(ierr);
        ierr = VecStrideNorm(Sol_upd, i, NORM_2, &sol_update_n2dof[i]); CHKERRQ(ierr);
    }

    // Store initial residual norms at iteration 0 for relative convergence checks
    if (it_number == 0) {
        for (i = 0; i < dof; i++) {
            user->norm0[i] = n2dof[i];
            if (user->norm0[i] < 1.0e-30) user->norm0[i] = 1.0;
        }
    }

    // -------------------------------------------------------------------------
    // Print header at iteration 0
    // -------------------------------------------------------------------------
    const PetscInt W_IT    = 2;
    const PetscInt W_FNORM = 11;
    const PetscInt W_N     = 11;
    const PetscInt W_R     = 9;
    const PetscInt W_S     = 9;

    if (it_number == 0) {
        char header_line[512];
        PetscInt nchar = 0;

        PetscSNPrintf(header_line, sizeof(header_line),
                      "    %*s | %*s | %*s | %*s | %*s | %*s | %*s | %*s | %*s | %*s | %*s",
                      (int)W_IT,    "it",
                      (int)W_FNORM, "fnorm",
                      (int)W_N,     "n0",
                      (int)W_R,     "r0",
                      (int)W_S,     "s0",
                      (int)W_N,     "n1",
                      (int)W_R,     "r1",
                      (int)W_S,     "s1",
                      (int)W_N,     "n2",
                      (int)W_R,     "r2",
                      (int)W_S,     "s2");

        PetscStrlen(header_line, (size_t*)&nchar);

        PetscPrintf(PETSC_COMM_WORLD, "    SNES DOF norms (2-norm) and ratios\n");
        PetscPrintf(PETSC_COMM_WORLD, "    %.*s\n", (int)nchar,
                    "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------");
        PetscPrintf(PETSC_COMM_WORLD, "%s\n", header_line);
        PetscPrintf(PETSC_COMM_WORLD, "    %.*s\n", (int)nchar,
                    "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------");
    }

    // -------------------------------------------------------------------------
    // Compute per-DOF relative residual and step size ratios
    // -------------------------------------------------------------------------
    PetscInt  nd_print = (dof < 3) ? dof : 3;
    PetscReal rel[3]   = {0.0, 0.0, 0.0};
    PetscReal step[3]  = {0.0, 0.0, 0.0};

    for (i = 0; i < nd_print; i++) {
        rel[i]  = (user->norm0[i] > 0.0) ? (n2dof[i] / user->norm0[i]) : 0.0;
        step[i] = (sol_n2dof[i]  > 0.0) ? (sol_update_n2dof[i] / sol_n2dof[i]) : 0.0;
    }

    // Print one line per Newton iteration
    PetscPrintf(PETSC_COMM_WORLD,
                "    %*d | %*.*e | %*.*e | %*.*e | %*.*e | %*.*e | %*.*e | %*.*e | %*.*e | %*.*e | %*.*e\n",
                (int)W_IT,    (int)it_number,
                (int)W_FNORM, 4, (double)fnorm,
                (int)W_N,     3, (double)((nd_print > 0) ? n2dof[0]  : 0.0),
                (int)W_R,     3, (double)((nd_print > 0) ? rel[0]    : 0.0),
                (int)W_S,     3, (double)((nd_print > 0) ? step[0]   : 0.0),
                (int)W_N,     3, (double)((nd_print > 1) ? n2dof[1]  : 0.0),
                (int)W_R,     3, (double)((nd_print > 1) ? rel[1]    : 0.0),
                (int)W_S,     3, (double)((nd_print > 1) ? step[1]   : 0.0),
                (int)W_N,     3, (double)((nd_print > 2) ? n2dof[2]  : 0.0),
                (int)W_R,     3, (double)((nd_print > 2) ? rel[2]    : 0.0),
                (int)W_S,     3, (double)((nd_print > 2) ? step[2]   : 0.0));

    // -------------------------------------------------------------------------
    // Retrieve solver tolerances
    // -------------------------------------------------------------------------
    PetscScalar atol, rtol, stol;
    PetscInt    maxit, maxf;
    ierr = SNESGetTolerances(snes, &atol, &rtol, &stol, &maxit, &maxf); CHKERRQ(ierr);

    *reason = SNES_CONVERGED_ITERATING;

    /* At iteration 0 only the absolute criterion is meaningful: rel[i]=1 by
     * construction (norm0 is set from this residual) and the solution update
     * is stale, so rtol/stol cannot be consulted.  But the absolute check
     * must fire here: if every DOF residual is already below atol there is
     * nothing for Newton to do, and forcing a step on a machine-noise
     * residual makes the KSP fail before iteration 1.  At dtmin that trapped
     * the timestep controller in an infinite reject/retry loop — the rtol
     * loosening below only runs at it >= 1, which was never reached. */
    if (it_number < 1) {
        PetscInt  nd0     = (dof < 3) ? dof : 3;
        PetscBool all_abs = PETSC_TRUE;
        for (i = 0; i < nd0; i++) {
            if (!(n2dof[i] < atol)) { all_abs = PETSC_FALSE; break; }
        }
        if (all_abs) {
            PetscPrintf(PETSC_COMM_WORLD,
                        "  SNES converged at it 0: all DOF residuals below atol\n");
            *reason = SNES_CONVERGED_FNORM_ABS;
        }
        PetscFunctionReturn(0);
    }

    // If a previous time step required reduction, loosen relative tolerance
    if (snes->prev_dt_red == 1) {
        rtol *= 10.0;
    }

    // -------------------------------------------------------------------------
    // Per-DOF convergence check
    // nd is clamped to actual dof count to avoid accessing uninitialized entries
    // -------------------------------------------------------------------------
    PetscInt  nd = (dof < 3) ? dof : 3;
    PetscBool conv_rel[3]  = {PETSC_FALSE, PETSC_FALSE, PETSC_FALSE};
    PetscBool conv_abs[3]  = {PETSC_FALSE, PETSC_FALSE, PETSC_FALSE};
    PetscBool conv_stol[3] = {PETSC_FALSE, PETSC_FALSE, PETSC_FALSE};
    PetscBool conv_dof[3]  = {PETSC_FALSE, PETSC_FALSE, PETSC_FALSE};

    for (i = 0; i < nd; i++) {
        conv_rel[i]  = (n2dof[i] <= rtol * user->norm0[i])          ? PETSC_TRUE : PETSC_FALSE;
        conv_abs[i]  = (n2dof[i] <  atol)                            ? PETSC_TRUE : PETSC_FALSE;
        conv_stol[i] = (sol_update_n2dof[i] <= stol * sol_n2dof[i]) ? PETSC_TRUE : PETSC_FALSE;
        conv_dof[i]  = (conv_rel[i] || conv_abs[i] || conv_stol[i]) ? PETSC_TRUE : PETSC_FALSE;
    }

    // All active DOFs must converge — for dof < 3, unset entries default to FALSE
    // so we only AND over the actual dof count
    PetscBool all_converged = PETSC_TRUE;
    for (i = 0; i < nd; i++) {
        if (!conv_dof[i]) { all_converged = PETSC_FALSE; break; }
    }

    if (all_converged) {
        PetscPrintf(PETSC_COMM_WORLD, "  SNES converged via DOF-wise criteria:\n");
        for (i = 0; i < nd; i++) {
            PetscPrintf(PETSC_COMM_WORLD, "    DOF %d:", i);
            if (conv_rel[i])  PetscPrintf(PETSC_COMM_WORLD, " REL (rtol)");
            if (conv_abs[i])  PetscPrintf(PETSC_COMM_WORLD, " ABS (atol)");
            if (conv_stol[i]) PetscPrintf(PETSC_COMM_WORLD, " STOL (stol)");
            PetscPrintf(PETSC_COMM_WORLD, "\n");
        }
        *reason = SNES_CONVERGED_FNORM_RELATIVE;
    }

    PetscFunctionReturn(0);
}
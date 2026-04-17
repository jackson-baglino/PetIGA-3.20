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
                      "    %*s | %*s | %*s | %*s | %*s | %*s | %*s | %*s | %*s | %*s | %*s | %*s | %*s | %*s",
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
                      (int)W_S,     "s2",
                      (int)W_N,     "n3",
                      (int)W_R,     "r3",
                      (int)W_S,     "s3");

        PetscStrlen(header_line, (size_t*)&nchar);

        PetscPrintf(PETSC_COMM_WORLD, "    SNES DOF norms (2-norm) and ratios\n");
        PetscPrintf(PETSC_COMM_WORLD, "    %.*s\n", (int)nchar,
                    "--------------------------------------------------------------------------------------------------------------------------------------");
        PetscPrintf(PETSC_COMM_WORLD, "%s\n", header_line);
        PetscPrintf(PETSC_COMM_WORLD, "    %.*s\n", (int)nchar,
                    "--------------------------------------------------------------------------------------------------------------------------------------");
    }

    // -------------------------------------------------------------------------
    // Compute per-DOF relative residual and step size ratios
    // -------------------------------------------------------------------------
    PetscInt  nd_print = (dof < 4) ? dof : 4;
    PetscReal rel[4]   = {0.0, 0.0, 0.0, 0.0};
    PetscReal step[4]  = {0.0, 0.0, 0.0, 0.0};

    for (i = 0; i < nd_print; i++) {
        rel[i]  = (user->norm0[i] > 0.0) ? (n2dof[i] / user->norm0[i]) : 0.0;
        step[i] = (sol_n2dof[i]  > 0.0) ? (sol_update_n2dof[i] / sol_n2dof[i]) : 0.0;
    }

    // Print one line per Newton iteration
    PetscPrintf(PETSC_COMM_WORLD,
                "    %*d | %*.*e | %*.*e | %*.*e | %*.*e | %*.*e | %*.*e | %*.*e | %*.*e | %*.*e | %*.*e | %*.*e | %*.*e | %*.*e\n",
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
                (int)W_S,     3, (double)((nd_print > 2) ? step[2]   : 0.0),
                (int)W_N,     3, (double)((nd_print > 3) ? n2dof[3]  : 0.0),
                (int)W_R,     3, (double)((nd_print > 3) ? rel[3]    : 0.0),
                (int)W_S,     3, (double)((nd_print > 3) ? step[3]   : 0.0));

    // -------------------------------------------------------------------------
    // Retrieve solver tolerances
    // -------------------------------------------------------------------------
    PetscScalar atol, rtol, stol;
    PetscInt    maxit, maxf;
    ierr = SNESGetTolerances(snes, &atol, &rtol, &stol, &maxit, &maxf); CHKERRQ(ierr);

    *reason = SNES_CONVERGED_ITERATING;

    // Enforce minimum of 3 Newton iterations before checking convergence
    if (it_number < 3) {
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
    PetscInt  nd = (dof < 4) ? dof : 4;
    PetscBool conv_rel[4]  = {PETSC_FALSE, PETSC_FALSE, PETSC_FALSE, PETSC_FALSE};
    PetscBool conv_abs[4]  = {PETSC_FALSE, PETSC_FALSE, PETSC_FALSE, PETSC_FALSE};
    PetscBool conv_stol[4] = {PETSC_FALSE, PETSC_FALSE, PETSC_FALSE, PETSC_FALSE};
    PetscBool conv_dof[4]  = {PETSC_FALSE, PETSC_FALSE, PETSC_FALSE, PETSC_FALSE};

    for (i = 0; i < nd; i++) {
        conv_rel[i]  = (n2dof[i] <= rtol * user->norm0[i]) ? PETSC_TRUE : PETSC_FALSE;
        conv_abs[i]  = (n2dof[i] <  atol)                  ? PETSC_TRUE : PETSC_FALSE;
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
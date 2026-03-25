#include "snes_convergence.h"
#include <petscstring.h>

PetscErrorCode SNESDOFConvergence(SNES snes, PetscInt it_number, PetscReal xnorm, 
  PetscReal gnorm, PetscReal fnorm, SNESConvergedReason *reason, void *cctx)
{
    PetscFunctionBegin;
    PetscErrorCode ierr;
    AppCtx *user = (AppCtx *)cctx;

    // Define vectors for residual, solution, and solution update
    Vec Res, Sol, Sol_upd;
    PetscInt  i, dof = user->iga->dof;
    PetscReal n2dof[dof], sol_n2dof[dof], sol_update_n2dof[dof];  /* component-wise 2-norms */
    PetscReal solv = 0.0, solupdv = 0.0;                          /* norms for DOF 2 when flag_tIC==1 */

    // Retrieve the residual vector from SNES and compute its component norms.
    ierr = SNESGetFunction(snes, &Res, 0, 0);CHKERRQ(ierr);
    for (i = 0; i < dof; i++) {ierr = VecStrideNorm(Res, i, NORM_2, &n2dof[i]);CHKERRQ(ierr);}

    ierr = SNESGetSolution(snes, &Sol);CHKERRQ(ierr);
    ierr = SNESGetSolutionUpdate(snes, &Sol_upd);CHKERRQ(ierr);

    /* Compute per-DOF norms of the solution and the solution update */
    for (i = 0; i < dof; i++) {
        ierr = VecStrideNorm(Sol,     i, NORM_2, &sol_n2dof[i]);CHKERRQ(ierr);
        ierr = VecStrideNorm(Sol_upd, i, NORM_2, &sol_update_n2dof[i]);CHKERRQ(ierr);
    }

    /* Optional special-case print helpers for DOF 2 */
    if (user->flag_tIC == 1 && dof > 2) {
        solv   = sol_n2dof[2];
        solupdv = sol_update_n2dof[2];
    }

    // Store initial residual norms at the first iteration for relative convergence checks.
    if (it_number == 0) {
        for (i = 0; i < dof; i++) {
            user->norm0[i] = n2dof[i];
            sol_n2dof[i] = sol_update_n2dof[i];
            if (user->norm0[i] < 1.0e-30) user->norm0[i] = 1.0;  /* prevent division by zero later */
        }
    }

    // Print iteration information and norm values (readable, aligned table)
    // Use identical fixed-width fields for header and data so columns line up.
    // (Assumes a monospaced terminal output.)

    // Column widths must match the data-row format below
    const PetscInt W_IT    = 2;
    const PetscInt W_FNORM = 11;
    const PetscInt W_N     = 11;
    const PetscInt W_R     = 9;
    const PetscInt W_S     = 9;

    // Build the header using the SAME field widths as the data row
    if (it_number == 0) {
        // Total printed width of the table line (excluding the trailing newline)
        // 4 spaces indent + (it) + 1 space + '|' separators and spaces
        // We compute this from the literal header we print.
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
        // Print full-width separator lines matching the header length
        PetscPrintf(PETSC_COMM_WORLD, "    %.*s\n", (int)nchar,
                    "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------");
        PetscPrintf(PETSC_COMM_WORLD, "%s\n", header_line);
        PetscPrintf(PETSC_COMM_WORLD, "    %.*s\n", (int)nchar,
                    "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------");
    }

    // Compute per-DOF relative residual and step ratios (limit to first 3 DOFs)
    PetscReal rel[3]  = {0.0, 0.0, 0.0};
    PetscReal step[3] = {0.0, 0.0, 0.0};
    PetscInt  nd_print = (dof < 3) ? dof : 3;

    for (i = 0; i < nd_print; i++) {
        rel[i] = (user->norm0[i] > 0.0) ? (n2dof[i] / user->norm0[i]) : 0.0;

        // step = ||dX|| / ||X|| (per DOF), with DOF-2 special-case for flag_tIC
        if (user->flag_tIC == 1 && i == 2) {
            step[i] = (solv > 0.0) ? (solupdv / solv) : 0.0;
        } else {
            step[i] = (sol_n2dof[i] > 0.0) ? (sol_update_n2dof[i] / sol_n2dof[i]) : 0.0;
        }
    }

    // One compact, aligned line per nonlinear iteration
    PetscPrintf(PETSC_COMM_WORLD,
                "    %*d | %*.*e | %*.*e | %*.*e | %*.*e | %*.*e | %*.*e | %*.*e | %*.*e | %*.*e | %*.*e\n",
                (int)W_IT,    (int)it_number,
                (int)W_FNORM, 4, (double)fnorm,
                (int)W_N,     3, (double)((nd_print > 0) ? n2dof[0] : 0.0),
                (int)W_R,     3, (double)((nd_print > 0) ? rel[0]  : 0.0),
                (int)W_S,     3, (double)((nd_print > 0) ? step[0] : 0.0),
                (int)W_N,     3, (double)((nd_print > 1) ? n2dof[1] : 0.0),
                (int)W_R,     3, (double)((nd_print > 1) ? rel[1]  : 0.0),
                (int)W_S,     3, (double)((nd_print > 1) ? step[1] : 0.0),
                (int)W_N,     3, (double)((nd_print > 2) ? n2dof[2] : 0.0),
                (int)W_R,     3, (double)((nd_print > 2) ? rel[2]  : 0.0),
                (int)W_S,     3, (double)((nd_print > 2) ? step[2] : 0.0));

    // Print a footer when we are about to return converged (i.e., last line before convergence report)
    // NOTE: we don't know convergence yet at this point, so we don't print the footer here.

    // Retrieve SNES solver tolerances (absolute, relative, step-size)
    PetscScalar atol, rtol, stol;
    PetscInt maxit, maxf;
    ierr = SNESGetTolerances(snes, &atol, &rtol, &stol, &maxit, &maxf);CHKERRQ(ierr);
    *reason = SNES_CONVERGED_ITERATING;

    if (it_number < 3) {
        PetscFunctionReturn(0);
    }

    // If a previous timestep required reduction, increase relative tolerance.
    if (snes->prev_dt_red == 1) {
        rtol *= 10.0;
    }

    PetscBool conv_rel[3], conv_abs[3], conv_stol[3], conv_dof[3];
    PetscInt  d, nd = (dof < 3) ? dof : 3;

    /* Evaluate convergence criteria per DOF */
    for (d = 0; d < nd; d++) {
        conv_rel[d]  = (n2dof[d] <= rtol * user->norm0[d]) ? PETSC_TRUE : PETSC_FALSE;
        conv_abs[d]  = (n2dof[d] < atol)                   ? PETSC_TRUE : PETSC_FALSE;
        conv_stol[d] = (sol_update_n2dof[d] <= stol * sol_n2dof[d]) ? PETSC_TRUE : PETSC_FALSE;

        /* SAME logic as your original OR condition */
        conv_dof[d] = (conv_rel[d] || conv_abs[d] || conv_stol[d]) ? PETSC_TRUE : PETSC_FALSE;
    }

    /* SAME logic as your original AND condition */
    if (nd == 3 && conv_dof[0] && conv_dof[1] && conv_dof[2]) {

    PetscPrintf(PETSC_COMM_WORLD,
                "  SNES converged via DOF-wise criteria:\n");

    for (d = 0; d < 3; d++) {
        PetscPrintf(PETSC_COMM_WORLD, "    DOF %d:", d);

        if (conv_rel[d])
        PetscPrintf(PETSC_COMM_WORLD, " REL (rtol)");
        if (conv_abs[d])
        PetscPrintf(PETSC_COMM_WORLD, " ABS (atol)");
        if (conv_stol[d])
        PetscPrintf(PETSC_COMM_WORLD, " STOL (stol)");

        PetscPrintf(PETSC_COMM_WORLD, "\n");
    }

    *reason = SNES_CONVERGED_FNORM_RELATIVE;
    }

    if (*reason == SNES_CONVERGED_ITERATING) {
        /* already iterating */
    }
    /* If we didn’t set a converged reason above, keep iterating */
    if (*reason == 0) *reason = SNES_CONVERGED_ITERATING;

    PetscFunctionReturn(0);  // Exit function safely.
}
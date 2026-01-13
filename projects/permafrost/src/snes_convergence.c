#include "snes_convergence.h"

PetscErrorCode SNESDOFConvergence(SNES snes, PetscInt it_number, PetscReal xnorm, 
  PetscReal gnorm, PetscReal fnorm, SNESConvergedReason *reason, void *cctx)
{
    PetscFunctionBegin;
    PetscErrorCode ierr;
    AppCtx *user = (AppCtx *)cctx;

    // Define vectors for residual, solution, and solution update
    Vec Res, Sol, Sol_upd;
    PetscInt i, dof=user->iga->dof;
    PetscScalar n2dof[dof], sol_n2dof[dof], sol_update_n2dof[dof];      // Norms of the residual vector components
    PetscScalar solv, solupdv;                                          // Norms of the solution and update vector

    // Retrieve the residual vector from SNES and compute its component norms.
    ierr = SNESGetFunction(snes, &Res, 0, 0);CHKERRQ(ierr);
    for (i = 0; i < dof; i++) {ierr = VecStrideNorm(Res, i, NORM_2, &n2dof[i]);CHKERRQ(ierr);}

    ierr = SNESGetSolution(snes, &Sol);CHKERRQ(ierr);
    ierr = SNESGetSolutionUpdate(snes, &Sol_upd);CHKERRQ(ierr);

    // For debugging: print the norm of the solution update vector.
    PetscReal upd_norm;
    ierr = VecNorm(Sol_upd, NORM_2, &upd_norm);CHKERRQ(ierr);

    for (i = 0; i < dof; i++) {
        ierr = VecStrideNorm(Sol, i, NORM_2, &sol_n2dof[i]);CHKERRQ(ierr);
        ierr = VecStrideNorm(Sol_upd, i, NORM_2, &sol_update_n2dof[i]);CHKERRQ(ierr);
    }

    // If temperature-dependent initial conditions are active, compute solution and update norms.
    if (user->flag_tIC == 1) {
        ierr = SNESGetSolution(snes, &Sol);CHKERRQ(ierr);
        ierr = SNESGetSolutionUpdate(snes, &Sol_upd);CHKERRQ(ierr);
        ierr = VecStrideNorm(Sol, 2, NORM_2, &solv);CHKERRQ(ierr);  // Norm of DOF 2 solution
        ierr = VecStrideNorm(Sol_upd, 2, NORM_2, &solupdv);CHKERRQ(ierr); // Norm of DOF 2 update
    }

    // Store initial residual norms at the first iteration for relative convergence checks.
    if (it_number == 0) {
        for (i = 0; i < dof; i++) {
            user->norm0[i] = n2dof[i];
            sol_n2dof[i] = sol_update_n2dof[i];
        }
        if (user->flag_tIC == 1) solupdv = solv;  // Initialize update norm to solution norm.
    }

    // Print iteration information and norm values for debugging.
    PetscPrintf(PETSC_COMM_WORLD, "    IT_NUMBER: %d ", it_number);
    PetscPrintf(PETSC_COMM_WORLD, "    fnorm: %.4e \n", fnorm);

    for (i = 0; i < dof; i++) {
        if (user->flag_tIC == 1 && i == 2) {
            PetscPrintf(PETSC_COMM_WORLD, "    x%d: %.2e s%d %.1e    ", i, n2dof[i], i, solupdv);
        } else {
            PetscPrintf(PETSC_COMM_WORLD, "    n%d: %.2e r%d %.1e    ", i, n2dof[i], i, n2dof[i] / user->norm0[i]);
        }
    }
    PetscPrintf(PETSC_COMM_WORLD, "\n");

    // Retrieve SNES solver tolerances (absolute, relative, step-size)
    PetscScalar atol, rtol, stol;
    PetscInt maxit, maxf;
    ierr = SNESGetTolerances(snes, &atol, &rtol, &stol, &maxit, &maxf);CHKERRQ(ierr);

    if (it_number < 3) {
        *reason = SNES_CONVERGED_ITERATING;
        PetscFunctionReturn(0);
    }

    // If a previous timestep required reduction, increase relative tolerance.
    if (snes->prev_dt_red == 1) {
        rtol *= 10.0;
    }

    // Convergence check based on flag_it0 setting
    if (user->flag_it0 == 1) {
        atol = 1.0e-12;  // Set absolute tolerance
    } else {
        atol = 1.0e-20;  // More strict absolute tolerance
    }

    PetscBool conv_rel[3], conv_abs[3], conv_stol[3], conv_dof[3];
    PetscInt  d;

    /* Evaluate convergence criteria per DOF */
    for (d = 0; d < 3; d++) {
        conv_rel[d]  = (n2dof[d] <= rtol * user->norm0[d]) ? PETSC_TRUE : PETSC_FALSE;
        conv_abs[d]  = (n2dof[d] < atol)                   ? PETSC_TRUE : PETSC_FALSE;
        conv_stol[d] = (sol_update_n2dof[d] <= stol * n2dof[d]) ? PETSC_TRUE : PETSC_FALSE;

        /* SAME logic as your original OR condition */
        conv_dof[d] = (conv_rel[d] || conv_abs[d] || conv_stol[d]) ? PETSC_TRUE : PETSC_FALSE;
    }

    /* SAME logic as your original AND condition */
    if (conv_dof[0] && conv_dof[1] && conv_dof[2]) {

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

    //  // Check for convergence using relative and absolute norms
    // if ((n2dof[0] <= rtol * user->norm0[0] || sol_update_n2dof[0] <= stol * sol_n2dof[0] || n2dof[0] < atol) &&
    //     (n2dof[1] <= rtol * user->norm0[1] || sol_update_n2dof[1] <= stol * sol_n2dof[1] || n2dof[1] < atol) &&
    //     (n2dof[2] <= rtol * user->norm0[2] || sol_update_n2dof[2] <= stol * sol_n2dof[2] || n2dof[2] < atol)) {
    //     *reason = SNES_CONVERGED_FNORM_RELATIVE;
    // }

    PetscFunctionReturn(0);  // Exit function safely.
}
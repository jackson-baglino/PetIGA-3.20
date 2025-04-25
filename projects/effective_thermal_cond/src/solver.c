#include "solver.h"

/*-----------------------------------------------------------
  Function: SetupAndSolve
  Assembles the system matrix and RHS vector, sets the initial
  condition for temperature, and solves the linear system using KSP.
  Returns: 0 on success, or a PETSc error code.
-----------------------------------------------------------*/
PetscErrorCode SetupAndSolve(AppCtx *user, IGA iga) {
  PetscErrorCode ierr;
  Mat A;
  Vec b;
  KSP ksp;

  PetscFunctionBegin;

  // Create system matrix and vectors
  ierr = IGACreateMat(iga, &A); CHKERRQ(ierr);
  ierr = IGACreateVec(iga, &user->T_sol); CHKERRQ(ierr);
  ierr = IGACreateVec(iga, &b); CHKERRQ(ierr);

  // Assemble system matrix and RHS vector
  ierr = IGASetFormSystem(iga, AssembleStiffnessMatrix, user); CHKERRQ(ierr);
  ierr = IGAComputeSystem(iga, A, b); CHKERRQ(ierr);

  // Set initial condition for temperature field
  ierr = ComputeInitialCondition(user->T_sol, user); CHKERRQ(ierr);

  // Solve the linear system using KSP
  ierr = IGACreateKSP(iga, &ksp); CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp, A, A); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);
  ierr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE); CHKERRQ(ierr);
  ierr = KSPSolve(ksp, b, user->T_sol); CHKERRQ(ierr);

  PetscPrintf(PETSC_COMM_WORLD, "KSP solve complete.\n\n");

  // Clean up
  ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
  ierr = MatDestroy(&A); CHKERRQ(ierr);
  ierr = VecDestroy(&b); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
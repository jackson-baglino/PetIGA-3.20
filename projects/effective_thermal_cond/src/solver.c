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
  PC pc;


  PetscFunctionBegin;

  // Create system matrix and vectors
  ierr = IGACreateMat(iga, &A); CHKERRQ(ierr);
  ierr = IGACreateVec(iga, &user->T_sol); CHKERRQ(ierr);
  ierr = IGACreateVec(iga, &b); CHKERRQ(ierr);

  // Assemble system matrix and RHS vector
  ierr = IGASetFormSystem(iga, AssembleStiffnessMatrix, user); CHKERRQ(ierr);
  ierr = IGAComputeSystem(iga, A, b); CHKERRQ(ierr);


  PetscScalar RHSsum;
  VecSum(b, &RHSsum);
  PetscPrintf(PETSC_COMM_WORLD,
            "Global RHS (sum of Neumann loads)  ΣF = %g\n", RHSsum);

  // Set initial condition for temperature field
  ierr = ComputeInitialCondition(user->T_sol, user); CHKERRQ(ierr);

 // Solve the linear system using KSP
  ierr = IGACreateKSP(iga, &ksp); CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp, A, A); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);
  // Set solver tolerances
  ierr = KSPSetTolerances(ksp, PETSC_SMALL, PETSC_SMALL, PETSC_DEFAULT, 4000); CHKERRQ(ierr); // rtol, abstol, dtol, maxits
  ierr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE); CHKERRQ(ierr);

  // Specify solver type (e.g., GMRES, CG, etc.)
  ierr = KSPSetType(ksp, KSPGMRES); CHKERRQ(ierr); // Replace KSPGMRES with desired solver type

  // Configure preconditioner
  ierr = KSPGetPC(ksp, &pc); CHKERRQ(ierr);
  ierr = PCSetType(pc, PCLU); CHKERRQ(ierr); // Use LU preconditioner

#if defined(PETSC_HAVE_MUMPS)
  if (size > 1) PetscCall(PCFactorSetMatSolverType(pc, MATSOLVERMUMPS));
#endif
  // ierr = PCFactorSetZeroPivot(pc, 1.0e-50); CHKERRQ(ierr);

  // Solve the system
  ierr = KSPSolve(ksp, b, user->T_sol); CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD, "KSP solve complete.\n\n");

  // Clean up
  ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
  ierr = MatDestroy(&A); CHKERRQ(ierr);
  ierr = VecDestroy(&b); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
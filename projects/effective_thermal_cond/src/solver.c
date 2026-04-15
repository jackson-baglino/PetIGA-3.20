#include "solver.h"
#include "assembly.h"

/*-----------------------------------------------------------------------------
  ComputeInitialCondition
  Zero the solution vector before each solve (warm start is handled by
  KSPSetInitialGuessNonzero — keeping T_sol non-zero between iterations).
-----------------------------------------------------------------------------*/
PetscErrorCode ComputeInitialCondition(Vec T, AppCtx *user)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  (void)user;
  ierr = VecSet(T, 0.0); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*-----------------------------------------------------------------------------
  CreateSolverObjects
  Allocate Mat A, Vec b, and KSP ksp once — before the per-file loop.

  KSP is configured with CG + GAMG (or HYPRE BoomerAMG when available).
  All settings can be overridden at runtime via PETSc command-line options.
-----------------------------------------------------------------------------*/
PetscErrorCode CreateSolverObjects(IGA iga, AppCtx *user,
                                   Mat *A, Vec *b, KSP *ksp)
{
  PetscErrorCode ierr;
  PC             pc;
  PetscInt       dim = user->dim;

  PetscFunctionBegin;

  /* Allocate matrix and vectors */
  ierr = IGACreateMat(iga, A); CHKERRQ(ierr);
  ierr = IGACreateVec(iga, &user->T_sol); CHKERRQ(ierr);
  ierr = IGACreateVec(iga, b); CHKERRQ(ierr);

  /* Convert BAIJ → AIJ so AMG can build its graph (PETSc 3.20 limitation) */
  {
    PetscBool isBAIJ = PETSC_FALSE;
    ierr = PetscObjectTypeCompare((PetscObject)*A, MATSEQBAIJ, &isBAIJ); CHKERRQ(ierr);
    if (!isBAIJ)
      ierr = PetscObjectTypeCompare((PetscObject)*A, MATMPIBAIJ, &isBAIJ); CHKERRQ(ierr);
    if (isBAIJ) {
      ierr = MatConvert(*A, MATAIJ, MAT_INPLACE_MATRIX, A); CHKERRQ(ierr);
    }
  }

  /* Mark symmetric so CG is valid */
  ierr = MatSetOption(*A, MAT_SYMMETRIC, PETSC_TRUE); CHKERRQ(ierr);

  /* Create and configure KSP */
  ierr = IGACreateKSP(iga, ksp); CHKERRQ(ierr);
  ierr = KSPSetType(*ksp, KSPCG); CHKERRQ(ierr);
  ierr = KSPGetPC(*ksp, &pc); CHKERRQ(ierr);

#if defined(PETSC_HAVE_HYPRE)
  ierr = PCSetType(pc, PCHYPRE); CHKERRQ(ierr);
  ierr = PCHYPRESetType(pc, "boomeramg"); CHKERRQ(ierr);
#else
  ierr = PCSetType(pc, PCGAMG); CHKERRQ(ierr);
#endif

  ierr = KSPSetTolerances(*ksp, 1e-9, PETSC_SMALL, PETSC_DEFAULT, 200); CHKERRQ(ierr);
  ierr = KSPSetInitialGuessNonzero(*ksp, PETSC_TRUE); CHKERRQ(ierr);

  /* Allow full runtime override */
  ierr = KSPSetFromOptions(*ksp); CHKERRQ(ierr);

  (void)dim; /* used implicitly through iga DOF count */
  PetscFunctionReturn(0);
}

/*-----------------------------------------------------------------------------
  Solve
  For the current ice field in user->ice:
    1. Zero and reassemble A and b.
    2. Build and attach constant-mode null space (one vector per DOF component).
    3. Reset KSP operators and solve.
    4. Enforce zero mean per component (parallel-safe via MPI_Allreduce).
-----------------------------------------------------------------------------*/
PetscErrorCode Solve(AppCtx *user, IGA iga, Mat A, Vec b, KSP ksp)
{
  PetscErrorCode ierr;
  PetscInt       dim = user->dim;

  PetscFunctionBegin;

  /* ---- Reassemble ---- */
  ierr = MatZeroEntries(A); CHKERRQ(ierr);
  ierr = VecZeroEntries(b); CHKERRQ(ierr);
  ierr = IGASetFormSystem(iga, AssembleStiffnessMatrix, user); CHKERRQ(ierr);
  ierr = IGAComputeSystem(iga, A, b); CHKERRQ(ierr);

  /* ---- Constant-mode null space (one vector per component) ----
   *
   * PARALLEL-SAFE: use VecGetOwnershipRange to access only the LOCAL portion
   * of the basis vectors.  Each rank sets the entries it owns.
   */
  {
    PetscInt  gsize, nnodes;
    Vec       basis[3];
    MatNullSpace nsp;

    ierr = VecGetSize(user->T_sol, &gsize); CHKERRQ(ierr);
    nnodes = gsize / dim;

    for (PetscInt m = 0; m < dim; ++m) {
      ierr = IGACreateVec(iga, &basis[m]); CHKERRQ(ierr);
      ierr = VecSet(basis[m], 0.0); CHKERRQ(ierr);

      /* Only touch the local portion */
      PetscInt     lo, hi;
      PetscScalar *y;
      ierr = VecGetOwnershipRange(basis[m], &lo, &hi); CHKERRQ(ierr);
      ierr = VecGetArray(basis[m], &y); CHKERRQ(ierr);
      for (PetscInt idx = 0; idx < hi - lo; ++idx) {
        if ((lo + idx) % dim == m) y[idx] = 1.0;
      }
      ierr = VecRestoreArray(basis[m], &y); CHKERRQ(ierr);
      ierr = VecAssemblyBegin(basis[m]); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(basis[m]); CHKERRQ(ierr);
      ierr = VecNormalize(basis[m], NULL); CHKERRQ(ierr);
    }

    ierr = MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_FALSE, dim, basis, &nsp); CHKERRQ(ierr);
    ierr = MatSetNullSpace(A, nsp); CHKERRQ(ierr);
    ierr = MatNullSpaceDestroy(&nsp); CHKERRQ(ierr);
    for (PetscInt m = 0; m < dim; ++m) {
      ierr = VecDestroy(&basis[m]); CHKERRQ(ierr);
    }
    (void)nnodes;
  }

  /* ---- Initial condition for T_sol (only on first call; later reuse) ---- */
  {
    PetscReal norm;
    ierr = VecNorm(user->T_sol, NORM_2, &norm); CHKERRQ(ierr);
    if (norm == 0.0) {
      ierr = ComputeInitialCondition(user->T_sol, user); CHKERRQ(ierr);
    }
  }

  /* ---- Solve ---- */
  ierr = KSPSetOperators(ksp, A, A); CHKERRQ(ierr);
  ierr = KSPSolve(ksp, b, user->T_sol); CHKERRQ(ierr);

  /* ---- Enforce zero mean per component (PARALLEL-SAFE) ----
   *
   * Accumulate local partial sums then reduce globally with MPI_Allreduce.
   * Only the LOCAL portion of the vector is accessed via VecGetArray.
   */
  {
    PetscInt         gsize, lo, hi, local_n;
    const PetscScalar *x;
    PetscScalar      *y;
    PetscReal         local_sum[3] = {0.0, 0.0, 0.0};
    PetscReal         global_sum[3];
    PetscInt          nnodes;

    ierr = VecGetSize(user->T_sol, &gsize); CHKERRQ(ierr);
    nnodes = gsize / dim;

    ierr = VecGetOwnershipRange(user->T_sol, &lo, &hi); CHKERRQ(ierr);
    local_n = hi - lo;

    ierr = VecGetArrayRead(user->T_sol, &x); CHKERRQ(ierr);
    for (PetscInt k = 0; k < local_n; ++k)
      local_sum[(lo + k) % dim] += PetscRealPart(x[k]);
    ierr = VecRestoreArrayRead(user->T_sol, &x); CHKERRQ(ierr);

    ierr = MPI_Allreduce(local_sum, global_sum, dim,
                         MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD);
    CHKERRQ(ierr);

    ierr = VecGetArray(user->T_sol, &y); CHKERRQ(ierr);
    for (PetscInt k = 0; k < local_n; ++k)
      y[k] -= global_sum[(lo + k) % dim] / (PetscReal)nnodes;
    ierr = VecRestoreArray(user->T_sol, &y); CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

/*-----------------------------------------------------------------------------
  DestroySolverObjects
  Free the matrix, RHS vector, and KSP after the file loop.
-----------------------------------------------------------------------------*/
PetscErrorCode DestroySolverObjects(Mat *A, Vec *b, KSP *ksp)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = KSPDestroy(ksp); CHKERRQ(ierr);
  ierr = MatDestroy(A);   CHKERRQ(ierr);
  ierr = VecDestroy(b);   CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

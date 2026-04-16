#include "solver.h"
#include "assembly.h"

/*-----------------------------------------------------------------------------
  ComputeInitialCondition
  Zero the solution vector (used on the very first call only).
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

  Default solver: PREONLY + LU (direct).  For parallel runs with more than
  one MPI process, use SuperLU_DIST automatically if available.

  The periodic homogenization system has one constant-mode null space per DOF
  component.  We handle this by pinning one node per component with
  MatZeroRowsColumns (see Solve), so the matrix is non-singular and direct LU
  works cleanly.  Iterative solvers (CG, GMRES) are unreliable for this class
  of problem and should be avoided.

  All settings can be overridden at runtime via PETSc command-line options,
  e.g.  -ksp_type preonly -pc_type lu
        -ksp_type preonly -pc_type lu -pc_factor_mat_solver_type superlu_dist
-----------------------------------------------------------------------------*/
PetscErrorCode CreateSolverObjects(IGA iga, AppCtx *user,
                                   Mat *A, Vec *b, KSP *ksp)
{
  PetscErrorCode ierr;
  PC             pc;
  PetscMPIInt    size;

  PetscFunctionBegin;

  /* Allocate matrix and vectors */
  ierr = IGACreateMat(iga, A); CHKERRQ(ierr);
  ierr = IGACreateVec(iga, &user->T_sol); CHKERRQ(ierr);
  ierr = IGACreateVec(iga, b); CHKERRQ(ierr);

  /* Create KSP and configure for direct solve */
  ierr = IGACreateKSP(iga, ksp); CHKERRQ(ierr);
  ierr = KSPSetType(*ksp, KSPPREONLY); CHKERRQ(ierr);
  ierr = KSPGetPC(*ksp, &pc); CHKERRQ(ierr);
  ierr = PCSetType(pc, PCLU); CHKERRQ(ierr);

  /* For parallel runs use SuperLU_DIST (available in this PETSc build) */
  MPI_Comm_size(PETSC_COMM_WORLD, &size);
  if (size > 1) {
    ierr = PCFactorSetMatSolverType(pc, MATSOLVERSUPERLU_DIST); CHKERRQ(ierr);
  }

  /* Allow full runtime override */
  ierr = KSPSetFromOptions(*ksp); CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

/*-----------------------------------------------------------------------------
  Solve
  For the current ice field in user->ice:
    1. Zero and reassemble A and b.
    2. Pin one DOF per component at node 0 (MatZeroRowsColumns) to eliminate
       the constant-mode null space and make the system non-singular.
    3. Solve with the direct KSP.
    4. Enforce zero mean per component (parallel-safe via MPI_Allreduce).

  Why DOF pinning instead of MatNullSpace?
  ─────────────────────────────────────────
  MatNullSpace works for iterative solvers but causes iterative methods (CG,
  GMRES) to converge to a numerically polluted solution for this problem class.
  DOF pinning makes the system uniquely solvable by any solver, including
  direct LU.  After pinning, t_m(node 0) = 0 for each component m; the
  zero-mean step below then centres the solution for good practice.
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

  /* ---- Pin one DOF per component to remove the constant-mode null space ----
   *
   * Global DOF layout (PetIGA interleaved): [c0_n0, c1_n0, c0_n1, c1_n1, ...]
   * Component m at node 0 has global index m.
   *
   * MatZeroRowsColumns zeros row m and column m, puts diag=1 on the diagonal,
   * and adjusts b.  Since T_sol = 0 at these DOFs, the adjustment to b is
   * zero and we obtain b[m] = 0 (homogeneous Dirichlet for the pinned DOF).
   *
   * This operation is collective — all ranks call it with the same global
   * row indices.
   */
  {
    PetscInt rows[3];  /* at most dim = 3 pinned DOFs */
    for (PetscInt m = 0; m < dim; m++) rows[m] = m;

    /* Zero T_sol at the pinned DOFs so the RHS adjustment is trivially zero */
    ierr = VecSet(user->T_sol, 0.0); CHKERRQ(ierr);

    ierr = MatZeroRowsColumns(A, dim, rows, 1.0, user->T_sol, b); CHKERRQ(ierr);
  }

  /* ---- Solve ---- */
  ierr = KSPSetOperators(ksp, A, A); CHKERRQ(ierr);
  ierr = KSPSolve(ksp, b, user->T_sol); CHKERRQ(ierr);

  /* ---- Enforce zero mean per component (PARALLEL-SAFE) ----
   *
   * The pinned-node solution has t_m(node 0) = 0 but a non-zero mean in
   * general.  Subtracting the mean is a gauge choice that does not affect
   * k_eff (which depends only on grad t), but it is conventional.
   *
   * Accumulate local partial sums then reduce globally with MPI_Allreduce.
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

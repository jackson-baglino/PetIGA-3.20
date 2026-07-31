#include "keff.h"
#include "material_properties.h"
#include "assembly.h"   /* Integration(), for the -keff_debug_phibar cross-check */

/* ---------------------------------------------------------------------------
 * keff.c -- lifecycle, options and guards for the in-line k_eff diagnostic.
 * See include/keff.h for the method and for why the corrector mesh is scalar.
 * ------------------------------------------------------------------------- */

/* Integrand returning 1, so IGAComputeScalar accumulates INT_Y dV. */
PetscErrorCode KeffVolumeIntegrand(IGAPoint point, const PetscScalar U[],
                                   PetscInt n, PetscScalar S[], void *ctx)
{
  PetscFunctionBegin;
  (void)U; (void)n; (void)ctx;
  if (!point->atboundary) S[0] += 1.0;
  PetscFunctionReturn(0);
}

/* ---------------------------------------------------------------------------
 * KeffParseOptions
 *
 * A self-contained options block rather than 40 more lines bolted onto the
 * ~360-line block in main(). -keff is read first so that a run without it pays
 * nothing.
 * ------------------------------------------------------------------------- */
static PetscErrorCode KeffParseOptions(KeffCtx *kc)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;

  PetscOptionsBegin(PETSC_COMM_WORLD, "",
                    "Effective thermal conductivity (periodic homogenization)", "");

  ierr = PetscOptionsBool("-keff",
      "Compute the effective thermal conductivity tensor in-line", "",
      kc->enabled, &kc->enabled, NULL); CHKERRQ(ierr);

  ierr = PetscOptionsInt("-keff_freq",
      "Sample k_eff every N accepted time steps", "",
      kc->freq, &kc->freq, NULL); CHKERRQ(ierr);

  ierr = PetscOptionsReal("-keff_t_interv",
      "Sample k_eff every this many simulated seconds (> 0 overrides -keff_freq)", "",
      kc->t_interv, &kc->t_interv, NULL); CHKERRQ(ierr);

  ierr = PetscOptionsBool("-keff_step0",
      "Always take a k_eff sample at step 0", "",
      kc->at_step0, &kc->at_step0, NULL); CHKERRQ(ierr);

  ierr = PetscOptionsBool("-keff_only",
      "Take one k_eff sample from the initial condition and exit (validation driver)", "",
      kc->only, &kc->only, NULL); CHKERRQ(ierr);

  ierr = PetscOptionsString("-keff_replay",
      "Compute k_eff for every sol_*.dat in this finished run directory, then exit", "",
      kc->replay_dir, kc->replay_dir, sizeof(kc->replay_dir), NULL); CHKERRQ(ierr);

  ierr = PetscOptionsString("-keff_replay_times",
      "Step->time map for replay (defaults to <replay dir>/SSA_evo.dat)", "",
      kc->replay_times, kc->replay_times, sizeof(kc->replay_times), NULL); CHKERRQ(ierr);

  ierr = PetscOptionsInt("-keff_replay_stride",
      "Process every Nth snapshot in replay mode", "",
      kc->replay_stride, &kc->replay_stride, NULL); CHKERRQ(ierr);

  ierr = PetscOptionsString("-keff_csv",
      "Output CSV path (default $folder/k_eff.csv)", "",
      kc->csv_path, kc->csv_path, sizeof(kc->csv_path), NULL); CHKERRQ(ierr);

  ierr = PetscOptionsBool("-keff_write_corrector",
      "Also write the corrector fields t_vec_%05d.dat and igakeff.dat", "",
      kc->write_corrector, &kc->write_corrector, NULL); CHKERRQ(ierr);

  ierr = PetscOptionsBool("-keff_pc_freeze",
      "Hold the corrector preconditioner across samples instead of rebuilding "
      "its coarse operators each time", "",
      kc->pc_freeze, &kc->pc_freeze, NULL); CHKERRQ(ierr);

  ierr = PetscOptionsInt("-keff_pc_refresh",
      "When frozen, force a full preconditioner rebuild every N samples", "",
      kc->pc_refresh, &kc->pc_refresh, NULL); CHKERRQ(ierr);

  ierr = PetscOptionsInt("-keff_max_its",
      "Corrector iterations above this force a preconditioner rebuild and one retry", "",
      kc->max_its, &kc->max_its, NULL); CHKERRQ(ierr);

  ierr = PetscOptionsBool("-keff_debug_phibar",
      "Cross-check the phi projection: mean ice fraction on the cloned mesh vs "
      "the solver mesh", "",
      kc->debug_phibar, &kc->debug_phibar, NULL); CHKERRQ(ierr);

  PetscOptionsEnd();

  PetscFunctionReturn(0);
}

/* ---------------------------------------------------------------------------
 * KeffCheckGuards
 *
 * Every one of these is a case where the homogenization would still produce a
 * number, and that number would be meaningless. Fail loudly instead.
 * ------------------------------------------------------------------------- */
static PetscErrorCode KeffCheckGuards(KeffCtx *kc)
{
  PetscErrorCode ierr;
  AppCtx        *app = kc->app;

  PetscFunctionBegin;

  if (kc->dim < 2)
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_INCOMP,
            "\n*** -keff requires -dim 2 or 3 (got %d) ***\n"
            "  The k_eff tensor is degenerate in 1D: with a single direction the\n"
            "  cell problem reduces to the harmonic mean, which is analytic and\n"
            "  needs no solve.\n", (int)kc->dim);

  /* A mapped (igakit) geometry does not tile space, so there is no cell to
   * upscale; its axes are not periodic either. Checked before periodicity so
   * the more specific message wins. */
  if (app->iga->geometry)
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_INCOMP,
            "\n*** -keff is incompatible with -geom_file ***\n"
            "  Periodic homogenization upscales a cell that TILES SPACE. An\n"
            "  igakit-mapped geometry (sediment bumps, axisymmetric r-z, ...)\n"
            "  does not tile, and the flux average would be divided by a cell\n"
            "  volume the deformed patch does not have.\n"
            "  Fix: use a uniform Cartesian mesh (-Nx/-Ny/-Lx/-Ly), or drop -keff.\n");

  if (app->axisym)
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_INCOMP,
            "\n*** -keff is incompatible with -axisym ***\n"
            "  The axisymmetric measure 2*pi*r dr dz is not a periodic-cell\n"
            "  measure, and the cell-problem assembly does not carry the r\n"
            "  weight that assembly.c applies to the phase-field residual.\n");

  /* Interrogate the AXES, not user->periodic: -geom_file bypasses the
   * IGAAxisSetPeriodic calls in main() entirely, so the flag does not imply a
   * periodic mesh. (The -geom_file guard above already catches that case, but
   * this stays correct if that guard is ever relaxed.) */
  for (PetscInt d = 0; d < kc->dim; d++) {
    IGAAxis   ax;
    PetscBool per;
    ierr = IGAGetAxis(app->iga, d, &ax); CHKERRQ(ierr);
    ierr = IGAAxisGetPeriodic(ax, &per); CHKERRQ(ierr);
    if (!per)
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_INCOMP,
              "\n*** -keff REQUIRES A FULLY PERIODIC DOMAIN ***\n"
              "  Axis %d is not periodic (run has -periodic %d).\n"
              "  This is periodic homogenization: the corrector t_m being\n"
              "  periodic on the cell boundary IS the definition of the upscaled\n"
              "  tensor. On a non-periodic window the boundary layer contaminates\n"
              "  the cell average, and the result is not an effective property of\n"
              "  anything -- it describes that particular box with those\n"
              "  particular walls.\n"
              "  Fix: run with -periodic 1, or drop -keff.\n",
              (int)d, (int)app->periodic);
  }

  PetscFunctionReturn(0);
}

/* ---------------------------------------------------------------------------
 * KeffCreate
 *
 * Called from main() right after the alph/mob allocation -- the first point at
 * which IGASetUp has run, -geom_file has corrected Nx/Ny/Nz, and the local
 * Gauss-point count is known. It is also early, so the large corrector matrix
 * allocation fails fast rather than six hours into a run.
 * ------------------------------------------------------------------------- */
PetscErrorCode KeffCreate(AppCtx *app)
{
  PetscErrorCode ierr;
  KeffCtx       *kc;

  PetscFunctionBegin;

  app->keff = NULL;

  ierr = PetscNew(&kc); CHKERRQ(ierr);
  kc->app           = app;
  kc->enabled       = PETSC_FALSE;
  kc->freq          = 1;
  kc->t_interv      = 0.0;
  kc->at_step0      = PETSC_TRUE;
  kc->only          = PETSC_FALSE;
  kc->replay_stride = 1;
  kc->pc_freeze     = PETSC_FALSE;
  kc->pc_refresh    = 20;
  kc->max_its       = 500;

  ierr = KeffParseOptions(kc); CHKERRQ(ierr);

  if (!kc->enabled) { ierr = PetscFree(kc); CHKERRQ(ierr); PetscFunctionReturn(0); }

  kc->dim = app->dim;
  ierr = KeffCheckGuards(kc); CHKERRQ(ierr);

  /* Scalar corrector mesh. IGAClone copies the axes INCLUDING periodicity, the
   * quadrature rules and the element partition, and deliberately does not copy
   * iga->form -- so the phase-field Dirichlet values for T and rho_v cannot
   * leak into the cell problem. */
  ierr = IGAClone(app->iga, 1, &kc->iga); CHKERRQ(ierr);
  ierr = IGASetMatType(kc->iga, MATAIJ); CHKERRQ(ierr);

  /* The cross-IGA contract: the clone must partition elements exactly as the
   * parent does, because kc->ice[] is indexed by the parent's Gauss-point
   * numbering and read back on the clone. */
  for (PetscInt d = 0; d < kc->dim; d++) {
    if (kc->iga->elem_width[d] != app->iga->elem_width[d] ||
        kc->iga->elem_start[d] != app->iga->elem_start[d])
      SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_PLIB,
              "k_eff: cloned IGA partitions axis %d differently from the solver "
              "IGA (start %d/%d, width %d/%d). The Gauss-point index mapping "
              "between the two would be invalid.", (int)d,
              (int)kc->iga->elem_start[d], (int)app->iga->elem_start[d],
              (int)kc->iga->elem_width[d], (int)app->iga->elem_width[d]);
  }

  ierr = IGACreateMat(kc->iga, &kc->A); CHKERRQ(ierr);

  /* Hand the matrix back its stock AIJ behaviour for CreateVecs and Duplicate.
   *
   * IGACreateMat composes an "IGA" object onto the matrix and overrides
   * MATOP_CREATE_VECS and MATOP_DUPLICATE with IGA-aware versions
   * (petigamat.c:385-391) that hard-error with "Matrix not generated from IGA"
   * if that composed object is absent. Algebraic multigrid builds a hierarchy of
   * coarse operators from this one, and those coarse matrices inherit the
   * overridden operations without inheriting the composed IGA -- so PCGAMG dies
   * during setup, on any mesh, every time.
   *
   * Nothing here needs the IGA-aware versions: the vectors are created directly
   * with IGACreateVec, and the matrix is only ever assembled by
   * IGAComputeMatrix (which uses MatSetValues) and consumed by KSP. Restoring
   * the defaults costs nothing and is what makes -keff_pc_type gamg usable.
   *
   * This is the whole of the "iterative solvers are unreliable for this problem"
   * story in the standalone project's README: not a numerical property of the
   * cell problem at all, but a PetIGA/GAMG plumbing incompatibility.
   *
   * MATOP_CREATE_VECS can simply be cleared -- MatCreateVecs falls back to a
   * layout-based default when the slot is empty. MATOP_DUPLICATE cannot:
   * clearing it yields "No method duplicate for Mat of type seqaij". PetIGA
   * stashes the original under "__IGA_MatDuplicate" (petigamat.c:390), so the
   * stock implementation is put back from there. */
  ierr = MatSetOperation(kc->A, MATOP_CREATE_VECS, NULL); CHKERRQ(ierr);
  {
    PetscErrorCode (*matduplicate)(Mat, MatDuplicateOption, Mat *) = NULL;
    ierr = PetscObjectQueryFunction((PetscObject)kc->A, "__IGA_MatDuplicate",
                                    &matduplicate); CHKERRQ(ierr);
    if (matduplicate) {
      ierr = MatSetOperation(kc->A, MATOP_DUPLICATE,
                             (PetscVoidFunction)matduplicate); CHKERRQ(ierr);
    }
  }
  /* Keep the sparsity pattern across MatZeroRowsColumns. Without this the
   * zeroed off-diagonals are REMOVED from the AIJ structure, which bumps
   * A->nonzerostate; PCSetUp then sees DIFFERENT_NONZERO_PATTERN and redoes
   * the full symbolic setup every sample, making AMG reuse a no-op. */
  ierr = MatSetOption(kc->A, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE); CHKERRQ(ierr);

  for (PetscInt m = 0; m < kc->dim; m++) {
    ierr = IGACreateVec(kc->iga, &kc->b[m]); CHKERRQ(ierr);
    ierr = IGACreateVec(kc->iga, &kc->T[m]); CHKERRQ(ierr);
    ierr = VecZeroEntries(kc->T[m]); CHKERRQ(ierr);
  }

  ierr = PetscMalloc1(app->ngp, &kc->ice); CHKERRQ(ierr);
  ierr = PetscMemzero(kc->ice, sizeof(PetscReal) * app->ngp); CHKERRQ(ierr);

  /* app->keff must be visible before KeffCheckGaussLayout, which reaches back
   * through it. Set here rather than at the end of the function. */
  app->keff = kc;

  /* Precondition for the flat Gauss-point index formula (see keff_field.c). */
  ierr = KeffCheckGaussLayout(app); CHKERRQ(ierr);

  ierr = KeffCreateSolver(app); CHKERRQ(ierr);

  /* Measure the cell volume by quadrature instead of assuming Lx*Ly[*Lz].
   * Costs one element loop at startup and is a genuine self-test: a mis-cloned
   * axis, a wrong quadrature rule or a partitioning bug shows up here as a
   * volume mismatch BEFORE any k_eff number is produced. T[0] is passed (rather
   * than NULL) only so IGAElementGetValues has a real array to read. */
  {
    PetscScalar vol = 0.0;
    ierr = IGAComputeScalar(kc->iga, kc->T[0], 1, &vol,
                            KeffVolumeIntegrand, kc); CHKERRQ(ierr);
    kc->vol = PetscRealPart(vol);
  }
  if (kc->vol <= 0.0)
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_PLIB,
            "k_eff: measured cell volume is %g, which is not positive.", (double)kc->vol);

  /* Default CSV path: alongside every other per-run output. monitoring.c uses
   * getenv("folder") for the same reason -- -output_path is parsed but never
   * read by anything that writes. */
  if (kc->csv_path[0] == '\0') {
    const char *dir = getenv("folder");
    if (dir) ierr = PetscSNPrintf(kc->csv_path, sizeof(kc->csv_path), "%s/k_eff.csv", dir);
    else     ierr = PetscSNPrintf(kc->csv_path, sizeof(kc->csv_path), "k_eff.csv");
    CHKERRQ(ierr);
  }
  if (kc->replay_dir[0] != '\0' && kc->replay_times[0] == '\0') {
    ierr = PetscSNPrintf(kc->replay_times, sizeof(kc->replay_times),
                         "%s/SSA_evo.dat", kc->replay_dir); CHKERRQ(ierr);
  }

  /* Summary. Print the analytic box next to the measured volume so a mismatch
   * is visible at run start rather than inferred from a wrong k_eff later. */
  {
    const PetscReal box = (kc->dim == 2) ? app->Lx * app->Ly
                                         : app->Lx * app->Ly * app->Lz;
    const char     *unit = (kc->dim == 2) ? "m^2" : "m^3";
    PetscInt        nrows;

    ierr = MatGetSize(kc->A, &nrows, NULL); CHKERRQ(ierr);

    ierr = PetscPrintf(PETSC_COMM_WORLD,
      "\n"
      "  ===============================================================================\n"
      "  >>> EFFECTIVE THERMAL CONDUCTIVITY (periodic homogenization)\n"
      "  ===============================================================================\n"
      "    corrector mesh    : dim %d, scalar (dof 1), %d unknowns x %d directions\n"
      "    cell volume (quad): %.12e %s\n"
      "    cell volume (box) : %.12e %s   [rel diff %.2e]\n"
      "    k_ice / k_air     : %g / %g W/m/K\n",
      (int)kc->dim, (int)nrows, (int)kc->dim,
      (double)kc->vol, unit, (double)box, unit,
      (double)(PetscAbsReal(kc->vol - box) / box),
      (double)app->thcond_ice, (double)app->thcond_air); CHKERRQ(ierr);

    if (kc->replay_dir[0] != '\0') {
      ierr = PetscPrintf(PETSC_COMM_WORLD,
        "    mode              : REPLAY of %s (stride %d)\n",
        kc->replay_dir, (int)kc->replay_stride); CHKERRQ(ierr);
    } else if (kc->only) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,
        "    mode              : -keff_only (one sample from the IC, then exit)\n"); CHKERRQ(ierr);
    } else if (kc->t_interv > 0.0) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,
        "    cadence           : every %g s of simulated time\n",
        (double)kc->t_interv); CHKERRQ(ierr);
    } else {
      ierr = PetscPrintf(PETSC_COMM_WORLD,
        "    cadence           : every %d accepted step(s)\n", (int)kc->freq); CHKERRQ(ierr);
    }

    ierr = PetscPrintf(PETSC_COMM_WORLD,
      "    output            : %s\n"
      "  ===============================================================================\n\n",
      kc->csv_path); CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

/* ---------------------------------------------------------------------------
 * KeffDebugPhiBar  (-keff_debug_phibar)
 *
 * The one assumption in this port with no prior test coverage anywhere: that
 * the SOLVER IGA and the CLONED corrector IGA number their local Gauss points
 * identically, so phi written by KeffProjectIce while walking the former is
 * read at the same physical location while assembling on the latter.
 *
 * If that fails, nothing crashes -- you get a plausible conductivity tensor for
 * a scrambled microstructure. So compare one scalar that both meshes can
 * compute independently: the mean ice fraction.
 *
 *   clone : (1/|Y|) INT phi dV, with phi read back out of kc->ice[]
 *   solver: TOT_ICE/|Y| from IGAComputeScalar(app->iga, U, ..., Integration),
 *           the same call Monitor makes (monitoring.c), which never touches
 *           kc->ice[] and shares no code with the projection
 *
 * A single mismatched point moves the mean, so agreement to round-off is strong
 * evidence the mapping is right.
 * ------------------------------------------------------------------------- */
PetscErrorCode KeffDebugPhiBar(AppCtx *app, Vec U)
{
  PetscErrorCode  ierr;
  KeffCtx        *kc = app->keff;
  PetscReal       mc[7] = {0,0,0,0,0,0,0};  /* from the clone, via kc->ice[] */
  PetscReal       ms[7] = {0,0,0,0,0,0,0};  /* from the solver mesh, via U   */
  PetscReal       worst = 0.0;
  PetscInt        worst_i = 0;
  const PetscReal tol = 1.0e-10;
  const char     *label2[7] = {"mean ice fraction",
                               "ice centroid x [m]", "ice centroid y [m]",
                               "ice <x^2> [m^2]",    "ice <y^2> [m^2]",
                               "", ""};
  const char     *label3[7] = {"mean ice fraction",
                               "ice centroid x [m]", "ice centroid y [m]",
                               "ice centroid z [m]",
                               "ice <x^2> [m^2]",    "ice <y^2> [m^2]",
                               "ice <z^2> [m^2]"};
  const char    **label = (kc && kc->dim == 3) ? label3 : label2;

  PetscFunctionBegin;
  if (!kc || !kc->debug_phibar) PetscFunctionReturn(0);

  ierr = KeffProjectIce(app, U); CHKERRQ(ierr);
  ierr = KeffIceMoments(app, PETSC_TRUE,  U, mc); CHKERRQ(ierr);
  ierr = KeffIceMoments(app, PETSC_FALSE, U, ms); CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD,
    "  [keff] Gauss-point projection check (clone vs solver mesh)\n"
    "         %-20s %22s %22s %11s\n",
    "quantity", "clone (kc->ice[])", "solver (U)", "rel diff"); CHKERRQ(ierr);

  for (PetscInt i = 0; i < 1 + 2 * kc->dim; i++) {
    /* Scale each residual by a fixed quantity of the right dimension, not by
     * the value itself: a centroid or moment legitimately near zero would
     * otherwise blow the ratio up. */
    const PetscReal scale = (i == 0)            ? PetscMax(PetscAbsReal(ms[0]), 1.0e-300)
                          : (i <= kc->dim)      ? app->Lx
                                                : app->Lx * app->Lx;
    const PetscReal rd    = PetscAbsReal(mc[i] - ms[i]) / scale;
    if (rd > worst) { worst = rd; worst_i = i; }
    ierr = PetscPrintf(PETSC_COMM_WORLD,
      "         %-20s %22.15e %22.15e %11.3e\n",
      label[i], (double)mc[i], (double)ms[i], (double)rd); CHKERRQ(ierr);
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD,
    "         worst = %.3e on '%s'  (tol %.1e)\n",
    (double)worst, label[worst_i], (double)tol); CHKERRQ(ierr);

  if (worst > tol)
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_PLIB,
            "k_eff: the cloned corrector mesh and the solver mesh disagree on "
            "'%s' by %.3e (tol %.1e). The Gauss-point index mapping between "
            "them is wrong, so every k(x) in the cell problem would be taken "
            "from the wrong location. See src/keff_field.c.",
            label[worst_i], (double)worst, (double)tol);

  PetscFunctionReturn(0);
}

/* ---------------------------------------------------------------------------
 * KeffDestroy
 *
 * Owns everything KeffCreate made. The standalone version split this ownership
 * (CreateSolverObjects allocated T_sol but main() destroyed it), which is an
 * easy leak once embedded.
 * ------------------------------------------------------------------------- */
PetscErrorCode KeffDestroy(AppCtx *app)
{
  PetscErrorCode ierr;
  KeffCtx       *kc = app->keff;

  PetscFunctionBegin;
  if (!kc) PetscFunctionReturn(0);

  ierr = KeffDestroySolver(app); CHKERRQ(ierr);
  if (kc->csv) { ierr = PetscViewerDestroy(&kc->csv); CHKERRQ(ierr); }
  for (PetscInt m = 0; m < kc->dim; m++) {
    if (kc->b[m]) { ierr = VecDestroy(&kc->b[m]); CHKERRQ(ierr); }
    if (kc->T[m]) { ierr = VecDestroy(&kc->T[m]); CHKERRQ(ierr); }
  }
  if (kc->A)   { ierr = MatDestroy(&kc->A);   CHKERRQ(ierr); }
  if (kc->iga) { ierr = IGADestroy(&kc->iga); CHKERRQ(ierr); }
  if (kc->ice) { ierr = PetscFree(kc->ice);   CHKERRQ(ierr); }

  ierr = PetscFree(kc); CHKERRQ(ierr);
  app->keff = NULL;

  PetscFunctionReturn(0);
}

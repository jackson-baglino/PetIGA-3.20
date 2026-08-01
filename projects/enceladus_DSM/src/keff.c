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

  /* The SOLVER mesh no longer has to be periodic -- the corrector builds its own
   * periodic mesh (see KeffCreate). But a warning is still owed when it is not,
   * because periodic homogenization computes the effective property of the
   * medium obtained by TILING the cell, and that presumes k(x) is periodic.
   *
   * Under -periodic 0 two things break that presumption, and neither is fixed by
   * giving the corrector its own mesh:
   *   - the packings are generated periodic, and the multi_grains IC applies
   *     minimum-image wrapping ONLY under -periodic 1 (initial_conditions.c),
   *     so grains that straddle the boundary are CUT rather than wrapped,
   *     leaving flat artificial faces at the walls from t = 0;
   *   - a zero-flux wall is a MIRROR, not a wrap, so the field evolves
   *     mirror-symmetric rather than periodic and the corrector's periodic
   *     condition sees a material mismatch across the seam.
   * The size of the resulting error scales with how much of the domain lies
   * within a grain radius of a wall.
   */
  {
    PetscInt nonper = 0;
    for (PetscInt d = 0; d < kc->dim; d++) {
      IGAAxis   ax;
      PetscBool per;
      ierr = IGAGetAxis(app->iga, d, &ax); CHKERRQ(ierr);
      ierr = IGAAxisGetPeriodic(ax, &per); CHKERRQ(ierr);
      if (!per) nonper++;
    }
    if (nonper) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,
        "\n"
        "  *** -keff WARNING: the solver mesh is NOT periodic (%d of %d axes) ***\n"
        "  The corrector still solves on its own periodic mesh, so k_eff is\n"
        "  computed -- but periodic homogenization upscales the medium obtained\n"
        "  by TILING this cell, and a non-periodic run does not produce one:\n"
        "    - grains straddling the boundary are cut, not wrapped, because the\n"
        "      multi_grains IC only applies minimum-image under -periodic 1;\n"
        "    - a zero-flux wall is a mirror, not a wrap, so the evolved field is\n"
        "      mirror-symmetric and k(x) does not match across the seam.\n"
        "  Treat the result as indicative. For a defensible tensor either run\n"
        "  with -periodic 1, or generate a packing whose grains lie fully inside\n"
        "  the box so nothing is cut.\n\n", (int)nonper, (int)kc->dim); CHKERRQ(ierr);
    }
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

  /* Scalar corrector mesh, ALWAYS PERIODIC regardless of the solver's own
   * boundary conditions.
   *
   * The two problems want different boundary conditions and there is no reason
   * they should share them. The phase field, temperature and vapour equations
   * are perfectly well posed with zero flux on every wall (-periodic 0); it is
   * only the CELL PROBLEM that needs periodicity, because a periodic corrector
   * is what makes the flux average an effective property rather than a property
   * of one box. PetIGA carries periodicity on the AXIS, i.e. on the knot vector
   * and DOF layout, so it applies to every field on an IGA at once -- which is
   * why the corrector needs its own mesh rather than a flag.
   *
   * Built from the solver IGA's element counts and extents rather than cloned,
   * so the periodicity can differ. The element PARTITION must still match
   * exactly, because kc->ice[] is filled walking the solver IGA and read back
   * walking this one; that is asserted immediately below and is the reason this
   * is safe. */
  {
    IGAAxis ax_src, ax_dst;
    ierr = IGACreate(PetscObjectComm((PetscObject)app->iga), &kc->iga); CHKERRQ(ierr);
    ierr = IGASetDim(kc->iga, kc->dim); CHKERRQ(ierr);
    ierr = IGASetDof(kc->iga, 1); CHKERRQ(ierr);
    for (PetscInt d = 0; d < kc->dim; d++) {
      PetscInt  p_d, N_d;
      PetscReal U0, U1;
      const PetscReal *Uknot;
      ierr = IGAGetAxis(app->iga, d, &ax_src); CHKERRQ(ierr);
      ierr = IGAAxisGetDegree(ax_src, &p_d); CHKERRQ(ierr);
      ierr = IGAAxisGetKnots(ax_src, &N_d, (PetscReal**)&Uknot); CHKERRQ(ierr);
      ierr = IGAAxisGetLimits(ax_src, &U0, &U1); CHKERRQ(ierr);
      ierr = IGAGetAxis(kc->iga, d, &ax_dst); CHKERRQ(ierr);
      ierr = IGAAxisSetPeriodic(ax_dst, PETSC_TRUE); CHKERRQ(ierr);
      ierr = IGAAxisSetDegree(ax_dst, p_d); CHKERRQ(ierr);
      ierr = IGAAxisInitUniform(ax_dst, app->iga->elem_sizes[d], U0, U1,
                                p_d - 1); CHKERRQ(ierr);
    }
    ierr = IGASetMatType(kc->iga, MATAIJ); CHKERRQ(ierr);
    ierr = IGASetUp(kc->iga); CHKERRQ(ierr);
  }

  /* The cross-IGA contract. kc->ice[] is FILLED walking the solver IGA and READ
   * walking this one, so the two must agree element for element. A periodic and
   * an open axis of the same element count have different NODE counts, and
   * PetIGA partitions from the node layout, so in parallel the element
   * decompositions can legitimately diverge -- this is not a formality. Fail
   * loudly rather than silently reading phi from the wrong place. */
  for (PetscInt d = 0; d < kc->dim; d++) {
    if (kc->iga->elem_width[d] != app->iga->elem_width[d] ||
        kc->iga->elem_start[d] != app->iga->elem_start[d])
      SETERRQ(PETSC_COMM_SELF, PETSC_ERR_PLIB,
              "k_eff: the periodic corrector mesh partitions axis %d differently "
              "from the solver mesh (start %d vs %d, width %d vs %d). This "
              "happens when the solver is NOT periodic: the two axes then carry "
              "different node counts and PetIGA can split the elements "
              "differently across ranks, which would make the Gauss-point index "
              "mapping read phi from the wrong location. Run this case on fewer "
              "ranks, or with -periodic 1 so both meshes share a layout.",
              (int)d, (int)kc->iga->elem_start[d], (int)app->iga->elem_start[d],
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
    /* Replaying writes alongside the run being replayed, which is where anyone
     * looking for it will look -- and where $folder usually is not pointing. */
    if (kc->replay_dir[0] != '\0')
      ierr = PetscSNPrintf(kc->csv_path, sizeof(kc->csv_path), "%s/k_eff.csv", kc->replay_dir);
    else if (dir)
      ierr = PetscSNPrintf(kc->csv_path, sizeof(kc->csv_path), "%s/k_eff.csv", dir);
    else
      ierr = PetscSNPrintf(kc->csv_path, sizeof(kc->csv_path), "k_eff.csv");
    CHKERRQ(ierr);
  }
  /* First scheduled time for the time-based cadence. Starting at t_interv
   * rather than 0 keeps step 0 the business of -keff_step0 alone. */
  kc->t_next = kc->t_interv;

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

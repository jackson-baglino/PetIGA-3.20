#ifndef KEFF_H
#define KEFF_H

#include "enceladus_types.h"

/* ---------------------------------------------------------------------------
 * keff.h -- in-line effective thermal conductivity by periodic homogenization.
 *
 * Merged in from projects/effective_thermal_cond (2026-07-31), which used to
 * run as a separate binary over sol_*.dat snapshots after a simulation had
 * finished. That coupled the k_eff sampling rate to the field-output rate, and
 * the two want very different cadences: a field snapshot at study resolution is
 * hundreds of MB (so you want a few hundred), while a k_eff sample is ~200
 * bytes (so you want as many as you can afford).
 *
 * THE METHOD. For each macroscopic direction m we solve the cell problem for
 * the corrector t_m on the periodic unit cell Y:
 *
 *     -div( k(x) grad t_m )  =  div( k(x) e_m )     in Y,   t_m  Y-periodic
 *
 * and then average the resulting microscopic flux:
 *
 *     k_eff[i][j] = (1/|Y|) * INT_Y  k(x) * ( d t_i / d x_j  +  delta_ij ) dV
 *
 * with k(x) = ThermalCond(phi(x)) taken from the live phase field.
 *
 * PERIODICITY IS NOT OPTIONAL. The corrector being periodic on the cell
 * boundary is what makes the cell average an effective property at all. On a
 * non-periodic window the boundary layer contaminates the average and the
 * result describes that particular box with those particular walls, not the
 * material. KeffCreate() therefore hard-errors on a non-periodic mesh; see the
 * guards there.
 *
 * WHY THE CORRECTOR MESH IS SCALAR (dof = 1). The cell-problem operator is the
 * same for every direction m: in the original blocked assembly the integrand
 * sum_g N1[i][g]*N1[j][g]*k never referenced m, so the identical scalar
 * Laplacian was being copied into all dim diagonal blocks. Solving dim scalar
 * systems against one shared operator instead of one dim-dof blocked system
 * cuts the matrix ~4x at dim=2 (9.6 GB -> 2.4 GB at Nx=Ny=2829) and hands GAMG
 * a textbook scalar Laplacian instead of a block system that is 3/4
 * structurally zero.
 * ------------------------------------------------------------------------- */

typedef struct KeffCtx {
  /* --- cadence ---------------------------------------------------------- */
  PetscBool enabled;         /* -keff */
  PetscInt  freq;            /* -keff_freq: sample every N accepted steps */
  PetscReal t_interv;        /* -keff_t_interv: seconds; > 0 overrides freq */
  PetscReal t_next;          /* next scheduled sample time */
  PetscBool at_step0;        /* -keff_step0 */
  PetscBool only;            /* -keff_only: one sample from the IC, then exit */
  char      replay_dir[PETSC_MAX_PATH_LEN];   /* -keff_replay */
  char      replay_times[PETSC_MAX_PATH_LEN]; /* -keff_replay_times */
  PetscInt  replay_stride;   /* -keff_replay_stride */

  /* --- discretisation (all DERIVED; never re-parsed from the options DB) -
   * The old standalone rebuilt an input IGA from -Nx/-p/-C/-periodic and
   * compared vector lengths to catch a mismatch. That check could not detect a
   * -periodic mismatch (its own comment said so), which is exactly the one that
   * silently produces wrong answers. Cloning the live IGA makes the entire
   * failure class unrepresentable. */
  IGA       iga;             /* IGAClone(app->iga, 1) -- SCALAR corrector mesh */
  PetscInt  dim;
  PetscReal vol;             /* cell volume MEASURED by quadrature, not Lx*Ly */

  /* --- linear algebra (created once, reused every sample) ---------------- */
  Mat       A;               /* scalar: INT k grad_Ni . grad_Nj */
  Vec       b[3];            /* one RHS per macroscopic direction */
  Vec       T[3];            /* one corrector per direction; also the initial
                              * guess for the next sample */

  /* --- microstructure ---------------------------------------------------- */
  PetscReal *ice;            /* phi at every local Gauss point; length app->ngp */
  PetscInt   cur_dir;        /* direction m the current form/scalar pass serves */

  /* --- output & bookkeeping ---------------------------------------------- */
  char        csv_path[PETSC_MAX_PATH_LEN];
  PetscViewer csv;           /* lazily opened, flushed per sample */
  PetscBool   write_corrector;
  PetscBool   debug_phibar;  /* -keff_debug_phibar: cross-IGA projection check */
  PetscInt    nsamples;

  AppCtx     *app;           /* back-pointer: thcond_*, iga, dim, ngp */
} KeffCtx;

/* Lifecycle. KeffCreate is a no-op unless -keff is set; both are safe to call
 * unconditionally. KeffDestroy owns everything KeffCreate made and nulls
 * app->keff. */
PetscErrorCode KeffCreate(AppCtx *app);
PetscErrorCode KeffDestroy(AppCtx *app);

/* Integrand returning 1.0, used once at create time to measure the cell volume
 * by quadrature rather than trusting Lx*Ly[*Lz]. */
PetscErrorCode KeffVolumeIntegrand(IGAPoint point, const PetscScalar U[],
                                   PetscInt n, PetscScalar S[], void *ctx);

/* --- keff_field.c ---------------------------------------------------------
 * THE Gauss-point indexing convention for this project is documented at the
 * top of src/keff_field.c: every per-quadrature-point array is indexed by
 *   point->index + point->count * point->parent->index
 * and never by a sequential counter. Read that comment before adding another
 * such array. */

/* Evaluate phi (dof 0 of the 3-dof solution) at every local quadrature point of
 * the solver IGA into kc->ice[], for the cell assembly to read on the clone. */
PetscErrorCode KeffProjectIce(AppCtx *app, Vec U);

/* Startup assertion: point->count uniform across elements of the clone, and
 * the clone's local Gauss-point population equals app->ngp. */
PetscErrorCode KeffCheckGaussLayout(AppCtx *app);

/* Ice moments: m[0] = mean ice fraction, m[1..dim] = ice centroid [m],
 * m[1+dim..2dim] = ice mean square per axis [m^2].
 * from_clone == TRUE reads phi from kc->ice[] on the corrector mesh;
 * FALSE reads it from the 3-dof solution U on the solver mesh.
 * The MOMENTS are the load-bearing part -- see the long comment in
 * keff_field.c for why the mean alone cannot detect a scrambled mapping. */
PetscErrorCode KeffIceMoments(AppCtx *app, PetscBool from_clone, Vec U, PetscReal m[7]);

/* -keff_debug_phibar: project phi onto the clone, then compare the ice mean and
 * centroid computed on each mesh. Hard-errors on disagreement. No-op unless the
 * flag is set. */
PetscErrorCode KeffDebugPhiBar(AppCtx *app, Vec U);

#endif /* KEFF_H */

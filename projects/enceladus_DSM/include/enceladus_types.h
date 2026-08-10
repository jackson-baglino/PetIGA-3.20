#ifndef ENCELADUS_TYPES_H
#define ENCELADUS_TYPES_H

#include <petsc.h>
#include "petiga.h"

#define SQ(x) ((x)*(x))
#define CU(x) ((x)*(x)*(x))

/* Field definitions for node data */
typedef struct {
  PetscScalar ice, tem, rhov;
} Field;
/* Application context structure */
typedef struct {
  IGA       iga;  // Isogeometric analysis (IGA) structure for managing geometry and basis functions
  SNES      snes; // Nonlinear solver handle (cached so Residual can call SNESSetFunctionDomainError)

  // Physical parameters related to phase field and thermodynamics
  PetscReal eps;  // Interface width parameter for phase field method
  PetscReal mob_sub;  // Mobility for ice phase evolution
  PetscReal Etai, Etam, Etaa;  // Surface energy terms: Sigma_i (ice-vapor), Etam (unused), Sigma_a (air-vapor side)
  PetscReal alph_sub;  // Substrate interaction coefficient
  PetscReal Lambd;  // Parameter related to thermal conductivity or latent heat (context-dependent)
  PetscReal beta_sub0, d0_sub0;  // Parameters related to phase change at the substrate

  // Thermophysical properties of different phases
  PetscReal thcond_ice, thcond_air;  // Thermal conductivities of ice and air
  PetscReal cp_ice, cp_air;  // Specific heat capacities of ice and air
  PetscReal rho_ice, rho_air;  // Densities of ice and air
  PetscReal dif_vap;  // Vapor diffusivity in air
  PetscReal lat_sub;  // Latent heat of sublimation
  PetscReal diff_sub;  // Diffusivity related to sublimation

  // Environmental conditions and threshold parameters
  PetscReal air_lim;  // Air phase fraction floor used in the vapor equation (avoids division/degeneracy as phi_a -> 0)

  // Initial and boundary condition parameters
  PetscReal T_melt;  // Melting temperature of ice
  PetscReal temp0, hum0;  // Initial temperature and humidity
  PetscReal grad_temp0[3];  // Initial temperature gradient in x, y, and z directions

  // Domain size and resolution
  PetscReal Lx, Ly, Lz;  // Physical domain dimensions in x, y, and z
  PetscInt  Nx, Ny, Nz;  // Number of elements in x, y, and z directions

  // Radius of curvature parameters (possibly for computing capillary effects)
  PetscReal RCice;  // Mean radius of curvature for ice grains
  PetscReal RCice_dev;  // Standard deviation of radius of curvature for ice grains

  // Per-grain radii for boundary-centered two-grain IC
  PetscReal RCice0, RCice1;   /* Outer ice radius of grain 0 / grain 1 (default: RCice) */

  // Sediment-grain "bump" geometry parameter (must match -geom_file's
  // build_geometry_sediment_grain.py R_sed; 0 => flat domain, no distortion)
  PetscReal geom_bump_R;

  // Multi-grain sediment "bump" geometry: the bottom edge of a -geom_file
  // bump geometry is the sum of n_sed_grains individual SedimentBump()
  // humps. Each bump has a center (sed_grain_x), half-width (sed_grain_R),
  // and peak height (sed_grain_h, defaults to sed_grain_R if not provided).
  // Must match build_geometry_multi_grain.py's SEDIMENT_GRAINS list.
  // n_sed_grains == 0 falls back to the single-bump geom_bump_R behavior.
#define MAX_SED_GRAINS 24
  PetscInt  n_sed_grains;
  PetscReal sed_grain_x[MAX_SED_GRAINS];
  PetscReal sed_grain_R[MAX_SED_GRAINS];
  PetscReal sed_grain_h[MAX_SED_GRAINS];  /* peak height; default = R */

  // Top-wall (ceiling) bumps — same C∞ shape, push DOWN from Ly.
  // Must match TOP_GRAINS in build_geometry_multi_grain.py.
  PetscInt  n_top_grains;
  PetscReal top_grain_x[MAX_SED_GRAINS];
  PetscReal top_grain_R[MAX_SED_GRAINS];
  PetscReal top_grain_h[MAX_SED_GRAINS];

  // Arrays storing geometry information for ice grains.
  //
  // Heap-allocated with capacity n_grain_max (-n_grain_max, default 2000).
  // These were fixed [200] arrays until 2026-07-27; granular packings from
  // preprocess/generate_packing.py routinely carry 400-500 grains at
  // 2 mm / 45 um mean radius, which silently overran them. Capacity is
  // allocated once in main() and freed at the end; indexing is unchanged.
  PetscInt  n_grain_max;   // Allocated capacity of the four arrays below
  PetscReal *cent[3];      // Ice grain centre coordinates: cent[dim][k]
  PetscReal *radius;       // Radii of individual ice grains (isotropic; default ax/ay)
  PetscReal *ice_grain_ax; /* ellipse semi-axis in x; defaults to radius[k] if -ice_grain_ax not set */
  PetscReal *ice_grain_ay; /* ellipse semi-axis in y; defaults to radius[k] if -ice_grain_ay not set */

  // Ice "shell" capping a floor bump at constant thickness, conformal to
  // the bump's own surface (true distance to the SedimentBumpField(x)
  // curve, not just a vertical offset -- see SedimentBumpFieldDeriv() and
  // the -ice_shell_x loop in FormInitialMultiGrains2D), windowed laterally
  // to [ice_shell_x[k]-ice_shell_R[k], ice_shell_x[k]+ice_shell_R[k]] so it
  // only covers the bump itself, not the whole floor. Added on top of the
  // ice_grain_* ellipses, not a replacement.
  PetscInt  n_ice_shells;
  PetscReal ice_shell_x[MAX_SED_GRAINS];
  PetscReal ice_shell_R[MAX_SED_GRAINS];
  PetscReal ice_shell_thickness[MAX_SED_GRAINS];

  // Flat ice layer encapsulating a floor bump: ice fills everything below
  // the ABSOLUTE height ice_flat_height[k] (not relative to the bump's own
  // surface, unlike ice_shell_*), windowed laterally to
  // [ice_flat_x[k]-ice_flat_R[k], ice_flat_x[k]+ice_flat_R[k]]. Gives a flat
  // (non-rounded) ice-air interface burying the bump, instead of a domed
  // cap or a conformal coating. Added on top of ice_grain_*/ice_shell_*.
  PetscInt  n_ice_flats;
  PetscReal ice_flat_x[MAX_SED_GRAINS];
  PetscReal ice_flat_R[MAX_SED_GRAINS];
  PetscReal ice_flat_height[MAX_SED_GRAINS];

  // Initial normal vector components (possibly for a structured interface)
  PetscReal norm0[3];  // Per-DOF initial residual norms for SNES convergence check

  // Flags for controlling different simulation options
  PetscInt  flag_tIC;        // IC geometry variant: 0=centered slab, 2=flat interface
  PetscInt  outp;            // output control flag
  PetscBool flag_Tdep;       // temperature-dependent material properties
  PetscBool decouple_phase_change;  // -decouple_phase_change 1: zero the ice_t-driven
                                     // source terms in R_tem (latent heat) and R_vap
                                     // (mass-balance) too, not just S_sub in R_ice --
                                     // isolates pure AC curvature relaxation from any
                                     // leakage into temperature/vapor (see assembly.c)

  // Numerical method and discretization parameters
  PetscInt p;  // Polynomial degree of basis functions (for IGA)
  PetscInt C;  // Continuity of basis functions
  PetscInt dim;  // Spatial dimension of the problem (2D or 3D)
  PetscInt dof;  // Degrees of freedom per node (ice, temperature, vapor)
  PetscInt periodic;  // Periodicity flag (0 = non-periodic, 1 = periodic boundaries)

  // Time stepping parameters
  PetscReal t_out;  // Output time interval
  PetscReal t_interv;  // Intermediate time step interval
  PetscReal t_IC;  // Total duration for initial condition phase

  // -t_out_log N: write N snapshots at LOG-SPACED times instead of the
  // every-N-steps (-outp) or time-uniform (-t_interv / -n_out) cadences.
  //
  // Power-law diagnostics need it. Fitting a growth exponent means reading a
  // slope in log-log, so the samples have to be spread over the decades, but
  // time-uniform output puts ~90% of them in the last decade and -outp puts
  // them wherever the adaptive dt happens to land. On a large mesh the
  // difference is not cosmetic: one snapshot of the 9792x2671 sintering grid
  // is ~630 MB, so the snapshot budget is ~40, and spending them uniformly
  // leaves the early neck completely unresolved.
  //
  // Times run from t_out_log_t0 (default: the first accepted step's t, which
  // is where the log axis can actually start) to t_final, inclusive. Schedule
  // is precomputed in enceladus_main.c; OutputMonitor consumes it in order.
  PetscInt   n_out_log;      // number of log-spaced snapshots; 0 = disabled
  PetscReal  t_out_log_t0;   // first scheduled time [s]; <= 0 = use first step's t
  PetscReal *t_out_log;      // ascending schedule, length n_out_log (PetscMalloc'd)
  PetscInt   i_out_log;      // index of the next unwritten entry

  // Counters for active ice grains
  PetscInt NCice;  // Number of ice grains
  PetscInt n_act;  // Number of currently active ice grains

  PetscInt npoints;

  /* Anti-trapping / stabilization multipliers (Moure & Fu 2024, K&P 2009).
   * xi_T scales thermal conduction and latent heat in R_tem.
   * xi_v scales vapor diffusion and the rho_ice part of the mass exchange
   *      source in R_vap — reduces spurious interface fluxes from the
   *      finite-width diffuse interface to the physical sharp-interface limit.
   * Defaults: xi_T = 1.0 (no effect), xi_v = 1e-3 (M&F sublimation value). */
  PetscReal xi_T;
  PetscReal xi_v;

  /* Phase-field bounds: simulation aborts if any phi leaves [phase_lo, phase_hi] */
  PetscReal phase_lo;     // lower bound for phi_ice, phi_air (default -0.25)
  PetscReal phase_hi;     // upper bound for phi_ice, phi_air (default  1.25)

  // NOTE: per-quadrature-point alph[]/mob[] arrays were removed 2026-07-31.
  // They were filled under -flag_Tdep and read by nothing: assembly.c uses the
  // scalars alph_sub/mob_sub everywhere. See the comment in monitoring.c's
  // flag_Tdep block. If spatially varying kinetics are ever wired into the
  // residual, reintroduce the array together with the code that reads it, and
  // index it as point->index + point->count * point->parent->index (the one
  // Gauss-point indexing convention in this project -- see keff_field.c).

  // Flag for reading input files
  PetscBool readFlag; // read initial field data from file

  // Output file path
  char output_path[PETSC_MAX_PATH_LEN];  // Path for output files
  char initial_cond[PETSC_MAX_PATH_LEN];  // Path for initial condition file
  char initial_PFgeom[PETSC_MAX_PATH_LEN];  // Path for initial geometry file
  char grains_file[PETSC_MAX_PATH_LEN];  // Grain list for -ic_type multi_grains_file

  // Capillary neck parameters
  PetscReal R1;  // Radius of capillary neck

  // Initial domain integrals (set at step 0, used for percentage reporting)
  PetscReal tot_ice_0;     // initial ∫ φ_i dΩ
  PetscReal tot_air_0;     // initial ∫ φ_a dΩ
  PetscReal tot_rhov_0;    // initial ∫ ρ_v φ_a dΩ (vapor mass in air phase)
  PetscReal tot_mass_0;    // initial ρ_ice·∫φ_i + ∫ρ_v·φ_a (total system mass)

  // Deferred bounds-rollback request — set by Monitor() when phase fields go
  // out of bounds, consumed by BoundsRollbackPreStep() before the next TSStep.
  // We can't call TSRollBack inside Monitor because ts->vec_sol is read-locked
  // there; the PreStep callback runs when the vector is writable.
  PetscBool bounds_violated;
  PetscReal bounds_new_dt;

  // Interface-CFL timestep limiter (InterfaceCFLMonitor in monitoring.c):
  // clamps the NEXT dt so max pointwise |d(phi)| per step stays below
  // cfl_dphimax, measured from the last accepted step's rate
  // ||phi^n - phi^{n-1}||_inf / dt. Lets dtmax be large (quiet phases
  // cruise) while fast interface events (grain collapse, neck formation)
  // throttle dt automatically instead of Gibbs-rippling the B-spline front.
  PetscBool flag_dtCFL;      // -dtCFL (default 1)
  PetscReal cfl_dphimax;     // -dtCFL_dphimax (default 0.2)
  Vec       cfl_U_prev;      // previous accepted solution (lazy-created)
  PetscReal cfl_t_prev;      // time of previous accepted step

  // Axisymmetric (r-z) mode: x = z (symmetry axis direction), y = r (radial),
  // axis on the y = 0 boundary. Every residual/Jacobian/monitor integrand is
  // weighted by the quadrature point's r-coordinate (2*pi*r dr dz measure);
  // the azimuthal curvature mode is generated by that weight automatically
  // (see docs/axisymmetric_plan.md section 1b). Grains must be centered on
  // the axis (ice_grain_cy = 0).
  PetscBool axisym;          // -axisym (default 0 = planar)

  // -ic_grain_union: build the multi_grains IC from the signed distance to the
  // UNION of the grains (sdf = min_k sdf_k, phi = 0.5-0.5*tanh(0.5*sdf/eps))
  // instead of SUMMING each grain's tanh profile.
  //
  // Why it matters for OVERLAPPING grains: the additive form's phi=0.5 contour
  // depends on eps at the neck, because both grains' tails contribute there.
  // (Away from the neck a lone grain crosses 0.5 at r=R for any eps, which is
  // why only the neck misbehaves.) With the union form phi=0.5 <=> sdf=0, i.e.
  // exactly the sharp union surface, for ANY eps -- so an eps series shares one
  // initial contour instead of one per eps (38.19/35.60/32.91 um measured in
  // the 2026-07-15 series).
  //
  // Caveat: the sharp union has a zero-radius concave crease at the neck, which
  // no diffuse interface can hold; AC dynamics rounds it to O(eps) within a few
  // steps. Contours coincide at t=0 and then diverge by O(eps) at the neck --
  // that residual IS the eps error being measured.
  //
  // Default 0 (additive) preserves every pre-2026-07-15 run. The two forms
  // agree for well-separated grains, where the far grain's tanh is negligible.
  PetscBool ic_grain_union;  // -ic_grain_union (default 0 = additive)

  // Persistent viewer for SSA_evo.dat. Opened ONCE on the first Monitor() call
  // and flushed every step, instead of re-opening a viewer per step (which
  // silently exhausted file descriptors ~step 15 -- the calls had no CHKERRQ --
  // truncating the file while outp.txt kept going). NULL until first use.
  PetscViewer ssa_view;

  // Number of quadrature points owned by this rank, sized exactly as the alph
  // and mob arrays are. Recorded here so the k_eff module can allocate its own
  // per-Gauss-point phi array without re-deriving elem_width x (p+1)^dim, and
  // so the cloned corrector IGA can assert it sees the same count.
  PetscInt ngp;

  // In-line effective thermal conductivity by periodic homogenization.
  // NULL unless -keff is set; owned by KeffCreate/KeffDestroy (see keff.h).
  // Declared as an incomplete type so enceladus_types.h stays independent of
  // keff.h, which includes it.
  struct KeffCtx *keff;

} AppCtx;/* Field definitions for node data */

#endif // ENCELADUS_TYPES_H
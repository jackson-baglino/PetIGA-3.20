
## 2026-07-01 — mob_sub sweep scripts + mesh/eps/kinetics decoupling

- Tightened K&P bounds: safety 0.5 → 0.25, eps halved to 2.4267e-7, Nx/Ny doubled to 170/210.
- Decoupled mesh refinement from kinetics: eps follows h=eps rule (~7.5 elements/interface); mob_sub and alph_sub frozen at safety=0.5 K&P values (6.5404e-10, 2.8363e6) via explicit opts overrides in geometry file.
- Added Studio mob_sub sweep script (scripts/Studio/run_mob_sweep.sh) for local sequential runs.
- Added HPC mob_sub sweep (scripts/HPC/submit_mob_sweep.sh): fans out one SLURM job per value via --wrap, each with its own subfolder under a shared parent; fixed bug where BATCH_OUT_DIR caused all jobs to write to the same folder (geom__exp naming collision).
- Saved memory: always push to GitHub before giving HPC run commands.

---

**Session ended:** 2026-07-01 14:32:36


---

**Session ended:** 2026-07-01 14:29:14


---

**Session ended:** 2026-07-01 14:27:08


---

**Session ended:** 2026-07-01 14:25:35


---

**Session ended:** 2026-07-01 14:24:28


---

**Session ended:** 2026-07-01 14:06:34


---

**Session ended:** 2026-07-01 14:02:43


---

**Session ended:** 2026-07-01 13:55:50


---

**Session ended:** 2026-07-01 13:54:50


---

**Session ended:** 2026-07-01 13:49:22


---

**Session ended:** 2026-07-01 13:44:43


---

**Session ended:** 2026-07-01 13:39:51


---

**Session ended:** 2026-07-01 13:29:14


---

**Session ended:** 2026-07-01 13:25:47


---

**Session ended:** 2026-07-01 12:22:34


---

**Session ended:** 2026-07-01 12:16:29


---

**Session ended:** 2026-07-01 11:59:31


---

**Session ended:** 2026-06-30 16:53:27


---

**Session ended:** 2026-06-30 16:41:48


---

**Session ended:** 2026-06-30 16:37:08


---

**Session ended:** 2026-06-30 16:36:26


---

**Session ended:** 2026-06-30 15:58:49


---

**Session ended:** 2026-06-30 15:57:31


---

**Session ended:** 2026-06-30 15:56:47


---

**Session ended:** 2026-06-30 15:56:42


---

**Session ended:** 2026-06-30 15:26:36


---

**Session ended:** 2026-06-29 15:54:09


---

**Session ended:** 2026-06-29 15:46:31


---

**Session ended:** 2026-06-29 15:35:53


---

**Session ended:** 2026-06-29 15:32:13


## 2026-06-29 — Switch to p=2/C1 basis, refine 2D two-grain mesh, overhaul parameter print

- Identified that all simulations were using linear basis (p=1, C=0) due to explicit settings in `inputs/solver.opts`; changed to p=2, C=1 (quadratic/C1) globally as the new default — both in `solver.opts` and the `permafrost2.c` code defaults.
- Updated `2D_two_ice_grains_boundary.opts`: increased mesh 1.5× (125→188, 154→231) to address tanh-profile kink artifact observed after ~100 steps in the 100yr run; eps updated to `2.5080e-07 m` (user-calibrated to achieve exactly 7.5 elements across the diffuse interface, slightly finer than the naive eps/1.5 scaling).
- Rewrote the simulation parameter print in `permafrost2.c`: removed the early scattered print blocks (scattered across before/after IGASetUp) and replaced with a single comprehensive block printed after all setup (IGASetUp + VI-bounds) so true Nx/Ny/p/C and derived kinetics are all available. Organized into sections: MESH & DISCRETIZATION, PHASE-FIELD INTERFACE, ENVIRONMENT & IC, TIME STEPPING, TRANSPORT & THERMOPHYSICAL PROPERTIES, PHASE-CHANGE KINETICS, SOLVER, BOUNDARY CONDITIONS.
- Ran smoke test (2-step IC verification): confirmed p=2/C1, Nx=188/Ny=231, eps=3.1099e-07 print correctly; output saved to `SimulationResults/permafrost/scratch/smoketest_2D_two_ice_grains_boundary_p2_refined1p5x/`.

---

**Session ended:** 2026-06-23 18:09:35


---

**Session ended:** 2026-06-23 17:35:57


---

**Session ended:** 2026-06-23 12:22:34


---

**Session ended:** 2026-06-23 12:11:01


---

**Session ended:** 2026-06-23 12:05:28


---

**Session ended:** 2026-06-23 11:58:05


---

**Session ended:** 2026-06-23 11:42:31


---

**Session ended:** 2026-06-23 11:35:49


---

**Session ended:** 2026-06-22 18:29:17


---

**Session ended:** 2026-06-22 18:26:55


---

**Session ended:** 2026-06-22 18:26:03


---

**Session ended:** 2026-06-22 18:15:42


---

**Session ended:** 2026-06-22 18:11:18


---

**Session ended:** 2026-06-22 17:45:27


---

**Session ended:** 2026-06-22 17:40:14


---

**Session ended:** 2026-06-22 17:38:18


---

**Session ended:** 2026-06-22 17:30:47


---

**Session ended:** 2026-06-22 16:39:59


---

**Session ended:** 2026-06-22 16:32:31


---

**Session ended:** 2026-06-22 16:17:14


---

**Session ended:** 2026-06-22 16:09:47


---

**Session ended:** 2026-06-22 16:09:38


---

**Session ended:** 2026-06-22 16:05:59


---

**Session ended:** 2026-06-22 16:04:42


---

**Session ended:** 2026-06-22 16:03:06


---

**Session ended:** 2026-06-22 16:00:54


---

**Session ended:** 2026-06-22 15:45:16


---

**Session ended:** 2026-06-22 15:40:56


---

**Session ended:** 2026-06-22 15:38:50


---

**Session ended:** 2026-06-22 15:17:16


---

**Session ended:** 2026-06-22 15:14:32


---

**Session ended:** 2026-06-22 15:13:39


---

**Session ended:** 2026-06-22 15:11:46


---

**Session ended:** 2026-06-22 11:11:16


---

**Session ended:** 2026-06-22 11:07:42


---

**Session ended:** 2026-06-22 11:04:13


---

**Session ended:** 2026-06-22 11:01:00


---

**Session ended:** 2026-06-22 10:54:29


---

**Session ended:** 2026-06-22 10:50:57


---

**Session ended:** 2026-06-22 10:49:24


---

**Session ended:** 2026-06-22 10:45:09


---

**Session ended:** 2026-06-22 10:15:10


---

**Session ended:** 2026-06-22 10:11:18


---

**Session ended:** 2026-06-22 10:10:30


---

**Session ended:** 2026-06-22 10:04:25


---

**Session ended:** 2026-06-22 08:59:47


---

**Session ended:** 2026-06-22 08:59:36


---

**Session ended:** 2026-06-22 08:56:55


---

**Session ended:** 2026-06-22 08:35:34


---

**Session ended:** 2026-06-22 08:33:24


---

**Session ended:** 2026-06-22 08:30:10


---

**Session ended:** 2026-06-22 08:27:48


---

**Session ended:** 2026-06-22 08:16:46


---

**Session ended:** 2026-06-22 08:15:37


---

**Session ended:** 2026-06-22 08:02:26


---

**Session ended:** 2026-06-22 07:53:18


---

**Session ended:** 2026-06-21 18:55:52


---

**Session ended:** 2026-06-21 17:42:56


---

**Session ended:** 2026-06-21 17:41:22


---

**Session ended:** 2026-06-21 17:25:20


---

**Session ended:** 2026-06-21 17:19:17


---

**Session ended:** 2026-06-21 17:15:29


---

**Session ended:** 2026-06-21 16:46:12


---

**Session ended:** 2026-06-21 16:43:31


---

**Session ended:** 2026-06-21 16:26:33


---

**Session ended:** 2026-06-21 15:58:08


---

**Session ended:** 2026-06-21 15:52:42


---

**Session ended:** 2026-06-21 15:32:25


---

**Session ended:** 2026-06-21 15:24:04


---

**Session ended:** 2026-06-21 14:13:09


---

**Session ended:** 2026-06-21 13:10:13


---

**Session ended:** 2026-06-21 13:05:06


---

**Session ended:** 2026-06-21 11:27:45


---

**Session ended:** 2026-06-21 10:20:58


---

**Session ended:** 2026-06-20 17:32:59


---

## 2026-06-20 — Mesh-convergence investigation on grain plateau; pivoted to pore-channel geometry

- 10-year run on 2D_two_ice_grains_boundary (half-size small grain) showed
  99% small-grain sublimation before plateauing -- visual check showed the
  grain's two opposing diffuse interfaces had merged into a degenerate
  blob (no flat phi_ice~1 core) by year 6.5, frozen for the remaining 3.5
  years. Diagnosed as a resolution-validity limit: once a grain's radius
  approaches the interface width, curvature degrades toward zero and the
  Gibbs-Thomson driving force has nothing left to act on.
- Built a 2x-refined geometry (Nx/Ny x2, eps /2, dx/eps held exactly
  constant at 0.70711) to test this. Result was the OPPOSITE of expected:
  the refined run plateaued earlier and at a much larger residual (-64%
  vs -99%), while staying visually well-resolved (clean core, no merged
  interface) throughout. The gap between original and 2x trajectories
  grows continuously from t=0, not just near the end -- revised diagnosis:
  the original (coarser) mesh was systematically overestimating curvature
  (Hessian-based, second-derivative-sensitive to discretization) as the
  grain shrank, artificially inflating the GT driving force throughout the
  run, not just failing once the grain got tiny. The original "near-total
  sublimation" result was not mesh-converged.
- Tried a 4x-refined geometry (same dx/eps ratio) plus -dtmax 5e5 (from
  1e5) to keep wall-clock down, since dtmax was being hit at step 105/day
  4.6 in the 2x run (essentially the whole run at the ceiling). User
  judged this was rabbit-holing -- mesh refinement alone isn't resolving
  the question productively, and called it (job left to finish in the
  background, not analyzed further). Decided to stop chasing full
  sublimation via mesh refinement.
- Discussed: this is a well-documented issue in the alloy/grain-growth
  phase-field literature ("small particle disappearance problem",
  "lattice/mesh pinning"), typically addressed via adaptive mesh
  refinement, explicit small-grain elimination rules, or multi-order-
  parameter grain tracking -- all bigger undertakings than today's time
  budget (conference presentation due Tuesday 2026-06-23) allows.
- Pivoted to a new deliverable: built 2D_pore_channel.opts, a tortuous
  pore-channel geometry through sediment-laden regolith (5 floor bumps +
  4 interleaved, differently-sized ceiling bumps, varied R/height for
  visual irregularity vs the uniform 5-bump production baseline), reusing
  the proven 2D_multi_grain_test.opts ice-grain layout (2 boundary
  attractor grains + 4 trough grains) and eps. Added --top-bumps to
  build_geometry_multi_grain.py (mirrors --bumps, was bottom-only before)
  to support this. Genuine LEFT/RIGHT side-wall intrusions (also
  requested) were explicitly deferred -- the current ruled-surface mesh
  construction only warps y as a function of x; doing both directions at
  once needs a true 4-edge Coons-patch blend, judged too risky to attempt
  for the first time under today's deadline.

---

## 2026-06-20 — Diagnosed d0_GT amplification failure; switched to tighter geometry + warmer temp

- Two 21-day sweeps on 2D_single_bump_two_grains at -20°C testing larger
  d0_GT (1e-7, 1e-8) to speed up Ostwald ripening: d0_GT=1e-7 caused
  I-A INTERF to jump 3.1e-12 -> 4.3e-12 in a few steps then plateau --
  visual check (igakit render of sol_00045/48/121) confirmed this is the
  diffuse profile itself fattening on the small grain, not faster mass
  transfer. d0_GT=1e-8 stayed smooth but moved INTERF only ~1.4% over 21
  days -- confirms no single d0_GT value threads the needle: too small is
  invisible, too large breaks the equilibrium profile (the GT forcing on
  S_sub overwhelms the gradient-energy/double-well restoring balance that
  holds the profile shape, since both act on similar wavelengths near the
  interface).
- Switched tactics per user: raised temperature to -5C instead (already
  had 30day_T-5_h1.00_GTamp/_GTphys.opts staged from earlier exploration,
  rationale: rho_vs(T) is exponential in 1/T, so GT's multiplicative
  correction produces a much larger absolute vapor-density difference at
  -5C than -20C without touching d0_GT at all). Ran both on
  2D_single_bump_two_grains/30day: GTphys (d0_GT=9.6e-10, physical)
  declined smoothly 3.081->2.978e-12 (-3.3%, monotonic, zero DIVERGED);
  GTamp (d0_GT=1e-8, 10x) oscillated mildly (3.081->3.095->3.048, ~1%
  swings, no instability) -- confirms -5C gives real headroom for the
  same 10x amplification that broke cleanly at -20C, and -5C alone is the
  single biggest lever found so far. Still far from full small-grain
  sublimation though.
- User's next idea: switch to the tighter 2D_two_ice_grains_boundary
  geometry (already existed, unused this whole session) -- per its own
  header it's purpose-built with a ~40-element vapor gap between a 9.375
  um and an 18.75 um grain on a near-square domain (4.12e-5 x 5.08e-5 m),
  vs. the 1.0e-4 m wide single_bump_two_grains domain used so far. IC-only
  local check confirms clean geometry (rendered, grains positioned as
  documented). Added 90day_T-5_h1.00_GTamp/_GTphys.opts (mirrors the
  existing 30-day versions, just t_final=90 days) for a longer-horizon
  test, per user's request to also try longer t_final.
- Not yet run: 2D_two_ice_grains_boundary x {30day,90day}_T-5_h1.00_{GTamp,GTphys}
  on the HPC. This combo (tighter grain spacing + warmer temp + longer
  time + already-validated-safe 10x GT amplification) stacks every lever
  found so far without revisiting the ones that failed (over-amplifying
  d0_GT alone, -20C alone).

---

**Session ended:** 2026-06-19 19:28:40


---

**Session ended:** 2026-06-19 19:17:59


---

**Session ended:** 2026-06-19 18:41:10


---

## 2026-06-19 — Confirmed bcgs fix on HPC, raised dtmax to 1.0e5

- job64418602 (2D_single_bump_two_grains, bcgs_fix) confirms the gmres ->
  bcgs switch fully resolved the dt stall: 0 DIVERGED events, 79 "Increase
  time step" vs 0 "Reduce time step" (vs the prior 325/272 near-1:1
  hunting cycle), dt cleanly hit the dtmax=1.0e4 ceiling 12 times. Full
  2-day run finished in 2.57 min wall-clock on 1 rank.
- Since dt was being limited purely by the ceiling (not any solver
  instability) once VI bounds + bcgs removed both prior failure modes,
  raised -dtmax from 1.0e4 to 1.0e5 (~1.16 days/step) in solver.opts --
  the NRmin/NRmax heuristic still self-limits based on Newton iteration
  count, so the higher ceiling only gets used during genuinely quiet
  stretches, which is most of a 21+ day Ostwald-ripening run.
- Next: rerun 2D_single_bump_two_grains (or 2D_single_bump_ice_cap) at
  21day_T-20_h0.95 on the HPC with the new dtmax to see ripening develop
  over a longer horizon, now markedly cheaper to run.

---

## 2026-06-19 — Fixed GT-induced dt stall: switched KSP from gmres to bcgs

- User reported job64416684 (2D_single_bump_two_grains, GT fix) confirmed
  Ostwald ripening works again, but dt stalls oscillating around 1e3
  instead of growing toward dtmax=1.0e4. Counts: 325 "Increase time step"
  vs 272 "Reduce time step" (near 1:1, a hunting/limit-cycle), versus the
  pre-GT baseline's 146 Increase vs only 15 Reduce -- confirms GT is the
  trigger, not a pre-existing issue.
- Reproduced locally with -snes_converged_reason -ksp_converged_reason:
  every "Reduce time step" event is `Linear solve did not converge due to
  DIVERGED_BREAKDOWN iterations 400` -> `Nonlinear solve did not converge
  due to DIVERGED_LINEAR_SOLVE iterations 0` -- GMRES breaking down on the
  very first Newton iteration once dt grows large, not a Newton-iteration-
  count (NRmax) issue and not a bound violation (zero WARN/bound-rollback
  events). Cause: the GT Jacobian's curvature chain-rule term divides by
  |grad_ice|^3 and |grad_ice|^5 in the diffuse band (|grad_ice| ~ 1/eps ~
  1e6), and at large dt the stabilizing shift*mass-matrix diagonal is too
  weak to keep the linear system well-conditioned for GMRES.
- Exact same failure mode (`GMRES DIVERGED_BREAKDOWN`) already has
  precedent in this project -- `2D_touching_grains_hires.opts` already
  overrides to `-ksp_type bcgs` per-file for it. Moved that fix to
  `solver.opts` globally instead, since GT is now active by default for
  every case, not just hires/touching-grain ones.
- Confirmed on HPC (job64418602) -- see entry above.

---

## 2026-06-19 — Restored Gibbs-Thomson curvature term (Ostwald ripening)

- After the VI-solver fix (below) eliminated the bound-violation artifact,
  user noticed the model no longer shows any curvature-driven sublimation
  or Ostwald ripening -- correctly diagnosed that what previously *looked*
  like curvature-driven ripening was actually phi going out of [0,1] near
  high-curvature points, eroding the phase field in a way that showed up as
  vapor. With that artifact gone, rho_vs(T) (flat-interface only, no
  curvature dependence) has no mechanism left to drive small-grain ->
  large-grain mass transfer.
- Found this term already existed in this codebase before the
  "rewrite/2phase-from-equations" branch's clean rewrite of assembly.c
  (commits 5108210/5d6424e) -- a working `Curvature()` function in
  material_properties.c (regularized kappa = -L/G + gHg/G^3 via
  `IGAPointFormHess`), a `-d0_GT` parameter, and full analytic Jacobian
  terms, all dropped when assembly.c was rewritten from the variational
  weak form. `docs/model_description.md` ​§4 already documents this exact
  failure mode ("Without GT the model freezes... no driving force for
  ripening").
- Ported `Curvature()` back into `material_properties.c`/`.h` verbatim.
  Wired `rhoI_vs_eff = rho_vs*(1 + d0_GT*kappa)` into the *current*
  (rewritten) `Residual_A1`/`Jacobian_A1` in `assembly.c` -- simpler than
  the old port since the current model's vapor equation sources from
  `ice_t` (not `S_sub` directly), so no `[vap,ice]` GT Jacobian block is
  needed, only `[ice,ice]`/`[ice,tem]`.
- `d0_GT` defaults to `user.d0_sub0` (the bare 1.0e-9 m capillary-length
  constant), NOT the locally-computed `d0_sub = d0_sub0/rho_rhovs` used to
  derive `alph_sub`/`mob_sub` just above it in `permafrost2.c` -- caught via
  a quick magnitude check before committing: `d0_sub` is rescaled by the
  huge `rho_ice/rho_vs(T)` ratio for the Karma-Plapp kinetic-coefficient
  asymptotics specifically, and using it for the standalone Kelvin equation
  would give a curvature correction to rho_vs of ~2e-8% for a 5um grain
  (utterly negligible) versus ~0.02% with the bare `d0_sub0` (close to the
  literature ice capillary length, ~9.6e-10 m). `-d0_GT 0` still disables it
  for diagnostics.
- Verified: clean build, IC-only sanity run shows `BOUNDS: phi_ice [0.0,
  1.0]` with `d0_GT (capillary length): 1.0000e-09 m [GT active]`. Full
  2-day `2D_single_bump_two_grains` validation run launched into
  `SimulationResults/permafrost/scratch/` (not /tmp, per standing
  instruction) to confirm Newton convergence and Ostwald ripening behavior
  before pushing further.

---

## 2026-06-19 — Bound-constrained (VI) Newton solve replaces dt/gate heuristics for phase bounds

- User reported that job64415277 (2D_single_bump_ice_cap, the v8 circle-fit
  encapsulation) saw phi_ice undershoot to -0.18 and -0.196 at steps 63/70,
  voiding those results -- much worse than the small (~-0.004), expected/
  tolerated overshoot seen by the end of the run. Traced every BOUNDS/WARN
  line in outp.txt: the violence was concentrated in steps 61-71 (the
  encapsulating ring's full extinction), where the existing rollback-and-
  retry (-phase_lo/-phase_hi, dt halving) fired repeatedly and still only
  got the violation down to -0.10 to -0.18 before accepting it -- confirmed
  via dense-evaluated (not just raw-DOF) fields that this was a genuine
  field-level undershoot, not a spline-coefficient artifact. dt-shrinking
  had diminishing returns: a signature of genuine AC-extinction stiffness,
  not just "dt was too big."
- Root fix: bound-constrained Newton solve. `permafrost2.c` now builds
  per-DOF bound vectors (Xl, Xu) via `IGACreateVec`/`VecStrideSet` and
  calls `SNESVISetVariableBounds(nonlin, Xl, Xu)` enforcing 0<=ice<=1
  exactly (temperature/vapor left unconstrained) -- since the field at any
  quadrature point is a convex combination of nearby DOFs, bounding the
  DOFs also bounds the field everywhere, with no retry/rollback needed.
  `-phase_lo`/`-phase_hi` and the rollback mechanism are left in
  `solver.opts` as a defensive backstop only.
- Tried `-snes_type vinewtonrsls` (reduced-space active-set) first: hit a
  PETSc-internal bug at the very first active-set change
  (`MatCreateSubMatrix_MPIBAIJ` "Nonconforming object sizes", with -n 6) --
  a PETSc/MPIBAIJ parallel limitation, not anything in this code. Switched
  to `-snes_type vinewtonssls` (semismooth reformulation, no active-set
  submatrix extraction): clean.
- Validated locally (`scripts/Studio/run_batch_tests.sh --tag
  vinewtonssls_validation`, full 2-day runs into
  `SimulationResults/permafrost/scratch/`): both
  `2D_single_bump_ice_cap:2day_T-20_h0.95` and
  `2D_single_bump_two_grains:2day_T-20_h0.95` (regression check on normal,
  non-extinction sintering dynamics) completed OK with zero WARNs and
  every BOUNDS line at machine-zero.
- Unexpected, important finding: with the fix, the ice-cap's encapsulating
  ring no longer fully sublimates within 2 days (I-A interfacial area only
  -3.7%, smooth and continuous) -- whereas the old, unconstrained run had
  it vanish completely (-23% interfacial area, but in a sudden cliff during
  exactly the steps 61-71 violation window, then completely flat for the
  remaining ~150,000 s of simulated time). Root cause: the AC residual's
  curvature/diffusion term (`3*mob_sub*eps*grad_N_dot_grad_ice` in
  assembly.c) uses the *raw*, unclamped ice gradient, not the bound-clamped
  copy used for material properties/reaction rate -- so the old run's
  -0.1 to -0.41 Newton-trial excursions created an artificially steep local
  gradient that injected spurious extra curvature-driven shrinkage. The
  old "full sublimation" was a numerical artifact, not real physics; the
  new, much slower rate is the physically trustworthy one.
- Restored `-dtmax` from 2.0e3 back to its pre-job64406550 default of
  1.0e4 in `solver.opts`: the original reason for capping it (dt growing
  into a full-sublimation AC singularity under the old unconstrained
  newtonls solve) is now structurally impossible under VI, regardless of
  dt. Added `inputs/experiment/21day_T-20_h0.95.opts` to test, on the HPC,
  whether the encapsulating ring fully sublimates given enough time at the
  now-correct (slower) rate.
- Note: an earlier autonomous commit on this branch (dadb7b0, already
  pushed) accidentally swept in unrelated pre-existing untracked files
  (`preprocess/igakit_complexgeocodes/random_bump_channel.*`) alongside
  the VI-solver fix -- flagged to the user, follow-up cleanup commit
  removes them from tracking (kept on disk).

## 2026-06-19 — Redesigned ice cap as concentric circle-fit encapsulation

- User clarified the actual physical intent behind the "ice cap" case,
  after every prior attempt (ellipse, conformal shell, flat layer) still
  looked wrong: the bump represents a sediment grain intruding from the
  domain boundary (modeled as a boundary bump, not a separate phase
  field, since direct phase-field modeling of sediment caused numerical
  issues months ago) and the ice around it should look like the grain is
  truly encapsulated, not coated/draped.
- User's own proposed design: fit a circle to the bump (through its two
  edges and peak), then use a second, concentric, larger circle (radius
  + thickness) as the ice region. Confirmed this is implementable with
  ZERO new C code -- it's just the existing -ice_grain_cx/cy/ax/ay circle
  mechanism with ax=ay and the center placed below the domain.
- Added `preprocess/comp_ice_encapsulation.py` implementing the sagitta/
  circular-segment fit formula. Updated
  `inputs/geometry/2D_single_bump_ice_cap.opts` (R_sed=5.0e-6,
  cy_sed=-3.0e-6, R_ice=7.0e-6 for the existing R=0.4e-5/H=0.2e-5 bump,
  thickness=0.2e-5). Verified locally: phi_ice stays in [0,1] exactly,
  uniform-thickness band wrapping the bump with no notch at the
  intrusion points, sediment grain itself renders as the correct
  no-data region below the domain boundary. Committed (0974838) and
  pushed.

---

**Session ended:** 2026-06-19 14:16:27

## 2026-06-19 — Generated first real movie; bounded highres disk cost

- User asked to actually run the new movie pipeline on job64410270's
  real output (5564 snapshots). Discovered plot_permafrost_highres.py's
  dense VTS output is ~70MB/file at n_per_elem=4 on this 608x122-element
  mesh -- converting all 5564 snapshots would be ~400GB, and even the
  ~1199 unique steps make_movie.py's default 600-frame, evenly-spaced-in-
  time sampling actually needs would be ~82GB, too close to the ~88GB
  free on disk at the time.
- Added `postprocess/select_movie_frames.py`: computes the same target
  times make_movie.py uses internally, finds the real snapshot(s)
  bracketing each one, and converts only that union via
  plot_permafrost_highres.py --steps -- run with matching --n-frames so
  the two scripts' target times line up exactly.
- Generated job64410270's first movie at --n-per-elem 2 (~20.5GB, safe
  margin) instead of the default 4: 600 frames, 1218x246, 20s @ 30fps,
  tightly cropped, with matching ice/vapor SVG colorbars. Cleaned up the
  ~21GB of intermediate highres VTS afterward (the movie itself is only
  ~105KB). Committed (a3fd664) and pushed.

---

**Session ended:** 2026-06-19 14:01:50

## 2026-06-19 — Fixed movie resolution, cropping, and added vector colorbars

- User reported `make_movie.py`'s frames were too low-res, not cropped
  tight to the geometry, and wanted colorbars as separate vector images
  for Inkscape editing rather than baked into the frames.
- Resolution now defaults to the reader's own native point grid (read
  from the data extent) instead of a fixed 1600px cap, with
  `--supersample`/`--max-width` controls.
- Root-caused the cropping bug: ParaView silently resets the camera (fit
  with padding) on the FIRST `Render()` of a newly-shown representation,
  undoing any camera settings applied earlier. Fixed by rendering once
  first to absorb that reset, then setting the explicit zero-margin
  camera afterward. Verified frames now fill edge to edge.
- Hid in-frame scalar bars; added standalone SVG colorbar export
  (`<out>_ice_colorbar.svg`/`_vapor_colorbar.svg`) reconstructed directly
  from the actual `lut.RGBPoints` so they're guaranteed to match the
  rendered frames. Also fixed an odd-pixel-dimension bug breaking
  libx264 encoding. Committed (9fdfbda) and pushed.

---

**Session ended:** 2026-06-19 13:42:10

## 2026-06-19 — ParaView movie script with linear-time playback

- Added `postprocess/make_movie.py` (run via pvpython) and
  `postprocess/make_cmocean_preset.py` (run via the regular venv, which
  has cmocean) to automate the user's manual ParaView workflow: IsoVolume
  of IcePhase in `[0.5,1.1]` colored by IcePhase (cmocean "ice" preset,
  exported to JSON once since pvpython's embedded Python lacks cmocean),
  plus a second IsoVolume in `[-0.1,0.5]` (air, since air=1-ice) colored
  by VaporDensity with an auto-computed (percentile-based, fixed for the
  whole movie) or manually set colorbar range.
- Addressed the non-linear-playback problem directly: adaptive dt makes
  raw per-snapshot frames dense early/sparse late, so a naive
  frame-per-snapshot movie crawls early and speeds up late. The script
  instead samples `--n-frames` evenly spaced in actual simulated time
  across `[t-start, t-end]`, with an optional `TemporalInterpolator` (on
  by default) so frames between sparse late snapshots blend smoothly.
  Documented in the script's docstring that this can't un-skip dense
  early dynamics under strictly linear time -- that's an inherent
  tradeoff of the request, not a bug; flagged a variable-speed
  two-segment mode as a follow-up if the user wants it later.
- Verified end-to-end against job64410270's `permafrost.pvd` (5564 real
  snapshots). Committed (5f65166) and pushed.
- Gave the user `submit_permafrost.sh` commands to rerun the three 2-day
  simple-case sanity sims (`2D_two_grains_flat`,
  `2D_single_bump_two_grains`, `2D_single_bump_ice_cap`, all still
  present and unchanged) with `2day_T-20_h0.95`.

---

**Session ended:** 2026-06-19 12:37:45

## 2026-06-19 — Diagnosed dt oscillation in job64410270; widened phase bounds

- User ran the 30-day retry (`job64410270`, `random_bumpy_floor_v2`) to let
  the small grains finish sublimating, then reported the timestep
  oscillating around a small value again. This time it wasn't a crash --
  the `dtmax`/`max_rej` fix from `job64406550` worked -- but a genuine
  limit cycle: `TOT_ICE`/`TOT_AIR`/`TEMP`/`TOT_RHOV`/`TOTAL_MASS` were all
  flat from step ~3000 onward (bulk system at steady state by t~9.8e5s,
  ~11.3 of 30 days), while one floor cap at `x=6.09e-5` had fully
  sublimated and settled at a static `phi_ice ~ -0.0004` -- never near the
  `-0.05` bound, but the *undamped* trial Newton step there swings much
  further before line search damps it, crossing `-0.05` at larger dt and
  tripping `SNESSetFunctionDomainError`. Every retry then trivially
  "converged" in 1 iteration (nothing else changing), so the growth
  heuristic immediately pushed dt back up past the same ceiling --
  oscillating ~2-6s forever instead of escaping.
- User chose to kill the job and resubmit fresh rather than let it ride
  out the plateau. (Turned out the job had already ended on its own by
  the time I asked, so no cancellation was needed.)
- Fix: widened `-phase_lo/-phase_hi` from `-0.05/1.05` to `-0.15/1.15` in
  `inputs/solver.opts` -- gives line search room to find a damped,
  in-bounds step at larger dt without changing the physical answer
  (material properties are still clamped to `[0,1]` downstream
  regardless, and every observed final accepted violation across both
  incidents has stayed under 0.3%). Committed (b735c6b) and pushed.

---

**Session ended:** 2026-06-18 19:32:52

## 2026-06-18 — Diagnosed job64406550 failure; raised dtmax/max_rej cap, more floor caps

- Investigated why `2D_random_bumpy_floor__...__job64406550` died at step
  716 of a 2-day run. `outp.txt` showed no explicit error (it just stops),
  but `TOT_ICE` flat while `I-A INTERF` steadily shrinks through steps
  650-716 pointed to Ostwald ripening, and rendering the field confirmed
  the floating ice grain at `cx=6.0e-5,cy=1.0e-5` fully sublimates away by
  step 690. dt had grown to 4-7k s during that quiet stretch; right as the
  grain finished disappearing, the AC reaction term went singular at that
  point, so no damped Newton step stayed within `[phase_lo,phase_hi]` --
  SNES failed before even reaching the convergence test. Counted exactly
  26 consecutive "Reduce time step" rejections after the last accepted
  step, one more than `-max_rej 25` -- PETSc's TS gave up and errored
  (to stderr, not the synced `outp.txt`, hence no visible error locally).
- Answered the user's question directly: tightening the SNES convergence
  criterion (rtol/atol/stol) would NOT have helped -- the failure happens
  before SNES ever reaches the convergence check, so it's a step-size/
  stiffness problem, not an accuracy problem.
- Fix in `inputs/solver.opts`: `-dtmax` 1.0e4 -> 2.0e3 (keeps dt from
  growing into the stiff regime at all), `-max_rej` 25 -> 50 (extra
  margin for the existing halving-recovery mechanism).
- Added 15 more floor-touching ice caps (cy=0, R=3.0-4.0e-6) along the
  bumpy floor per request, greedily spaced with clearance from each other
  and the existing 7 grains -- `2D_random_bumpy_floor.opts` now has 22
  ice grains total. Verified locally with `-t_final 0`. Committed
  (e17371c) and pushed.

---

**Session ended:** 2026-06-18 19:25:27


---

**Session ended:** 2026-06-18 17:16:18

## 2026-06-18 — Reverted ice cap; new full-coverage random bumpy floor geometry

- Neither the conformal shell nor the flat-height layer matched what the
  user wanted for the ice cap case — reverted
  `inputs/geometry/2D_single_bump_ice_cap.opts` to a plain `-ice_grain_*`
  ellipse (same style as the other simple cases). The `-ice_shell_*`/
  `-ice_flat_*` C-code mechanisms stay available but unused.
- Picked back up the earlier-shelved "random bumpy floor" idea: user asked
  for a new, larger geometry with the entire bottom bumpy and ice grains
  placed within the domain. Confirmed via AskUserQuestion: domain ~2x
  wider (Lx=2.0e-4 m, Ly unchanged at 4.0e-5 m), 5-8 ice grains scattered
  independently of bump positions (not one-per-bump).
- Added `generate_random_bumps()`/`validate_bump_field()` to
  `preprocess/build_geometry_multi_grain.py`: sequential left-to-right walk
  (each bump's center pulled toward the previous one by a randomized
  overlap fraction) guarantees full coverage with no flat gaps, unlike
  evenly-spaced bumps with troughs between them. Also exposed `--Lx/--Ly`
  on the CLI (previously hardcoded globals). New
  `inputs/geometry/2D_random_bumpy_floor.opts`: 24 bumps (MAX_SED_GRAINS
  cap) spanning the full 2.0e-4 m width, R in [5e-6,7e-6] m, same dx/eps
  as production (608x122 mesh). 7 ice grains of varying size (2 boundary,
  5 interior — mix of floor-touching and floating) scattered across the
  domain.
- Verified locally with `-t_final 0`: full floor coverage, no gaps, all 7
  grains at expected positions, floor grains nucleate naturally into the
  local bump shape. Committed (83f86df) and pushed.

---

**Session ended:** 2026-06-18 16:56:22

## 2026-06-18 — Flat ice layer to fully encapsulate the bump (case 3 final)

- Bumped `-ice_shell_thickness` from 0.1e-5 to 0.2e-5 (= bump's own height)
  on user feedback that the v5 distance-to-surface shell still read as too
  thin. Same shape/mechanism, just a bigger thickness.
- User then clarified the actual desired picture: a sediment grain fully
  ENCAPSULATED in ice with a FLAT (non-rounded) ice-air interface, not a
  domed conformal coating — suggested approximating it with a giant circle
  poking up from below the domain. Recognized that idea is mathematically
  equivalent to a flat height threshold, but a literal giant-radius ellipse
  would also blow up the existing tc_k=tc*sqrt(ax*ay) interface-sharpness
  scaling (ultra-sharp, non-physical transition) — implemented a dedicated
  primitive instead.
- Added `-ice_flat_x/-ice_flat_R/-ice_flat_height` (`include/NASA_types.h`,
  `src/permafrost2.c` parsing, `FormInitialMultiGrains2D` in
  `src/initial_conditions.c`): ice fills everything below an absolute
  height (unlike `-ice_shell_thickness`, which is relative to the bump's
  own surface), windowed laterally the same way as `-ice_shell_*`, using
  the same `tc=1/(sqrt(2)*eps)` scaling as everywhere else (no artificial
  sharpening). Updated `inputs/geometry/2D_single_bump_ice_cap.opts`
  (R=0.4e-5 matching the bump, height=0.4e-5 = 2x the bump's height).
- Verified locally: the bump is now fully buried, ice reaches the floor at
  both edges (no gap at the intrusion points), and the top is flat with
  rounded corners only at the lateral edges — matches the requested
  picture. Committed (91e973a) and pushed. The `-ice_shell_*` conformal-
  coating mechanism from the previous session stays in the C code,
  unused by this opts file now but available if a true coating is wanted
  elsewhere.

---

**Session ended:** 2026-06-18 16:32:44


---

**Session ended:** 2026-06-18 16:18:13

## 2026-06-18 — Iterated case-3 ice cap to a true distance-to-surface shell

- User wanted the bump (a sediment grain intruding into the domain) coated
  in roughly-constant-thickness ice, all the way down to where the grain
  meets the floor — not a cap sitting on top.
- v3 (single ellipse, tight-fit + thickness=1/2 bump height): pinches to a
  point at the bump's edges, leaving no ice near y=0 right at the
  intrusion points — looked like a frozen droplet on top of the bump, per
  explicit feedback.
- v4 (constant-thickness band, `-ice_shell_x/-R/-thickness`, with R fixed
  to match the bump's own R=0.4e-5 exactly): much closer, but the lateral
  and vertical windows are independent tanh functions multiplied together,
  so right at the bump's edge both sit at half-strength and their product
  (0.25) falls under the ice=0.5 contour — a small residual notch right at
  the intrusion point. Asked the user whether this was good enough or
  worth fixing properly; they chose to fix it.
- v5 (final): replaced the two-window product with an actual distance-to-
  the-floor-curve calculation in `FormInitialMultiGrains2D`
  (`src/initial_conditions.c`): perpendicular distance to the local
  tangent line inside the shell's footprint (using a new
  `SedimentBumpFieldDeriv()`, the analytic derivative of
  `SedimentBumpField()`), and Euclidean distance to the footprint's fixed
  endpoint outside it (naturally rounding the coating's lateral ends). The
  two pieces agree exactly at the segment boundary (slope and bump height
  both vanish there), so the result is continuous with no seam and no
  corner notch — verified locally: the ice=0.5 contour now reaches y=0 at
  both of the bump's edges. Committed (8ba2018, dc3db33) and pushed.
- Also fixed an HPC "Stale file handle" build race: 3 jobs submitted
  together all ran `make clean && make all` concurrently against the same
  shared `obj/`. Pointed the user at the existing `SKIP_COMPILE=1` escape
  hatch in `run_permafrost.sh` (build once on the login node, then
  `--export=ALL,SKIP_COMPILE=1` on each submission).

---

**Session ended:** 2026-06-18 16:03:14


---

**Session ended:** 2026-06-18 15:52:53

## 2026-06-18 — Constant-thickness ice-shell primitive for the bump-cap case

- Diagnosed an HPC "Stale file handle" build failure: 3 jobs submitted
  together each ran `make clean && make all` concurrently against the same
  shared `obj/` directory and raced. `run_permafrost.sh` already has a
  documented `SKIP_COMPILE=1` escape hatch for this; pointed the user at
  building once on the login node then resubmitting with
  `--export=ALL,SKIP_COMPILE=1`.
- All 3 sanity-case simulations ran and looked good. Follow-up request: make
  the case-3 "ice cap" a true constant-thickness shell wrapping the bump,
  not the ellipse approximation from the previous session.
- Added a new IC primitive (`-ice_shell_x/-ice_shell_R/-ice_shell_thickness`)
  to `FormInitialMultiGrains2D` (`src/initial_conditions.c`): a smooth band
  of fixed vertical thickness sitting directly on `SedimentBumpField(x)` —
  i.e. it literally follows the bump's own floor curve — windowed laterally
  via a second tanh so it only covers the bump, not the whole floor. New
  `AppCtx` fields in `include/NASA_types.h`, options parsing in
  `src/permafrost2.c` (mirrors the existing `-sed_grain_x` array pattern).
  Added on top of the existing `-ice_grain_*` ellipse summation, not a
  replacement for it elsewhere.
- Updated `inputs/geometry/2D_single_bump_ice_cap.opts` to use the new
  shell (R=0.5e-5, thickness=0.1e-5) in place of the old ellipse grain.
  Verified locally with `-t_final 0`: the ice now visibly wraps the bump at
  constant thickness, with small rounded bulges at the lateral edges
  (expected smoothing artifact of the two combined tanh windows). Committed
  (f748487) and pushed.

---

**Session ended:** 2026-06-18 14:04:32


---

**Session ended:** 2026-06-18 13:59:46

## 2026-06-18 — 3 simple sanity-case geometries; fix GrevilleAbscissae bug

- Diagnosed why the HPC run still used 44 CPUs after the eps/mesh resize: the
  resize commits (75a414e..f846c33) were never pushed to GitHub, so the
  HPC-side `git pull` had nothing new to fetch. Pushed them (commit range now
  on `origin/rewrite/2phase-from-equations`).
- Planned (plan mode) and shelved a full-coverage random-bumpy-floor geometry
  design for later; the full design is in the conversation transcript if
  needed again (the plan file itself was reused/overwritten for the plan
  below, per this session's plan-mode workflow).
- Per a request to validate simpler cases before more geometry complexity,
  planned and built 3 new sanity-case geometries, none touching the active
  production `multi_grain_test.dat`/`2D_multi_grain_test.opts` (live HPC run
  in flight): (1) `2D_two_grains_flat.opts` — flat domain, no `-geom_file`,
  4 ice grains (1 large/1 medium boundary + 2 small bottom); (2)
  `2D_single_bump_two_grains.opts` — one floor bump + 2 boundary grains; (3)
  `2D_single_bump_ice_cap.opts` — same bump + 1 boundary grain + an
  approximate ice-cap ellipse draped over the bump.
- Found and fixed a real latent bug while building case 1: `GrevilleAbscissae()`
  in `src/initial_conditions.c` returned un-normalized parametric values.
  This happened to work for `-geom_file` meshes (whose knots already span
  [0,1] by the igakit builder's convention) but silently broke for the
  default (`-Nx/-Ny`, no `-geom_file`) axis, whose knots span `[0,Lx]/[0,Ly]`
  directly — `FormInitialMultiGrains2D` multiplied by Lx/Ly a second time,
  collapsing every sample point into a tiny corner near (0,0) and leaving
  ice phase exactly 0 everywhere. `multi_grains` had only ever been run with
  `-geom_file` before, so this was never exercised. Fixed by normalizing by
  the knot span in `GrevilleAbscissae()` (commit 69dd1fa).
- Added a `--bumps "cx,R,h;..."` CLI override to
  `preprocess/build_geometry_multi_grain.py` so a one-off bump layout (the
  shared single-bump mesh for cases 2/3) doesn't require hand-editing the
  script's hardcoded `SEDIMENT_GRAINS` list; default behavior unchanged.
- Verified all 3 cases locally with `-t_final 0` + `plot_permafrost_highres.py`;
  the ice-cap ellipse needed one retuning pass (ax/ay/cy) after the first
  attempt floated visibly above the flat floor away from the bump. Committed
  (8fabe59).

---

**Session ended:** 2026-06-18 13:16:00


---

**Session ended:** 2026-06-17 18:49:48


---

**Session ended:** 2026-06-17 18:45:22


---

**Session ended:** 2026-06-17 18:40:48

## 2026-06-17 — Diagnose HPC crash, resize mesh via comp_eps.py, rescale eps

- Diagnosed job64387928 (eps_model75, 600x240, eps=1.8659e-07) crash:
  `TSStep failed due to DIVERGED_STEP_REJECTED` -> MPI_Abort (exit 137).
  Last accepted step (STEP 69, t~1.11e4 s) showed phase bounds widening
  to [-0.0255, 1.0255] just before the failure cascade as dt grew past
  ~1.3e3 s; every retry showed identical it=0 fnorm (SNES "converging"
  trivially) while the bounds/rollback check kept rejecting anyway,
  geometrically shrinking dt until `-max_rej 25` was exhausted. Root
  cause: `eps_model = 0.75*eps` (introduced to fix interface growth)
  made the AC reaction term too stiff at the old, manually-shrunk
  eps=1.8659e-07, on top of pre-existing bump-flank mesh shear.
- Re-derived eps and mesh size from `comp_eps.py` for the actual domain
  (Lx=1.0e-4, Ly=4.0e-5, Rave=3.5e-6, T0=-20, --n 4): eps=4.6648e-07
  (full Kaempfer & Plapp bound, not manually shrunk), Nx=304, Ny=122
  (down from 600x240) -- confirmed the earlier 600x240 mesh size was
  driven by a ParaView control-point-rendering artifact, not a real
  resolution requirement (see prior session's plot_permafrost_highres.py).
- Added `--C` continuity option to `build_geometry_multi_grain.py`
  (previously hardcoded to C=P-1); used to verify a P=2/C0 mesh, then
  rebuilt the production mesh at the intended P=2/C=1.
- Ran a local `-t_final 0` sanity check on the new 304x122 mesh (commit
  75a414e), converted with both `plotpermafrost.py` and
  `plot_permafrost_highres.py`, confirmed clean IC. Discovered
  `monitoring.c` reads the output directory from a `folder` env var
  (set by `run_permafrost.sh`'s `export folder`), not `-output_path`
  directly -- not a bug, just missed when invoking the binary by hand.
- User measured ~9.5 elements across the 1%-99% diffuse band on this
  mesh in ParaView (target ~7.5); rescaled eps by 7.5/9.5 to
  eps=3.68274e-07 (commit 423a051), still below comp_eps.py's ceiling.
- Archived the mesh/opts as `inputs/geometry/multi_grain/compEps_304x122/`
  via `build_geometry_multi_grain.py --variant`.
- Handed off the HPC submission command
  (`./scripts/HPC/submit_permafrost.sh 2D_multi_grain_test 2day_T-20_h0.95 compEps304x122`)
  for a full t_final=172800 (2-day) run to check whether the new mesh/eps
  avoids the dt-rejection-cascade crash and keeps the interface width
  stable over time.

---

**Session ended:** 2026-06-17 18:34:04


---

**Session ended:** 2026-06-17 18:27:20


---

**Session ended:** 2026-06-17 18:03:09


---

**Session ended:** 2026-06-17 17:46:12


---

**Session ended:** 2026-06-17 17:35:35


---

**Session ended:** 2026-06-17 17:24:06


---

**Session ended:** 2026-06-17 17:18:08


---

**Session ended:** 2026-06-17 16:48:15


---

**Session ended:** 2026-06-17 16:46:21


---

**Session ended:** 2026-06-17 16:33:40


---

**Session ended:** 2026-06-17 16:26:53


---

**Session ended:** 2026-06-17 16:13:28


---

**Session ended:** 2026-06-17 16:12:40


---

**Session ended:** 2026-06-17 16:06:44


---

**Session ended:** 2026-06-17 15:34:34


---

**Session ended:** 2026-06-17 15:20:47


---

**Session ended:** 2026-06-17 12:41:03


---

**Session ended:** 2026-06-17 11:56:44


## 2026-06-17 — Fix dt stall: lower SNES convergence guard; submit BL Phase 1 job

- Diagnosed root cause of dt barrier at ~2.4e-02 s in job64369939 (600×240 5-bump):
  - BT line search at Newton step 0 backtracks due to transient phi < phase_lo near trough
    grain floor contacts; lands at ||G|| ~ 2e-17 (machine precision) with lambda very small.
  - `SNESDOFConvergence` guard `it_number < 3` prevented convergence declaration at it=1,
    forcing a 2nd Newton KSP solve.
  - At dt >= 3.162e-02 s, the TSALPHA Jacobian J = M/(alpha_m*dt) + alpha_f² dF/dX1 is
    insufficiently regularised: KSP (GMRES, max_it=500) diverges at the 2nd solve.
  - `SNESCheckKSPSolve` exits `SNESSolve_NEWTONLS` with SNES_DIVERGED_LINEAR_SOLVE before
    `SNESConverged` is called; TSAdaptCheckStage rejects → cascade back to dt=2.371e-02.
- Fix: changed guard to `it_number < 1` in `src/snes_convergence.c`. The DOF convergence
  criteria (rtol/atol/stol) still cannot fire at it=0 (rel[i]=1 by construction), so no
  spurious early convergence is possible.
- Secondary hardening: increased `ksp_max_it 500 → 1000` and `sub_pc_factor_levels 2 → 3`
  in `inputs/solver.opts` to better handle intermediate-dt ill-conditioned Jacobians.
- Committed both changes together: 4d33a89.

---

**Session ended:** 2026-06-17 09:51:10


---

**Session ended:** 2026-06-17 09:42:33


## 2026-06-17 — Phase 1 BL y-refinement implemented

- Diagnosed root cause of dt stall across multiple 12-bump runs: touching-support bumps create
  continuous distortion zones → constant phase-bound violations → dt can't grow past ~1 s.
- Confirmed fix: 5 isolated bumps with 1.2e-5 m gaps (proven stable in job64345061). Reverted
  Lx=1.0e-4 domain to 5-bump geometry.
- Analysed colleague's igakit mesh scripts (MeshGenerate.py, test2_blayer_ref.py,
  nsch_axisym_solidtopextra_cosine.py) — key technique: non-uniform knot refinement via
  surf.refine() or direct Uy construction to densify elements near boundary.
- Implemented Phase 1 BL refinement in preprocess/build_geometry_multi_grain.py:
  - Replaced open_uniform_knots(Ny, P) with a two-zone non-uniform Uy built directly.
  - BL zone [0, 8 µm]: 96 elements (h_BL = 8.33e-8 m, 2× finer than bulk).
  - Bulk zone [8–40 µm]: 192 elements (h_bulk = 1.667e-7 m, unchanged).
  - eps/h_BL = 5.6 in grain zone (up from 2.8), eps unchanged at thermodynamic bound.
- Regenerated multi_grain_test.dat (602×290 control points; +694 kB).
- Updated DOF_GRID comment in 2D_multi_grain_test.opts to 602 290.
- Committed and pushed: 555f173.

---

**Session ended:** 2026-06-17 09:23:20


---

**Session ended:** 2026-06-17 09:13:11


---

**Session ended:** 2026-06-16 22:42:04


---

**Session ended:** 2026-06-16 22:23:33


---

**Session ended:** 2026-06-16 22:14:48


---

**Session ended:** 2026-06-16 19:10:35


---

**Session ended:** 2026-06-16 19:08:17


---

**Session ended:** 2026-06-16 19:03:29


---

**Session ended:** 2026-06-16 17:34:04


## 2026-06-16 — Add non-symmetric ceiling bumps to channel geometry

- Added `TOP_GRAINS` (6 ceiling bumps of varied size/position, offset from floor bumps) to `build_geometry_multi_grain.py`; added `_top_bump_field()` and two-sided y-mapping in `build_surface()`.
- Added `TopBumpField()` in `initial_conditions.c`; updated y-mapping in both IC functions to `y = y_bot + v*(Ly - top_bump - y_bot)`.
- Added `n_top_grains`/`top_grain_x/R/h` to `NASA_types.h` and `-top_grain_x/-top_grain_R/-top_grain_h` CLI to `permafrost2.c`.
- Updated `2D_multi_grain_test.opts` with ceiling bump parameters; regenerated `multi_grain_test.dat`.
- Build clean; pushed to remote.

---

**Session ended:** 2026-06-16 17:31:50


---

**Session ended:** 2026-06-16 17:26:42


## 2026-06-16 — Complex channel geometry: Lx=1.0e-4, elliptical ice grains

- Added `ice_grain_ax[200]` / `ice_grain_ay[200]` fields to `AppCtx` in `NASA_types.h`; bumped `MAX_SED_GRAINS` 16→24.
- Added `-ice_grain_ax` / `-ice_grain_ay` CLI options in `permafrost2.c`; default to `radius[k]` for backward compatibility.
- Rewrote grain distance loop in `FormInitialMultiGrains2D` to use elliptical normalized distance `sqrt((dx/ax)^2+(dy/ay)^2)` with `tc_k = tc*sqrt(ax*ay)` scaling; circular grains (ax=ay=R) produce identical output to the old code.
- Scaled domain: `Lx 4.0e-5 → 1.0e-4 m`, `Nx 240 → 600` (h unchanged at 1.667e-7 m); 12 alternating-height bumps (was 5 uniform); 8 ice grains (was 6): 2 circular boundary + 6 elliptical trough puddles (ax=4.5e-6, ay=2.5e-6).
- Regenerated `inputs/geometry/multi_grain_test.dat` (602×242 control points).
- Build clean (no errors).

---

**Session ended:** 2026-06-16 16:53:57


---

**Session ended:** 2026-06-16 16:46:58


---

**Session ended:** 2026-06-16 16:38:51


---

**Session ended:** 2026-06-16 16:34:16


---

**Session ended:** 2026-06-16 15:48:48


---

**Session ended:** 2026-06-16 15:02:59


---

**Session ended:** 2026-06-16 14:54:55


---

**Session ended:** 2026-06-16 13:56:18


---

**Session ended:** 2026-06-16 12:33:15


---

**Session ended:** 2026-06-16 12:13:37


---

**Session ended:** 2026-06-16 12:11:53


---

**Session ended:** 2026-06-16 08:07:29


---

**Session ended:** 2026-06-16 07:59:09


---

**Session ended:** 2026-06-16 07:44:33


---

**Session ended:** 2026-06-16 07:33:59


---

**Session ended:** 2026-06-15 22:21:26


---

**Session ended:** 2026-06-15 22:19:15


---

**Session ended:** 2026-06-15 21:42:28


---

**Session ended:** 2026-06-15 19:09:23


---

**Session ended:** 2026-06-15 19:05:51


---

**Session ended:** 2026-06-15 17:02:34


---

**Session ended:** 2026-06-15 17:01:52


## 2026-06-15 — Mesh rebuild, basis revert, HPC prep

- Diagnosed `-geom_file` + `-Nx/-Ny` conflict: `-geom_file` overrides `-Nx/-Ny` entirely; user was setting -Nx 488 which was ignored because the `.dat` file has 240x240 baked in.
- Rebuilt P=1 mesh at 480x480 (doubled from 240) with eps halved to 9.72e-08.
- Reverted to P=2/C1 basis functions (P=1 interfaces are inherently jagged regardless of resolution).
- Resized P=2 mesh to 122x122 elements (comp_eps.py `--n 4` minimum) with eps=4.6648e-07; regenerated `multi_grain_test.dat`.
- Fixed HPC `run_permafrost.sh` `compute_optimal_nprocs()` to read `# DOF_GRID:` comment and `-dof` from solver.opts (same fix applied earlier to Studio version; without it NPROCS collapses to 1 for -geom_file meshes).
- Added `dt_ramp_test.opts` experiment file (20 s run for dt-ramp diagnostics).

---

**Session ended:** 2026-06-15 16:52:40


---

**Session ended:** 2026-06-15 16:33:13


---

**Session ended:** 2026-06-15 16:30:19


---

**Session ended:** 2026-06-15 16:05:58


---

**Session ended:** 2026-06-15 15:52:19


---

**Session ended:** 2026-06-15 15:48:55


---

**Session ended:** 2026-06-15 15:42:58


---

**Session ended:** 2026-06-15 15:16:44


---

**Session ended:** 2026-06-15 15:13:48


---

**Session ended:** 2026-06-15 15:11:17


---

**Session ended:** 2026-06-15 14:48:06


---

**Session ended:** 2026-06-15 14:43:31


---

**Session ended:** 2026-06-15 14:43:16


---

**Session ended:** 2026-06-15 14:42:59


---

**Session ended:** 2026-06-15 14:42:51


---

**Session ended:** 2026-06-15 14:42:29


---

**Session ended:** 2026-06-15 14:42:22


---

**Session ended:** 2026-06-15 14:42:12


---

**Session ended:** 2026-06-15 14:41:51


---

**Session ended:** 2026-06-15 14:41:35


---

**Session ended:** 2026-06-15 14:41:23


---

**Session ended:** 2026-06-15 14:40:59


---

**Session ended:** 2026-06-15 14:40:42


---

**Session ended:** 2026-06-15 14:40:28


---

**Session ended:** 2026-06-15 14:40:08


---

**Session ended:** 2026-06-15 14:28:05


---

**Session ended:** 2026-06-15 14:23:27


---

**Session ended:** 2026-06-15 14:23:10


---

**Session ended:** 2026-06-15 14:22:48


---

**Session ended:** 2026-06-15 14:22:38


---

**Session ended:** 2026-06-15 14:22:19


---

**Session ended:** 2026-06-15 14:22:01


---

**Session ended:** 2026-06-15 14:21:46


---

**Session ended:** 2026-06-15 14:21:43


---

**Session ended:** 2026-06-15 14:21:09


---

**Session ended:** 2026-06-15 14:20:49


---

**Session ended:** 2026-06-15 14:20:30


---

**Session ended:** 2026-06-15 14:20:13


---

**Session ended:** 2026-06-15 14:19:58


---

**Session ended:** 2026-06-15 14:19:42


---

**Session ended:** 2026-06-15 14:19:23


---

**Session ended:** 2026-06-15 14:19:04


---

**Session ended:** 2026-06-15 14:18:50


---

**Session ended:** 2026-06-15 14:18:37


---

**Session ended:** 2026-06-15 14:18:35


---

**Session ended:** 2026-06-15 14:17:57


---

**Session ended:** 2026-06-15 14:17:37


---

**Session ended:** 2026-06-15 14:17:22


---

**Session ended:** 2026-06-15 14:17:13


---

**Session ended:** 2026-06-15 14:17:11


---

**Session ended:** 2026-06-15 14:16:44


---

**Session ended:** 2026-06-15 14:16:25


---

**Session ended:** 2026-06-15 14:16:08


---

**Session ended:** 2026-06-15 14:15:57


---

**Session ended:** 2026-06-15 14:15:17


---

**Session ended:** 2026-06-15 14:14:59


---

**Session ended:** 2026-06-15 14:14:44


---

**Session ended:** 2026-06-15 14:14:29


---

**Session ended:** 2026-06-15 14:14:12


---

**Session ended:** 2026-06-15 14:13:51


---

**Session ended:** 2026-06-15 14:13:35


---

**Session ended:** 2026-06-15 14:13:20


---

**Session ended:** 2026-06-15 14:12:59


---

**Session ended:** 2026-06-15 14:12:41


---

**Session ended:** 2026-06-15 14:12:22


---

**Session ended:** 2026-06-15 14:12:03


---

**Session ended:** 2026-06-15 14:11:48


---

**Session ended:** 2026-06-15 14:11:33


---

**Session ended:** 2026-06-15 14:11:17


---

**Session ended:** 2026-06-15 14:11:04


---

**Session ended:** 2026-06-15 14:10:04


---

**Session ended:** 2026-06-15 14:09:11


---

**Session ended:** 2026-06-15 14:08:50


---

**Session ended:** 2026-06-15 14:08:30


---

**Session ended:** 2026-06-15 14:08:18


---

**Session ended:** 2026-06-15 14:08:06


---

**Session ended:** 2026-06-15 14:07:57


---

**Session ended:** 2026-06-15 14:07:50


---

**Session ended:** 2026-06-15 14:07:44


---

**Session ended:** 2026-06-15 14:07:38


---

**Session ended:** 2026-06-15 14:07:35


---

**Session ended:** 2026-06-15 14:07:25


---

**Session ended:** 2026-06-15 14:07:19


---

**Session ended:** 2026-06-15 14:07:10


---

**Session ended:** 2026-06-15 14:06:52


---

**Session ended:** 2026-06-15 14:06:37


---

**Session ended:** 2026-06-15 14:06:28


---

**Session ended:** 2026-06-15 14:06:24


---

**Session ended:** 2026-06-15 14:06:21


---

**Session ended:** 2026-06-15 14:05:38


---

**Session ended:** 2026-06-15 13:57:05

## 2026-06-15 — Greville-abscissa IC mapping, p=2/C1 geometry, eps/mesh fix, dt start

- Added a degree-agnostic `GrevilleAbscissae()` helper in
  `src/initial_conditions.c` and switched `FormInitialMultiGrains2D` to map
  DOFs to physical (x,y) via Greville abscissae (works for any (P, C^{P-1})).
- Rewrote `preprocess/build_geometry_multi_grain.py` to build the multi-grain
  geometry directly as degree-(2,2)/C1 NURBS via Greville-abscissa control
  points (no more igakit `cad.ruled`/`refine_surface`); regenerated
  `inputs/geometry/multi_grain_test.dat`.
- Worked through eps/mesh-resolution relationship: eps is a fixed physical
  bound from `comp_eps.py` (4.6648e-07 for this domain/T), independent of
  mesh; two earlier attempts to rescale eps with mesh size (2.9155e-07,
  3.5569e-07) were superseded. Final: eps fixed at 4.6648e-07, mesh refined
  to 240x240 (n~8 elements across the diffuse interface) to fix a jagged
  ice-air interface seen at 160x160.
- Changed `-delt_t` in `inputs/geometry/2D_multi_grain_test.opts` from
  1.0e-8 (unknown origin, ~51 ramp-up steps/~8.5 min on the 240x240 mesh)
  to 1.0e-4 (historical minimum starting dt); validated clean via smoke
  test with t_final=5e-4.
- Launched `mesh240_2day` full run (240x240, degree(2,2)/C1, eps=4.6648e-07,
  old delt_t=1e-8) — left running to completion; future runs benefit from
  the new delt_t=1e-4 default.

---

**Session ended:** 2026-06-15 13:48:00


---

**Session ended:** 2026-06-15 13:42:20


---

**Session ended:** 2026-06-15 13:40:40


---

**Session ended:** 2026-06-15 13:21:43


---

**Session ended:** 2026-06-15 13:13:15


---

**Session ended:** 2026-06-15 13:07:38


---

**Session ended:** 2026-06-15 12:45:01


---

**Session ended:** 2026-06-15 12:40:39


---

**Session ended:** 2026-06-15 12:31:48


---

**Session ended:** 2026-06-15 12:30:47


---

**Session ended:** 2026-06-15 12:19:44


---

**Session ended:** 2026-06-15 12:18:38


---

**Session ended:** 2026-06-15 12:09:33


---

**Session ended:** 2026-06-15 11:58:43

## 2026-06-15 — Validate refined multi-grain mesh (160x160) over 2 days

- Re-ran `2D_multi_grain_test` + `2day_T-20_h0.95` (tag `multigrain_refine_2day`)
  with the 160x160 mesh and enlarged LHS ice grain (R=6.0e-6) from commit
  322d880.
- 202 steps to t=1.74e5s, Total mass Delta=-0.001%, no [ABORT]/NaN.
- Snapshots (step 141 and final step 202, via fixed plot2D_snapshot.py)
  show smooth circular ice-air interfaces -- the jaggedness reported on the
  old 80x80/R=0.9375e-6 setup is resolved. Small grains (LHS boundary +
  2 trough nuclei) are nearly fully sublimated by t~1.74e5s via Ostwald
  ripening, leaving only the large RHS grain (R=1.5e-5).
- Also reverted 2day_T-20_h0.95.opts to every-step (-outp 1) output for
  future runs, and expanded run_permafrost.sh's header comment with
  instructions for adding new experiments/geometries.

---

**Session ended:** 2026-06-15 11:53:07


---

**Session ended:** 2026-06-15 11:50:28


---

**Session ended:** 2026-06-15 11:45:02


---

**Session ended:** 2026-06-15 11:33:23


---

**Session ended:** 2026-06-15 11:32:01

## 2026-06-15 — Generalize sediment bumps and ice IC to N grains (multi-grain geometry)

- Added `-sed_grain_x`/`-sed_grain_R` arrays (summed via new
  `SedimentBumpField()`) so the bottom-edge bump geometry supports N
  sediment grains instead of a single `-geom_bump_R` hump. `n_sed_grains==0`
  preserves the original single-bump behavior (`two_ice_grains_boundary`
  unchanged).
- Added `-ice_grain_cx`/`-ice_grain_cy`/`-ice_grain_R` arrays and new
  `FormInitialMultiGrains2D` (`-ic_type multi_grains`), placing N ice grains
  via summed/clamped tanh distance profiles, reusing the `cent`/`radius`
  AppCtx fields.
- New `preprocess/build_geometry_multi_grain.py` generates a 3-sediment-grain
  geometry (`inputs/geometry/multi_grain_test.dat`) with 4 ice grains
  (2 boundary grains + 2 nucleating in the inter-grain troughs), config in
  `inputs/geometry/2D_multi_grain_test.opts`.
- Smoke test (`smoke_short`, t_final=1e-6): no `[ABORT]` bound violations,
  grains placed at expected coordinates, mass conserved to ~0.000%.
- Ran a 2-day validation run (`2day_T-20_h0.95`, tag `multigrain_2day`,
  135 steps to t=1.74e5s): Total mass Δ=+0.0031%, no `[ABORT]` violations,
  timestep history clean (geometric growth with a small dip/recovery near
  t~1e4s, as seen on the single-bump geometry). Multi-grain mechanism
  validated end-to-end.

---

**Session ended:** 2026-06-15 11:17:36


---

**Session ended:** 2026-06-15 11:04:51

## 2026-06-15 — Verify fixes with fresh 2-day sediment-grain run (sedgeom_recheck)

- Re-ran `2D_sediment_grain_test` + `2day_T-20_h0.95` (tag `sedgeom_recheck`)
  with the updated `monitoring.c`/plotting scripts.
- `mass.png` now matches `outp.txt` exactly: Total mass Δ = -4.87e-04%,
  Ice Δ ≈ -4.9e-04% (effectively conserved), Vapor Δ = +5.05% (but vapor
  mass is ~6.4e-13 kg/m, negligible in absolute terms).
- `timestep.png` is a clean monotonic geometric dt growth (3e-8s → 6.05e3s,
  133 steps, t up to 1.74e5s) -- both post-processing bugs confirmed fixed.

---

**Session ended:** 2026-06-15 11:00:46

## 2026-06-15 — Fix plot_timestep.py and plot_mass.py post-processing bugs

- `plot_timestep.py`'s `_DOMAIN_NFIELDS = 10` was stale (left over from the
  old 4-phase 9-column header). The current domain-integral rows have 8
  fields after STEP, but SNES/Newton iteration rows coincidentally also have
  10 fields, so the regex was matching those instead -- producing a chaotic
  dt-vs-t plot (t spanning ~1e-20 to ~1e-8 s, "steps: 3") that looked like
  the simulation moved backwards in time. Fixed to 8; replotting
  `sedgeom_icfix_2day` now gives a correct monotonic geometric dt growth
  from ~3e-8 s to ~6e3 s over 133 steps, t up to ~1.7e5 s.
- `plot_mass.py` previously recomputed mass from `sol_*.dat` snapshots using
  `dV = (bounding-box area)/(N control points)`, which overcounts the area
  for the non-rectangular sediment-grain bump geometry (bounding box
  includes area under the bump) and ignores non-uniform cell volumes. This
  produced a spurious **+4.77% mass "gain"**, contradicting `outp.txt`'s
  IGA-quadrature `TOTAL_MASS` (-0.000%). There was also a separate time-axis
  bug: sparse `sol_*.dat` output (43 files) was positionally paired with the
  first 43 rows of `SSA_evo.dat` (134 rows), truncating the plotted time
  range to ~0.0053s instead of the real ~1.7e5s.
- `monitoring.c` now appends `tot_air`, `tot_rhov`, `tot_mass` (already
  computed via `IGAComputeScalar`) to every row of `SSA_evo.dat`.
  `plot_mass.py` now prefers these per-step quadrature-correct values
  directly (matching `outp.txt` exactly) for new runs; for older runs
  without these columns it falls back to the `sol_*.dat` approach but
  matches each snapshot to its `SSA_evo.dat` row by step number (parsed from
  the filename) instead of positional slicing.
- Confirmed for `sedgeom_icfix_2day` (old SSA format, fallback path): fixing
  the time axis alone did NOT change the +4.77% figure, confirming the dV
  approximation (not the time axis) is the source of the discrepancy --
  `outp.txt`'s -0.000% is the trustworthy number. A fresh run with the
  updated `monitoring.c` is needed to get a corrected `mass.png` for this
  geometry.

---

**Session ended:** 2026-06-15 10:45:59

## 2026-06-15 — Fix non-circular IC grain on sediment-grain geometry; resolves ~3% mass drift

- Confirmed via direct derivation that `R_vap`'s `(rho_ice - rhov)*ice_t`
  term already contains the product-rule `rho_v*air_t` piece (since
  `air_t = -ice_t`) -- no missing term in `src/assembly.c`.
- Root-caused the ~3% TOTAL_MASS increase seen in the prior 2-day
  sediment-grain run (commit `051a4cc`): `FormInitialTwoIceGrainsBoundary2D`
  placed grain centers/circles using index-ratio coordinates
  `(Lx*i/mx, Ly*j/my)`, which only equals the true physical (x,y) for a
  flat rectangular domain. On the `-geom_file` bump geometry, the real
  map is `y_phys = bump(x) + v*(Ly - bump(x))`, so the larger grain
  (R1=1.875e-5, which extends into the bump's x-range) came out visibly
  non-circular -- an unphysical IC that the solver then had to relax away
  over the first ~100 steps, injecting the mass drift.
- Fixed `src/initial_conditions.c`: added `SedimentBump()` (duplicate of
  `build_geometry_sediment_grain.py`'s `_bump()`) and a new
  `-geom_bump_R` option (`include/NASA_types.h`, `src/permafrost2.c`,
  default 0 = flat domain, no behavior change for other geometries). The
  IC now computes each node's true physical (x,y) via the bump map before
  building the tanh distance fields (commit `00fd2dd`). Set
  `-geom_bump_R 1.0e-5` in `inputs/geometry/2D_sediment_grain_test.opts`
  to match the geometry's `R_sed`.
- Verified: smoke test (`smoke_short`) shows both grains now circular
  (clipped only at domain edges as intended) with clean bounds
  `phi_ice/phi_air in [0,1]` at step 0. Re-ran the full 2-day case
  (`2day_T-20_h0.95`, tag `sedgeom_icfix_2day`): 133 steps to
  t=1.739e5s, TOT_ICE -0.000%, TOT_AIR +0.000%, **TOTAL_MASS -0.000%**
  (down from +~3% with the old IC), TOT_RHOV +5.05% (consistent with
  prior runs). A few transient bounds excursions around steps 88-111
  (phi_air up to ~1.37) recovered via the existing rollback/dt-halving
  logic; final bounds clean at [0,1]. The IC fix alone resolved the mass
  drift -- no further mesh/basis-order tuning needed.

---

**Session ended:** 2026-06-15 10:37:09


---

**Session ended:** 2026-06-15 10:30:34


---

**Session ended:** 2026-06-15 10:28:42

## 2026-06-15 — 2-day run on sediment-grain geometry; fix stale opts comment

- Fixed a stale comment in `inputs/geometry/2D_sediment_grain_test.opts`
  claiming degree (2,1)/~83x80 elements from arc-length refinement --
  the current geometry (from the smooth-bump rewrite, commit `c4098ad`)
  is degree (1,1), C0, 80x80 elements (commit `9d8e0c5`).
- Explained why grain 1 (the larger ice grain, RCice1=1.875e-5) looks
  non-circular in the geometry smoke test: `FormInitialTwoIceGrainsBoundary2D`
  places grains using index-ratio coordinates `(Lx*i/mx, Ly*j/my)`, treating
  the parametric grid as uniform physical space. But the actual map is
  `y_phys = bump(x) + v*(Ly - bump(x))`, which compresses y near the bump
  (x in [1.0e-5, 3.0e-5]). Grain 1 extends from x=Lx inward to x~2.125e-5,
  overlapping the bump region, so part of its "circle" gets squashed --
  grain 0 (smaller, R0=9.375e-6) stays just outside the bump region and
  looks circular. Purely an IC/smoke-test artifact, not a physics issue.
- Ran a 2-day (-20C/95% humidity) simulation on the sediment-grain geometry
  (`2D_sediment_grain_test` + `2day_T-20_h0.95`, tag `sedgeom_2day`):
  completed 141 steps to t=1.728e5s (full target), 44 VTK snapshots.
  TOT_ICE +0.006%, TOT_AIR -0.005%, TOTAL_MASS +0.006%, TOT_RHOV +5.06%
  (consistent with prior rectangular-domain runs). Bounds stayed close to
  [0,1] (phi_ice in [-0.0019, 1.0044], phi_air in [-0.0044, 1.0019]) for
  the whole run. A transient bounds excursion around steps 117-120
  (phi_air down to ~-0.44, likely AC-interface overshoot interacting with
  the distorted mesh near the bump) was recovered automatically by the
  existing rollback/dt-halving logic, and the run finished cleanly at
  dtmax=8.625e3.

---


---

**Session ended:** 2026-06-15 10:19:56

## 2026-06-15 — Fix geometry fold with smooth-bump ruled-graph patch

- Verified numerically that the prior Coons-patch geometry had
  ~0.16% of its parameter domain with det(J)<0, localized at the two
  reflex corners where the sharp semicircular bite met the bottom
  edge -- and that degree elevation (-p 2/-C 1) cannot fix this, since
  elevation reproduces the identical geometric map/Jacobian field with
  more control points.
- Replaced the line-arc-line + `cad.coons` construction with a smooth
  C-infinity "bump" function y=g(x) for the bottom boundary (zero, with
  all derivatives vanishing, outside |x-Lx/2|<R_sed -- no cusps/reflex
  corners). Both bottom and top curves share the same u<->x
  parametrization, so `cad.ruled(bottom, top)` gives vertical
  v-isolines and a Jacobian that is positive everywhere (checked: min
  detJ > 0 over the whole parameter domain). Result is degree (1,1)
  with C0 interior knots, matching `solver.opts`' `-p 1 -C 0` (commit
  `c4098ad`).
- Re-ran the `2D_sediment_grain_test` + `smoke_short` smoke test: now
  completes (13 steps), TOT_ICE/TOTAL_MASS drift 0.000%, bounds clean
  (phi_ice/phi_air in [0,1]). Multi-patch geometry is NOT needed for
  this case.
- Avoided rerunning the doomed `-p 2/-C 1` smoke test (would have
  reproduced the same fold for the reason above).

---

## 2026-06-15 — IGARead integration smoke test; geometry folds at reflex corners

- Confirmed `-p 1 -C 0` is intentional/correct (was raised from p=2/C=1
  only during earlier debugging); restored it in `inputs/solver.opts`
  (commit `85265d1`, after an earlier mistaken commit/revert cycle
  `1bbd825`/`b7429d5` while waiting for confirmation).
- `src/permafrost2.c`: added `-geom_file` option that loads an
  igakit-generated `.dat` geometry via `IGARead`, overriding
  `-p/-C/-Nx/-Ny/-Nz` axis setup. Per-axis degree is now read back from
  the IGA (`IGAAxisGetDegree`) to size the `alph`/`mob` quadrature-point
  arrays, since a loaded geometry's degree can differ from `-p`.
- `preprocess/build_geometry_sediment_grain.py`: switched to
  arc-length-based knot refinement (vs. uniform-in-parameter) and added
  a legacy-VTK structured-grid export
  (`preprocess/sediment_grain_geometry.vtk`) for ParaView inspection.
- Smoke test (`2D_sediment_grain_test` + `smoke_short`, commit `85e50f0`):
  IGARead + IC + assembly setup all ran correctly, but the solve failed
  with "Non-positive det(Jacobian)" -- the single-Coons-patch geometry
  folds near the two reflex corners where the semicircular bite meets
  the bottom edge. This is a geometry-construction issue (a non-convex
  domain can't be validly mapped by one tensor-product patch), not a
  solver bug. Likely fix is multi-patch geometry (new infrastructure,
  not yet started) or a different geometric approach -- discussed with
  user, decision deferred to a future session.

---

## 2026-06-14 — 2-day Ostwald run completes; start igakit geometry exploration

- Ran the 2-day -20C/95% humidity Ostwald-ripening case
  (`2D_two_ice_grains_boundary` + `2day_T-20_h0.95` + `ostwald_v1`) via
  `run_permafrost.sh`, with time-equally-spaced output
  (`-outp 0 -n_out 100`, new `inputs/experiment/2day_T-20_h0.95.opts`).
  Completed cleanly in 134 steps, final t=1.738e5s (target 1.728e5s):
  TOT_ICE/TOTAL_MASS drift -0.000%, TOT_RHOV +5.15% (consistent with the
  earlier 1-day run, not a concern), dt grew up to 5.3e3s. Only 39
  snapshots were actually written (not 100) because adaptive dt grew
  past t_interv near the end -- acceptable for now, flagged as a future
  refinement if denser late-time output is needed.
- Found and fixed an uncommitted, previously-flagged change to
  `inputs/solver.opts` (`-p 2/-C 1` -> `-p 1/-C 0`) that had never been
  confirmed with the user. First mistakenly committed it with an
  unverified justification (`1bbd825`); reverted that part back to
  `-p 2/-C 1` in a follow-up commit (`b7429d5`) pending clarification.
- Started exploring igakit's NURBS/CAD toolbox (`cad.line`, `cad.circle`,
  `cad.join`, `cad.coons`, `igakit.io.PetIGA`) for building non-square
  domains (rectangle with a semicircular "sediment grain" bite). Wrote
  a standalone prototype `preprocess/build_geometry_sediment_grain.py`
  that constructs such a geometry via line/arc/line `join()` + `coons()`,
  plots it, and writes a PetIGA `.dat` geometry file. The boundary shape
  is correct but the element grid is badly distorted near the arc
  (uniform parameter-space knot refinement doesn't give uniform physical
  element sizes across the line/arc/line transition) -- next step is
  arc-length-based knot insertion. No changes to `permafrost2.c` yet;
  that integration (switching from `IGAAxisInitUniform` to `IGARead`,
  plus BCs on the curved sediment edge) is a separate follow-up.

---

## 2026-06-14 — Fix latent-heat coefficient: rho_ice instead of mixture rho

- Local 1-day test run (via `run_permafrost.sh`, 2-grain boundary IC) of
  the S_sub->ice_t fix showed TEMP integral dropping from -4.189e-08 to
  -5.086e-08 (~4.3 degC domain-average cooling) within just t~30s, and
  adaptive dt stalling near 1-2s -- both signs of a spurious energy source.
- Root cause: R_tem's latent-heat term used the local mixture density
  `rho(phi_i)` (919 in ice, 1.34 in air) for `-rho*lat_sub*ice_t`, while
  R_vap's mass-exchange term correctly uses the constant `rho_ice`
  (rho_SE per Moure & Fu 2024 eq.9). At the diffuse interface these
  differ by up to ~700x, injecting far more "latent heat" than the mass
  actually exchanged warranted. Since thermal diffusion equilibrates
  this ~40um domain in ~1e-4s, any real net energy injection should show
  up almost immediately as a near-uniform T shift -- the 4 degC/30s
  drift was that shift.
- Fixed `src/assembly.c`: R_tem and Jacobian_A1 [tem,ice] now use
  `rho_ice * lat_sub * ice_t` (commit `fbde305`).
- Re-ran the same local test: TEMP held at exactly -4.189e-08 (0.000%
  drift) through t=12725s (14.7% of the 1-day target), dt grew properly
  to ~4217s (toward dtmax=1e4, vs stuck at ~1-2s before). TOT_ICE/TOTAL_MASS
  drifted only -0.10%/-0.11% (vs -41% by t=9938s in the old v3 run).
  I-A_INTERF rose to a peak (+35%) around t~400s then began slowly
  decreasing -- consistent with the AC interface profile overshooting
  then relaxing to its equilibrium width. Bounds clean (phi_ice/phi_air
  in [0,1]).
- Noted (not yet resolved): `inputs/solver.opts` has an uncommitted
  change from `-p 2 -C 1` to `-p 1 -C 0` with a stale comment still
  referencing "C^1 for p=2" -- flagged to user, not yet addressed.

---

## 2026-06-14 — Fix mass-conservation bug: source T/vapor with d(phi_i)/dt, not S_sub

- Diagnosed the v3 60-day run's 41% mass loss (TOT_ICE/TOTAL_MASS both
  -41.331% by step 206/t=9938s) as a missing-term conservation bug, not
  a humidity/domain/mobility issue: R_tem and R_vap were sourced only by
  the kinetic term S_sub, while the Allen-Cahn double-well/curvature part
  of dphi_i/dt was completely unbalanced in the T/vapor equations.
- Cross-checked against Moure & Fu (2024) SI, Sec. 1.2 "Sublimation"
  (the exact two-phase equivalence model this code implements), eqs.(6)-(7):
  T and vapor must be sourced by the FULL dphi_i/dt (ice_t), not just S_sub.
  Confirmed algebraically that with this change,
  rho_ice*dphi_i/dt + d(phi_a*rho_v)/dt = -div(flux), which integrates to
  zero with no-flux BCs -> exact mass conservation independent of mob_sub.
- `src/assembly.c`: Residual_A1 R_tem/R_vap now use `ice_t` instead of
  `S_sub`; Jacobian_A1 [tem,*] and [vap,*] blocks rewritten to match
  (J[a][1][b][2] dropped, now zero). [ice,*] block (R_ice/S_sub) unchanged.
  Updated module docstring to describe the new sourcing.
- `src/permafrost2.c`: removed a stray `10x` factor in
  `user.alph_sub = 10 * lambda_sub / tau_sub` -> `lambda_sub / tau_sub`,
  matching eq.(9) exactly (mob_sub formula itself was already correct).
- Verified `mob_sub = eps/(3*tau_sub)` is correctly derived from the
  Gibbs-Thomson parameters `beta_sub0`/`d0_sub0` per eq.(9) -- answers
  user's question about mobility definition; no change needed there.
- Ran a verification diagnostic (`-t_final 10000`, 2-grain boundary IC):
  TOT_ICE and TOTAL_MASS stayed EXACTLY unchanged (6.921e-10 / 6.361e-7,
  0.000% drift) through t~42.7s (98 steps), vs. measurable drift by a
  similar point in the old code. Bounds stayed in [0,1]. Adaptive dt
  stalled near 1-1.3s (much smaller than v3's dt~100), so a full 60-day
  run with this fix will be substantially more expensive -- flagged as a
  follow-up to investigate (NR iteration counts / dtmax tuning).

---

## 2026-06-14 — Trim two-grain boundary domain, relaunch 60-day run

- User asked to shrink the 190x190 domain (excess empty space): reduce
  Lx/Nx while keeping ~40 elements of vapor gap between the two grains
  (to ensure coarsening is via vapor-diffusion Ostwald ripening, not
  direct contact), and reduce Ly/Ny proportionally with ~20-element
  margins above/below the larger grain.
- Resized `inputs/geometry/2D_two_ice_grains_boundary.opts`: Lx=4.123140e-5
  (Nx=125), Ly=5.079708e-5 (Ny=154), dx unchanged (3.2985e-7). Mesh
  36100 -> 19250 nodes (~47% smaller).
- Verified with a 28-step run: TOT_ICE stays exactly 6.921e-10 (-0.000%),
  bounds in [0,1] -- matches the un-eroded K2P-fix baseline.
- Relaunched the 60-day practice run (v3) on this geometry via
  `./scripts/Studio/run_permafrost.sh 2D_two_ice_grains_boundary
  60day_T-20_h0.95 practice60day_v3`.

---

**Session ended:** 2026-06-14 15:28:48

## 2026-06-14 — Fix ice double-well coefficient (K2P) causing grain erosion

- User reported the v2 (190x190, eps/R~0.025) run's interface started at
  ~10.5 elements wide and collapsed around step 60, eroding both grains.
- Root cause: `K2P = 3*mob_sub*(Sigma_i+Sigma_a+2*Lambda)/Sigma_i` in
  `(K2P/eps)*f1(phi_i)` (src/assembly.c). Compared against
  dry_snow_metamorphism's 3-phase ice equation with met=0 (our true
  2-phase case): ETA=Etaa*Etai, fmet=0, fice=Etai*f1(ice), so the whole
  expression collapses to `mob*3/eps*f1(ice)` -- Sigma_i/Sigma_a/Lambda
  cancel completely with no third phase. With current defaults
  (Sigma_i=0.109, Sigma_a=0.132, Lambda=1.0) the old K2P was ~20.6x too
  large, making the PDE's equilibrium interface ~4.5x narrower than the
  eps-based IC -> collapse onto an unresolvable width.
- Fixed: `K2P = 3*mob_sub` in both Residual_A1 and Jacobian_A1 (commit
  `728ca51`). Removed now-unused Sigma_i/Sigma_a/Lambda locals.
- Verified with a 62-step run: TOT_ICE and I-A INTERF stay essentially
  constant (no decay) through t=1.69s, phi_ice/phi_air in [0,1].
- Killed the in-flight v2 60-day run (invalid under the old K2P) and
  cleaned up /tmp/pf_k2p_check*.

---

**Session ended:** 2026-06-14 15:25:51


---

**Session ended:** 2026-06-14 15:19:20


---

**Session ended:** 2026-06-14 15:01:50


---

**Session ended:** 2026-06-14 14:40:53


---

**Session ended:** 2026-06-14 13:23:42

## 2026-06-14 — Resize two-grain boundary IC interface/grain ratio

- User flagged the diffuse interface as "way too diffuse" in the
  practice run. Diagnosis: with eps=4.6648e-07 and RCice0/1=2e-6/4e-6,
  eps/R ~ 0.23 — the interface width (~1.32e-6 m) was comparable to the
  grain radii, so the grains themselves looked like diffuse blobs.
- Resized `inputs/geometry/2D_two_ice_grains_boundary.opts` to mirror
  `2D_touching_grains.opts`: same eps/dx, RCice1=1.875e-5 (that file's
  RCice), RCice0=half of RCice1, domain Lx=Ly=6.25e-5, Nx=Ny=190 (via
  comp_eps.py). Now eps/R ~ 0.025 — a sharp interface relative to grain
  size, matching the established touching-grains setup.
- Verified with a short run: phi_ice/phi_air stay in [0,1],
  TOT_ICE/TOT_AIR/TOTAL_MASS conserved to ~0.000%.

---

**Session ended:** 2026-06-14 01:02:38


---

**Session ended:** 2026-06-14 00:47:16

## 2026-06-14 — Fix disk-fill crash, resolve mesh/eps, run 60-day practice sim

- Fixed `src/snes_convergence.c`: an extra `%*.*e` placeholder in the
  per-iteration print (11 placeholders for 10 supplied values) caused
  undefined-behavior padding-space output that filled `/tmp` to 322GB
  and crashed a long run. Removed the extra placeholder (commit `824ac39`).
- Diagnosed a phase-field bounds violation on the original 64x32/1mm
  mesh: `eps=9.0e-7` was far smaller than `dx≈1.56e-5`, leaving the
  diffuse interface unresolved. Used `preprocess/comp_eps.py` to size a
  new domain (Lx=2.6389e-5, Ly=1.5833e-5, Nx=81, Ny=49, eps=4.6648e-07,
  RCice0=2e-6, RCice1=4e-6) with ~4 elements across the interface.
- Added `inputs/geometry/2D_two_ice_grains_boundary.opts` and
  `inputs/experiment/60day_T-20_h0.95.opts` for a 60-day, -20C/95% RH
  practice run via `run_permafrost.sh`.
- Updated `inputs/solver.opts` to `-dof 3` and removed the orphaned
  `-difvap_pen`/`-k_pen` vapor-penalty block left over from the dropped
  sediment phase.
- Ran the full 60-day simulation via `./scripts/Studio/run_permafrost.sh
  2D_two_ice_grains_boundary 60day_T-20_h0.95 practice60day`. With
  solver.opts' NRmin=5/NRmax=15, dt grew steadily from 1e-8 toward
  dtmax=1e4 with no oscillation (previously observed oscillation was a
  CLI-default NRmin/NRmax issue, resolved by using the proper opts file).
  Both ice grains fully sublimated by ~4.4 days (consistent with the
  Kelvin/Gibbs-Thomson effect on micron-scale grains under 95% RH).
- Ported `postprocess/plot_mass.py` from the old 4-DOF (ice/T/vapor/sed)
  format to 3-DOF (ice/T/vapor): removed sediment mass/series/plots and
  the `-t_sed_freeze` annotation machinery; `sol_*.dat` now reshapes to
  `(-1, 3)` and `phi_a = 1 - phi_i`. Verified `mass.png` and per-phase
  plots generate correctly for the completed run.

---

**Session ended:** 2026-06-14 00:34:07


---

**Session ended:** 2026-06-14 00:16:47


---

**Session ended:** 2026-06-14 00:15:25


---

**Session ended:** 2026-06-14 00:15:17


---

**Session ended:** 2026-06-13 17:06:39

## 2026-06-13 — Complete 3-DOF hard fork (drop sediment phase)

- Finished the 9-commit hard fork from the 4-DOF (ice/temperature/vapor/sediment)
  system to a 3-DOF (ice phi_i, temperature T, vapor density rho_v) system, with
  phi_a = 1 - phi_i computed algebraically.
- Commit 8 (`permafrost2.c`): dof=3, IGASetFieldName trimmed to
  phaseice/temperature/vap_density, replaced the gamma_im/gamma_iv/gamma_mv
  Kim-Steinbach derivation with direct -Sigma_i/-Sigma_a CLI params feeding
  F_dub(phi_i) = C*phi_i^2(1-phi_i)^2, C=(Sigma_i+Sigma_a)/2+Lambda. Removed all
  sediment CLI options/defaults (-NCsed, -RCsed*, -grain_sep, -x_slab_frac,
  -mob_sed, -thcond_sed/-cp_sed/-rho_sed, -difvap_pen/-k_pen, -t_relax) and the
  Phi_sed0 snapshot loop. -ic_type trimmed to
  two_ice_grains_boundary (new default) | ice_slab | single_ice.
- Commit 9 (`grain_initialization.c`/`.h`): removed entirely. All functions
  (InitialSedGrains, InitialSedGrainsGravity, InitializeSedimentFromInputSolution,
  InitialIceGrains, etc.) were unreferenced anywhere in the codebase, and the
  sediment-grain helpers referenced AppCtx fields already removed in commit 1
  (no longer compiled). Dropped the include from NASA_main.h.
- Build verified clean (only a pre-existing unrelated format-string warning in
  snes_convergence.c). Smoke-tested
  `-ic_type two_ice_grains_boundary -dim 2 -Nx 64 -Ny 32 -Lx 1e-3 -Ly 0.5e-3
  -RCice0 1e-4 -RCice1 2e-4 -t_final 4e-4`: bounds stay within [0,1] for
  phi_ice/phi_air, SNES convergence table shows exactly 3 DOF columns, and
  TOT_ICE/TOT_AIR/TOTAL_MASS are conserved to within rounding over 4 steps.
- Follow-up (out of scope, flagged in the plan): postprocess/*.py and
  scripts/*.py still assume the old 4-field (ndof=4, "sediment" field) output
  format and will need updating before plotting new results.

---

**Session ended:** 2026-06-11 15:06:11


---

**Session ended:** 2026-06-11 15:02:38


---

**Session ended:** 2026-06-11 14:43:01


---

**Session ended:** 2026-06-11 14:18:51


---

**Session ended:** 2026-06-11 14:18:41


---

**Session ended:** 2026-06-11 14:18:36


---

**Session ended:** 2026-06-11 14:18:33


---

**Session ended:** 2026-06-11 14:18:29


---

**Session ended:** 2026-06-11 14:18:20


---

**Session ended:** 2026-06-11 14:18:14


---

**Session ended:** 2026-06-11 14:18:04


---

**Session ended:** 2026-06-11 14:18:01


---

**Session ended:** 2026-06-11 14:17:53


---

**Session ended:** 2026-06-11 14:12:08


---

**Session ended:** 2026-06-11 14:09:26


---

**Session ended:** 2026-06-11 13:23:15


---

**Session ended:** 2026-06-11 13:21:05


---

**Session ended:** 2026-06-11 13:18:39


---

**Session ended:** 2026-06-11 13:07:56


---

**Session ended:** 2026-06-11 12:49:04


---

**Session ended:** 2026-06-11 12:48:41


---

**Session ended:** 2026-06-11 12:01:01


---

**Session ended:** 2026-06-04 10:58:15


---

**Session ended:** 2026-06-02 12:14:53


---

**Session ended:** 2026-06-02 11:11:09


---

**Session ended:** 2026-06-01 11:50:53


---

**Session ended:** 2026-06-01 11:37:25


---

**Session ended:** 2026-06-01 08:14:29


---

**Session ended:** 2026-05-31 07:06:39


---

**Session ended:** 2026-05-31 07:03:18


---

**Session ended:** 2026-05-31 06:48:08


---

**Session ended:** 2026-05-30 17:16:06


---

**Session ended:** 2026-05-30 11:27:19


---

**Session ended:** 2026-05-30 11:06:22


---

**Session ended:** 2026-05-30 11:05:37


---

**Session ended:** 2026-05-30 10:04:43


---

**Session ended:** 2026-05-29 15:39:52


---

**Session ended:** 2026-05-29 15:08:52


---

**Session ended:** 2026-05-29 14:45:55


---

**Session ended:** 2026-05-29 14:08:36


---

**Session ended:** 2026-05-29 13:48:52


---

**Session ended:** 2026-05-29 13:40:57


---

**Session ended:** 2026-05-29 13:23:15


---

**Session ended:** 2026-05-28 14:30:06

## 2026-05-28 — Diagnose advisor_T-20_nopen freeze; vap_src reformulation + GT re-enable + atol

- Diagnosed why the 2026-05-27 advisor_T-20_nopen HPC batch froze: every 2D case did a brief Allen-Cahn relaxation burst (2D_touching_grains steps ~123-129) then went bit-for-bit static for 99.4% of the 30 simulated days. Root cause: Gibbs-Thomson disabled (d0_GT=0) → no curvature driving force → genuine steady state at saturation once AC relaxes. Loose snes_atol=1e-4 amplified it (vapor residual ~rhov/dt fell below atol at large dt). Confirmed penalties were genuinely off; the static-vapor look matched the penalty case only because the end state is identical.
- Reformulated vap_src from -2*rho_ice*air*ice_t to -rho_ice*sub_src (Moure & Fu 2024): driven by the phase-change term only, so mass-neutral Allen-Cahn curvature motion no longer launders non-physical "mass loss" into spurious vapor (was ~0.6% TOTAL_MASS drift). Reverted Jacobian [vap,*] to the sub_src-derived rows; removed ice_t read.
- Re-enabled Gibbs-Thomson: restored Hessian reads, Curvature() call, rhoI_vs_eff in sub_src, [ice,ice] + [vap,ice] GT Jacobian blocks, N2 shape funs, rhoI_vs_eff/d_rhovs_eff_dtem in [ice,*]/[tem,*]/[vap,*]. Gated on d0_GT (0=off). This is the LSW driver that prevents the freeze.
- Tightened snes_atol 1e-4 → 1e-6 so weak cold-T vapor dynamics are resolved, not skipped at large dt.
- Added experiment file 30day_T-20_h1.00_nopen_GT.opts (penalties off + d0_GT=9.6e-10) for the advisor-slides rerun.
- Verified: clean build (no warnings), smoke test 2D_touching_grains T=-20 penalties-off d0_GT=9.6e-10 t_final=50s — clean exit, no NaN/divergence, mass conserved to display precision, dt grows normally under tighter atol.
- Updated model_description.md (§3.4 V_src=-rho_ice*S_sub, §4 GT active, §12 table, §14 LaTeX) and branch log §28.

---

**Session ended:** 2026-05-28 14:15:43


---

**Session ended:** 2026-05-27 06:40:48


---

**Session ended:** 2026-05-27 06:32:03

## 2026-05-27 — Penalty-off diagnostic + docs sync + HPC full-suite submitter

- Added `inputs/experiment/30day_T-20_h1.00_nopenalty.opts` (k_pen=0, difvap_pen=1.0) to diagnose model behaviour with both vapor-equation penalties disabled. Result: cleaner phase bounds and vapor field than any configuration with penalties active.
- Wrote `spurious_ice_sed_air_branch_log.md` §27 explaining why penalty-off is now the cleanest configuration — independent fixes (mass-conserving Stefan source §26, diffusivity-penalty direction §25, latent heat paired with `S_sub`, constant `rhov_eq` anchor, NRmin=5 dt unlock) collectively removed every job the penalties were originally compensating for.
- Brought `docs/model_description.md` up to date with the current code: §3.3 latent heat pairing with `S_sub` (not `∂φ_a/∂t`), §3.4 vapor equation with `V_src = -2·ρ_ice·φ_a·∂φ_i/∂t` and constant anchor `rhov_eq = h₀·ρ_vs(T)`, §4 GT-disabled clarification, §9 solver param refresh (snes_atol=1e-4, NRmin=5, ksp_rtol=1e-8), §12 parameter table update.
- Added `docs/model_description.md` §14 — self-contained LaTeX source block for the model equations (free energy, ice/sed AC, temperature, vapor, S_sub, D_eff, G_pen) plus §14a pointwise mass-balance derivation, paste-ready for manuscript prep.
- Added "Status note (2026-05-27)" banners to `ice_mass_loss_analysis.md`, `spurious_air_bug_analysis.md`, and `material_parameters.md` so readers of these historical analyses are routed to the current branch-log §27 and model_description.md.
- Added `scripts/HPC/submit_full_suite.sh` — thin wrapper around `submit_batch.sh` that mirrors `scripts/Studio/run_batch_tests.sh::DEFAULT_TESTS` byte-for-byte and supports `--skip-1d`, `--skip-hires`, `--dry-run`, and pass-through sbatch flags. Lets the same curated 21-test suite run on the cluster with one command.

---


---

**Session ended:** 2026-05-26 20:23:20


---

**Session ended:** 2026-05-26 13:54:23


---

**Session ended:** 2026-05-26 13:45:05


---

**Session ended:** 2026-05-26 13:32:50


---

**Session ended:** 2026-05-26 13:29:57


---

**Session ended:** 2026-05-26 13:07:19


---

**Session ended:** 2026-05-26 07:35:53


---

**Session ended:** 2026-05-26 07:14:23


---

**Session ended:** 2026-05-26 07:07:15


---

**Session ended:** 2026-05-25 21:57:56


---

**Session ended:** 2026-05-25 21:54:42


---

**Session ended:** 2026-05-25 18:21:14


---

**Session ended:** 2026-05-25 17:13:16


---

**Session ended:** 2026-05-25 16:38:10


---

**Session ended:** 2026-05-25 16:16:58


---

**Session ended:** 2026-05-25 15:57:27


---

**Session ended:** 2026-05-25 15:50:11


---

**Session ended:** 2026-05-25 15:07:16


---

**Session ended:** 2026-05-25 15:02:12


---

**Session ended:** 2026-05-25 14:59:06


---

**Session ended:** 2026-05-25 14:50:20


---

**Session ended:** 2026-05-25 14:41:52


---

**Session ended:** 2026-05-25 13:56:51


---

**Session ended:** 2026-05-25 13:53:23


---

**Session ended:** 2026-05-25 13:38:55


---

**Session ended:** 2026-05-25 13:34:44


---

**Session ended:** 2026-05-25 13:28:50


---

**Session ended:** 2026-05-25 13:12:07


---

**Session ended:** 2026-05-25 13:00:27


---

**Session ended:** 2026-05-25 12:46:20


---

**Session ended:** 2026-05-25 12:30:41


---

**Session ended:** 2026-05-25 12:15:17


---

**Session ended:** 2026-05-25 12:12:43


---

**Session ended:** 2026-05-25 11:56:17


---

**Session ended:** 2026-05-25 11:47:26


---

**Session ended:** 2026-05-25 11:40:16


---

**Session ended:** 2026-05-25 11:36:04


---

**Session ended:** 2026-05-25 11:20:51


---

**Session ended:** 2026-05-25 11:17:58


---

**Session ended:** 2026-05-25 11:15:42


---

**Session ended:** 2026-05-25 10:50:09


---

**Session ended:** 2026-05-25 10:18:37


---

**Session ended:** 2026-05-25 09:55:51


---

**Session ended:** 2026-05-25 08:18:57


---

**Session ended:** 2026-05-22 14:03:37


---

**Session ended:** 2026-05-22 13:59:49


---

**Session ended:** 2026-05-20 23:07:19


---

**Session ended:** 2026-05-20 21:23:52


---

**Session ended:** 2026-05-20 15:02:43


---

**Session ended:** 2026-05-20 13:01:56


---

**Session ended:** 2026-05-20 12:36:10


---

**Session ended:** 2026-05-20 11:58:41


---

**Session ended:** 2026-05-20 10:30:25


---

**Session ended:** 2026-05-20 09:44:45


---

**Session ended:** 2026-05-20 09:35:50


---

**Session ended:** 2026-05-20 08:37:02


---

**Session ended:** 2026-05-20 08:28:35


---

**Session ended:** 2026-05-20 07:35:04


---

**Session ended:** 2026-05-20 07:10:00


---

**Session ended:** 2026-05-19 21:50:42


---

**Session ended:** 2026-05-19 16:07:04


---

**Session ended:** 2026-05-19 14:28:45


---

**Session ended:** 2026-05-19 12:36:45


---

**Session ended:** 2026-05-19 12:18:17


---

**Session ended:** 2026-05-19 09:41:27


---

**Session ended:** 2026-05-19 08:59:01


---

**Session ended:** 2026-05-19 08:49:38


---

**Session ended:** 2026-05-19 08:19:24


---

**Session ended:** 2026-05-18 20:15:07


---

**Session ended:** 2026-05-18 19:04:48


---

**Session ended:** 2026-05-18 18:39:15


---

**Session ended:** 2026-05-18 18:21:00


---

**Session ended:** 2026-05-18 18:00:46


---

**Session ended:** 2026-05-18 17:34:04


---

**Session ended:** 2026-05-18 17:09:57


---

**Session ended:** 2026-05-18 16:54:13


---

**Session ended:** 2026-05-18 14:34:16


---

**Session ended:** 2026-05-18 14:33:37


---

**Session ended:** 2026-05-18 14:11:32


---

**Session ended:** 2026-05-18 11:53:02


---

**Session ended:** 2026-05-18 11:32:37


---

**Session ended:** 2026-05-18 11:16:45


---

**Session ended:** 2026-05-18 11:11:38


---

**Session ended:** 2026-05-18 10:15:16


---

**Session ended:** 2026-05-18 09:58:17


---

**Session ended:** 2026-05-18 09:45:37


---

**Session ended:** 2026-05-17 17:24:01


---

**Session ended:** 2026-05-17 17:17:52


---

**Session ended:** 2026-05-17 15:50:44


---

**Session ended:** 2026-05-17 10:51:13


---

**Session ended:** 2026-05-16 22:49:54


---

**Session ended:** 2026-05-16 22:40:53


---

**Session ended:** 2026-05-16 19:36:20


---

**Session ended:** 2026-05-16 19:12:18


---

**Session ended:** 2026-05-16 18:12:54


---

**Session ended:** 2026-05-16 18:06:41


---

**Session ended:** 2026-05-16 15:18:35


---

**Session ended:** 2026-05-16 14:57:44


---

**Session ended:** 2026-05-16 14:49:41


---

**Session ended:** 2026-05-16 14:30:12


---

**Session ended:** 2026-05-16 14:18:55


---

**Session ended:** 2026-05-16 14:07:08


---

**Session ended:** 2026-05-16 13:49:07


---

**Session ended:** 2026-05-16 13:34:36


---

**Session ended:** 2026-05-16 13:15:49


---

**Session ended:** 2026-05-16 11:42:42


---

**Session ended:** 2026-05-16 10:48:56


---

**Session ended:** 2026-05-16 10:10:47


---

**Session ended:** 2026-05-16 09:40:49


---

**Session ended:** 2026-05-16 09:14:45


---

**Session ended:** 2026-05-16 07:47:04


---

**Session ended:** 2026-05-16 07:41:22


---

**Session ended:** 2026-05-16 07:34:36


---

**Session ended:** 2026-05-15 20:25:41


---

**Session ended:** 2026-05-15 19:53:35


---

**Session ended:** 2026-05-15 19:42:56


---

**Session ended:** 2026-05-15 19:08:44


---

**Session ended:** 2026-05-15 19:00:51


---

**Session ended:** 2026-05-15 18:57:38


---

**Session ended:** 2026-05-15 18:12:00


---

**Session ended:** 2026-05-15 17:54:02


---

**Session ended:** 2026-05-15 17:47:16


---

**Session ended:** 2026-05-15 17:46:10


---

**Session ended:** 2026-05-15 17:26:11


---

**Session ended:** 2026-05-15 13:39:38


---

**Session ended:** 2026-05-15 13:25:47


## 2026-05-20 — Documented Gibbs-Thomson + vapor-transport fixes in branch log + model_description

- Added §24 (Gibbs-Thomson curvature dependence with analytical Jacobian), §25 (vapor diffusivity penalty direction fix — was on phi_ice, should have been on phi_ice+phi_sed), and §26 (vap_src reformulation: pair with sub_src + restore localized k_pen) to `docs/spurious_ice_sed_air_branch_log.md`.
- Updated §23c and the §23 arc summary to remove the now-stale claim that "k_pen=0 settled as the operating configuration" — that was a diagnostic config; the penalty was later re-engaged at moderate strength in §26.
- Extended the trajectory summary at the bottom of the branch log with a 5th stage covering the GT enablement and the cascade of vapor-transport fixes it forced.
- Updated `docs/model_description.md`: rewrote §3.4 (vapor density equation) and §4 (sublimation kinetics) to reflect the current model — GT-corrected rhoI_vs_eff in sub_src, vap_src paired with sub_src, diffusivity penalty switching on g(ice+sed), shifted PenaltyWeight band [0.90, 1.00]. Updated parameter summary table values (k_pen=1e3, α_pen=1e-8, ξ_v=1, ξ_T=1, t_sed_freeze=10, Λ=1e4) and added d0_GT row. Added forward-pointer at the top of §13 noting that the post-bug-cascade model fixes are documented in the branch log.

---

## 2026-05-19 — Documented penalty + Ostwald-ripening arc in branch log

- Added §23 to `docs/spurious_ice_sed_air_branch_log.md` covering the full arc of changes since the previous log: identifying that `k_pen` suppressed Gibbs-Thomson, the `PenaltyWeight()` localization experiment, discovering the localized penalty was still acting as a directional mass sink (manufacturing a counterfeit GT signal), settling on `k_pen = 0`, observing both Ostwald ripening and necking at -5°C on the narrow-gap geometry, the `grain_sep` 5→20 µm widening for clean LSW isolation, and the HPC regression sweep.
- Extended the trajectory summary at the bottom of the branch log with a 4th stage ("Vapor-penalty refinement and Ostwald-ripening diagnostics") so the document can carry the narrative arc for slide preparation.
- No code or opts changes today.

---

## 2026-05-18 — Localized vapor penalty + warmer-T Ostwald ripening experiment

- Diagnosed why disabling `k_pen` left a numerical artifact: with no penalty anywhere, tiny ice motions near solid boundaries deposit cumulative vapor mass inside sed (vap_src = -ρ_ice·ice_t), and `difvap_pen=1e-8` is too small to redistribute it. This produced the negative ρ_v ~ -1e-1 the user observed in ParaView inside sediment, and explains most of the -60% TOT_RHOV drift (integral pulled down by negative values in solid, not real bulk depletion).
- Implemented "option B" for the penalty: added `PenaltyWeight()` in `material_properties.{h,c}` — a SmoothHeaviside of `(ice+sed - 0.85) / 0.15` clamped to [0,1]. Zero at the diffuse ice-air interface (ice+sed ≈ 0.5) so Gibbs-Thomson can emerge; unity in solid interior (ice+sed → 1) so rhov stays pinned to rhov_eq.
- Replaced `g_phiiphis` / `dg_phiiphis` with `g_pen` / `dg_pen` in the penalty terms of `Residual_A1` and `Jacobian_A1` only (the difvap_pen weighting still uses the original `g_phia`).
- Restored `-k_pen 1.0e3` in `solver.opts` (was 0); strength is moderate but penalty is localized so it does not suppress GT.
- Created `inputs/experiment/30day_T-5_h1.00.opts` for the warmer-T diagnostic — at -5°C, rho_v_sat is ~12× larger than at -20°C, giving the system a much bigger vapor budget for grain-to-grain transfer in the same 30-day window.

---

## 2026-05-15 — Ice mass loss analysis, material parameters, and monitoring improvements

- Diagnosed ice shrinkage as Allen-Cahn curvature coarsening (not sublimation): v_AC/v_phys ≈ 2.5×10⁸; proved via mass balance (ΔTOT_RHOV ≈ ρᵥₛ × ΔTOT_AIR, not ρᵢ × |ΔTOT_ICE|).
- Created `docs/ice_mass_loss_analysis.md`: quantitative root-cause analysis, parameter recommendations, negative φᵢ at sed-air interface explanation.
- Created `docs/material_parameters.md`: surface energies (γ_iv, γ_is, γ_sv), thermal properties, and grain geometry for water ice and lunar regolith with 14 literature citations.
- Extended `AppCtx` in `include/NASA_types.h` with 4 initial-integral fields (tot_ice_0, tot_air_0, tot_sed_0, tot_rhov_0).
- Updated `src/monitoring.c` to store initial integrals at step 0 and print percentage-change row after every domain-integral output.
- Added dashed semi-opaque reference lines at initial mass values to `postprocess/plot_mass.py`.
- Reduced `d0_sub0` (1e-9 → 1e-11) and `beta_sub0` in `src/permafrost2.c` to suppress Allen-Cahn coarsening.
- Updated `inputs/universal.opts`: tightened solver tolerances (snes_atol 1e-12, ksp_atol 1e-12), raised k_pen (1e5→1e7), difvap_pen (1e-4→1e-8), Lambda (1e2→1e3), t_sed_freeze (1→300).
- Updated `inputs/tests/test_2D_IceSedPair.opts`: humidity set to 1.00 for AC-coarsening diagnostic.

---

**Session ended:** 2026-05-15 12:19:43


---

**Session ended:** 2026-05-15 11:55:48


---

**Session ended:** 2026-05-15 11:24:23


---

**Session ended:** 2026-05-15 10:41:31


---

**Session ended:** 2026-05-15 10:32:44


---

**Session ended:** 2026-05-15 09:42:06


## 2026-05-15 — deep root cause analysis of spurious air bug

- Expanded docs/model_description.md §13 from 4 root causes to 8, with full quantitative analysis.
- PRIMARY BUG identified: vapor penalty term in old assembly.c was missing the xi_v factor (`k_pen * g * drhov` instead of `xi_v * k_pen * g * drhov`). With xi_v=1e-3 and k_pen=1e9, this made the penalty 1000× too large, rendering the vapor equation essentially algebraic.
- SECONDARY BUG: old IC set rhov = h0 * rho_vs everywhere, including inside solid. Correct formula is rhov = rho_vs * (h0 * phi_air + (1 - phi_air)), which gives rhov = rho_vs inside solid. With h0=0.5 and k_pen=1e9, the initial disequilibrium inside ice created a penalty residual 50× larger than the time derivative.
- Documented cascade failure: IC bug → enormous t=0 residual → FD Jacobian errors of O(k_pen * eps_machine) = O(10) → Newton diverges or converges to wrong root → phi_i < 0 → spurious air accumulates over time.
- Added quantitative table comparing all 12 changed parameters; added per-bug magnitude analysis and cascade walkthrough.

---

**Session ended:** 2026-05-15 09:24:34


---

**Session ended:** 2026-05-15 08:48:30


## 2026-05-15 — model description document + t_sed_freeze tuning

- Created docs/model_description.md: comprehensive model description covering governing equations, three-phase free energy, Allen-Cahn dynamics, vapor penalty scheme, sublimation kinetics, material properties, IGA/B-spline discretisation, generalized-α time integration, Newton-Krylov nonlinear solver, per-DOF convergence test, GMRES+ASM+ILU(2) linear solver, adaptive time stepping, and initial conditions. Includes analysis of why k_pen=1e9+FD Jacobian caused spurious air vs. the corrected k_pen=1e5+analytical Jacobian.
- User changed t_sed_freeze 300 → 1 in universal.opts to more rapidly enter 2-phase mode.

---

## 2026-05-15 — tighten solver tolerances and improve preconditioner

- Switched preconditioner from bjacobi+ILU(1) to ASM+ILU(2) with pc_asm_overlap=1 in universal.opts; additive Schwarz with node overlap improves conditioning at MPI block boundaries for the 4-DOF coupled B-spline system.
- Tightened SNES atol 1e-8→1e-10 and raised snes_max_it 10→15 to prevent false convergence and handle stiff first timesteps.
- Tightened KSP rtol 1e-5→1e-6 and atol 1e-8→1e-10 for consistency with SNES and to preserve Newton's quadratic convergence rate.
- Reduced GMRES restart 500→200; the original restart=max_it meant unrestarted GMRES storing ~500 Krylov vectors; 200 sufficient for well-preconditioned problem.
- Raised NRmax default 5→8 in universal.opts (consolidates per-file overrides from ~15 test files).
- Added direct LU override (ksp_type preonly + pc_type lu) to all 10 1D test opts files; N≤760 DOFs so LU is trivially cheap and exact.
- Removed redundant -NRmax 8 override blocks from all 21 test opts files.

---

**Session ended:** 2026-05-14 16:18:20


---

**Session ended:** 2026-05-14 16:12:24


---

**Session ended:** 2026-05-14 16:08:39


---

**Session ended:** 2026-05-14 15:55:02


## 2026-05-14 — add single_ice and ice_sed_pair initial conditions

- Added FormInitialSingleIceGrain2D/1D: single pure ice grain centred in the domain (no sediment). Square domain Lx=Ly=6.25e-5 m, standard Nx=Ny=96, hi-res Nx=Ny=192.
- Added FormInitialIceSedPair2D/1D: one ice grain at (Lx/2, Ly/4) and one sediment grain at (Lx/2, 3Ly/4), evenly spaced with equal wall clearance and grain-surface gap. Same 2:1 domain as TouchingGrainPair; standard 95×190, hi-res 190×380.
- Both ic_types wired into permafrost2.c dispatch (1D and 2D); error message updated; -ic_type option string updated.
- Created 8 new opts files: test_{1D,2D}_SingleIceGrain{,_hires}.opts and test_{1D,2D}_IceSedPair{,_hires}.opts.
- Added single_ice→SingleIceGrain and ice_sed_pair→IceSedPair to subfolder mapping in both run scripts.
- Compiles clean (zero warnings).

---

## 2026-05-14 — add plot_mass.py phase mass tracking

- Created postprocess/plot_mass.py: integrates ice, sediment, and vapor masses from sol_*.dat snapshots using first-order quadrature over the uniform mesh; plots all four curves (ice, sed, vap, total) with a vertical dashed line at t_sed_freeze marking the 3-phase → 2-phase model switch.
- Time axis mapped from SSA_evo.dat; time units auto-selected by run length; mass units dimension-aware (kg m⁻² / kg m⁻¹ / kg).
- Updated postprocess/run_postprocess.sh to call plot_mass.py automatically when igasol.dat and sol_*.dat are present.

---

**Session ended:** 2026-05-14 15:43:21


## 2026-05-14 — opts eps/mesh cleanup, EnclosedGrainPair split, auto folder naming

- Set eps=7.12e-07 in all 6 standard-resolution test opts files (mesh sizes unchanged).
- Set Nx/Ny = 2× standard in all 8 hi-res opts files (eps=3.56e-07 unchanged); updated header comments.
- Split 4 generic EnclosedGrainPair opts files into 8: 4 TouchingGrainPair (grain_sep=0) + 4 SeparatedGrainPair (grain_sep=5.0e-6); each has standard and hires companion.
- Both run scripts (Studio + HPC): output folder auto-derived from ic_type → subfolder (IceSlab, EnclosedGrainPair, etc.) and opts basename + timestamp; second arg is now an optional tag.
- HPC run script: removed all post-processing calls; results are staged for local post-processing after rsync.
- submit_permafrost.sh: title arg made optional; SLURM --job-name set from opts basename.
- Created postprocess/run_postprocess.sh: self-contained post-processing runner copied into every run folder; handles VTK, 1D profiles, scalars, and timestep diagnostic.

---

**Session ended:** 2026-05-14 15:16:18


---

**Session ended:** 2026-05-14 14:38:27


---

**Session ended:** 2026-05-14 14:12:28


---

**Session ended:** 2026-05-14 14:10:02


---

**Session ended:** 2026-05-13 14:12:58


---

**Session ended:** 2026-05-13 13:23:54


---

**Session ended:** 2026-05-13 12:56:01


## 2026-05-13 — Centralize penalty params in universal.opts; clean all test opts files

- Moved `t_sed_freeze` (300 s), `difvap_pen` (1e-4), `k_pen` (1e5), and `Lambda` (1e2) from individual test opts files into `universal.opts`, using values from the working `test_1D_IceSlab.opts`.
- Rewrote all 8 remaining test opts files (`test1–test5`, `test_1D_*`, `test_2D_IceSlab`) with consistent header blocks and section labels.
- Added `-flag_avenue 1` to all files that were missing it.
- Removed duplicate SNES/KSP/PC settings (already in `universal.opts`) from test1/test2/test4.
- `test_2D_IceSlab.opts` retains `t_sed_freeze 60` as an explicit override (geometry needs shorter freeze than default 300 s).
- `test1_IceCap.opts` retains `-dof 3` override (no sediment DOF for ice-cap geometry).

---

**Session ended:** 2026-05-13 12:19:11


---

**Session ended:** 2026-05-13 12:05:04


---

**Session ended:** 2026-05-13 11:45:35


## 2026-05-13 — Analytical Jacobian for Avenue 1 + speedups

- Replaced `IGAFormIJacobianFD` with fully analytical `Jacobian_A1` in `src/assembly.c`.
- Jacobian covers all three Avenue 1 modes: relaxation, 3-phase pre-freeze, 2-phase frozen-sediment.
- PETSc Jacobian checker confirms accuracy: `||J - Jfd||/||J|| = 1.9e-9` (relaxation), `1e-7–9e-6` (active physics).
- Hoisted 8 loop-invariant scalars out of `Residual_A1`'s per-node loop.
- Disabled Avenues 2 and 3 (bodies replaced with no-op stubs; dispatcher raises error with clear message).
- Switched `permafrost2.c` to use analytical `Jacobian` instead of `IGAFormIJacobianFD`.
- Compiles clean (zero warnings); Newton converges in 1–3 iterations per step.

---

**Session ended:** 2026-05-13 11:07:12


---

**Session ended:** 2026-05-13 10:34:57


## 2026-05-13 — Add slab_and_grains 2D initial condition for ice migration

- Added `FormInitialSlabAndGrains2D` to `src/initial_conditions.c`: solid ice slab on the right (`x_slab_frac` fraction of Lx), random non-overlapping ice + sediment grains on the left, with grains placed using the same algorithm as `FormInitialRandomPackedPermafrost2D` but restricted to the grain region.
- Used `PetscMalloc1`/`PetscFree` for grain arrays (no VLAs, per code style).
- Added `x_slab_frac` field to `AppCtx` in `include/NASA_types.h` (default 0.175).
- Registered `-x_slab_frac` option in `src/permafrost2.c`; added `"slab_and_grains"` to IC dispatch.
- Declared function in `include/initial_conditions.h`.
- Created `inputs/tests/test5_SlabAndGrains.opts` (500×250 mesh, 15 ice + 12 sed grains, eps=5e-7, 10-day run).
- Verified with `mpiexec -n 4`: IC sets up correctly, phase bounds [0,1], SNES converges in 1-2 iterations.

---

**Session ended:** 2026-05-13 08:13:34


---

**Session ended:** 2026-05-13 08:07:31


---

**Session ended:** 2026-05-12 17:02:20


---

**Session ended:** 2026-05-12 16:47:50


---

**Session ended:** 2026-05-12 16:39:04


---

**Session ended:** 2026-05-12 16:25:48


---

**Session ended:** 2026-05-12 16:15:27


---

**Session ended:** 2026-05-12 16:11:04


## 2026-05-12 — Merge feature branch, delete stale branches, modernize HPC script

- Merged `feature/refactor-modular-parallel-clean` into `main` (already done at session resume).
- Committed HPC run script modernization: project-root auto-resolution from script location, `compute_optimal_nprocs()` for auto-sizing MPI ranks from Nx/Ny/Nz, `run_1d_plotting()` for automatic 1D post-processing, `srun`/`mpiexec` dispatch.
- Deleted local branch `feature/refactor-modular-parallel-clean` (merged).
- Deleted stale remote branches: `feature/homogenization`, `feature/permafrost-io-system`, `safe_branch`, `feature/refactor-modular-parallel-clean`.
- Repository now has a single branch: `main`.

---

**Session ended:** 2026-05-12 16:01:29


---

**Session ended:** 2026-05-12 15:55:52


---

**Session ended:** 2026-05-12 15:44:56


---

**Session ended:** 2026-05-12 15:01:48


## 2026-05-12 — Fix rhov initialization to saturated vapor density in all IC functions

- Updated all IC functions in `src/initial_conditions.c` to use the phase-weighted formula `rhov = rho_vs * (hum0 * phi_air + (1 - phi_air))` so vapor density equals `rho_vs` (saturated) inside ice/sediment and `hum0 * rho_vs` in air.
- Previously all functions set `rhov = hum0 * rho_vs` uniformly, under-saturating solid regions at t=0 and potentially stiffening the first Newton step.
- `FormInitialIceSlab2D` was already correct and left unchanged.
- Functions fixed: FormInitialLayeredPermafrost2D, FormInitialFlatSedIceCap2D, FormInitialEnclosedPermafrost2D, FormInitialContactSedPermafrost2D, FormInitialRandomEnclosedPermafrost2D, FormInitialRandomPackedPermafrost2D, FormInitialCondition2D, FormInitialCondition3D, FormLayeredInitialCondition2D, FormInitialCondition1D, FormInitialEnclosed1D, InitializeFromInputSolution, FormIC_grain_ana.
- Compiled cleanly.

---

**Session ended:** 2026-05-12 14:34:26


---

**Session ended:** 2026-05-11 12:11:26


---

**Session ended:** 2026-05-11 11:51:10


---

**Session ended:** 2026-05-11 11:10:29


---

**Session ended:** 2026-05-11 10:41:17


---

**Session ended:** 2026-05-08 08:17:50


---

**Session ended:** 2026-05-07 08:17:46


---

**Session ended:** 2026-05-06 16:18:58


---

**Session ended:** 2026-05-06 14:28:51


---

**Session ended:** 2026-05-06 14:18:20


---

**Session ended:** 2026-05-06 14:08:39


---

**Session ended:** 2026-05-06 14:01:55


---

**Session ended:** 2026-05-06 13:57:32


---

**Session ended:** 2026-05-06 13:54:42


---

**Session ended:** 2026-05-06 12:19:58


---

**Session ended:** 2026-05-06 12:19:53


---

**Session ended:** 2026-05-06 12:16:47


---

**Session ended:** 2026-05-06 12:16:16


---

**Session ended:** 2026-05-06 12:16:14


---

**Session ended:** 2026-05-06 11:21:10


---

**Session ended:** 2026-05-06 11:04:31


---

**Session ended:** 2026-05-06 10:40:05


## 2026-05-06 — tune_sed_freeze.py sweep script + n_relax smoke test

- Added `scripts/tune_sed_freeze.py`: 2-D parameter sweep over `t_sed_freeze × k_pen` to find the minimum sediment-freeze delay that avoids spurious air generation at the ice-sediment boundary.
- Metrics: `dt_final`, `dt_growth_frac`, `ice_sed_chg_pct` (from `relax_monitor.dat`), `air_growth_pct` (from `outp.txt`), `max_snes_iters` (nonlinear only); outputs ASCII table, CSV, and 6-panel heatmap.
- Smoke-tested (2 runs) and ran full 5×5 sweep (25 runs) — all PASS on 1D slab test (expected: 1D geometry has no ice-sediment contact, so metrics don't differentiate; real discrimination requires 2D grain-pair case).
- Set `n_relax 12` in `test_1D_IceSlab.opts` (was `0`).
- Fixed `run_permafrost.sh`: `make clean && make all` to prevent stale object files.

---

**Session ended:** 2026-05-05 16:38:34


---

**Session ended:** 2026-05-05 16:27:51


## 2026-05-05 — n_relax relaxation phase + Permafrost_output fix

- Fixed missing `sol*.dat` output: `test_1D_IceSlab.opts` had `-Permafrost_output 0`, which prevented `OutputMonitor` from registering; restored to 1 and corrected the misleading "VTK" comment in all three opts files.
- Added `n_relax` parameter: AC-only relaxation steps before full physics starts. Phase fields evolve under 3-phase Allen-Cahn (no sublimation, no T/vapor coupling) for first N steps, then full physics resumes.
- Replaced defunct `flag_tIC == 1` / `nsteps_IC` mechanism with dedicated `flag_relax` (PetscBool) + `n_relax` (PetscInt) fields in AppCtx.
- `flag_tIC` is now geometry-only (0=centered slab, 2=flat interface).
- Sediment freeze transition (`t_sed_freeze`) is deferred until after relaxation completes, even if `t_sed_freeze = 0`.
- Updated `test_1D_IceSlab.opts` with `-flag_tIC 0` and `-n_relax 0` (disabled by default).

---

**Session ended:** 2026-05-05 16:03:17


---

**Session ended:** 2026-05-05 15:57:25


---

**Session ended:** 2026-05-05 15:54:56


---

**Session ended:** 2026-05-05 15:48:48


---

**Session ended:** 2026-05-05 15:45:38


---

**Session ended:** 2026-05-05 15:37:18


---

**Session ended:** 2026-05-05 14:53:44


---

**Session ended:** 2026-05-05 14:52:19


## 2026-05-05 — Avenue 1 restructure, D_pen fix, opts cleanup

- Fixed `difvap_pen`: changed from absolute diffusivity to dimensionless factor so `D_pen = difvap_pen * difvap`; default changed from `3e-5` to `1e-5`; updated in both `Residual_A1` and `Residual_A2`.
- Restructured `Residual_A1` into explicit THREE-PHASE and TWO-PHASE blocks; removed the nested `flag_2ph_ice` patch in favor of a clean top-level branch; thermal and vapor residuals are now clearly shared by both formulations.
- Simplified sediment freeze logic: `t_sed_freeze <= 0` starts immediately in 2-phase; `flag_sed_mode` removed from opts parsing.
- Removed `k_pen` sentinel (now has explicit default `1e7`); kept `k_sed_pen` sentinel for Avenue 2.
- Cleaned up `test_1D_IceSlab.opts`: removed `flag_sed_mode` and `k_sed_pen`; updated `difvap_pen` to `1e-5` (factor); reorganised sections.

---

**Session ended:** 2026-05-04 11:51:14


---

**Session ended:** 2026-05-04 11:25:20


---

**Session ended:** 2026-05-04 10:52:22


---

**Session ended:** 2026-05-04 10:41:15


---

**Session ended:** 2026-05-04 10:38:15


---

## 2026-05-04 — Physics fixes and postprocessing rewrite

- `assembly.c`: Simplified A1 vapor residual time term from `N0[a] * (air_eff * rhov_t + air_t * rhov)` to `N0[a] * rhov_t` — removes the air-phase weighting on the vapor mass time derivative to reduce spurious coupling with phase evolution.
- `initial_conditions.c`: Fixed 2D ice-slab rhov IC to blend saturation density inside ice/sed regions (`rho_vs`) with undersaturated air (`hum0 * rho_vs`), improving physical consistency at initialization.
- `material_properties.c` + `.h`: Added `SmoothHeavisidePoly()` utility (extracted from inline Mobility code). Changed `Mobility()` to use linear interpolation (hi=ice, hs=sed, ha=air) instead of cubic Hermite.
- `test_2D_IceSlab.opts`: Added `-difvap_pen 1e-5` and `-k_pen 1e7` as defaults so the sweep script's opts file is self-contained.
- `postprocess/plot1D_profiles.py`: Major rewrite — twin y-axes for thermal plots (T + ρ_v on one panel), first/last comparison always produced (no flag needed), removed GIF/animation support, removed cmocean colormap variables (just checks availability), cleaner helper structure.

---

**Session ended:** 2026-04-30 15:03:44


---

## 2026-04-30 — Rewrite tune_difvap_pen.py as independent 2D penalty sweep

- Replaced 1D `difvap_pen`-only sweep (with derived `k_pen = difvap_pen/eps²`) with an independent 2D grid sweep over `difvap_pen × k_pen` using `itertools.product`.
- Only `-difvap_pen` and `-k_pen` are passed as CLI overrides to the binary; `t_final`, `t_sed_freeze`, and all physics flags come from the opts file (`test_2D_IceSlab.opts`) unchanged.
- Removed two-run A/B structure (flat-saturated + sublimation); replaced with a single run per grid point.
- Ported metrics from `tune_sed_params.py`: `delta_tot_air` (primary spurious-air metric), `ice_sed_drop_pct`, `max_snes_iters`, phase bounds check.
- Added `_get_t_sed_freeze(opts_file)` helper to parse freeze time from opts so the script can locate the freeze event in monitor output.
- New outputs: `sweep_penalty2d.csv` and `sweep_penalty2d.png` (4-panel 2D heatmap on log-scale axes).
- Default sweep: 5×5 = 25 runs; timeout 600 s/run.

---

**Session ended:** 2026-04-30 13:49:28


---

## 2026-04-29 — Add inline comments to all .opts files

- Annotated all 17 `.opts` files in `inputs/tests/` with inline `#` comments on every parameter line.
- Comments explain the purpose/units of each flag: domain params (`-dim`, `-dof`, `-Nx`, `-Lx`), element continuity (`-p`, `-C`), time stepping, model flags (`-flag_sed_mode`, `-flag_avenue`, `-flag_2ph_ice`, `-flag_tIC`, `-flag_BC_Tfix`), SNES/KSP/PC solver settings, and grain geometry params (`-NCsed`, `-RCsed`, etc.).
- No parameter values changed — comments only.

---

**Session ended:** 2026-04-29 15:47:34


## 2026-04-29 — Three-avenue modular residual refactor + unused flag cleanup

- Rewrote `assembly.c` with three distinct, self-contained residual functions (`Residual_A1`, `Residual_A2`, `Residual_A3`) plus a `Residual` dispatcher.
  - **Avenue 1**: Allen-Cahn + penalty vapor; after `t_sed_freeze`, sediment RHS = 0 (frozen by inertia); optional 2-phase ice switch via `flag_2ph_ice`.
  - **Avenue 2** (default): Allen-Cahn + penalty vapor; after `t_sed_freeze`, sediment Allan-Cahn + restoring penalty `k_sed*(phi_s - phi_s0)`; optional 2-phase ice via `flag_2ph_ice`.
  - **Avenue 3**: Cahn-Hilliard biharmonic form for ice and sediment (uses `IGAPointFormHess` + `N2` shape function Hessians; requires p≥2, C≥1); standard vapor diffusion, no penalty parameters.
- Added `ChemPot_dsed` static helper — computes `∂fi/∂sed`, `∂fs/∂sed`, `∂fa/∂sed` analytically for the CH gradient term in A3.
- Added `flag_avenue` (1/2/3) and `flag_2ph_ice` (0/1) to `AppCtx` (NASA_types.h) and registered as CLI options in `permafrost2.c`.
- Removed unused flags/fields: `flag_xiT`, `flag_it0`, `mob_sed_arr` from `AppCtx` and all source files; removed `flag_sedgrav` CLI option (was a no-op local variable).
- Simplified `flag_sed_mode`: removed mode 2 (was the old Avenue 1 experiment); only -1/0/1 remain.
- Updated `monitoring.c`: removed `mob_sed_arr` assignment, removed `flag_it0` reset, restored freeze trigger to mode-1-only check.
- Updated `assembly.h` to declare all three avenue functions.
- Updated `test_2D_IceSlab.opts` to use `flag_avenue 1`, `flag_sed_mode 1`, `flag_2ph_ice 1` (equivalent to prior `flag_sed_mode 2` experiment).
- Build: clean, zero errors.

---

**Session ended:** 2026-04-29 10:54:06


---

**Session ended:** 2026-04-29 10:07:10


---

**Session ended:** 2026-04-28 17:14:42


---

**Session ended:** 2026-04-28 17:04:45


---

**Session ended:** 2026-04-28 16:40:27


---

**Session ended:** 2026-04-28 16:12:21


---

**Session ended:** 2026-04-28 16:04:47


---

**Session ended:** 2026-04-28 15:50:35


---

**Session ended:** 2026-04-28 15:43:47


---

**Session ended:** 2026-04-28 15:41:14


## 2026-04-28 — Add phase-field out-of-bounds detection for difvap sweep

- `monitoring.c`: added a quadrature-point loop each monitor step that computes global min/max of `phi_ice`, `phi_sed`, and `phi_air` via `MPI_Allreduce`, then prints a `BOUNDS:` line to stdout/`outp.txt`.
- `tune_difvap_pen.py`: added `_parse_bounds_violation()` that scans the captured output for any reported phase value outside `[-0.25, 1.25]`; wired `pass_bounds` into the per-run `PASS` logic and added a `bounds` column to the summary table.
- Rebuilt binary cleanly (no new errors).

---

**Session ended:** 2026-04-28 15:26:17


---

**Session ended:** 2026-04-28 15:15:56


---

**Session ended:** 2026-04-28 15:07:14


---

**Session ended:** 2026-04-28 14:55:47


---

**Session ended:** 2026-04-28 14:47:27


---

**Session ended:** 2026-04-28 13:18:23


---

**Session ended:** 2026-04-28 11:51:27


---

**Session ended:** 2026-04-23 18:23:33


---

**Session ended:** 2026-04-23 18:11:34


---

**Session ended:** 2026-04-23 17:33:56


---

**Session ended:** 2026-04-23 17:20:21


---

## 2026-04-23 — Propagate stability-triggered sediment freeze to all opts files

- Added `-sed_freeze_tol 1.0e-3` to all 17 `.opts` files in `inputs/tests/` to make the default explicit.
- Raised `-nsteps_sed` from 10/45 → 200 on long production runs with ice-sed contact (1D/2D IceSlab, EnclosedGrainPair, IceCap, CapillaryBridge, ContactSedPair, T05, T11, T12); quick 10-20 step regression tests and no-sediment flat-interface tests keep `nsteps_sed 10`.
- Changed `test3_relax.opts` to `flag_sed_mode -1` (never freeze), consistent with its purpose of studying unrestricted 3-phase relaxation.

---

**Session ended:** 2026-04-23 17:11:54


---

**Session ended:** 2026-04-23 16:47:08


---

## 2026-04-23 — Fix phase-variable mismatch: always-3-phase ice equation with gated sediment penalty

- Root cause identified: the old code changed **both** the ice equation and added a sediment penalty when `flag_sed_frozen` fired. This created an equilibrium discontinuity — the 3-phase ice equilibrium is incompatible with the 2-phase ice equation, causing spurious air and negative ice at the ice-sediment interface.
- Fix in `assembly.c`: widened `else if (flag_sed_frozen == 0)` to a plain `else` so the 3-phase ice equation is always used. The sediment penalty `R_sed += k_sed*(sed-sed0)*N0[a]` is now gated on `flag_sed_frozen` inside the 3-phase branch — it activates when sediment needs pinning but never changes the ice residual form.
- Deprecated 2-phase ice branch preserved as a commented block with explanation of why it must not be restored.
- Fixed stale comment in `permafrost2.c`: default `k_sed_pen` is `1e-7/eps²`, not `1e-3/eps²`.
- Updated `test_2D_IceSlab.opts`: changed `flag_sed_mode 0` (always-frozen, no relaxation) to `flag_sed_mode 1` (10-step free relaxation before penalty activates).
- This change also addresses vanishing sediment and grain sintering in 3-phase mode: once the penalty activates, each grain is pinned near its initial shape (`sed0`), preventing Gibbs-Thomson shrinkage and grain-grain merging.

---

**Session ended:** 2026-04-23 14:58:52


## 2026-04-23 — flag_sed_mode CLI option and 2-phase sediment drift fix

- Added `-flag_sed_mode` to `AppCtx` and registered as a PETSc CLI option: `-1` = always 3-phase, `0` = always 2-phase, `1` (default) = switch after `nsteps_sed` steps.
- Changed default `nsteps_sed` from 0 to 10 so mode 1 runs 10 relaxation steps before freezing sediment.
- Fixed sediment drift bug in 2-phase mode (`assembly.c`): the 2-phase `R_sed` branch was missing the `k_sed*(sed-sed0)*N0[a]` restoring penalty, allowing numerical drift through SNES tolerance.
- Updated `monitoring.c` freeze logic to only trigger when `flag_sed_mode == 1` (previously always triggered when `nsteps_sed > 0`, which prevented mode -1 from working).
- Added `-nsteps_sed 10` and `-flag_sed_mode 1` to all 17 `.opts` files in `inputs/tests/`.

---

**Session ended:** 2026-04-23 12:38:28


---

**Session ended:** 2026-04-23 12:37:15


---

**Session ended:** 2026-04-23 12:22:25


---

**Session ended:** 2026-04-23 12:17:22


---

**Session ended:** 2026-04-23 12:12:31


---

**Session ended:** 2026-04-23 11:30:16


---

**Session ended:** 2026-04-23 11:25:51


---

**Session ended:** 2026-04-22 21:50:40


## 2026-04-22 — Penalty parameter tuning framework

- Exposed `difvap_pen`, `k_pen`, and `k_sed_pen` as PETSc CLI options so parameter sweeps no longer require recompilation; defaults unchanged.
- Added `scripts/tune_vapor_penalty.py`: sweeps `difvap_pen` and `k_pen` independently; measures rhov drift, sublimation rate, SNES iters, in-ice rhov error; outputs CSV + PNG.
- Added `scripts/tune_sed_penalty.py`: sweeps `k_sed_pen` prefactor over 7 decades; measures sediment drift, spatial shape error, SNES iters; outputs CSV + PNG.

---

**Session ended:** 2026-04-22 18:04:13


---

**Session ended:** 2026-04-22 16:17:23

## 2026-04-22 — 2D ice slab IC and improved 1D plotting

- Added `FormInitialIceSlab2D` in `src/initial_conditions.c`: 2D equivalent of the 1D centered-slab IC (flag_tIC=0), using tanh profiles in x uniform across Ly.
- Registered `ic_type = "ice_slab"` in `src/permafrost2.c` and declared in `include/initial_conditions.h`.
- Created `inputs/tests/test_2D_IceSlab.opts`: thin 2D domain (114×23 nodes, Lx=1e-4 m, Ly=2e-5 m) matching the 1D slab parameters.
- Updated `postprocess/plot1D_profiles.py`: added dashed black sum-of-phases line (φ_i+φ_s+φ_a) to per-step phase plots as a partition-of-unity diagnostic.
- Added `plot_thermal_steps()` in plotting script: saves per-step T(x) and ρ_v(x) figures on individual subplots to `thermal_steps/` directory, running alongside phase images by default.

---

**Session ended:** 2026-04-22 10:54:16


## 2026-04-22 — Interface relaxation analysis for test3_EnclosedGrainPair

- Created `inputs/tests/test3_relax.opts`: copy of test3 with `nsteps_sed=0` and `t_final=100000` to run full 3-phase relaxation without sediment freezing.
- Ran simulation (125 steps, t=0–100,409 s) to generate `relax_monitor.dat` with uninterrupted sediment evolution.
- Analysis finding: `AssemblePhase_reinit.dat` IC is already fully relaxed — per-step changes in `ice_air_interf` and `ice_sed_interf` are <1e-7 from step 0, well below the 1e-3 threshold. Sediment freezing at `nsteps_sed=60` was overly conservative; sediment degrades numerically past ~step 100 (ice-sed interface drops to 6% of initial by step 124 without freezing).
- Updated `nsteps_sed` in `test3_EnclosedGrainPair.opts`: 60 → 10 (safe conservative buffer, no performance penalty).
- Enhanced `scripts/plot_relax.py` with `--slideshow`, `--outdir`, `--fmt` flags; generates 5 separate 16:9 high-DPI PNG files for slideshow use (combined, normalized, rate, early-zoom, physical-time axis); switched to `constrained_layout` to fix `twiny` compatibility.
- Saved 5 slideshow plots to `SimulationResults/.../test3_relax_.../relax_plots/`.

---

**Session ended:** 2026-04-22 10:10:26


---

**Session ended:** 2026-04-21 13:52:54


---

**Session ended:** 2026-04-21 13:31:32


---

**Session ended:** 2026-04-21 07:35:59


---

**Session ended:** 2026-04-20 16:07:48


## 2026-04-20 — Add contact_sed IC and T18/T19 tests; Mobility refactor; 17/17 pass

- Replaced `MobilitySub` with `Mobility` (Hermite-interpolated blending of ice/sed/air mobilities); added `mob_air` to AppCtx.
- Renamed parameter `met` → `sed` and `Fwat` → `Fsed` throughout headers and material_properties.c.
- Removed dead `mob_sed` variable and commented-out `flag_Tdep` mobility block from assembly.c.
- Added `FormInitialContactSedPermafrost2D`: two sediment grains just touching, each with a concentric ice shell; fixed SNES divergence at contact point by switching from additive sum to min-distance tanh union.
- Changed `gamma_im` 0.33→0.07; set `mob_sed = mob_air = mob_sub` so all phases share equal mobility by default.
- Added `test4_ContactSedPair.opts` (long run) and `test_T18_contact_sed_quick.opts` (10-step smoke run).
- Added T18 (2D contact-sed smoke) and T19 (2D IC accuracy with lens-area corrected expected ice volume); all 17/17 tests pass.

---

**Session ended:** 2026-04-19 14:24:40


## 2026-04-19 — Extend test suite to T11–T17 (15/15 pass)

- Added T11 (sublimation deceleration in finite domain), T12 (deposition at hum=1.5), T13 (temperature field consistency), T16 (mass conservation), T17 (Dirichlet temperature BC fix).
- Created 3 new opts files: `test_T11_sublim_rate.opts` (1000 steps), `test_T12_deposition.opts` (hum=1.5), `test_T17_temp_bc.opts` (flag_BC_Tfix=1).
- T11/T13/T16 share a single 1000-step sublimation run; T12 and T17 each have dedicated runs.
- Diagnosed `temp` monitor column as ∫T dx [°C·m], not domain-average — fixed T13/T17 criteria accordingly.
- Discovered and documented finite-domain vapour depletion: sublimation rate drops 21× over 0.1 s as accumulated vapour quenches the driving force.
- Updated `test/TEST_ANALYSIS.md` with full descriptions for T11–T17, updated cross-test analysis and parameter tables.
- All 15/15 tests pass.

---

**Session ended:** 2026-04-18 08:12:40


---

**Session ended:** 2026-04-18 08:05:57


## 2026-04-18 — Systematic test suite (T01–T10), all 10 pass

- Rewrote `preprocess/comp_eps.py`: added `compute_params()`, `rho_vs_sat()`, `patch_opts()`, and full CLI (`--patch`, `--quiet`); eps and Nx now computed from domain geometry and Kaempfer & Plapp bounds.
- Patched all 5 test opts files with `eps=9.3295e-07`, `Nx=152`.
- Fixed `test/run_tests.py`: corrected `_parse_monitor()` regex to reject 14-column SNES DOF rows, fixed `rho_vs()` to match C ASHRAE polynomial, replaced fragile `min(arr, default=)` with `_safe_min()`, removed `-output_path` arg.
- Fixed `src/permafrost2.c`: added `PetscOptionsInt("-flag_tIC", …)` so the flag is actually read from opts files.
- Fixed `src/initial_conditions.c`: `flag_tIC==2` (flat interface) now sets `n_actsed=0`; previously placed a large sediment slab filling the air region causing DIVERGED_STEP_REJECTED.
- Redesigned T07 (Bergeron): changed metric from undetectable ice-volume change to `tot_rhov` ratio at t=0; raised `grad_temp0` to `1e5 K/m` (ΔT=10 K, ratio=1.27). All 10 tests pass.

---

**Session ended:** 2026-04-17 19:44:58


---

**Session ended:** 2026-04-17 15:52:15


## 2026-04-17 — Phase-dependent mobility + 4-DOF SNES convergence fix

- Updated `MobilitySub` in `material_properties.c` from stub returning constants to `M(φ) = m₀·φ·(1−φ)` for both ice and sediment.
- In `assembly.c`: compute `mob_eff` / `mob_sed_eff` after mobility selection; split `C3` into `C3_ice` / `C3_sed` so each phase equation uses its own effective mobility in both the Laplacian and driving-force terms.
- Fixed `norm0[3]` → `norm0[4]` in `NASA_types.h` to prevent out-of-bounds writes during SNES convergence checking.
- Expanded all 3-DOF caps in `snes_convergence.c` to 4 (arrays, header, print format) so sediment (DOF 3) residuals are shown during Newton iterations.
- Extended `Integration` in `assembly.c` to return `S[6] = sed`; updated `Monitor` in `monitoring.c` to extract `tot_sed` and display a `TOT_SED` column in the time-step table.

---

**Session ended:** 2026-04-17 14:02:23


---

**Session ended:** 2026-04-17 11:55:47


## 2026-04-16 — Post-processing updates and bug fixes

- Updated `postprocess/analyze_interface.py`: read sediment DOF (sol[:,3]), fix air = 1−φ_i−φ_s, add sediment volume fraction to metrics and 4×2 panel plot
- Updated `postprocess/plot2D_snapshot.py`: expand `_plot_cuts` from 2×3 to 2×4 to include sediment cross-sections
- Updated `postprocess/plotpermafrost.py`: add SedPhase DOF export, fix file-exists path bug, add CLI flags
- Updated `scripts/plotpermafrost.py`: fix `❌ Error processing` bug caused by erroneous `np.newaxis` in air-phase computation
- Set up `.claude/ACTIVITY_LOG.md` and Stop hook for session logging; added activity log instruction to `CLAUDE.md`

---

**Session ended:** 2026-04-16 21:55:05


---

**Session ended:** 2026-04-16 21:51:32


---

**Session ended:** 2026-04-16 21:45:52


---

**Session ended:** 2026-04-16 21:42:34

# Activity Log — Permafrost Project

Newest entries appear at the top. Each session is separated by `---`.

---

## 2026-04-16 — Configure activity log

Set up `.claude/ACTIVITY_LOG.md` and a Stop hook to automatically timestamp each session end. Added instruction to `CLAUDE.md` requiring a summary entry before stopping.

---
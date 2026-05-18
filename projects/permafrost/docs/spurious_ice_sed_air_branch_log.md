# Branch log: `fix/spurious-ice-sed-air-penalty`

This document narrates the changes made on the branch `fix/spurious-ice-sed-air-penalty`,
ordered from most significant (architectural / model-level) to least significant
(cleanups). For each entry: *why* the change was made and *what happened* after.

Branch base: `c8ed4915` on `main` (the commit immediately before this branch started).

The original motivation for the branch was a specific failure: with the 2-phase
ice formulation under frozen sediment, **spurious ice was forming at sed-air
interfaces** — symmetric to the original spurious-air-at-ice-sed problem the
2-phase switch was meant to solve. Diagnosing and fixing that opened a cascade
of related model, solver, and workflow improvements summarized below.

---

## 1 — Use the 3-phase ice equation under frozen sediment   (`da57e1e`)

**Why.** The original code, after `t = t_sed_freeze`, switched the ice residual to
a 2-phase Kim-Steinbach form pinned to the frozen sed boundary. That form was
derived assuming the frozen sed always borders ice; at sed-air interfaces it
introduced a structural forcing that drove ice to grow where ice should be
zero. With asymmetric grains and large `dt`, this could trigger
`phi_ice → +1.6` and `phi_air → −0.6` excursions and crash the run.

**What happened.** Keeping the 3-phase ice equation when sediment is frozen
(only the sed equation degenerates to `sed_t = 0`) eliminated the
catastrophic sed-air spurious-ice mode. The 3-phase ice equation has the
correct symmetry — `ice = 0` is a stable equilibrium against a frozen sed
boundary just as it is against air. A mild residual `~ 10⁻³` spurious-air
floor at ice-sed remains (a known feature of the 3-phase form's slight
interface-width mismatch), but that's the trade-off baseline rather than a
catastrophic failure.

This is the single largest correctness fix on the branch.

---

## 2 — Three-file opts architecture: solver / geometry / experiment   (`0971afd`, `377b89c`)

**Why.** The pre-existing pattern (`universal.opts` + per-test `tests/*.opts`)
conflated three independent axes: numerical solver settings, geometric setup,
and environmental experiment parameters. Sweeping one axis required
duplicating everything else. Geometry-specific tweaks (`-mob_sub` reductions,
`-ksp_type bcgs`) were trapped inside specific tests rather than living with
the geometry they applied to.

**What happened.** Now a run is composed:
```
mpiexec ./permafrost \
    -options_file inputs/solver.opts \
    -options_file inputs/geometry/<geom>.opts \
    -options_file inputs/experiment/<exp>.opts
```
- `inputs/solver.opts` (renamed from `universal.opts`) holds 100% numerical settings.
- `inputs/geometry/` has 21 files covering all the canonical 1D and 2D geometries × resolution variants.
- `inputs/experiment/` has the environmental parameters (`t_final`, `t_sed_freeze`, `T`, `humidity`, `∇T`).

Sweep one axis at a time: trivial. Geometry-specific tuning (`-mob_sub
2e-10`, `-ksp_type bcgs` for stiff 2D hires) lives in the geometry file
itself.

Per-test `tests/test_1D_*.opts` and `tests/test_2D_*.opts` files were deleted;
the older legacy `test{1,2,4,5}_*.opts` (IceCap, CapillaryBridge,
ContactSedPair, SlabAndGrains) are still present awaiting a port.

---

## 3 — Mass-balance fix: drop `xi_v` from the vapor source term   (`3d3ab8c`)

**Why.** `Residual_A1` had
`vap_src = xi_v · ρ_ice · air_t`
with `xi_v = 1e-4`. That meant when ice grew by `Δm` of mass, vapor only lost
`xi_v·Δm = 10⁻⁴·Δm` — the other 99.99% appeared from nowhere. Symptom: in
the early SeparatedGrainPair runs, `TOT_ICE` drifted **+11%** over 75 simulated
days with no boundary inflow that could account for it.

**What happened.** Removing the `xi_v` (the Stefan-condition closure has no
business with a time-scaling factor) and updating the matching Jacobian rows
restored mass balance: `TOT_ICE` drift dropped from +11% to near zero. The
`xi_v` factor still scales the *regularization* terms (diffusion + penalty)
but not the conservative mass-exchange source. This is the second-largest
correctness fix.

---

## 4 — Decouple the vapor source from sediment motion   (`d6ca749`)

**Why.** After fix #3, `vap_src = ρ_ice · air_t = −ρ_ice · (ice_t + sed_t)`. The
`sed_t` term was generating spurious vapor whenever sed moved (e.g., during
the pre-pin AC relaxation window). Sediment doesn't sublime — sed motion
shouldn't drive vapor exchange.

**What happened.** Changed to `vap_src = −ρ_ice · ice_t` (only ice motion
contributes). In the post-pin regime (`sed_t = 0`) the two forms are
identical, but during pre-pin the new form correctly attributes vapor
exchange only to ice motion. Also dropped the corresponding `[vap, sed]`
shift contribution from the Jacobian. The "negative vapor accumulation
near the sediment" artifact that had been visible in ParaView disappeared.

---

## 5 — Bounds-rollback safety net + TSPreStep deferral   (`5bd50bd`, `b9e501b`)

**Why.** When a stiff transient produced a Newton-converged-but-physically-bad
state (e.g., `phi_ice > 1.1` somewhere), the original code aborted the entire
run. We wanted to instead roll the step back, halve `dt`, and retry — letting
the integrator self-heal across topology changes.

The first attempt called `TSRollBack(ts)` directly from the monitor. That
crashed with
`Vector ... was locked for read-only access in TSMonitor()`
on the second rollback in a row.

**What happened.** Final design:
- `Monitor()` detects an out-of-bounds state and sets a flag (`user->bounds_violated`)
  plus the proposed new `dt` (`user->bounds_new_dt`). No `vec_sol` mutation.
- A new `BoundsRollbackPreStep()` callback runs *before each TSStep*, where
  `ts->vec_sol` is writable. If the flag is set, it does the `TSRollBack` +
  `TSSetTimeStep` and clears the flag.

The two-tier safety net (this + #6 below) now works reliably; the
`Vec locked` crash is gone.

---

## 6 — SNES domain-error catch during Newton line search   (`0ae0e18`)

**Why.** Catching bad states *after* the step (via the rollback) only works if
the post-step solution is recognizable as bad. Some unstable spatial modes
look fine globally (no individual `phi` outside bounds) but represent an
unstable Newton basin that next-step dynamics amplify. The earlier
`test_2D_TouchingGrainPair` crash at t ≈ 19 days was exactly this — `dt` was
cut 125× and the instability persisted.

**What happened.** `Residual_A1` now checks the trial Newton iterate against
`[phase_lo, phase_hi]` before clamping. If any `phi` is outside, it calls
`SNESSetFunctionDomainError(user->snes)` — telling SNES the current line-search
trial is in an invalid region. Line search treats the evaluation as having
infinite residual, backtracks, and Newton retries with a smaller update.
Catches the instability *during* the iteration, before it ever becomes a
"converged" solution.

Plumbing: `AppCtx` gained a `SNES snes` field set in `permafrost2.c` right
after `TSGetSNES`.

---

## 7 — Restore physical xi_v and xi_T (= 1)   (`1557a36`, `2b226ab`)

**Why.** Old defaults were `xi_v = 1e-4`, `xi_T = 1e-2`, which artificially
throttled vapor diffusion 10,000× and thermal diffusion 100× below physical.
That pushed the model deeply into a diffusion-limited regime, which dry-snow
metamorphism is *not* — physically it sits in the mixed kinetics/diffusion
regime with `Sherwood ~ 0.01–0.1`.

**What happened.** Set both `xi_v = xi_T = 1` in `solver.opts`. Effective
penalty `xi_v · k_pen` increases 10⁴× and SNES has to work harder at
interfaces, but the model now lives in the right kinetic regime.
Surface kinetics (`alph_sub`) sets the rate — not the numerical throttle.
Also exposed both as CLI options for per-test override (e.g., dial down to
`1e-1` if Newton starts struggling on a particular config).

---

## 8 — Constant phase-independent mobility   (`e85e4ad`)

**Why.** `Mobility()` returned `mob_sub·ice + mob_sed·sed + mob_air·air`. With
`mob_sed = 0` and `mob_air = mob_sub`, this collapsed to
`mob_sub · (1 − sed)` — full mob in pure ice/air, but **0.5× at ice-sed
contacts**. That slowdown was originally protective (preventing AC from
dissolving sed during a 3-phase pre-pin window). With `t_sed_freeze = 0`
(sed pinned from the very first step), there is no pre-pin AC window — so
the slowdown only artificially reduced the ice AC rate at ice-sed contacts,
slowing legitimate sintering near sediment grains.

**What happened.** `Mobility()` now returns `mob_sub` everywhere. Ice AC runs
at its proper rate at ice-sed contacts; sintering near sediment is no longer
half-speed. The `(ice, sed)` arguments are retained in the signature for
caller compatibility but go unused.

---

## 9 — Sediment AC mobility inert by default (`mob_sed = 0`)   (`093d1a0`)

**Why.** There was a latent bug: `permafrost2.c` unconditionally set
`user.mob_sed = user.mob_sub` after the CLI option block, silently
overriding any `-mob_sed` value and contradicting the inline comment
("Sediment is inert by default"). Harmless while `mob_sub` was small, but
after restoring physical `d0_sub0/beta_sub0` (#10) the resulting
`mob_sub ~ 2e-9` propagated to sed and dissolved 7.6% of `TOT_SED` in 300 s.

**What happened.** Removed the unconditional override. `mob_sed` now defaults
to 0 (set by `PetscMemzero`), which means the local mobility
`mob_sub·ice + 0·sed + mob_air·air` vanishes inside the sed grain — sed AC
can't dissolve the interior. `-mob_sed <value>` still works as a CLI override
for tests that want a softer sed phase.

---

## 10 — Restore physical Gibbs-Thomson parameters   (`0de386c`)

**Why.** `d0_sub0` and `beta_sub0` had been reduced (1e-9 → 1e-11, 1.4e5 → 1e-3)
earlier in the campaign to suppress AC coarsening artifacts. That gave
`mob_sub ≈ 6e-11` — 40× slower than the physically motivated value.
Suppressed the artifacts but also suppressed the kinetics: Ostwald ripening
between a 12 µm and 25 µm grain barely moved over 100 simulated days.

**What happened.** Restored `d0_sub0 = 1e-9`, `beta_sub0 = 1.4e5`. `mob_sub`
jumped 40× to `~ 1.9e-9`, `alph_sub` decreased ~3× to `3.2e8`. Sintering and
ripening dynamics became visibly active on day-to-week timescales. The
sed-dissolution problem (#9) emerged at the same time and was fixed by
`mob_sed = 0`.

---

## 11 — `t_sed_freeze = 0` default   (`e698260`)

**Why.** With faster AC dynamics (#10), the pre-pin 3-phase relaxation window
let sed AC eat into the grain interior before pinning. A short `t_sed_freeze`
(10s) was tried but the cleanest result came from `t_sed_freeze = 0` — pin
sediment from t = 0 entirely.

**What happened.** With sed pinned from step 0, `TOT_SED` drift dropped to
`-0.001%` (essentially machine noise) on the SeparatedGrainPair test, vs
`-7.6%` with a 300s pre-pin window. The ice equation still operates in
3-phase mode against the pinned sed profile (per #1).

---

## 12 — Tightened phase-field bounds: `[-0.05, 1.05]`   (`5bd50bd`, `d203f02`)

**Why.** Original `[-0.25, 1.25]` was very permissive. In one of the early
TouchingGrainPair crashes, `phi_ice` peaked at 1.17 for ~500 monitor steps,
slipping below the 1.25 threshold and never triggering the rollback. The
spatial-mode instability propagated unchecked and the run accumulated
`−12.5%` ice mass before finally failing.

**What happened.** Tightening to `[-0.05, 1.05]` makes the rollback (#5) and
the SNES domain error (#6) fire on much smaller excursions, catching the
instability while it's still recoverable. Combined with the rollback +
domain-error mechanisms, the integrator now self-heals through topology
events that previously crashed it.

---

## 13 — `dtmax` cap   (`5bd50bd`, then `1e4 → 1e5` in `d203f02`)

**Why.** Default behavior set `dtmax = 0.5 · t_interv ≈ days`. With `dt`
growing 30%/step under healthy Newton convergence, `dt` reached `~10⁵ s`
quickly and the integrator skipped through topology transitions (merging
grains, etc.) in a single step. Newton "converged" but to the wrong
post-merge basin.

**What happened.** Started at `1e4`, raised to `1e5` after the rollback +
SNES safety nets were in place. The cap is now a soft ceiling that the
integrator approaches but rarely binds against; the safety nets handle the
stiff transients.

---

## 14 — `TOTAL_MASS` column + per-phase mass plots   (`8bab5ab`, `c29a091`, `1b2ec9a`)

**Why.** Monitoring `TOT_ICE`, `TOT_AIR`, `TOT_RHOV` independently was
misleading: a `-60%` drift on vapor sounded catastrophic but represented
`~10⁻¹³ kg` against a total system mass of `~10⁻⁶ kg`. Hard to evaluate
whether the model was actually conserving.

**What happened.** Monitor table gained a `TOTAL_MASS` column =
`ρ_ice · TOT_ICE + ρ_sed · TOT_SED + TOT_RHOV` plus its percent drift line.
Same in `plot_mass.py`, which now also produces a `<run>/mass_plots/`
subdirectory with one figure per phase (`total.png`, `ice.png`,
`sediment.png`, `vapor.png`) so subtle drifts visible at the dominant-phase
scale on the combined plot don't get crushed.

Downstream parsers (`plot_timestep.py`, `plotpermafrost.py`,
`tune_sed_freeze.py`) updated to handle the new 11-field monitor row.

---

## 15 — Insulating BCs as the explicit default   (`b910426`)

**Why.** `flag_BC_Tfix` defaulted to `PETSC_TRUE` in source (Dirichlet
temperature at boundaries) while `universal.opts` overrode it to 0. A new
opts file that forgot `-flag_BC_Tfix 0` would silently get Dirichlet T,
hiding mass loss as a BC artifact.

**What happened.** Source default flipped to `PETSC_FALSE` (insulating).
Added a startup banner that prints the resolved BC state per DOF, so the
question "what BCs ran?" has an obvious answer in `outp.txt`. The
documentation of the natural-Neumann mechanism (no `IGASetBoundaryValue` +
`if (pnt->atboundary) return 0` in the residual = zero-flux) lives in the
commit message.

---

## 16 — Single sediment grain IC + test   (`988e0ba`)

**Why.** With `mob_sed = 0` and `t_sed_freeze = 0`, a single sediment grain
in air should be exactly stationary. Useful as a direct verification of
the inertness fix.

**What happened.** Added `FormInitialSingleSedGrain1D` / `2D` (mirror of
`FormInitialSingleIceGrain*` with ice/sed swapped) and a matching
`-ic_type single_sed` dispatch. Confirmed: sed stays put.

---

## 17 — `mob_sub` and `xi_v` / `xi_T` exposed as CLI options   (`2b226ab`, `d203f02`)

**Why.** Tests with stiffer geometries (touching grains) needed to dial down
`mob_sub` without recompiling. `xi_v` was hard-coded too.

**What happened.** Both now overridable per-test. Used in the touching/
separated grain geometry files (`-mob_sub 2.0e-10`) to trade kinetics speed
for stability through merge events. The `-mob_sub` override hook runs after
the physical `mob_sub` computation so the formula stays correct as the
default.

---

## 18 — Linear solver swap to BCGS for hires merge configs   (`d203f02`)

**Why.** The hires touching/separated grain pairs hit repeated KSP
`DIVERGED_BREAKDOWN` with GMRES — the Jacobian becomes severely
ill-conditioned during the merge transition and the Krylov basis goes
unstable.

**What happened.** `-ksp_type bcgs` (Bi-CGSTAB) added to the hires geometry
files for touching and separated grain pairs. BCGS is more robust against
breakdowns; the simulations completed where GMRES had failed.

---

## 19 — Adaptive time stepping: `snes_max_it` raised to 25   (`b096fea`)

**Why.** With `xi_v = 1`, the effective penalty `xi_v · k_pen = 10⁵`
stiffens the Jacobian at interfaces. On fine 1D grids with stiff geometry
(touching/separated pairs), Newton sometimes needs 20+ iterations and was
hitting `max_it = 15` then thrashing through dt-reduction → retry loops.

**What happened.** Raised to 25 in `solver.opts`. Gives Newton enough
headroom without much wall-clock cost in the common case where it converges
in 3-5 iterations.

---

## 20 — Batch test runner   (`6538242`, then iterated)

**Why.** Manual one-at-a-time invocation of 20+ tests was painful and
inconsistent. Wanted a single command that compiles once, runs a curated
list, captures per-test output in a shared parent directory, and emits a
summary.

**What happened.** `scripts/Studio/run_batch_tests.sh` does exactly that.
Iterated several times to: (a) call `plotpermafrost.py` directly instead
of via a wrapper that required staged `postprocess/` (`1b2ec9a`); (b)
flag suspicious successes where `t_final < 60 s` as `OK?` instead of
plain `OK` (`5255bd7`); (c) adopt the new three-file composition format
(`0971afd`).

---

## 21 — Removed `np.clip()` on `AirPhase` in `.vts` output   (`dcb75c2`)

**Why.** `plotpermafrost.py` had been clipping `AirPhase = 1 − ice − sed`
to `[0, 1]` before write. That hid small negative excursions (the
`~−10⁻³` ice-sed residual) from ParaView even though the simulation knew
about them. The mismatch between BOUNDS output and ParaView visualizations
caused a lot of confusion.

**What happened.** Removed the clip. ParaView now shows the actual field
values; users set the colormap range manually if they want a cleaner
display.

---

## 22 — Smaller fixes and configuration tweaks

- `k_pen 1e7 → 1e5` and `Lambda 1e3 → 1e4` in solver defaults (`53f2f43`).
  Relaxed the vapor penalty so it doesn't dominate the residual; raised
  the triple-junction penalty so binary interfaces stay sharp.
- Unified all canonical test `t_final` to 1 day (`22f80f4`) — earlier
  values (10 days, 100 days) made smoke tests slow to iterate.
- Per-test debug `t_final = 1e-2` typos fixed in `b096fea`.
- `CLAUDE.md` commit-policy guidance (`3283217`) — explicit authorization
  for autonomous commits, overriding the harness default.

---

## Summary trajectory

The model went through three identifiable stages on this branch:

1. **Diagnose-and-fix** (commits 1–10): Each issue revealed the next.
   Spurious ice at sed-air → 3-phase ice eq fix → mass non-conservation
   visible → xi_v on vap_src fix → spurious vapor near sed → vap_src
   decoupled from sed_t. By #10 the model was both physically reasonable
   *and* numerically stable on standard configs.

2. **Stability hardening** (commits 11–19): The remaining failures were
   stiff-Jacobian / topology-change cases — touching grains during
   sintering, hires merging events. Two safety nets (bounds rollback in
   PreStep + SNES domain error in line search) caught these, supported
   by physical `xi_v=1`, constant mobility, BCGS for the hardest hires
   cases, and tighter phase bounds.

3. **Workflow / diagnostics** (commits 20–22 + the three-file refactor):
   Better tooling for visualization (per-phase mass plots, unclipped
   AirPhase), batch orchestration, and a clean three-axis opts
   architecture (solver × geometry × experiment).

The starting question — *why is spurious ice forming at sed-air interfaces?* —
turned out to be downstream of a fundamental model-formulation choice in the
frozen-sed branch of the residual. Once that was understood, the rest fell
out as a series of independent improvements to model fidelity, numerical
stability, and developer ergonomics.

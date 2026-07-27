# Branch log: `fix/spurious-ice-sed-air-penalty`

> ⚠️ **HISTORICAL (banner 2026-07-22).** Sediment-era analysis (pre-2026-06-13 fork) of the removed three-phase model. Kept as reference for Effort 2 (`studies/icy_regolith/explicit_sediment_phase/`), not a description of the current two-phase code.


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

```bash
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

## 23 — Vapor-penalty refinement and Ostwald-ripening diagnostics   (`e40ab95`, `cf70046`, `cc50c23`, `b5e3e6a`)

This is a self-contained second-phase of the branch — once stability was in
hand (§1–22) the model was numerically robust but had never demonstrated
clean Ostwald ripening. Diagnosing why led to a series of changes that
fundamentally rewrote the role of the `k_pen` term.

### 23a — Identifying that `k_pen` was suppressing Gibbs-Thomson   (`e40ab95`)

**Why.** Running the asymmetric `2D_separated_grains` test (12 µm + 25 µm
grains, 5 µm gap) at −20 °C for 30 days with `k_pen = 1e5` showed **no
Ostwald ripening**: contour plots of `phi_ice = 0.5` at t = 0 and t = final
were identical. Around the same time, ParaView showed an "outline of nearly-
saturated vapor" hugging both grain surfaces — a strong hint that something
was pinning rhov to its flat-interface equilibrium right where the curvature
correction should have been expressed.

The penalty form is
`vap_pen = xi_v · k_pen · g(ice+sed) · (rhov − rhov_eq)`
where `g(ice+sed)` is the standard smooth Heaviside. At the diffuse ice-air
interface mid-thickness (`ice+sed ≈ 0.5`), `g ≈ 0.25`, so the effective
penalty is `0.25 · k_pen = 2.5×10⁴`. That's strong enough at the *interface*
— exactly where the curvature-dependent rhov_eq needs to emerge — to drive
rhov back to its flat-interface value and erase the Gibbs-Thomson signal.

The original justification for the penalty had been *mass conservation* for
geometries with sediment cores enclosed in concentric ice shells, where
vapor in the shell interior was prone to drift. That motivation no longer
held with the post-§4 vapor equation (vap_src no longer coupled to sed
motion).

**What happened.** As a diagnostic, set `-k_pen 0` in `solver.opts`. Result
on the −20 °C run: still no measurable ripening (`TOT_ICE` +0.001 %, no
visible grain motion), but rhov inside the sediment interior drifted to
order −10⁻¹ over 30 days — large in magnitude but isolated to a region with
no physical effect (sub_src `∝ ice²·air² = 0` in pure solid). The diagnostic
confirmed that the penalty *was* what suppressed GT, but exposed a cosmetic
artifact in solid that needed addressing if the penalty was to stay
disabled.

### 23b — Localized penalty via `PenaltyWeight()`   (`cf70046`)

**Why.** Goal: keep some penalty so rhov stays near saturation in solid
interiors (no cosmetic drift), but turn it *off* at the diffuse ice-air
interface so GT can emerge.

**What happened.** Added `PenaltyWeight(phi)` in `material_properties.{h,c}`
— a `SmoothHeavisidePoly` of the shifted variable `(ice+sed − 0.85)/0.15`,
clamped to `[0, 1]`. This gives `g_pen = 0` for `ice+sed ≤ 0.85` (so zero
through the entire diffuse interface, where `ice+sed ≤ 0.5` at midthickness
and ramps to ~0.85 at the inner edge) and `g_pen = 1` for `ice+sed ≥ 1`
(deep solid). The smooth-Heaviside form keeps `dg_pen/dphi` continuous.

`Residual_A1` and `Jacobian_A1` now use `g_pen` / `dg_pen` for the *penalty*
terms only; the existing `g_phia = SmoothHeavisidePoly(ice)` weighting for
`difvap_pen` is untouched. `k_pen` restored to `1e3` (much weaker than the
old `1e5`; the localization does the work).

**Result on the −5 °C / 30-day run.** Sediment rhov well-behaved (no drift).
TOT_ICE drifted **−0.518 %** — initially appeared to be GT-driven sublimation
finally working. TOT_RHOV still drifted **−66 %**, no asymmetric grain
motion visible.

### 23c — Realizing the localized penalty was still a mass sink   (`b5e3e6a`)

**Why.** Re-running with `-k_pen 0` (PenaltyWeight infrastructure still in
place but the penalty contribution is zero) was supposed to be the
*reference*: it should give the cosmetic rhov drift in sed and nothing else.
Surprise result: TOT_ICE drift collapsed from −0.518 % to **+0.001 %**.

That meant the −0.518 % ice loss in the localized-penalty run was *fake* —
the penalty active at the inner edge of the interface (`ice+sed ∈ [0.85, 1]`)
was draining vapor from a region where `vap_src = −ρ_ice · ice_t > 0` was
adding it, creating a local `rhov < rhoI_vs` deficit that the bulk
`sub_src = alph_sub · ice² · air² · (rhov − rhoI_vs) / ρ_ice` then
"corrected" by sublimating ice. The penalty was, in effect, manufacturing a
counterfeit Gibbs-Thomson signal by acting as a directional mass sink.

**What happened.** Adopted `-k_pen 0` as a diagnostic configuration to
isolate the "raw" model physics from any penalty artifact. The
`PenaltyWeight` infrastructure stays in the codebase; it is later
re-enabled at moderate strength in §26 once the diffusion-direction
fix (§25) and `vap_src` reformulation (§26) exposed a real numerical
need for a deep-solid sink that the penalty correctly provides.

### 23d — Wide-separation geometry for clean LSW   (`cc50c23`)

**Why.** With `k_pen = 0` at −5 °C on the existing `2D_separated_grains`
geometry (`grain_sep = 5 µm`), running 30 days produced **both Ostwald
ripening and a visible neck** between the two grains. The 5 µm gap is
roughly twice the diffuse interface width (`~2·2eps ≈ 3 µm`) — the two
grains' interfaces overlap, so vapor transport between them happens partly
through the overlapping diffuse layer (Kuczynski-style sintering) rather
than purely through bulk air diffusion (LSW Ostwald ripening). The
combined-mode result is good evidence the model can produce both
phenomena, but isolating LSW for diagnostic comparison needs the grains
further apart.

**What happened.** Bumped `grain_sep` from `5e-6` to `2e-5` (4× wider, ~7×
the interface width). Domain `Ly = 125 µm` still fits the new geometry
(grain centers at ±28.5 µm). With this spacing the two grains' diffuse
interfaces are well-separated, so any mass transfer must traverse the bulk
air gap — pure LSW.

### 23e — Warmer-T diagnostic experiment   (`cf70046`)

**Why.** Even with the penalty fixed, the −20 °C runs show essentially no
ripening over 30 days. The bottleneck is the vapor budget: at −20 °C,
`rho_v_sat ≈ 1 × 10⁻³ kg/m³`, giving a total vapor mass in the domain of
~3.7 × 10⁻¹² kg/m. Even fully redistributed, that mass would only shift
TOT_ICE by ~10⁻⁶ relative — below the integral-monitor's detection
threshold for asymmetric ripening over 30 days.

**What happened.** Added `inputs/experiment/30day_T-5_h1.00.opts` (same as
the −20 °C variant but with `-temp -5.0`). At −5 °C, Clausius-Clapeyron
gives `rho_v_sat` ≈ 12× larger — much bigger vapor budget for mass
transfer in the same 30-day window. This is the experiment where the
penalty-free run finally produced visible Ostwald ripening on the
narrow-gap geometry, motivating §23d for clean LSW.

### 23f — HPC regression sweep   (`cc50c23`)

**Why.** With `k_pen` removed and `PenaltyWeight` infrastructure added,
the other geometries (single grains, ice-sed pair, touching grains, 1D
variants) need to be re-validated to confirm nothing else broke.

**What happened.** Added `scripts/HPC/submit_regression.sh`, which fires
nine geometries × `1day_T-20_h1.00` plus the wide-sep
`2D_separated_grains × 30day_T-5_h1.00` diagnostic, in parallel via
`sbatch`. Wraps `run_permafrost.sh` directly (the legacy
`submit_permafrost.sh` still expects the pre-§2 single-opts-file format).
The sweep covers `1D_*` (ice_sed_pair, ice_slab, separated_grains,
single_ice, touching_grains) and `2D_{single_ice, single_sed,
ice_sed_pair, touching_grains}` plus the LSW diagnostic.

### Summary of the §23 arc

Before §23: model was numerically stable but Ostwald ripening was being
secretly suppressed by the `k_pen` penalty, which had been included for a
mass-conservation reason that had since been fixed at a different layer.

After §23: the penalty is temporarily disabled (`k_pen = 0`) as a
diagnostic and ParaView confirms that at −5 °C on the narrow-gap
geometry the model produces **both** Ostwald ripening (small grain
shrinks toward the larger one) **and** necking (diffuse interfaces meet
and form a connecting bridge). Clean LSW isolation requires the wide-gap
geometry. The `PenaltyWeight` machinery is retained; §26 re-enables it
at a tighter band once §24–25 are in place.

---

## 24 — Explicit Gibbs-Thomson curvature dependence in `sub_src`   (`73aca7e`, `2c80a78`)

**Why.** After §23 the model could ripen at narrow grain separation (5 µm
gap), but the mechanism turned out to be sintering through overlapping
diffuse interfaces (Kuczynski-style) rather than bulk-air LSW Ostwald
ripening. At 20 µm separation no measurable mass transfer occurred even
over 90 days at −5 °C. The reason: the existing `sub_src` closure used
`(rhov − rhoI_vs)` with a flat-interface `rhoI_vs(T)`, so the local
equilibrium at a curved ice surface had no curvature dependence. Without
it, two grains of different curvature see *the same* equilibrium vapor
density and there is nothing to drive vapor from one to the other.

**What happened.** Added an explicit Gibbs-Thomson correction:

```
rhoI_vs_eff(x) = rhoI_vs(T) * (1 + d0_GT * kappa(x))
sub_src        = alph_sub * ice^2 * air^2 * (rhov - rhoI_vs_eff) / rho_ice
```

where `kappa = -div(grad_ice / |grad_ice|)` is computed from the phase
field's gradient and Hessian, and `d0_GT = γ_iv·v_m / (R_g·T)` is the
capillary length (~9.6×10⁻¹⁰ m for ice at −5 °C).

Implementation details:

- `include/NASA_types.h`: `AppCtx.d0_GT` field (default 0 = disabled,
  recovers flat-interface behavior bit-for-bit).
- `src/permafrost2.c`: `-d0_GT` CLI option + startup printout
  identifying whether GT is active.
- `include/material_properties.h`, `src/material_properties.c`:
  `Curvature()` function that returns `kappa` plus optional
  `dkappa_dg[]` (length `dim`) and `dkappa_dH[]` (length `dim*dim`)
  partials for the analytical Jacobian chain rule. Regularization
  `|grad|^2 + (0.01/eps)^2` keeps `kappa` bounded in bulk regions
  where `|grad_ice| → 0`. `dim == 1` short-circuits to `kappa = 0`
  (curvature is identically zero in 1D).
- `src/assembly.c::Residual_A1`: reads the Hessian via
  `IGAPointFormHess` (which uses `p->shape[2]`; the IGA's default
  `order` already matches the polynomial degree p=2 so no extra
  setup is required), computes `kappa`, and uses `rhoI_vs_eff` in
  `sub_src`.
- `src/assembly.c::Jacobian_A1`: reads `N2` (second-derivative shape
  functions) and `Curvature`'s partials, updates the existing
  `[ice, *]` sub_src rows to use `rhoI_vs_eff` and `d_rhovs_eff_dtem`,
  and adds a new analytical `J[ice, ice]` block coupling `phi_ice`
  through `N1[b][l] * dkappa/dg_l + N2[b][k][l] * dkappa/dH_{kl}`.
  (After §26 the same chain is wired into `J[vap, ice]` with the
  opposite sign because `vap_src = -rho_ice * sub_src` couples
  symmetrically.)

A new pair of experiment files exercises GT in an A/B sweep:
`inputs/experiment/30day_T-5_h1.00_GTphys.opts` (`-d0_GT 9.6e-10`,
physical magnitude) and `30day_T-5_h1.00_GTamp.opts` (`-d0_GT 1.0e-8`,
~10× amplified for an unambiguous diagnostic signal).

The Jacobian is analytical (not FD) because finite-difference Jacobians
have historically caused Newton convergence issues on this model — the
`d0_GT * kappa` correction is ~10⁻⁴ relative to `rhoI_vs`, well below
the FD step size noise.

---

## 25 — Vapor diffusivity penalty direction fix   (`436c1ac`)

**Why.** Enabling §24 with `d0_GT` non-zero quickly showed that the
diffuse-interface vapor field had a strange pattern in ParaView: a
"ring" of saturated vapor pinned tightly to each grain that didn't
diffuse outward across the bulk-air gap. Investigation revealed a
long-standing bug in the diffusivity penalty: the switch was on
`g(phi_ice)`, so the *penalty* (1e-8 × physical) applied in **air** and
*full physical* diffusivity applied in **ice** — the opposite of what
the `inputs/solver.opts` comment claimed, and the opposite of what is
physical.

Numerically, effective bulk-air vapor diffusion was ~3×10⁻¹³ m²/s
versus the physical ~3×10⁻⁵ m²/s — 8 orders of magnitude too slow. So
even with GT producing a curvature-dependent local equilibrium at each
grain's interface, the vapor field had no bulk-air channel to carry
mass from one grain to the other. LSW was throttled at the most basic
transport step.

**What happened.** In `Residual_A1` and `Jacobian_A1`, the diffusivity
weighting was changed from `g(phi_ice)` to `g(phi_ice + phi_sed)`
(renamed `g_phia → g_solid` for clarity), and the linear-combination
coefficient ordering was swapped so the penalty `D_pen` now multiplies
the solid-side weight and the full physical `difvap` multiplies the
air-side weight:

```
difvap_eff = D_pen * g(ice + sed) + difvap_raw * (1 - g(ice + sed))
```

In pure air this gives the full physical `difvap`; in pure ice or pure
sed it gives `D_pen = difvap_pen * difvap`; the diffuse interface gets
a smooth ramp.

The Jacobian's `d_difvap_eff_dice` and `d_difvap_eff_dtem` derivatives
were updated to match, and a new `d_difvap_eff_dsed` term was added to
the `[vap, sed]` block (previously zero because `g(phi_ice)` had no
sed dependence). The `inputs/solver.opts` comment matches the code's
behavior again.

---

## 26 — `vap_src` Stefan closure: pair with `sub_src`, restore localized penalty   (`8477df0`)

**Why.** The §25 fix gave vapor a fast diffusion channel in bulk air,
but exposed a second issue: `vap_src = −ρ_ice · ice_t` over-counts the
Stefan condition. `ice_t` includes both `sub_src` (real mass exchange
with vapor) and AC interface motion (mass-neutral rearrangement). At
ice-sed boundaries the AC contribution is nonzero — ice can move along
the boundary without sublimating — but it does not produce vapor
physically. The old `vap_src` formula injected vapor at those boundaries.
With §25's correct diffusivity, the inner edge of every ice-sed
interface had nowhere to put that spurious vapor: `g_solid ≈ 1`,
`air_eff` clamped to `air_lim = 10⁻⁶`, effective diffusion 3×10⁻¹⁹ m²/s.
Runs on `2D_separated_grains` blew up within seconds; `phi_ice` and
`phi_air` left bounds, `rhov` reached ~10⁴ kg/m³ inside the solid.

**What happened.** Two changes that together restore stability:

1. `vap_src = -rho_ice * sub_src` (replacing `-rho_ice * ice_t`).
   Algebraically `-alph_sub * ice² · air² * (rhov - rhoI_vs_eff)`. This
   pairs the vapor source directly with the sublimation closure, so:
   - mass exchange between ice and vapor is exactly equal and opposite
     by construction (pointwise mass balance);
   - the `ice² · air²` factor naturally localizes `vap_src` to the
     ice-air diffuse band and zeroes it at ice-sed boundaries (where
     `air = 0`) and in bulk solid (`ice` or `air` = 0).

   The Jacobian's `[vap, *]` block lost the old `shift * rho_ice * Na_Nb`
   contribution (no `ice_t` dependence anymore) and gained sub_src-
   derived `dloc_dice`, `dloc_dsed`, `d_rhovs_eff_dtem`, and (via GT)
   `dkappa/dg`, `dkappa/dH` terms — mirroring the `[ice, *]` block
   with opposite signs because `vap_src = -rho_ice * sub_src`.

2. `PenaltyWeight` tightened (`PENALTY_PHI_LO`: 0.85 → 0.90) and
   `-k_pen` restored from `0` to `1.0e3`. The new `vap_src`
   localization fixes the ice-sed boundary issue, but a smaller
   residual issue remains: the inner edge of the ice-air diffuse
   interface (`ice + sed ∈ [0.9, 1.0]`) has small but nonzero
   sub_src-driven source and very slow effective diffusion (it's
   solid-side), so `rhov` could still drift to large values over long
   integration windows. The penalty acts as a sink there, pinning
   `rhov` to `rhov_eq`. With the ramp pushed to `ice + sed > 0.9` the
   entire diffuse interface midpoint is *outside* the penalty support
   so the GT signal still emerges; the penalty bites only in the
   deep-solid edge of the diffuse band where there's no physics
   that depends on a free `rhov`.

Validates on the 30-day −5 °C wide-separation `2D_separated_grains`
run: no blowup, clean phase bounds, finite-`d0_GT` runs progress
normally.

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

## 27 — Diagnostic: both vapor-equation penalties OFF gives the cleanest results so far

**Setup.** `2D_separated_grains` (20 µm separation), 30 days at T = −20 °C,
h = 1.00, run with the new experiment file
`inputs/experiment/30day_T-20_h1.00_nopenalty.opts` which overrides:

```
-k_pen 0.0          # equilibrium penalty OFF
-difvap_pen 1.0     # D_pen = 1.0 * difvap_raw, so the diffusivity penalty
                    # linear combination collapses to full physical D_v
                    # everywhere, regardless of g_pen.
```

This was the first run with **no artificial stabilisers on the vapor field at
all** — vapor sees full physical diffusion in air, in the diffuse interface,
and inside solid; the rhov → rhov_eq equilibrium-sink term is zero
everywhere.

**Result.** Cleaner phase bounds, cleaner vapor field, no spurious migration,
no halos, no inner-band rhov drift. Better behaviour than every prior
configuration with either penalty active.

**Why this is the case.** Several things have changed between the
configurations that *needed* the penalties (§13–§23) and the current state of
the residual that does not:

1. **Mass-conserving `vap_src` (§26).** The original justification for
   `k_pen` was that `vap_src = -ρ_ice · ice_t` over-injected vapor at
   ice-sed boundaries (where AC motion is mass-neutral). The penalty was
   needed as a sink. Replacing `vap_src` with `-2·ρ_ice·air·ice_t` localised
   the source to the ice-air diffuse band and zeroed it at ice-sed
   boundaries. Once `vap_src` is no longer over-counting the Stefan
   condition, there is no spurious vapor for the penalty to sink.

2. **Diffusivity-penalty direction fix (§25).** When the penalty was
   inverted, bulk air had `D_eff ≈ 3×10⁻¹³` m²/s — eight orders of
   magnitude too slow. Vapor had no way to escape the diffuse interface, so
   it accumulated. The penalty looked load-bearing because it was draining
   an accumulation the inverted-diffusivity bug was causing. Direction
   fixed, bulk air now carries vapor at full physical `D_v`, accumulation
   never builds up.

3. **Latent heat paired with `sub_src` (this branch).** Pairing latent heat
   with `air_t = -(ice_t + sed_t)` conflated AC interface motion (mass-
   neutral) with real ice ↔ vapor exchange. That generated/absorbed
   spurious latent heat during the initial relaxation, distorted T, and
   shifted `ρ_vs(T)` and `rhov_eq`. With latent heat paired directly with
   `sub_src` instead, T is no longer perturbed by AC interface motion;
   `rhov_eq` and `ρ_vs` stay quiet at the diffuse interface, and the
   penalty has nothing to chase.

4. **NRmin=5 unlocked dt growth.** With dt stuck at ~0.013 s, the
   adaptive scheme was integrating very early-time transients forever.
   Many of those transients are exactly what the penalty was tuned to
   suppress (sub-millisecond AC equilibration overshoots). Once dt can
   grow into the seconds-to-hours range where the system is genuinely
   quasi-static, there are no fast transients left for the penalty to act
   on.

5. **Constant-anchor `rhov_eq = hum0 · ρ_vs(T)`.** The earlier formulation
   `rhov_eq = (ice+sed)·ρ_vs + air·rhov` was self-referential — the
   penalty pulled `rhov` toward itself in air, and toward `ρ_vs` in solid.
   That introduced a one-way drain from air-side to solid-side along the
   diffuse band. Replacing with a constant anchor `hum0·ρ_vs(T)` removed
   that drain, but also removed most of the penalty's purpose: there's no
   strong physical reason to pin `rhov` toward the initial humidity inside
   solid; the only motivation was numerical sink-of-last-resort for the
   inner-band accumulation that no longer happens (item 1).

**What this means.** Items 1–5 are independently-derived fixes for distinct
underlying bugs. The two vapor penalties were each calibrated to mask one or
more of those bugs. With all the bugs now fixed in their own right, the
penalties have nothing left to do — and turning them on adds artificial
weight to a vapor equation that is otherwise behaving correctly, which
manifests as suppressed Gibbs-Thomson signal in the diffuse interface
(see §23) and minor parasitic drains (§26).

**Recommendation, pending broader confirmation.** Adopt `k_pen = 0`,
`difvap_pen = 1.0` as the default in `solver.opts`. Keep the `PenaltyWeight`
infrastructure in the code as dead code (zero cost: the penalty terms
multiply by `k_pen = 0`) — both for clean revert in case of a not-yet-seen
failure mode, and because the same machinery is the right place to plug a
future *physical* surface-kinetics term if one is needed. Validate on the
two-grain `2D_touching_grains` sintering test and the `1D_separated_grains`
narrow-gap case before promoting it; those exercise different parts of the
residual than the wide-gap test that exposed the win.

---

## 28 — The penalty-off run froze: GT was off, and `vap_src` laundered AC mass loss

**Symptom.** The `2026-05-27 advisor_T-20_nopen` HPC batch (penalties off, the
§27 recommendation) was run on all seven advisor-slide geometries. Every 2D
case did a brief Allen-Cahn relaxation burst (`2D_touching_grains`: steps
≈123–129, t ≈ 3000–17000 s) — the touching grains forming a neck and the
interfaces settling — and then **froze bit-for-bit** for the remaining 99.4 %
of the 30 simulated days. `TOT_ICE`, `TOT_RHOV`, `TEMP` all constant to 7
significant figures; dt pinned at `dtmax = 1e4`.

**Diagnosis.** Two independent causes:

1. **Gibbs-Thomson was disabled (`d0_GT = 0`).** §24 added GT; §26→ later
   iterations (`9d9fdf4`) stripped it back out as a stability measure. With
   GT off, every grain shares one equilibrium vapor density `rhoI_vs(T)`, so
   there is *no curvature driving force* for vapor-mediated mass transfer.
   The early "coarsening" the run showed was **not** Ostwald ripening — it was
   one-time AC curvature relaxation of the initial condition, which is a
   gradient flow that stops at its local minimum. At saturation (h = 1.00)
   with GT off, `sub_src ≈ 0` once vapor equilibrates, so the frozen state is
   a genuine steady state of the model as configured. (It *looked* like the
   "no-diffusion-in-solid + active rhov_eq penalty" case because the end
   state — static vapor at uniform saturation — is identical, but the cause is
   the opposite: total absence of a driving force, not pinning.)

2. **Loose `snes_atol = 1e-4` amplified it.** At T = −20 °C the vapor field is
   O(`rho_vs`) ≈ 8.5e-4. Once dt grew to `dtmax`, the vapor residual
   (~`rhov/dt`) fell below `atol`, so SNES declared "converged" without ever
   resolving the slow dynamics — the snap to the static state was instant.

**Concern raised about `vap_src`.** Separately, the `vap_src =
-2·rho_ice·air·ice_t` form (from `f2c797b`) couples *all* ice motion to vapor,
including the mass-neutral AC curvature relaxation. During the burst this
laundered AC "mass loss" into spurious vapor — a measured ~0.6 % `TOTAL_MASS`
drift. The AC mass change is a numerical artifact and should **not** be
conserved by manufacturing vapor.

**Fix (`70999f5`).** Three coupled changes:

1. **`vap_src = -rho_ice · sub_src`** (Moure & Fu 2024). Driven by the
   phase-change term only, so AC curvature motion produces no vapor and the
   physical ice↔vapor exchange is exactly equal-and-opposite pointwise (the
   `rho_ice` cancels the `1/rho_ice` in `sub_src`). This reverts the
   `f2c797b` change and restores the §26 intent — but the §26 mass-leak worry
   that motivated `2·air·ice_t` is now explicitly retired: that "leak" was
   non-physical AC mass that should not be conserved. The Jacobian `[vap,*]`
   block was reverted to the `sub_src`-derived rows + GT block; the `ice_t`
   read in `Jacobian_A1` is gone.

2. **Re-enabled Gibbs-Thomson.** Restored the Hessian reads, `Curvature()`
   call, `rhoI_vs_eff = rhoI_vs·(1 + d0_GT·kappa)` in `sub_src`, the
   `[ice,ice]` GT Jacobian block, `N2` shape funs, and `rhoI_vs_eff` /
   `d_rhovs_eff_dtem` throughout the `[ice,*]` / `[tem,*]` / `[vap,*]` blocks.
   Still gated on `d0_GT` (0 = off). This is the LSW driver that prevents the
   freeze. New experiment file `30day_T-20_h1.00_nopen_GT.opts` sets
   `d0_GT = 9.6e-10` with penalties off.

3. **`snes_atol 1e-4 → 1e-6`** so weak cold-T vapor dynamics are resolved
   instead of skipped at large dt.

The thermal `[tem,*]` Jacobian keeps its modified-Newton approximation (no
explicit GT-curvature term and the empirical `[tem,tem]` sign from `4a82338`);
the GT correction there is a ~1e-4 relative perturbation of an already-small
latent-heat term, so it doesn't affect convergence. The residual is fully
GT-consistent.

Smoke-tested `2D_touching_grains`, T = −20, penalties off, `d0_GT = 9.6e-10`,
`t_final = 50 s`: clean exit, no NaN/divergence, `TOTAL_MASS` conserved to
display precision, dt grows normally under the tighter atol. Long-run
validation (does ripening actually proceed over 30 days?) is the next HPC step.

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

4. **Vapor-penalty refinement and Ostwald-ripening diagnostics** (§23,
   commits `e40ab95`, `cf70046`, `cc50c23`, `b5e3e6a`): once stability
   was in hand, the question shifted from *does the model run?* to
   *does it produce the right physics?* The `k_pen` penalty — included
   under an obsolete mass-conservation justification — turned out to be
   suppressing Gibbs-Thomson at the diffuse ice-air interface, killing
   Ostwald ripening. A localized-penalty experiment (`PenaltyWeight()`)
   exposed a subtler failure mode: even confined to the inner interface,
   the penalty was still manufacturing a counterfeit GT signal by
   draining vapor where `vap_src` was adding it. The penalty was
   *temporarily* disabled (`k_pen = 0`) as a diagnostic. ParaView
   confirmed both ripening and necking on the narrow-gap −5 °C run.
   The geometry was widened (`grain_sep` 5 → 20 µm) to isolate LSW.

5. **Enabling Ostwald ripening as a first-class physics: GT closure
   and the cascade of vapor-transport fixes it forced** (§24–26,
   commits `73aca7e`, `2c80a78`, `436c1ac`, `8477df0`): the 90-day
   wide-separation diagnostic showed no measurable LSW signal —
   the model lacked explicit curvature dependence in the local
   equilibrium. Adding it (§24) gave grains of different curvature
   different `rhoI_vs_eff` values, the missing driving force for
   bulk-air mass transfer. Turning it on immediately exposed two
   pre-existing model bugs that the earlier configuration had been
   masking: (a) the vapor diffusivity penalty was inverted — it
   applied in air rather than solid, throttling bulk-air diffusion
   by 8 orders of magnitude (§25); (b) `vap_src = -ρ_ice · ice_t`
   over-counted the Stefan condition by including mass-neutral AC
   interface motion, producing spurious vapor at ice-sed boundaries
   (§26). Fix (a) was a one-line direction swap. Fix (b) replaced
   `ice_t` with `sub_src` directly, so `vap_src` is paired
   one-for-one with `sub_src` and inherits its `ice² · air²`
   localizer. The `PenaltyWeight` infrastructure was then re-engaged
   at moderate strength (`k_pen = 1e3`) with the ramp pushed deeper
   into solid (`PENALTY_PHI_LO = 0.90`) to provide a numerical sink
   in the inner-band diffusion-starved zone without intruding on
   the GT-active midpoint of the diffuse interface.

6. **Both vapor penalties retired** (§27): once the bugs the penalties were
   masking were fixed independently (mass-conserving `vap_src`, correctly
   directed diffusivity, latent heat paired with `sub_src`, dt-growth
   unlocked, constant `rhov_eq` anchor), the penalties themselves became
   counter-productive — suppressing real GT signal and adding parasitic
   drains. The `2D_separated_grains` 30-day −20 °C diagnostic with
   `k_pen = 0` and `difvap_pen = 1.0` produced cleaner results than any
   prior configuration. The penalty infrastructure (`PenaltyWeight()`,
   the residual/Jacobian wiring) stays in the code as dead code at zero
   cost — useful both as a clean revert path and as the natural home for
   any future *physical* surface-kinetics term.

The starting question — *why is spurious ice forming at sed-air interfaces?* —
turned out to be downstream of a fundamental model-formulation choice in the
frozen-sed branch of the residual. Once that was understood, the rest fell
out as a series of independent improvements to model fidelity, numerical
stability, developer ergonomics, the central physics phenomenon of dry-snow
metamorphism (Ostwald ripening), a cascade of formerly-hidden vapor-transport
bugs that turning on real Ostwald ripening exposed, and finally the discovery
that the two vapor-equation penalties — originally added to fight those
hidden bugs — were no longer needed once the bugs themselves were fixed.

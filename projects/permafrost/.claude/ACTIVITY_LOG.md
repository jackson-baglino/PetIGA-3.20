
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
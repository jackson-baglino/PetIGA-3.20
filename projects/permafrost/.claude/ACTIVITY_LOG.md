
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
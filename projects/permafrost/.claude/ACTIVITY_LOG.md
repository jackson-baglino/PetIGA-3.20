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
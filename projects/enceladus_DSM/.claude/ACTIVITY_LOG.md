# Activity log — enceladus_DSM

Newest entries first.

## 2026-07-31 — Merge effective_thermal_cond in; k_eff becomes an in-line diagnostic

- Merged `projects/effective_thermal_cond` into this solver. k_eff is now
  computed during the time loop (`-keff 1`) on a cadence independent of field
  output (`-keff_freq` / `-keff_t_interv`), because a snapshot is hundreds of MB
  while a k_eff sample is ~200 bytes — the old coupling capped k_eff(t) at the
  snapshot count. `-keff_replay <dir>` covers finished runs and replaces the
  standalone binary, which is now marked superseded.
- New: `include/keff.h`, `src/keff{,_field,_cell,_solve,_sample}.c`. Guards
  hard-error on non-periodic meshes, `-geom_file`, `-axisym` and `dim < 2`.
- Renamed `NASA_types.h`/`NASA_main.h` → `enceladus_types.h`/`enceladus_main.h`;
  deleted dead `env_helper.c/h`; raised `TARGET_DOFS_PER_CORE` 50k → 80k.
- Deleted the write-only `alph[]`/`mob[]` Gauss-point arrays and standardised on
  one indexing convention, `point->index + point->count*point->parent->index`,
  documented at the top of `src/keff_field.c`. The read side is a PetIGA form
  callback, so a sequential counter was never an option there.
- Fixed `-ic_type ice_slab`, which built a circular blob of radius `Ly` rather
  than a slab and self-overlapped on a periodic cell.
- **The "iterative solvers are unreliable for this problem" claim was a bug**,
  not numerics: `IGACreateMat` overrides `MATOP_CREATE_VECS`/`MATOP_DUPLICATE`
  with IGA-aware versions, AMG's coarse operators inherit them without the
  composed IGA, and `PCGAMG` died in setup on every mesh. CG+GAMG now agrees
  with LU to all ten printed digits and is 17× faster at 322k unknowns.
- Verification suite in `studies/snow_thermal/verification/` — analytic gate,
  solver comparison, IC resolution study, figures, raw logs, README.
- Two measurement notes that matter for the study: on the production mesh rule
  the IC's volume-fraction error is `≈118·eps`, so the discretisation floor on
  the safety 0.5 vs 1.0 comparison is **0.024%**, not the ~0.2% first estimated;
  and the diffuse band (9.2·eps = 9.2 µm at safety 0.5, 18.4 µm at safety 1.0)
  exceeds the median pore throat (4.0–12.6 µm) in most shakedown packings, so
  throat closure is the mechanism to watch, biasing k_eff high at low porosity.

---

## 2026-07-27 — Project created from lunar_regolith_DSM

- Created `enceladus_DSM` as a verbatim copy of `lunar_regolith_DSM` (itself
  the renamed `sublimation_pf`), so both projects start from the same
  extensively-tested two-phase solver and diverge only as the two papers
  require.
- The previous `dry_snow_metamorphism` solver was **discarded**, not merged:
  it had a wrong latent-heat density (mixture `rho` instead of `rho_ice`),
  stale `xi_v`/`xi_T` (1e-5/1e-4 vs 1e-3/1.0), Libbrecht kinetics on by
  default, `eps` sized at 273.15 K on a 2.8x-too-coarse mesh, and `p=1/C=0`
  which is incompatible with `effective_thermal_cond`. Its final state is
  preserved at tag `archive/dry_snow_metamorphism-legacy`.
- The three-phase (ice/sediment/air) formulation is gone from both projects.
- `studies/snow_thermal/` moved here from the lunar project; it is this
  project's effective-thermal-conductivity study.

---

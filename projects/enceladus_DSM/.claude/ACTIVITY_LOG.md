# Activity log — enceladus_DSM

Newest entries first.

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

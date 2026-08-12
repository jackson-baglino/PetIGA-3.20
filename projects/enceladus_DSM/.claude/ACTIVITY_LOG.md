# Activity log — enceladus_DSM

Newest entries first.

## 2026-08-12 (later) — `--axisym` correction; the arms agree after all

- The fine arm's `neck_width.csv` had been extracted **without `--axisym`**.
  These are axisymmetric r–z runs, so the measured chord is the neck *radius*
  and the flag doubles it; every fine-arm width was half its true value.
- Corrected, the fine arm spans 9.3 → 51.1 µm (not 4.7 → 25.6). This withdraws
  the entry below: the arms differ by **1.9× at 32.8 µm, 1.2× by 51 µm**, not
  36–150×, and the fine arm *does* reach the experiment's starting width
  (15.05 h), so nothing needs extrapolating.
- The fine arm is above its own 31.6 µm floor over the compared range — the
  resolved arm the pair was built to produce. a = **0.283 ± 0.001**, inside the
  Demmenie band; the coarse arm (below its 60.1 µm floor throughout) gives
  0.229. Molaro's data gives 0.204 ± 0.053, so model and experiment do *not*
  agree on the exponent once the resolved arm is used.
- Rate offset vs Molaro: 152× (coarse), 156× (fine) — agreeing to 3 % across a
  3.6× change in ε, which is what marks it as the α_c choice, not a mesh effect.
- Rewrote `compare_mesh_pair.py` down to a single axes as asked: both arms, the
  Molaro points, a dotted power-law fit per arm, each simulation clock shifted
  so t = 0 lands where its own neck reaches 32.81 µm. No mechanism guide lines.
  Dropped the now-subsumed `meshpair_both_*` and `meshpair_aligned_*` sets.

---

## 2026-08-12 — Fine arm arrives: the mesh pair does not converge (WITHDRAWN — see above)

- Extracted the fine arm's neck curve (`neck_width.py`, 57 snapshots) and
  compared both arms against the Molaro −20 °C data.
- **The coarse mesh is not adequate.** At matched neck width the arms differ by
  36–150× in elapsed time, and re-zeroing at a shared neck does not collapse
  them — so it is not the coarse arm merely being ahead on one trajectory.
- The exponent is not converged either: equal over each arm's *own* range
  (0.19–0.21 vs 0.21–0.24), but 0.157–0.238 (coarse) vs 0.261–0.285 (fine) over
  the range they share. Refining ε moves the exponent toward 1/3; the fine arm
  lands inside the Demmenie band, the coarse arm does not.
- This supersedes 2026-08-11's headline: the coarse arm's agreement with the
  Molaro exponent, and the ~152× rate factor, are mesh properties. README
  updated where it previously said the opposite.
- Neither arm resolves its own neck (coarse 53.5 µm vs a 60.1 µm floor; fine
  25.6 µm vs 31.6 µm). The pair was expected to give one resolved arm, gave none.
- The requested alignment on the experiment's first neck (32.81 µm) works for
  the coarse arm (−7.74 h) and is impossible for the fine one — it tops out
  7.25 µm below that width, extrapolating to ~325 h (4× its run).
- New `analysis/compare_mesh_pair.py` (3-panel comparison + rate-ratio CSV).
- Fixed `--anchor-neck` in `fit_neck_growth.py`: it interpolated the anchor
  after the `t > 0` filter had dropped the t = 0 sample, which for a pre-necked
  experiment is the sample carrying the starting width — so it silently dropped
  the very dataset it was written to align. `Series` now keeps raw arrays.

---

## 2026-08-11 — mesh_pair coarse arm: exponent fits vs the Molaro −20 °C data

- Fitted the tangent-contact coarse run (79 h, α_c = 1e-3, ε = 0.87 µm) from
  `batch_2026-08-11__17.24.16_mesh_pair` against `molaro2019_fig11_T-20.csv`.
  a = 0.19–0.21 on all three fit forms, vs a = 0.18–0.25 for the data — inside
  the data's CI, at m ≈ 5, not Demmenie's m = 3.
- The tangent IC is the reason: it gives the run a physical clock zero, so the
  three forms agree instead of spanning 2x as they did on the pre-necked runs
  in `results/molaro_prenecked/` (a = 0.09–0.12 there).
- Recorded the caveats with the numbers: the whole curve sits below this arm's
  resolution floor √(12ε/R) = 0.346 (default protocol → no fittable window), so
  the fine arm settles it; the absolute rate is ~155x slow at α_c = 1e-3; the
  first ~15 samples are diffuse-interface relaxation off tangent contact.
- New `analysis/run_mesh_pair_fits.sh` (picks up the fine arm automatically) and
  `results/mesh_pair/` with three fit windows, both figures, and a README.
- Added `analysis/plot_neck_linear.py`: neck **width** vs time on linear axes,
  `w = A·t + C`, Molaro overlaid. A = 0.290 ± 0.025 µm/h over the shared width
  range (R² = 0.96); a line does *not* fit the full 79 h (R² = 0.77).
- Best result of the day: stretching the Molaro clock by one factor S = 152
  (elapsed time at equal neck size) puts the data on the model curve within its
  error bars — model vs experiment is a pure rate offset over 32.8–53.5 µm.
- The requested `w = A(t − t0) + C` is over-parameterized (t0 and C trade off
  exactly). Resolved it by anchoring: both clocks re-zeroed at a common neck
  width (32.81 µm; the model needs 7.74 h to get there from tangent contact),
  which fixes the origin from the data and leaves the one-parameter pinned form
  `w = A·t' + w_anchor`. A = 0.349 ± 0.027 µm/h model vs 52.91 ± 12.57 Molaro.
- That cross-validates the rate factor: pinned ratio 151.5× vs model-free
  elapsed-time ratio 151.8×, two independent routes to the same number. The
  free-slope 175× was an artifact of the data's truncated fit window.

---

## 2026-08-10 — Sintering growth exponent vs Demmenie et al. (2025)

- New target: Demmenie, Woutersen & Bonn, *J. Phys. Chem. Lett.* 16(8)
  2104–2109 (2025), doi:10.1021/acs.jpclett.5c00050 — two ~1 mm ice spheres at
  −3 °C in a box at ice-saturation, 2.5 h, giving r ~ t^alpha with alpha =
  0.29/0.33/0.30/0.26 (±0.01), i.e. Kuczynski m = 3, evaporation–condensation.
  That is this model's own transport route, and alpha is a *shape* observable,
  so it is independent of the ~50 % rate deficit against Molaro.
- Pulled `neck_width.py`, `plot_neck_convergence.py`,
  `calibrate_neck_geometry.py` and the Molaro validation CSVs out of the
  scratch quarantine.
- Added `postprocess/fit_neck_growth.py`. Reports three fit forms side by side
  (`d_free` = Demmenie's `C(t+t0)^a` with t0 free, `d_fixed`, Kuczynski
  `r^m − r0^m = Kt`) plus the local log-slope, because the exponent turns out
  to be mostly a property of the protocol: the same model curve gives a = 0.087
  over the solver's dense output, 0.078 at the data's sample times, 0.110 with
  t0 released.
- **Found `neck_width.py` could never read a snapshot.** Its local VTS parser
  called a `decode` helper it never imported; a bare `except Exception` then
  reported every NameError as "empty/corrupt file, likely an incomplete
  transfer" and the script exited 0 with a header-only CSV. Broke when pplib
  absorbed the decoder (a09f873) while the script sat in scratch/. Now
  delegates to `pplib.read_vts`, narrows the except, refuses to write an empty
  CSV.
- Added `-t_out_log N` (log-spaced snapshot cadence). Neither existing cadence
  spreads samples across decades, and one snapshot of the 9792×2671 sintering
  mesh is ~630 MB, so the budget is ~40 files. t0 defaults to `delt_t`, not
  `dtmin`: entries below the first accepted step are silently subtracted from
  the budget (`-t_out_log 8` gave 7 files before the fix).
- Stage-1 reanalysis (`studies/sinter_exponent/`, no compute): **the Molaro
  geometry cannot resolve an exponent at all.** The neck fillet `rho ~ r²/(2R)`
  is only resolved above `r/R >= sqrt(12·eps/R)` = 0.22 there, and the whole
  dataset spans 0.19–0.38 — a quarter decade sitting on the floor, over which
  1/3 and 1/7 are indistinguishable. The eps 0.60 and 0.86 µm arms have no
  fittable window whatsoever. Fix is bigger grains, not finer mesh: the floor
  falls only as sqrt(eps).
- Also: the Molaro *data* sits at a flat local slope of 0.185, between 1/7 and
  1/5, so it does not show Demmenie's exponent either — consistent with their
  claim that the 1/3-to-1/7 literature scatter tracks humidity control. And the
  model/data gap is far smaller on the r0-subtracting Kuczynski form (m = 5.32
  vs 4.04) than on the naive one (11.5 vs 5.4).
- Added five opts files: strict/coarse Demmenie arms (1 mm spheres, exactly
  tangent, `-ic_grain_union 1` — the additive IC would open with a spurious
  eps-dependent bridge), their shared −3 °C experiment file, and the Molaro
  pair restarted from tangent contact with an 8 h experiment. IC verified at
  the full pilot mesh: TOT_ICE = 1.047e-9 = the analytic two-sphere volume.
- Not yet run. Pilot batch is `studies/sinter_exponent/pilot_batch.txt`
  (~$5); production strict arm ~$60–110, gated on the pilot showing the neck
  climbs past r/R0 ≈ 0.15.

---

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

# Effort 3 — Snow thermal properties: detailed plan

Paper 2. Quantify how dry-snow metamorphism (DSM) reshapes icy granular packings
and how that changes the **effective (homogenized) thermal conductivity** k_eff.
Branch: `feat/effective-thermal-conductivity`. Shares the `sublimation_pf` master
solver; k_eff comes from the **existing** `projects/effective_thermal_cond`
solver (periodic homogenization) — no new estimator is written.

See `README.md` for the one-paragraph overview; this file is the execution plan.

---

## Stage A — Harvest reusable assets from `dry_snow_metamorphism`

`projects/dry_snow_metamorphism` is a sibling in the same monorepo. It is stale
(pre-rewrite) as a *solver*, but holds assets Effort 3 needs. Port these into the
shared tooling, then retire the directory (Stage F).

| Asset | Source | Destination | Note |
|---|---|---|---|
| 36 mm-scale 3D packings | `dry_snow_metamorphism/inputs/grains__phi=*` | `studies/snow_thermal/packings/` | φ 0.24–0.30 × 1/1.5/2/3 mm × seeds 7/21/22; each has `metadata.json` (achieved porosity), `grains.dat/.opts` |
| Packing index | `.../inputs/manifest.csv` | `studies/snow_thermal/packings/` | achieved-vs-target porosity table |
| FCC packing generator | `.../preprocess/FCC_grain_packing.py` | `preprocess/` | unique; not in sublimation_pf |
| Tangency circle packing | `.../preprocess/circle_packing_touching.py` | `preprocess/` | unique |
| Memory/rank planner | `.../preprocess/petsc_resource_planner.py` | `preprocess/` | estimates mem/ranks/mem-per-cpu; nothing equivalent exists |
| Band-clamping node planner | `.../scripts/HPC/batch_dsm.sh:113-200` | fold into `scripts/lib/alloc.sh` + submit scripts | clamps tasks/node into [28,32], adjusts node count both ways — more robust than the current fixed-32 |
| Molaro digitized curves | `.../postprocess/MolaroSims/wpd_datasets{20,5}.csv` | `postprocess/` (if not already) | validation reference |
| Clean options parsing | `.../src/options_helper.c` (77 lines) | consider adopting in `src/` | one clean PetscOptionsBegin/End vs the inlined parsing in permafrost2.c — optional refactor |

**Do NOT port** DSM's ξ_v/ξ_T (its `CLAUDE.md` says 1.0 but source sets
1e-5/1e-4 — inconsistent). Effort 3 uses `sublimation_pf`'s ξ_T=1, ξ_v=1e-3.

**Compatibility check:** DSM packings are 3-DOF-compatible geometries, but their
`.opts` predate the current option set (they carry sediment-era + `.env`-era
flags). Run each through `scripts/check_ic_types.sh` and the option scan; convert
`-ic_type`/removed options as needed. The packings are 3D (`-dim 3`) — confirm
the 3D IC path (`multi_grains` / grain-file readers) works in the current solver.

## Stage B — DSM runs on packings (large domains)

- These are mm-scale, up to 1713×1713×52 (~1.5×10⁸ DOF) — far larger than
  anything the solver has run (its largest is 1.35 mm × 0.27 mm, all 2D). Expect
  to raise the target: `TARGET_DOFS_PER_CORE=80000 ./scripts/HPC/submit_*` (the
  env override added in `scripts/lib/alloc.sh`), and lean on `--half-cores`.
- **eps per packing/temperature**: recompute with `comp_eps.py` for each packing's
  grain size and the run temperature. The new `-eps_valid_temp` guard (now in the
  solver) will abort a mismatched `-temp`, so emit it in the generated `.opts`.
- Gate B: a short run on the smallest (1 mm) packing that passes mass
  conservation (`postprocess/plot_mass.py`) before committing HPC time to the
  3 mm cases. Run on HPC, not locally.

## Stage C — Effective conductivity via `effective_thermal_cond` (the handoff)

The solver exists (`effective_thermal_cond/src/effective_k_ice_homog.c`): it
solves the periodic cell problem `-div(k∇t_m)=div(k·e_m)` per direction and
integrates the Voigt k_eff tensor. The work is streamlining the snow→k_eff
handoff (currently `field_init.c:224` reads `<init_dir>/igasol.dat` +
`sol_%05d.dat`, and it still uses the `.env→.opts` `gen_opts.py` path):

1. Point k_eff directly at a `sublimation_pf` run directory. The snow solution is
   3-DOF (ice, T, vapor); k_eff needs only the **ice** field — make the DOF-offset
   explicit in `field_init.c` (confirm it reads component 0, not a 1-field file).
2. Share domain/eps metadata so k_eff's cell matches the snow mesh without a
   hand-copied opts; retire `gen_opts.py`'s `.env` dependence in favour of the
   snow run's staged `.opts`.
3. Keep k_eff a **separate solver** (different equation, different build). Add a
   thin adapter in `studies/snow_thermal/analysis/` that locates snapshots and
   launches k_eff over a time series → k_eff(t).

## Stage D — Validation gate (before any science claim)

Reproduce analytic k_eff bounds at known porosity **before** trusting
microstructure results:
- Wiener bounds (arithmetic/harmonic means) — loose.
- Hashin–Shtrikman bounds — tight, the real test.
The 36 harvested packings have recorded porosities; run k_eff on the *initial*
(pre-DSM) microstructures and confirm k_eff falls within HS bounds for each φ.
Only then report k_eff(t) evolution under DSM.

## Stage E — Relaxation / safety-factor testing

Systematically test how far `comp_eps.py --safety` (default 0.5) can be relaxed
for these mm-scale domains, and which of the four K&P bounds actually binds
(reuse the recorded s100/s075/s050/s025 sweep pattern). Larger safety = coarser
mesh = cheaper; find the loosest value that preserves k_eff to tolerance.

## Stage F — Retire `dry_snow_metamorphism`

After harvest: tag its final state (`git tag archive/dry_snow_metamorphism-final`),
replace the directory with a README stub pointing to `studies/snow_thermal/` and
the harvested tooling. Nothing lost (tag + history), no stale second codebase.

---

## Open questions to settle before Stage B

- Which temperature(s)? The validated kinetics are at −20 °C; lunar/snow science
  may want a series. Each temperature needs its own eps/mesh per packing.
- 2D or 3D k_eff? The packings are 3D; k_eff supports it but cost is high.
  Possibly start 2D slices for the validation gate, then 3D for production.
- Which porosity × seed subset is the production matrix vs. exploratory?

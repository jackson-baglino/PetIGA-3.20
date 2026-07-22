# Study: Snow thermal properties (Paper 2)

**Scientific question.** How does dry-snow metamorphism (DSM) reshape icy
granular packings over time, and how does the resulting microstructural
evolution change the **effective (homogenized) bulk thermal conductivity** of
the pack? Quantify k_eff(t) across porosity / temperature / packing.

**This is a separate paper** from the icy-regolith study (`../icy_regolith/`),
sharing the same master solver.

**Pipeline.**
1. Generate / harvest mm-scale granular packings → `packings/`.
   Source: the 36 pre-generated 3D packings in
   `../../../dry_snow_metamorphism/inputs/grains__phi=*` (φ 0.24–0.30 ×
   1/1.5/2/3 mm × seeds 7/21/22), plus `FCC_grain_packing.py` /
   `circle_packing_touching.py`.
2. Run DSM to evolve the microstructure → opts in `opts/`. Use the master
   model's ξ_v = 1e-3 / ξ_T = 1 (NOT dry_snow_metamorphism's stale values).
3. Compute k_eff on snapshots with the **existing** homogenization solver at
   `../../../effective_thermal_cond` (periodic cell problem; reads
   `igasol.dat` + `sol_NNNNN.dat`). Effort 3 is streamlining that handoff, not
   writing a new estimator. → `analysis/`.

**Parameter regime.** Few-mm domains (up to 3 mm, ~1.5×10⁸ DOF in 3D) — much
larger than anything the master model has run (its largest is 1.35 × 0.27 mm,
all 2D). Expect to raise the DOF/core target (env var
`TARGET_DOFS_PER_CORE=80000`) and relax the `comp_eps.py --safety` factor;
testing how far is part of the work.

**Status (2026-07-21).** Not started. Branch `feat/effective-thermal-conductivity`.
First gate: reproduce analytic k_eff bounds (Wiener / Hashin–Shtrikman) at known
porosity before trusting microstructure results.

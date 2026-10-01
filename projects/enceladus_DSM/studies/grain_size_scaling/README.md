# Grain size: what it does to the mesh, the interface and the cost

`grain_size_scaling.py` sweeps the mean grain radius from 0.05 µm to 500 µm through
the campaign's own sizing code (`preprocess/comp_eps.py` `compute_eps`,
`generate_study_opts.derived_vn`/`tau_sub_of`). The settings are −20 °C,
α_c = 1e-3, safety 0.5 and L = 40 R_ave. At 50 µm it reproduces the
production run exactly: eps = 1 µm, Nx = 2829, τ_sub = 7822 s, ≈ $3.9 per run
(measured $3.85). Outputs are `grain_size_scaling.png/.pdf` and
`grain_size_scaling.csv`.

![grain size scaling](grain_size_scaling.png)

## The mesh is scale-free

The campaign sizes eps from R_feat = R_ave/25, so eps = R_ave/50, and the box is
L = 40 R_ave. Both rules are ratios, so the element count does not depend on
grain size: **Nx = 2829 and 24.0 M unknowns at every R_ave**. The B-KINETIC
bound binds across the whole range. The heat and vapour ceilings
(D·β_HK ≈ 59 / 139 µm) would only take over above R_ave ≈ 1.5 mm. Everything
physical shrinks with the grains:

| R_ave | L | eps | element h | visible band (1–99 %) |
|---|---|---|---|---|
| 500 µm | 20 mm | 10 µm | 7.1 µm | 92 µm |
| 50 µm (campaign) | 2 mm | 1 µm | 0.71 µm | 9.2 µm |
| 1 µm | 40 µm | 20 nm | 14 nm | 0.18 µm |
| **0.1 µm** | **4 µm** | **2 nm** | **1.4 nm** | **18 nm** |
| 0.05 µm | 2 µm | 1 nm | 0.7 nm | 9 nm |

Keeping today's absolute mesh spacing (eps = 1 µm) instead is not an option in
either direction:
- Below R_ave = 20 µm, eps would exceed R/20 and the grains would not be
  resolved.
- At 500 µm the mesh would need 2.4 × 10⁹ unknowns, about $10⁵ per run.

## The cost is scale-free too, while the run stays attachment-limited

All the model's timescales scale as R²: τ_sub ∝ eps², and the curvature-driven
interface speed ∝ 1/R over a distance ∝ R. So the time step cap and the time
to reach the same sintering state shrink together, and the **step count stays
at ~300**. The simulated time is the thing that changes. Our 30-day state at
50 µm corresponds to:

| R_ave | simulated time | dtmax | steps | $/run (−20 °C) |
|---|---|---|---|---|
| 500 µm | ~28 yr | 9.8 d | ~1100 | ~$11 |
| 50 µm | 30 d | 2.4 h | 374 | $3.9 |
| 10 µm | 22.6 h | 5.6 min | 310 | $3.3 |
| **0.1 µm** | **7.6 s** | **0.034 s** | **294** | **$3.2** |

The departure above ~50 µm is vapour transport. The crossover is
L* = D_v·β_HK = 139 µm at α_c = 1e-3, and the series estimate
t ∝ R²(1 + R/L*) grows faster than R² past it. The equivalent time is that
estimate, normalized to our 30 d at 50 µm. The step count assumes the run
spends most steps at dtmax, as the campaign runs do (79%).

## At 0.1 µm the computation would be fine; the physics would not be

The model would run at R_ave = 0.1 µm for the same ~$3, but four of its
assumptions break at that size:

1. **Vapour is no longer a continuum.**
   - The air mean free path at −20 °C and 1 atm is ~57 nm, about the median
     pore throat (≈ 0.49 R = 49 nm). The Knudsen number is ~1, so vapour
     diffusion is transitional, not the continuum D_v the model uses.
   - A 100 nm pore holds only ~28 vapour molecules at saturation, so a
     continuous ρ_v field means little.
2. **Gas conduction is suppressed.** At Kn ~ 1, air in the pores conducts several
   times less than bulk air (the Knudsen effect), and k_eff uses bulk k_air.
   k_eff is otherwise scale-free (it depends only on geometry and k_ice/k_air),
   so this is the one change at fixed geometry.
3. **Surface diffusion would dominate sintering.**
   - Herring's scaling has the time for vapour transport (evaporation–condensation)
     ∝ R² and for surface diffusion ∝ R⁴. Surface diffusion therefore gains a
     factor ~(50/0.1)² ≈ 2.5 × 10⁵ in relative speed going from 50 µm to 0.1 µm.
   - The Molaro validation already puts the missing surface-diffusion share near
     50% at ~100 µm grains (`studies/molaro_2019/`). At 0.1 µm the vapour-only
     model would be missing nearly the whole mechanism.
4. **The mesh approaches molecular scale.** h = 1.4 nm is about 4 water molecules,
   and eps = 2 nm is twice the capillary length d₀ = 1.0 nm. Gibbs–Thomson is
   still linear for the grains (d₀/R = 1%), but not for fine necks: a 10 nm neck
   radius gives d₀κ ≈ 10%.

The grey band in the figure (R_ave ≲ 1 µm) marks where points 1–4 begin to
matter. Point 3 is already significant at the campaign size. The blue band
marks R > L*, where the run becomes vapour-diffusion-limited and temperature
would no longer act as a pure time rescaling.

**Bottom line.** Changing grain size changes no mesh and almost no cost. It
rescales simulated time as R² (30 d → 7.6 s at 0.1 µm). Smaller grains would
need a different physical model (Knudsen transport and gas conduction, surface
diffusion), not a different mesh.

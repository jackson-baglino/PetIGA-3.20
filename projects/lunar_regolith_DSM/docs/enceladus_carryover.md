# What was carried over from `enceladus_DSM`, and why

2026-10-01. `enceladus_DSM` was split from this project on 2026-07-27 and the
two share `src/` apart from the main file. Several results that bear directly
on the wedge and contact-angle work were produced and recorded over there —
most of them on a *lunar* run, the 2026-08-07 `wedge_bc` batch. This file
records what was brought back so the history lives with the project it
describes. The enceladus activity log
(`../enceladus_DSM/.claude/ACTIVITY_LOG.md`, entries 2026-08-12 to 2026-08-14
and 2026-09-02) is the original source.

## Files copied

| here | from enceladus | what it is |
|---|---|---|
| `docs/curvature_driven_growth.md` | same path | derivation of the interface condition the solver obeys |
| `docs/figures/gt_balance.{csv,png}` | same path | the Gibbs–Thomson balance measured on the wedge run |
| `postprocess/gt_balance.py` | same path | tests `σ = d0·κ + β·v_n` on the φ = 0.5 contour of any run |
| `scripts/paraview_macros/plot_rhovsI.py` | same path | **replaces** the older lunar copy; adds curvature and GT fields |
| `scripts/paraview_macros/plot_growth_rate.py` | same path | `v_n` as a ParaView field, split into curvature and phase-change parts |
| `scripts/paraview_macros/verify_curvature.py` | same path | analytic gate for the curvature computation |

References to `enceladus_main.c` in the copied files were retargeted to
`lunar_main.c`. The lunar `plot_rhovsI.py` was the 2-commit ancestor of the
enceladus one, so nothing lunar-specific was lost by replacing it.

Not copied, because they do not touch this project: the k_eff homogenisation
campaign, the Molaro and sintering-exponent studies, the packing generator,
and the manuscript figure tooling.

## Findings that matter for the velocity study

**1. Gibbs–Thomson is emergent, not imposed (2026-08-13).** Capillarity is
carried by the Allen–Cahn term, not by ρ_vs. Projecting the residual onto the
translation mode gives

```
v_n = −3·M·ε·κ + (ε·α_sub / (5·ρ_ice))·(ρ_v − ρ_vs)
σ   = d0·κ + β·v_n,   d0 = 15·M·ρ_ice/(α_sub·ρ_vs),   β = 5·ρ_ice/(ε·α_sub·ρ_vs)
```

On the wedge run the median `|σ − d0·κ|/|d0·κ|` falls from 1.000 at t = 0 to
0.103 at 90 d: the interface *relaxes onto* the GT line. Do **not** add an
explicit `d0·κ` term to `RhoVS_I`; it would double-count. The lever is
`-d0_sub0`.

**2. d0 is exact; β runs 22 % high.** d0 from the run's own M and α_sub is
1.0166e-9 m, against the requested 1.0166e-9 and the physical
γ·V_m/(R·T) = 1.0168e-9. β is 22 % above the requested `-beta_sub0`, and the
gap equals the thin-interface correction
`5·a2·ε·(1/diff_sub + 1/dif_vap)·ρ_ice/ρ_vs` exactly. This is the same 1.22
factor `wedge_gt_velocity.py` found independently, and it is why
`studies/contact_angle/velocity_plan_2026-10-01/` uses β_eff = 1.22·β_sub0.
`a2` enters β only, never d0. The offset is additive (8.7e5 s/m at this ε),
not a fixed percentage: +11 % at α_c = 1e-3, +110 % at 1e-2.

**3. v_n is a small residual of two large terms.** On the wedge run the
curvature and phase-change contributions are ~1e-11 m/s each and cancel by a
median of 88.5 %, leaving v_n ~ 5e-13 m/s (concave +5.2e-13, convex −5.0e-13).
Any velocity comparison therefore needs σ and κ to better than 1 % each;
`wedge_gt_velocity.py` makes the same point about σ.

**4. Why a growing interface reads undersaturated.** At the concave meniscus
d0·κ = −5.05e-6 and the measured σ = −4.71e-6: undersaturated against a flat
surface, but *super*saturated by +3.4e-7 against its own depressed
equilibrium. This is the picture behind the zero crossings at σ∞ = d0·χ in
batches B and D.

**5. Curvature must be computed as −div(n), normalising first (2026-08-12).**
Two earlier forms were wrong in ways that looked like physics: a Tikhonov
regularisation left a residue ∝ φ'' that flipped sign across φ = 0.5 and
reached −1e5 /m against a real 4.5e3; and differencing discretised second
derivatives amplified error by `|1 − 2φ|/(ε·κ)` at the band edges ("bunny
ears", 36 % peak-to-trough). The macro now builds the unit normal pointwise
and differentiates that. Curvature does **not** change sign across one
interface. Read curvature on the φ = 0.5 contour, not as a band median.

**6. β is set by α_c.** β_sub ∝ 1/α_c. The September contact-angle runs here
used `-beta_sub0 3.9704e6` (α_c = 2e-3 at −20 °C); the enceladus k_eff campaign
settled on a constant α_c = 1e-3. From 2026-10-01 the velocity study runs
every batch at α_c = 1e-3 (lunar) and 1e-2 (laboratory). The literature band
is 1e-3 < α_c < 1e-1.

**7. ParaView macro trap.** `plot_rhovsI` runs a Programmable Filter that
star-imports over `max`/`min`/`sum` in `__main__`, which macros share. Macros
that run after it must bind those from `builtins`.

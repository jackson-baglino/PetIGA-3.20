# Manuscript plan — figures and sections (k_eff part)

Written 2026-10-06 against the draft `Literature/BaglinoJackson_DSM_manuscript_v2.1.pdf`
(34 pp., 8 figures). Figures are built now on the data in hand and rebuilt when
3f and the redo stage land. Claims and their limits: `TODO.md`, "Manuscript
argument". Run inventory: `RUN_TABLE.md`.

## Figure plan (8 in the draft → 8 main + supplement)

Figures 6 and 7 share one normalization since 2026-10-08: each run's value at the reference age θ_r = 30 (subscript r), not the opening frame. Figure 7 is the four-panel version (k_eff(t) by porosity + state law), chosen by the user 2026-10-08; the three-panel `figB_closure` is no longer a manuscript figure.
Numbering: fixed 2026-10-06, then the Demmenie figure inserted as 4 on 2026-10-07 (everything after it moved up one). All main-text figures are collected flat in `LOC/ManuscriptFigures/` by `figures/collect_manuscript_figures.py`.
Original note: numbering fixed 2026-10-06; the manuscript folders carry these numbers. Gallery and collapse swapped the same day (user): show what the model does to the aggregates, then the temperature collapse, then the weaker porosity result. NOTE: the iCloud folders and files for 4 and 5 still carry the OLD numbers until the rename there is done (access was denied from the session).

| new | content | status | replaces | section |
|---|---|---|---|---|
| 1 | (a) sintering mechanisms schematic; (b) `aggregate_strip`: close-up of a few dozen grains of the master run at 0/10/19/30 d. Decided 2026-10-06: one figure, two panels; the `setting_panel` sample is not used | **built**: `figures/fig1_mechanisms_aggregate.py` → `Figure1__SinteringMechanisms/Figure1_mechanisms_aggregate.*` | Fig. 1 | 1 |
| 2 | **Homogenization method** on the master run: (a) cell, (b) zoom with mesh and φ = 0.01/0.5/0.99, (c) corrector t_x, (d) local heat flux per unit mean gradient | **built** with the corrector replay: `figures/fig_Figure2_homogenization_method.py` → `Figure2__HomogenizationMethod/` | new | 2.3 |
| 3 | Molaro grain-pair validation: (a) sections, (b) neck width, (c)(d) grain shrinkage D/D₀ of the large and small grain (`Figure3_molaro_validation`, 170 × 181 mm; `molaro_combined` is the version without shrinkage). Note for the text: the −20 °C wall humidity was set from the measured large-grain shrinkage, so that panel is a calibration, not a test; the small grain and the −5 °C pair are the independent checks. Add the Demmenie-conditions exponent when that run lands | built: `postprocess/plot_molaro_validation.py`; Demmenie run queued | Fig. 2 | 3.1 |
| 4 | **Neck growth at saturation against t^(1/3)** (Demmenie conditions): laid out as the companion of Fig. 3: (a) the pair at 0 and 100 h, ice over the supersaturation, instants 1–2; (b) neck width with the free and the one-third fit. Two panels only (user, 2026-10-07); the exponent-vs-window and grain-size plots are diagnostics in the run folder. Result to state plainly: a = 0.21 after relaxation (0.21–0.23 by window), below 0.26–0.33 | **built**: `studies/molaro_2019/demmenie/fig_demmenie.py` → `Figure4_saturated_neck_growth.*` | new | 3.1 |
| 5 | **What happens to the aggregates**: five porosities × four instants, qualitative colour bars | built (sample): `Figure5_gallery` | Fig. 7 | 3.2 |
| 6 | **k_eff evolution and its temperature collapse**: (a) four snapshots of one packing; (b) k/k₀ vs days, five temperatures; (c) the same vs θ = t/τ_sub; (d) SSA/SSA₀ vs θ | **built**: `figures/fig_keff_collapse.py` → `Figure6__KeffEvolution/Figure6_keff_collapse.*` | Figs. 3, 4, 5 | 3.2 |
|   | *Alternate (`Figure6_keff_collapse_alt`, `fig_keff_collapse.py --other-seeds`)*: the other four packings of φ 0.325 in (c) and (d) as thin grey lines, each its own bundle of five collapsed temperatures. Number for the text: k_eff at matched SSA agrees across temperature within 0.27% over all 100 pairs (median 0.07%; 0.17% at φ 0.325; 0.24% for φ ≤ 0.375) | built 2026-10-08; NOT used: the user keeps the original without grey lines (2026-10-08) | | |
| 7 | **The state law**: (a) k vs SSA, well-connected packings, one power law; (b) level vs porosity; (c) exponent vs porosity, grey above φ 0.40 | built (sample): `Figure7_state_law` or `Figure7_state_law_alt`; message below | Figs. 7, 8 | 3.2 / 4.2 |
| 8 | **The clock extrapolated**: time to a given sintering age vs temperature and grain radius, with Enceladus conditions and Choukroun's 180 K point | built (sample): `Figure8_timescales`; needs the surface-diffusion check | new | 4.3 |

## Table plan (2 main + 3 supporting; decided 2026-10-06)

All in `tables/tables.tex`, one command each; every table fits a 6.5 in text
block at 11 pt.

| table | command | where | place it right after |
|---|---|---|---|
| 1 | `\TableKeffParameters` (compact, 11 rows) | §3.2, set-up paragraph | the sentence that states the matrix: five porosities × five temperatures × five packings, 2 mm periodic cells |
| 2 | `\TableKeffTemperature` | §3.2, before Fig. 6 | the definition of τ_sub and θ = t/τ_sub, and the statement that 30 d spans 43 to 1256 τ_sub |
| S1 | `\TableMaterialProperties` | SI; cited from §2.1/2.2 | the governing equations, where the constitutive relations are named |
| S2 | `\TableKeffParametersFull` | SI; cited from §2.4 and Table 1's note | the numerical-method paragraph (discretization, time stepping, k_eff sampling) |
| S3 | `\TableGrainPairParameters` | SI; cited from §3.1 | the description of the two-grain set-up, before Fig. 3 |

**Supplement / thesis** (all data in hand):
- S1 domain-size convergence (`rve_convergence/`): the rise is size-independent.
- S2 interface width (`eps_sensitivity/`): late rate robust, early transient not.
- S3 time step and sampling (`dtmax_check/`, `compare/kinks/`).
- S4 run table (`RUN_TABLE.md`).
- S5 anisotropy and k_xy (`rve_anisotropy/`, `rve_convergence/kxy_vs_L`).
- S6 packing construction: the deposition storyboard already made
  (`postprocess/make_deposition_movie.py`).

## The message of each results figure

- **Fig. 6:** temperature changes how fast a deposit matures, not the path it
  takes. (b) fans out, (c) and (d) do not. SSA(θ) is included because it shows
  the microstructure itself collapses; k follows from it.
- **Fig. 5:** orientation, no analysis. Shows necks forming, pores rounding and
  closing, and the loss of a connected ice network at high porosity.
- **Fig. 7:** *conductivity is a power law of SSA with one exponent (−0.80)
  for every well-connected packing; porosity sets the level.* Built on SSA, not
  θ, because the SSA law is tight (0.1–0.2% per packing) and its exponent does
  not move with the fitting window, whereas SSA(θ) and k(θ) are only roughly
  power laws (exponents drift, −0.088 → −0.095 for SSA between θ ≥ 30 and
  θ ≥ 100). The reference state (θ = 30) only marks the end of the
  width-dependent transient.
- **Fig. 8:** vapour sintering matures a deposit only where it is warm and the
  grains are fine. Read as a bound on the vapour route.

**SSA against age is logarithmic after the relaxation** (Jackson's reading of
Fig. 6d, confirmed 2026-10-06): SSA/SSA₀ = A − B·log₁₀θ with B = 0.133–0.136
for φ ≤ 0.375 (A = 0.95–0.98), rms 0.10–0.19% per packing for θ ≥ 30, half the
scatter of a power law. k against θ is not logarithmic (it bends upward),
which follows from combining this with the SSA power law of Fig. 7.

**The master simulation** is φ = 0.325, seed 1702, −20 °C. Every single-run
illustration (Fig. 6, the methods figure, the φ 0.325 column of the gallery)
uses it.

**Manuscript folders (2026-10-06).** `Manuscript/Figures/` holds
Figure1__SinteringMechanisms … Figure8__SinteringTimescales as numbered above,
and FigureS1–S5 for the supplement. `Manuscript/Tables/tables.tex` is the one
tables file: `\input` it in the preamble, then `\TableMaterialProperties`,
`\TableKeffParameters`, `\TableKeffTemperature`, `\TableGrainPairParameters`
where each table goes. Source: `studies/keff_sintering/tables/tables.tex`.

Chain of the argument: Fig. 6 gives SSA = SSA(θ), one curve for all T.
Fig. 7 gives k = k(SSA; φ). Fig. 8 turns θ into years for a given T and R.

## Section by section

**Key points / abstract.** Rewrite. The draft's first key point (an empirical
temperature dependence of the kinetic parameters, "without retuning") is no
longer what the model does: α_c is constant (0.1 in the Molaro runs, a stated
tuning choice; 1e-3 in the packings). The second ("SSA … for both varying
temperatures and varying porosities") must be narrowed: exact across
temperature; across porosity only the exponent, and only for φ ≤ 0.375.

**1 Introduction.** Keep the structure. Add: the near-subsurface framing; the
sintering-age idea up front; Molaro (2019) and Choukroun (2020) as the two
studies this follows.

**2 Methods.** 2.3: fill the missing citations ("?"); state the tensor
conductivity law (prior art Nicoli, Plapp & Henry 2011); add packing
construction, the acceptance tests, and the definition of τ_sub and θ.

**3.1 Validation.** The text says 3D; the runs are axisymmetric. Add the
Demmenie-conditions result when it lands.

**3.2 Effective thermal conductivity. Replace entirely.** Everything in the
draft is from the old setup (1142² mesh, 98% humidity, 28 d, nested packings
at φ 0.24–0.30). New order: set-up and run table → Fig. 5 → Fig. 6 → Fig. 7.
Also correct: the off-diagonal components are NOT orders of magnitude smaller;
k_xy is up to ~12% of k_xx for one packing, zero in the seed mean, and shrinks
with domain size.

**4.1** becomes "Temperature sets the clock": τ_sub, why the collapse holds
(attachment-limited, α_c fixed, same mesh), what would break it.
**4.2** "SSA as the state variable": compare with Calonne et al.; scope to
well-connected packings; the level–porosity relation; what 2D does to
connectivity.
**4.3** "Implications": Fig. 8; near subsurface; Choukroun's activation energy
(24 kJ/mol) against ours (48 kJ/mol) and what that means for which mechanism
leads where; limits (2D, air-filled pores, vapour only, unresolved necks).

**5 Conclusion.** Rewrite around the three claims.

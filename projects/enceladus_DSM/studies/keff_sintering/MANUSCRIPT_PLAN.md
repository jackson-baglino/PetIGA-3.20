# Manuscript plan — figures and sections (k_eff part)

Written 2026-10-06 against the draft `Literature/BaglinoJackson_DSM_manuscript_v2.1.pdf`
(34 pp., 8 figures). Figures are built now on the data in hand and rebuilt when
3f and the redo stage land. Claims and their limits: `TODO.md`, "Manuscript
argument". Run inventory: `RUN_TABLE.md`.

## Figure plan (8 in the draft → 6 main + supplement)

| new | content | status | replaces | section |
|---|---|---|---|---|
| 1 | Sintering mechanisms schematic. Suggest adding a panel for the setting: a near-subsurface plume deposit with sealed pores, which is what the periodic cell models | exists; setting panel to draw | Fig. 1 | 1 |
| 2 | Molaro grain-pair validation; add the Demmenie-conditions exponent when that run lands | exists; Demmenie run queued | Fig. 2 | 3.1 |
| 3 | **k_eff evolution and its temperature collapse**: (a) four snapshots of one packing; (b) k/k₀ vs days, five temperatures; (c) the same vs θ = t/τ_sub; (d) SSA/SSA₀ vs θ | **built**: `figures/fig3_keff_collapse.py` → `Figure3__KeffEvolution/keff_collapse_snapshots.*` | Figs. 3, 4, 5 | 3.2 |
| 4 | **What happens to the aggregates**: five porosities × four instants, qualitative colour bars | built (sample): `figD_gallery` | Fig. 6 | 3.2 |
| 5 | **The state law**: (a) k vs SSA, well-connected packings, one power law; (b) level vs porosity; (c) exponent vs porosity, grey above φ 0.40 | built (sample): `figB_closure`; message below | Figs. 7, 8 | 3.2 / 4.2 |
| 6 | **The clock extrapolated**: time to a given sintering age vs temperature and grain radius, with Enceladus conditions and Choukroun's 180 K point | built (sample): `figC_timescales`; needs the surface-diffusion check | new | 4.3 |

**Supplement / thesis** (all data in hand):
- S1 domain-size convergence (`rve_convergence/`): the rise is size-independent.
- S2 interface width (`eps_sensitivity/`): late rate robust, early transient not.
- S3 time step and sampling (`dtmax_check/`, `compare/kinks/`).
- S4 run table (`RUN_TABLE.md`).
- S5 anisotropy and k_xy (`rve_anisotropy/`, `rve_convergence/kxy_vs_L`).
- S6 packing construction: the deposition storyboard already made
  (`postprocess/make_deposition_movie.py`).

**Suggested additions to the early sections**
- **Methods figure (new, §2.3):** the homogenization on one packing. Panels: the
  periodic cell; a zoom on one contact showing the diffuse interface and the
  mesh; the two corrector fields. The solver can write the correctors
  (`-keff_write_corrector`); one local replay gives the data. This would make
  §2.3, currently equations only, concrete.
- **Setting panel in Fig. 1** (above).

## The message of each results figure

- **Fig. 3:** temperature changes how fast a deposit matures, not the path it
  takes. (b) fans out, (c) and (d) do not. SSA(θ) is included because it shows
  the microstructure itself collapses; k follows from it.
- **Fig. 4:** orientation, no analysis. Shows necks forming, pores rounding and
  closing, and the loss of a connected ice network at high porosity.
- **Fig. 5:** *conductivity is a power law of SSA with one exponent (−0.80)
  for every well-connected packing; porosity sets the level.* Built on SSA, not
  θ, because the SSA law is tight (0.1–0.2% per packing) and its exponent does
  not move with the fitting window, whereas SSA(θ) and k(θ) are only roughly
  power laws (exponents drift, −0.088 → −0.095 for SSA between θ ≥ 30 and
  θ ≥ 100). The reference state (θ = 30) only marks the end of the
  width-dependent transient.
- **Fig. 6:** vapour sintering matures a deposit only where it is warm and the
  grains are fine. Read as a bound on the vapour route.

Chain of the argument: Fig. 3 gives SSA = SSA(θ), one curve for all T.
Fig. 5 gives k = k(SSA; φ). Fig. 6 turns θ into years for a given T and R.

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
at φ 0.24–0.30). New order: set-up and run table → Fig. 3 → Fig. 4 → Fig. 5.
Also correct: the off-diagonal components are NOT orders of magnitude smaller;
k_xy is up to ~12% of k_xx for one packing, zero in the seed mean, and shrinks
with domain size.

**4.1** becomes "Temperature sets the clock": τ_sub, why the collapse holds
(attachment-limited, α_c fixed, same mesh), what would break it.
**4.2** "SSA as the state variable": compare with Calonne et al.; scope to
well-connected packings; the level–porosity relation; what 2D does to
connectivity.
**4.3** "Implications": Fig. 6; near subsurface; Choukroun's activation energy
(24 kJ/mol) against ours (48 kJ/mol) and what that means for which mechanism
leads where; limits (2D, air-filled pores, vapour only, unresolved necks).

**5 Conclusion.** Rewrite around the three claims.

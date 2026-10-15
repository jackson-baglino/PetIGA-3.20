# bouvetSnowEquiTemperatureMetamorphism2022 — Bouvet, Calonne, Flin & Geindreau 2022, JAMES 14, e2022MS002998 (doi 10.1029/2022MS002998)

PDF p.1 = HAL cover; journal page "k of 18" = PDF p.(k+1).

**Summary.** Applies the Bretin et al. (2019) volume-conserving mean-curvature-flow phase-field model ("Snow3D") to 3D μCT images of snow to simulate equi-temperature (isothermal) metamorphism at −2 °C. **Kinetics-limited: no vapor diffusion, no settling.** Physical time is obtained by fitting simulated to measured SSA (Flin et al. 2004 series), giving α = (9.8 ± 0.7)×10^-4; validated on Hagenmuller et al. 2019 series (would need α = 13.1×10^-4). Then computes k_eff, vapor diffusivity and permeability tensors by homogenization (Auriault et al. 2009) on simulated images of four samples over 70–80 days.

**Findings relevant to the manuscript**
- Model is kinetics-limited: p.5 (4 of 18) "Such a model implies that we assume kinetic-limited metamorphism: vapor transport in the pore space is not described."
- α rescales time only: Eq. 8 p.6 (5 of 18) t = t̃ dx²/(α v_kin d0); "The condensation coefficient α is needed to determine the physical time." Abstract p.2: "Relating the numerical time ... to the real physical time involves the condensation coefficient α, a poorly constrained parameter".
- α value and validity: p.7 (6 of 18) "The resulting condensation coefficient is α = (9.8 ± 0.7) × 10^−4 ... the calibrated model can only be used to simulate ETM at this temperature."
- Literature range: p.3 (2 of 18) "Numerous values can be found in the literature, usually ranging from 10^−4 to 10^−1 (see, e.g., Libbrecht & Rickerby, 2013)." Discussion p.16 (15 of 18): low supersaturation (0.1–0.01%) gives low α; L&R 2013 α ≈ 0.06 at σ = 0.1%, −3 °C; Sei & Gonda α ≈ 0.15; Beckmann & Lacmann 0.14 at high σ.
- **Isothermal k_eff evolution was computed**: p.12 (11 of 18) I17 "conductivity increases from 0.078 to 0.081 W m^−1 K^−1"; p.13 Grad3 0.37 → 0.42. p.15 (14 of 18): "For the thermal conductivity, the density changes seem to drive most of the evolution and the impact of the microstructural changes is not visible in our ETM simulation." Density artifacts 1–3% (p.14, 16).
- Anisotropy enhanced under ETM for anisotropic starts (Grad3 1.44 → 1.64) p.11.
- SSA "decreases exponentially ... as classically reported for ETM" p.11 (10 of 18).
- Off-diagonal tensor terms negligible (p.9, 8 of 18).
- γ = 0.109 J m^-2 (Table 1, Libbrecht 2005), same as manuscript Table 1.

**Claims it supports / contradicts / qualifies**
- Draft L276 "Bouvet et al. (2022) used a phase-field model to predict microstructural evolution for dry snow under equitemperature conditions": supported — but should say it is curvature-driven, kinetics-limited, no vapor diffusion.
- Draft L276 "In both these works [Bouvet 2022, Braun 2024], it was required to tune the condensation coefficient ... consistent with experimental μCT data for a given temperature": supported for Bouvet (p.7).
- Draft L439 / Table 4 note "α_c = 0.1 ... at the upper end of the range found in the literature (Braun 2024; Bouvet 2022)": partially. Bouvet gives the range 10^-4–10^-1 (p.3), so 0.1 is its upper end, but Bouvet also quotes α ≈ 0.14–0.15 measured at high supersaturation (p.16). "Upper end of the commonly quoted range" is safer.
- Draft L501 (attributed to Bouvet 2023): the time-rescaling-through-α statement belongs to **this** paper (Eq. 8, p.6) — wrong-citation fix: cite bouvetSnowEquiTemperatureMetamorphism2022. But note: their α rescales time *exactly* because the model is kinetics-limited (no vapor diffusion) — consistent with the outline's correction that in our model α also enters Λ.
- Draft L279 / outline gap paragraph "How metamorphism influences the macroscopic properties in the absence of compaction is still unquantified": **contradicted in part** — Bouvet 2022 computed k_eff (and D, K) on isothermally metamorphosed (simulated, no settling) snow and found small k_eff change (0.078→0.081) that they attribute mostly to density artifacts. The outline's proposed narrowing ("what is not quantified is the isothermal case") must therefore acknowledge Bouvet 2022: they did isothermal, constant-density (nominally), but kinetics-limited, −2 °C only, and found microstructure effect "not visible". Suggested gap: isothermal sintering *from a granular packing with point contacts*, across temperature and porosity, with vapor diffusion — and a conductivity–SSA relation.
- Outline L1050 "SSA against log time in isothermal snow" — Bouvet says SSA decreases "exponentially" (p.11) in their simulations; form differs; do not cite Bouvet for log-time.
- Outline L1215–1216 "\Note Bouvet et al. fit α_c = (9.8 ± 0.7)×10^-4 to tomography at −2 °C": supported (p.7).

**Candidate citations / claims**
- Add Bouvet 2022 to the Intro "what is known about conductivity" paragraph: the only prior isothermal-metamorphism k_eff computation in the library; it found microstructural changes under ETM did not visibly change k_eff at fixed density (p.15). This is the strongest contrast with our result (k_eff +17% at fixed porosity) and must be discussed (difference: rounded snow vs. sintering point contacts; necks).
- Use Bouvet 2022 p.6 Eq. 8 for "α rescales time in a kinetics-limited model".

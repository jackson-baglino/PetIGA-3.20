# braunRigorousApproachSpecific2024 — Braun, Fourteau & Löwe (2024), The Cryosphere 18, 1653–1668, doi:10.5194/tc-18-1653-2024 (printed)

**Summary.** Steady heat and vapour diffusion is solved by FE (Elmer) on µCT images from two Snowbreeder temperature-gradient series (Pinzer et al. 2012): series 1 at −8.1 °C mean and 47 K/m for 384 h; series 2 at −7.6 °C and 55 K/m for 665 h. Both start as rounded grains (φ_i 0.31 / 0.28, SSA_V 20 / 24 mm⁻¹). A Hertz–Knudsen interface law is used with one constant α, and latent heat is neglected. The interface velocities feed the exact SSA evolution equation of Krol & Löwe 2018, and α is fitted by matching the measured SSA.

**Findings relevant to the manuscript**
- p.4: "we rely on the simplifying assumption of a single and constant α value. It should thus rather be understood as an effective condensation coefficient."
- p.8: "αbest = 10−1 for series 1 and αbest = 10−2.25 for series 2". p.12 Fig. 7: series 1 needs α = 10⁻¹.⁵ after t = 160 h, attributed (p.11) to "the predominant emergence of facets on down-facing surfaces".
- p.12: the values "lie in the commonly found range of 10−3 < α < 10−1 (Libbrecht and Rickerby, 2013), which is also used by Kaempfer and Plapp (2009)". They are compared with Bouvet 2022's α ≈ 9.8 × 10⁻⁴ from isothermal metamorphism. "It is well known that α is difficult to measure experimentally."
- p.13: "On facets, one expects that α is significantly reduced by orders of magnitude"; latent heat is neglected in the interface condition.
- p.11: "Our diffusion model is essentially physically equivalent to Kaempfer and Plapp (2009) in the steady-state limit".

**Claims supported / qualified**
- Outline L434–438 ("commonly found range" citing Libbrecht & Rickerby; K&P use it; Braun's values inside it): **supported** (p.12). Caveat: Libbrecht & Rickerby 2013 do not state this range themselves (they give α = A·exp(−σ0/σ), A ≈ 1; see my notes). Cite Braun for "commonly found", and K&P 2009 Table I for the range used in phase-field DSM.
- Outline L469–490 (Braun: steady heat and vapour on µCT under a gradient, −8.1 and −7.6 °C, α 10⁻¹ → 10⁻¹.⁵ and 10⁻²·²⁵, facets, "effective"): **supported**. The exact phrase is "effective condensation coefficient" (p.4), not "effective kinetic coefficient". Fix the quote at L486.
- Outline L1215–1217 (Braun find 10⁻³–10⁻¹ changing as surfaces facet, an effective surface-averaged value): **supported**, with "fitted values 10⁻²·²⁵ to 10⁻¹, within the commonly cited 10⁻³–10⁻¹" being more precise.
- Draft L276 (Braun "modeled the steady state thermal and vapor diffusion equations"): **supported** but incomplete, as the outline notes.
- Draft L501 ("α only changes the rate"): Braun treats α as changing the SSA trajectory (different α give different SSA decay curves, Figs. 4 and 7) under a temperature gradient. It does not support the claim that α only rescales time in general. Combined with Legagneux 2005 [47], this claim needs qualifying to the attachment-limited isothermal case.

**Candidate additions.** None new. Braun p.4 "effective condensation coefficient" is the quote to use for treating α_c as effective.

# braunRigorousApproachSpecific2024 — Braun, Fourteau & Löwe 2024, The Cryosphere 18:1653–1668 (doi 10.5194/tc-18-1653-2024)

**Summary.** Steady (Laplace) heat and vapor diffusion solved by FE (Elmer) on μCT images of snow under a temperature gradient, with a Hertz–Knudsen + Gibbs–Thomson Robin condition at the interface; interface velocities fed into the exact surface-area equation ṡ = 2s⟨v_n H⟩ (Krol & Löwe 2018) to predict SSA. Data: Pinzer et al. (2012) TGM series 1 (T̄ = −8.1 °C, 47 K m^-1, 384 h) and series 2 (−7.6 °C, 55 K m^-1, 665 h). α fitted by RMSE. Latent heat neglected in the interface condition.

**Findings relevant to the manuscript**
- Model: steady diffusion, p.3 (1655) "we employ the common assumption of a small particle Péclet number ... and consider stationary heat and mass diffusion equations"; latent heat neglected: p.3 "the latent heat during the sublimation and deposition is neglected for reduced model complexity."
- α best: p.8 (1660) "αbest = 10^−1 for series 1 and αbest = 10^−2.25 for series 2"; transition p.11 (1663) "α = 10^−1.5 for the final stage" after t ≈ 160 h in series 1, linked to facets on down-facing surfaces.
- Range: p.12 (1664) "the obtained values of 10^−1 and 10^−1.5 for series 1 and of 10^−2.25 for series 2 lie in the commonly found range of 10^−3 < α < 10^−1 (Libbrecht and Rickerby, 2013), which is also used by Kaempfer and Plapp (2009)." Jafari et al. 2020 translates to α ≈ 5×10^-7.
- Effective parameter: p.12 "Using a single, constant value of α ... must be therefore understood as an effective kinetic coefficient."; p.4 (1656) "It should thus rather be understood as an effective condensation coefficient."
- Isothermal case has one timescale: p.12 (1664) "the temperature gradient case considered here is governed by two different timescales instead of only one in the isothermal case (Bouvet et al., 2022)". Useful for Section 5.2 (temperature as a rescaling of time).
- Isothermal rare in nature: p.1 (1653) "gradient-free (i.e., isothermal) conditions exist at most in deep polar firn."
- REV: p.6 (1658) cites Calonne 2011 "representativeness for sample side lengths of between 2.5 and 5 mm" (for conductivity).
- α ∈ [0,1] as sticking probability; α ≈ 1 → diffusion dominated (p.3).

**Claims it supports / contradicts / qualifies**
- Outline L434–438 (α range 10^-3–10^-1 "commonly found range", citing L&R 2013, Kaempfer & Plapp use it, Braun's own values inside): **supported** (p.12).
- Outline L469–491 (Braun description: steady heat and vapor diffusion on μCT under gradient, two series −8.1 and −7.6 °C, exact SSA equation; α 10^-1 then 10^-1.5, and 10^-2.25; "effective kinetic coefficient"): **supported** in every detail.
- Draft L276 "Braun et al. (2024) modeled the steady state thermal and vapor diffusion equations to predict the rates of ice growth": supported but incomplete (as outline says).
- Draft L439 α_c = 0.1 "upper end of the range found in the literature (Braun 2024; Bouvet 2022)": supported by Braun (upper bound 10^-1 of the "commonly found range"; Braun's own best value for series 1 is 10^-1).
- Outline L1215–1218 "Braun et al. find 10^-3–10^-1, changing as surfaces facet, and call it an effective, surface-averaged value": supported (p.12–13). Strictly their fits are 10^-2.25–10^-1; 10^-3–10^-1 is the literature range they cite.
- Outline Section 5.2 temperature rescaling: Braun p.12 supports that the isothermal problem has a single timescale (Bouvet 2022's dimensional analysis), unlike TGM.

**Candidate citations / claims**
- Cite Braun p.12 for "in the isothermal case the dynamics are governed by a single timescale" in Section 5.2 (with Bouvet 2022).
- Note for Methods: Braun neglect latent heat; our model includes it — one clause in the model-lineage paragraph if wanted.

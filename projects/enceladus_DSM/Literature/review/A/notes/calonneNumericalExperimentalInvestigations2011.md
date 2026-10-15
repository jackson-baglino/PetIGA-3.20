# calonneNumericalExperimentalInvestigations2011 — Calonne, Flin, Morin, Lesaffre, Rolland du Roscoat & Geindreau 2011, GRL 38, L23501 (doi 10.1029/2011GL049234)

**Summary.** Full k_eff tensor computed by periodic homogenization (GeoDict) on 30 μCT images spanning all seasonal snow types (PP, DF, RG, FC, DH, MF; 103–544 kg m^-3), conduction in ice and air only, k at 271 K (k_air 0.024, k_ice 2.107) and some at 203 K; compared to needle-probe and flux/gradient measurements (cold-room TG experiment 43 K m^-1, 300 kg m^-3, 3 weeks).

**Findings relevant to the manuscript**
- Cell problem exactly as in the draft (eqs. 1–11, p.2), with "Without loss of generality, the medium is assumed to be Ω-periodic" (p.2, para 9) — this is where the draft's (mis-typed) phrase comes from. Separation of scales cites Auriault 1991 and Auriault et al. 2009 (p.2, para 8); "Using the homogenization method for periodic structures, it can be shown that the corresponding macroscopic description is [Auriault et al., 2009]".
- Regression p.3: k_eff = 2.5×10^-6 ρ² − 1.23×10^-4 ρ + 0.024 (Eq. 12), "The correlation coefficient for the 30 samples is 0.985" (p.4).
- Anisotropy p.4 (para 15): "The ratio kz/kxy can be as low as 0.7" (RG deep in snowpack) and "as high as 1.5" (FC, DH); "This anisotropy induces a slight scatter around the overall correlation between keff and r".
- Off-diagonal: p.3 "Non-diagonal terms of keff are more than 100 times lower than diagonal terms".
- **REV size** p.3 (para 11): "REVs for thermal conductivity estimation range from about 2.5 mm (PP) to 5.5 mm (MF)."
- Conduction in air vital: ignoring it underestimates kz by ×10 for PP, up to ×2 for denser snow (p.5).
- k_eff increases with decreasing T, following k_ice (p.4).
- Conductive vs. apparent: p.5 "we recommend purely conductive effects (i.e., conduction through ice and interstitial air) to be considered separately from non-conductive processes. At present, the only method able to unambiguously characterize the purely conductive effective thermal conductivity is numerical modeling based on microstructure 3D images."
- Intro motivation p.1: "Understanding heat and mass transfer through the snowpack is critical to assess the climatic and hydrological role of snow in the cryosphere." No mention of permafrost or glacier stability.
- k_eff regression "follows closely the regression curve proposed by Yen [1981] ... corroborates the use of this relationship in several snowpack models (e.g., Crocus)" (p.5).

**Claims it supports / contradicts / qualifies**
- Draft L256 "the effective thermal conductivity of snowpacks and firn determines how heat moves through the icy porous medium, affecting the stability of glacial ice, snowpacks, and permafrost (Calonne 2011)": **partially** — heat transfer through snow matters for "the climatic and hydrological role of snow in the cryosphere" (p.1); no firn, glacier, or permafrost content (samples are seasonal snow). Outline's removal of the permafrost clause is right.
- Draft L353 periodicity "as described in Calonne 2011 and Auriault 1991": supported (Calonne uses "without loss of generality"). Outline's \tC ("periodicity is an assumption, not a convenience") is a matter of wording; Calonne and Auriault both frame it as without loss of generality for the *form* of the law (see Auriault note). Keep "assume", fix "less" → "loss".
- Draft L362 macroscopic description and K_eff definition (Eqs. 9–11): supported (p.2, Eqs. 5–11) — draft equations are Calonne's. Draft Eq. (11d) has t where Calonne Eq. (10) has ∇t → typo in draft confirmed.
- Outline L519–528 (\Note Calonne 2011: full tensor, 30 samples, Eq. 12, r = 0.985, kz/kxy 0.7 to 1.5, Fig. 1, Sect. 3.2): **supported**.
- **Outline L852–856 "Domain size L/R̄ = 40, at the 1–2 mm representative volume of Calonne et al. (2011)": contradicted.** Calonne's REVs for k are "about 2.5 mm (PP) to 5.5 mm (MF)" (p.3). Our 2 mm domain is below their smallest REV. Defensible statement: the REV in grain units — our domain-size study (Fig. S1) — not Calonne's mm values; or "Calonne et al. (2011) found REVs of 2.5–5.5 mm for snow, whose grains are larger than ours" (snow grain sizes not tabulated here, so ratio claim needs care).
- Outline L1276–1278 (Discussion 5.3 \Note): supported.
- Outline L565–568 / Limitations: Calonne supports the conduction-only definition of k_eff (p.5).
- Draft L471 / outline: our off-diagonal ~12% on one packing vs Calonne's >100× lower for snow — a difference worth one clause (2D packings are less isotropic).

**Candidate citations / claims**
- Cite Calonne 2011 p.5 for "the purely conductive k_eff is what pore-scale computation yields; vapor/latent effects must be treated separately" (Section 2.3 / Limitations).
- Cite Calonne 2011 (with Bouvet 2022) in the Intro "what is known about conductivity" paragraph: density explains most (r = 0.985), anisotropy the scatter.

**Verified on rendered page (PDF p.3):** 'REVs for thermal conductivity estimation range from about 2.5 mm (PP) to 5.5 mm (MF)' — units are mm on the page.

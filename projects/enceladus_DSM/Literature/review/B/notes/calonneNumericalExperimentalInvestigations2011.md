# calonneNumericalExperimentalInvestigations2011 — Calonne, Flin, Morin, Lesaffre, Rolland du Roscoat & Geindreau (2011), GRL 38, L23501, doi:10.1029/2011GL049234 (printed)

**Summary.** Full k_eff tensor computed by periodic homogenization (Geodict) on 30 µCT images spanning all seasonal snow types (PP, DF, RG, FC, DH, MF; 103–544 kg/m³). Conduction only, through ice and air, steady state. k at 271 K is k_air = 0.024 and k_ice = 2.107 W/m/K (the paper misprints the second as "kair"); a subset was also run at 203 K and with k_air = 0. The results are compared with needle-probe and flux/gradient measurements from a 43 K/m slab experiment.

**Findings relevant to the manuscript (all rendered, pp.2–5)**
- p.2 Eqs. 1–11: the microscopic problem and the periodic homogenization, with k_eff = (1/|Ω|)[∫_air k_air(∇t_air + I) + ∫_ice k_ice(∇t_ice + I)] and Ω-periodic vectors t_ice and t_air with zero mean (Eq. 11), following Auriault et al. 2009.
- p.3 [11]: "REVs for thermal conductivity estimation range from about 2.5 mm (PP) to 5.5 mm (MF). Calculations performed on a lower volume than the REV were disregarded."
- p.4 Eq. 12: k_eff = 2.5×10⁻⁶ρ² − 1.23×10⁻⁴ρ + 0.024, with "correlation coefficient for the 30 samples is 0.985" and residual s.d. 0.025 W/m/K. It agrees with Yen 1981, k = 2.22362(ρ_snow/ρ_water)^1.885.
- p.4 [15]: "The ratio kz/kxy can be as low as 0.7" (deep rounded grains) and "as high as 1.5" (FC, DH); "This anisotropy induces a slight scatter around the overall correlation".
- p.4 [17]: k_z rises as T falls, mostly through k_ice. p.5: neglecting air conduction underestimates k_z ×10 for PP and up to ×2 for denser snow (Table 1).
- p.5 [21]: "we recommend purely conductive effects ... to be considered separately from nonconductive processes. At present, the only method able to unambiguously characterize the purely conductive effective thermal conductivity is numerical modeling based on microstructure 3D images." Also: "the effective thermal conductivity of snow is strongly correlated with snow density".
- p.1 [2]: "Understanding heat and mass transfer through the snowpack is critical to assess the climatic and hydrological role of snow in the cryosphere." Permafrost is not mentioned.

**Claims supported / qualified**
- Draft L256 (bulk k of snowpacks and firn "affecting the stability of glacial ice, snowpacks, and permafrost"): **partially / unsupported** beyond snowpack. The paper is about seasonal snow only (no firn or permafrost); p.1 supports the "climatic and hydrological role of snow". This agrees with outline L292.
- Draft L362 (macroscopic description as in Auriault 1991 and Calonne 2011): **supported** (p.2 Eqs. 5–11). The vector naming t follows Calonne.
- Outline L515–520 (30 samples, Eq. 12, r = 0.985, k_z/k_xy 0.7–1.5, Fig. 1, Sect. 3.2): **supported**.
- Outline L853 ("Domain size L/R̄ = 40, at the 1–2 mm representative volume of Calonne et al. (2011)"): **contradicted**. Calonne reports k_eff REVs of ~2.5 mm (PP) to 5.5 mm (MF) (p.3 [11]). Our 2-mm domain is below their smallest REV, and the comparison is in physical size for natural 3-D snow, not in L/R̄. Fix: state the REV in grain radii from our own domain-size test (L/R̄ = 40), and cite Calonne only as "snow REVs for k_eff are 2.5–5.5 mm, i.e. tens of grain sizes". The project memory's "≥2 mm" also does not match the paper.
- Outline L297 (k_eff "sets how fast heat crosses them and so the temperature profile"): Calonne 2011 supports the snow case in general terms only (p.1). There is no temperature-profile statement.
- Outline L1355 (Calonne 2019 on T-dependent conductivities): Calonne 2011 p.4 [17] (k_eff rises as T falls, via k_ice) also supports this and is closer to our use of constant k_ice and k_air.
- Steady-state conduction only, matching our approach (p.2 "Neglecting convection, phase change and under steady state conditions").

**Candidate additions.** p.5 [21] supports our framing: separating conduction from vapour and latent effects, with conduction computed from microstructure. Also report the 271 K values (k_ice 2.107, k_air 0.024) in comparison with our Table 1 (2.29 / 0.02, from K&P).

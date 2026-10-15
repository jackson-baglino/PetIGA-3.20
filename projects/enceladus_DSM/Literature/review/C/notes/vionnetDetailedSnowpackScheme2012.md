# vionnetDetailedSnowpackScheme2012 — Vionnet, Brun, Morin, Boone, Faroux, Le Moigne, Martin & Willemet, Geosci. Model Dev. 5, 773–791 (2012)
DOI (PDF p.1): 10.5194/gmd-5-773-2012. Read in full (19 pp).

## Summary
Model description of the detailed multilayer snowpack scheme Crocus in SURFEX v7.2: layering rules, empirical metamorphism laws (dendricity, sphericity, grain size; dry and wet), compaction, wind drift, albedo, energy balance, heat conduction, melt, water flow; evaluation at Col de Porte (18 seasons), distributed Alps run, coupled Dome C run.

## Findings relevant to the manuscript
- Snow conductivity in Crocus is a density-only law: "The snow effective thermal conductivity, k ... follows the expression of Yen (1981): k = k_ice (ρ/ρ_w)^1.88" (Eq. 18, p.9–10/781–782).
- Conductivity sets the temperature profile: wind-drift densification parameterization "turned out to be necessary in polar environments to reproduce correctly the snow thermal conductivity and, therefore, the snow temperature profile" (p.8/780).
- Metamorphism is empirical: dry metamorphism laws in temperature and temperature gradient (Table 1, p.5/777), Arrhenius-like factor e^{−6000/T}; "mostly based on empirical fits to experimental data" (p.6/778). Sintering enters only indirectly ("decrease with time due to sintering (increase of s)", p.7/779).
- Heat diffusion solved within the snowpack (p.9/781).

## Claims supported / contradicted / qualified
- Outline L297–299 (Intro para 2 addition: k_eff "is an input of the detailed snowpack models used to follow their evolution (Lehning 2002; Vionnet 2012)"): supported — Crocus computes heat conduction with k from Eq. 18 (p.9–10).
- Outline Section 5.3 "Density is not enough" (L1273–1281): Crocus is a direct example of the density-only parameterization that cannot change k while density stays fixed — useful contrast (k depends on ρ only, Eq. 18).

## Candidate additions
- Section 5.3: "Detailed snowpack models parameterize conductivity by density alone (e.g., Crocus, k = k_ice(ρ/ρ_w)^1.88; Vionnet et al., 2012, Eq. 18), so they cannot represent the rise in conductivity of a sintering deposit at fixed density." Closely matching.

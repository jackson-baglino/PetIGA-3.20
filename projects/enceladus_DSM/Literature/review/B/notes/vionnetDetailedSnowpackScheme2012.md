# vionnetDetailedSnowpackScheme2012 — Vionnet et al. (2012), GMD 5, 773–791

**Summary.** Model-description paper for the detailed 1-D multilayer snowpack scheme Crocus inside SURFEX v7.2 (snowfall, layering, empirical metamorphism laws for dendricity/sphericity/grain size, compaction, drift, radiation, heat conduction, melt). Evaluated at Col de Porte (10 seasons), a distributed alpine run, and coupled run at Dome C. No microstructure physics; metamorphism is empirical (Brun et al. 1992; Marbouty 1980).

**Findings relevant to the manuscript**
- Snow thermal conductivity is a density-only parameterization: PDF p.9 (journal p.781), Sect. 3.8: "The snow effective thermal conductivity, k, ... follows the expression of Yen (1981): k = k_ice (ρ/ρ_w)^1.88" (Eq. 18).
- Heat diffusion in the snowpack solved with this k (Sect. 3.8, p.9): "The heat diffusion within the snow cover is computed by SNOWCROSOLVT using the implicit backward-difference integration scheme of ISBA-ES".
- p.8 (journal 780): drift parameterization "turned out to be necessary in polar environments to reproduce correctly the snow thermal conductivity and, therefore, the snow temperature profile" — conductivity controls the modeled temperature profile.
- Dry metamorphism laws (Table 1, p.5) are empirical functions of T and gradient G (G≤5, 5–15, >15 K/m); grain size static under weak gradient for non-dendritic snow.

**Claims supported / qualified**
- Outline L294–300 (Intro ¶2 add): "conductivity ... is an input of the detailed snowpack models used to follow their evolution (Lehning et al., 2002; Vionnet et al., 2012)" — **supported** for Vionnet (Eq. 18, p.9). Note it enters as a function of density only, which also supports outline §5.3 ¶6 "density-based formula cannot describe a sintering deposit" as the contrast (Crocus would predict no change at fixed density).
- "sets the temperature profile" (Calonne et al. claim) — Vionnet p.8 is consistent but is a model statement.

**Candidate additions.** In §5.3 "Density is not enough": operational detailed snow models (Crocus) parameterize k from density alone (Vionnet et al. 2012, Eq. 18, after Yen 1981), so a sintering deposit at constant porosity would be predicted to keep constant k — the paper's 17 % rise at fixed porosity is the counterexample within the simulated range. Close match; recommended.

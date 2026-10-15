# vionnetDetailedSnowpackScheme2012 — Vionnet et al. 2012, Geosci. Model Dev. 5:773–791 (doi 10.5194/gmd-5-773-2012)
Read (model description; heat/conductivity/metamorphism sections closely). Eq. 18 rendered p.781.
- Sect 3.8, Eq. 18: "The snow effective thermal conductivity, k ... follows the expression of Yen (1981): k = k_ice (ρ/ρ_w)^1.88" — density only.
- Metamorphism in Crocus is via grain variables (dendricity, sphericity, size / SSA in later versions) but conductivity does not use them.
Claims:
- Outline L297–299 (conductivity "is an input of the detailed snowpack models (Lehning 2002; Vionnet 2012)"): SUPPORTED for Vionnet (Eq. 18).
- Strong support for Sec 5 "Density is not enough": the operational Crocus scheme sets k from density alone (Yen 1981), so sintering at fixed density changes nothing in it. Candidate sentence: "Crocus computes k from density alone (Vionnet et al., 2012, Eq. 18)".

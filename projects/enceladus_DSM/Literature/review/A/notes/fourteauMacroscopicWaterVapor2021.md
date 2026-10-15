# fourteauMacroscopicWaterVapor2021 — Fourteau, Domine & Hagenmuller 2021, The Cryosphere 15:389–406 (doi 10.5194/tc-15-389-2021)

**Summary.** Theory and steady FE simulations (Elmer) of macroscopic water-vapor diffusion in snow under a temperature gradient (idealized disconnected spheres; six μCT samples DF–MF, 125–380 kg m^-3; mean 258 K; 5–200 K m^-1), with Hertz–Knudsen boundary conditions for α = 0 (inert) to 1, infinitely fast kinetics, and a nonlinear α = exp(−σ0/σ). Shows effective vapor diffusivity in snow is below free-air diffusivity for any α; refutes "hand-to-hand" enhancement. Latent heat not included in the temperature field.

**Findings relevant to the manuscript**
- Slow/fast transition p.10 (398): "The transition between the fast- and slow-surface-kinetics regimes occurs for values of α around 10^−3." Table 1: D_norm at α = 10^-3 is 0.65–0.85; α = 0.1 within ~2% of infinitely fast.
- Curvature effects neglected under TG (p.5): "In the presence of a large-enough thermal gradient, the dependence of the saturation concentration on the local curvature of the ice surface can be neglected" (Colbeck 1983).
- Nonlinear α makes D_eff gradient-dependent, "not intrinsic" (p.11, 399).
- Constants (Table D1, Eq. 11): Clausius–Clapeyron with ΔH_s = 51×10^3 J mol^-1, P0 = 611 Pa at 273.15 K — i.e., the temperature dependence of saturation vapor density is an Arrhenius-like factor with ≈51 kJ mol^-1.
- α depends on local saturation and crystallography (Saito; Libbrecht & Rickerby 2013), p.8.

**Claims it supports / contradicts / qualifies**
- Outline L1357: "your .bib has two [Fourteau 2021]; the point about apparent conductivity is the first" — **correct**: this paper (fourteauMacroscopicWaterVapor2021) is about vapor diffusivity, not conductivity; it is not the source for the apparent-conductivity point.
- Outline L1258–1259 (clock activation energy ~50 kJ mol^-1, "that of the saturation vapor density"): consistent with Fourteau's ΔH_s = 51 kJ mol^-1 in the Clausius–Clapeyron form (and Choukroun's "vapor (~51 kJ/mol)"). τ_sub ∝ 1/ρ_vs also has a weaker √T and d0 ∝ 1/T dependence, so "about 50" is right.
- Outline Section 5.2 (Λ ≈ 2.8, diffusion and attachment comparable at α_c = 10^-3): Fourteau's transition near α ≈ 10^-3 (for mm-scale snow at 258 K) is consistent in kind.

**Candidate citations / claims**
- Optional: cite for "α ≈ 10^-3 marks the transition between attachment- and diffusion-limited vapor transport in snow (Fourteau et al., 2021)" in Section 2.2/5.2 — closely matching statement; note their grains are larger (mm), so the transition α scales with grain size (our Λ).

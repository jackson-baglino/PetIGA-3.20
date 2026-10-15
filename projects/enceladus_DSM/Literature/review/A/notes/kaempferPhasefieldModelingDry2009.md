# kaempferPhasefieldModelingDry2009 — Kaempfer & Plapp 2009, Phys. Rev. E 79, 031502 (doi 10.1103/PhysRevE.79.031502)

Table I and Eqs. 37–43 verified on rendered PDF pages 4 and 8 (text layer scrambles them).

**Summary.** Phase-field model of dry snow metamorphism (phase-field eq. with double-well and driving force u − u_eq(T), heat equation with latent heat, one-sided vapor diffusion), built by analogy with alloy-solidification models; parameters related to the sharp-interface problem (Gibbs–Thomson with kinetic coefficient β) via the **sharp-interface limit**; time-scaling ξ of diffusion and source terms (speed-up up to 10^4); validated in 1D against a pseudo-analytical solution and in 2D against Stehle's bubble-migration experiments; qualitative 2D/3D runs on μCT snow under 190 K m^-1. Aimed at TG metamorphism ("temperature-gradient driven", conclusion p.15).

**Findings relevant to the manuscript**
- Saturation vapor pressure Eq. 2 (p.3): P_vs(T) = exp(Σ_{j=0}^{4} K_j T^{j−1} + K_5 ln T), coefficients from Wexler (1977) (Table I: −0.5865×10^4, 0.2224×10^2, 0.1375×10^-1, −0.3403×10^-4, 0.2697×10^-7, 0.6918). Eq. 1: x_s = 0.62 P_vs/(P_a − P_vs); ρ_vs = ρ_a x_s (Eq. 3).
- Table I values (p.4, rendered): k_i, k_a = 2.29, 0.02 W m^-1 K^-1 [Engineering Toolbox]; ρ_i, ρ_a = 918.9, 1.341 kg m^-3; D_v^STP = 2.178×10^-5 m² s^-1, exponent 1.81; L_sg = 2.60×10^9 J m^-3; γ = 1.09×10^-1 J m^-2; d0 = 1.3×10^-9 m; "10^−3 < α < 10^−1" [Libbrecht 2006 report]; 3×10^4 < β < 3×10^6 s m^-1; evaluated at T0 = −10 °C (263 K).
- Gibbs–Thomson with kinetics Eq. 11: ρ_v = ρ_vs(T)(1 + d0 K + β v_n); β = (ρ_i/ρ_vs) √(2πm/kT)/α (Eq. 16).
- Sharp- vs thin-interface p.3: "the intrinsic properties of the ice-vapor interface (in particular its slow kinetics) and the experimental conditions ... justify the use of the simpler sharp-interface limit". Thin-interface correction only sketched (Eqs. 39–41, p.8): β' = a1(τ/(Wλ) − a2 W/α_T − a3 W/D_v), a1 = (5/8)√2, "a3 depends on the local interface temperature"; conditions W ≪ (κ/C)(ρ_vs/ρ_i)β0 etc. (Eq. 43) → W < 5×10^-7 to 5×10^-5 m (Eq. 44); also W < 10^-6–10^-4 m (Eq. 46).
- Validation: 1D quantitative; 2D bubble migration "within the order of magnitude of the experimental values" (p.13); 3D/2D snow "only qualitative" (p.14).
- Isotropic constant β and d0 assumed; faceting would need anisotropic β and a thin-interface limit (p.15).
- Neglected: "grain boundary diffusion or the anisotropies created by the crystal structure" (p.3); "The influence of grain boundaries, volume diffusion in ice, viscous flow, or plastic deformations is neglected" (p.15).
- Heat flow concentrates along favoured paths in 2D snow (p.14).

**Claims it supports / contradicts / qualifies**
- Draft L276 "Kaempfer & Plapp (2009) were the first to attempt to use a phase-field model to capture the dynamics of DSM, by minimizing the thermodynamic free energy, coupled with temperature and mass balance equations": supported (Bouvet 2022 p.3 also calls them "pioneers"). Their target was TG metamorphism.
- Draft L339 "Our model reproduces, under the conditions outlined by K&P, the kinetics as defined by the Gibbs–Thomson condition": supported for the sharp-interface limit (K&P Eqs. 11, 37–38, 42–46).
- **Outline L712–718 (thin-interface term Δβ = a1 a2 ε (D_T^-1 + D_v^-1) ρ_ice/ρ_vs; cite K&P 2009, Moure & Fu 2024)**: partially. K&P give only the dimensional form with *separate* constants a2 (thermal) and a3 (vapor, T-dependent) and explicitly do not use it ("a complete analysis ... beyond the scope", p.8); they work in the sharp-interface limit. The single-a2 form and a1 = 5, a2 = 0.1581 (Table 1 note) are not from K&P (K&P's a1 = 5√2/8 ≈ 0.884 for their tanh-type model); check Moure & Fu for the source. Cite K&P for the form only.
- Draft L494 "quantitatively validate their model in 1D and 2D by comparing them to a pseudoanalytical solution and bubble migration experiments ... also qualitatively validate ... μCT images": 1D quantitative ✓; 2D bubble migration agreement only "within the order of magnitude" — "quantitatively" overstates for 2D. μCT runs "qualitative" ✓.
- Manuscript Table 1 (k_ice 2.29, k_air 0.02, ρ 919/1.341, D_v0 2.178×10^-5, exponent 1.81, γ 0.109, m 2.99×10^-26) matches K&P Table I exactly — cite K&P as the source; K&P evaluate them at −10 °C.
- Table 1 note "P_vs (ASHRAE polynomial)" — K&P (and the draft Eq. 3 form) use Wexler 1977; supports changing the note to Wexler (1977).
- Draft Eq. (3) "k5 + ln(T+T0)": K&P Eq. 2 has K5 ln T → typo confirmed again.
- Outline Section 5.1 / Methods: time scaling ξ (K&P Sec. II D) scales diffusion *and* source terms together, valid while quasi-steady — supports the manuscript's ξ_v = 10^-3 (Tables 3, 4) with a citation to K&P; outline L663–664 asks to explain ξ_v: cite K&P Sec. II D.
- Braun 2024's statement that K&P use α in 10^-3–10^-1: confirmed (Table I).

**Candidate citations / claims**
- Methods: "ξ_v reduces vapor diffusivity and source together, which leaves interface velocities unchanged while the vapor field is quasi-steady (Kaempfer & Plapp, 2009, Sect. II D)."
- Table 1: cite K&P 2009 for the material constants.

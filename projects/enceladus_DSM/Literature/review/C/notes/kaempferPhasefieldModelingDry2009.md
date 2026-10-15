# kaempferPhasefieldModelingDry2009 — Kaempfer & Plapp, Phys. Rev. E 79, 031502 (2009)
DOI (PDF p.1): 10.1103/PhysRevE.79.031502. Read in full (17 pp). Table I rendered and checked (p.4).

## Summary
Phase-field model of dry snow metamorphism: phase field φ ∈ [−1, 1] (double well), asymmetric heat diffusion, one-sided vapor diffusion, Gibbs–Thomson + kinetic boundary condition; parameters chosen in the sharp-interface limit (justified by slow kinetics, W bounds Eq. 44/46), plus a time-scaling procedure (diffusion and source terms scaled together by κ). Validated in 1D against a pseudo-analytical solution (T = 261/260 K, 200 K m⁻¹) and in 2D against Stehle's (1965) cylindrical-hole migration experiments; qualitative 2D/3D runs on µCT snow (268 kg m⁻³, 190 K m⁻¹). Finite differences, PETSc.

## Findings relevant to the manuscript
- Table I (p.4, rendered; T0 = −10 °C): d0 = 1.3×10⁻⁹ m; α "10⁻³ < α < 10⁻¹" (ref [31], Libbrecht 2006, a departmental source); β range 3×10⁴–3×10⁶ s m⁻¹, 5.5×10⁵ in simulations; γ = 1.09×10⁻¹ J m⁻² (ref [26], snowcrystals.com); κi, κa = 2.29, 0.02 W m⁻¹ K⁻¹; Ci, Ca = 1.8×10⁶, 1.4×10³ J m⁻³ K⁻¹; ρi, ρa = 918.9, 1.341 kg m⁻³; Lsg = 2.60×10⁹ J m⁻³; Dv(STP) = 2.178×10⁻⁵ m² s⁻¹, exponent 1.81 (Massman 1998); Kj = {−0.5865×10⁴, 0.2224×10², 0.1375×10⁻¹, −0.3403×10⁻⁴, 0.2697×10⁻⁷, 0.6918} (Wexler 1977). Typo in their table: "2.99×20⁻²⁶ kg".
- Pvs(T) = exp(Σ_{j=0}^{4} K_j T^{j−1} + K_5 ln T) (Eq. 2, p.3) — a product K5·lnT, confirming the draft Eq. 3 "+ k5 + ln(T+T0)" is a typo.
- ρvs = ρa·0.62·Pvs/(Pa − Pvs) (Eqs. 1, 3, p.3, p.5): matches draft Eq. 2.
- Gibbs–Thomson/kinetic condition ρv = ρvs(T)(1 + d0 K + β vn) (Eq. 11, p.5); β = (1/α)(ρi/ρvs)√(2πm/kT) (Eq. 16, p.6). Draft Eqs. (GibbsThomson) and β(α) match.
- "the influence of interface curvature on the water vapor pressure ... is very small for snow grains ... of the order of 0.1 mm. Nevertheless, the difference of local curvatures ... provides a driving force ... that cannot be neglected" (p.5–6).
- Time scaling: "both fluxes are scaled by the same factor ... the relative magnitude of the different transport mechanisms is preserved" (p.10); valid for κ ≳ 10⁻⁴ in 1D tests (p.12); accelerates by up to four orders of magnitude (p.15).
- Validation: 1D pseudo-analytical (p.11–12); 2D: "For all the four experiments, our simulated velocities are within the order of magnitude of the experimental values" (Stehle, p.13). 3D µCT runs "only qualitative" (p.14).
- Neglected: "The influence of grain boundaries, volume diffusion in ice, viscous flow, or plastic deformations is neglected" (p.15); isotropic d0 and β (p.15).
- Sharp-interface requirement W ≪ 1/K and W bounds 5×10⁻⁷–5×10⁻⁵ m depending on β (Eq. 44, p.9).

## Claims supported / contradicted / qualified
- Draft L276 "were the first to attempt to use a phase-field model to capture the dynamics of DSM, by minimizing the thermodynamic free energy, coupled with temperature and mass balance equations": supported in substance (p.2: "no numerical model operates at similar length scales ... while at the same time containing more complete fundamental physics"; Eq. 33 "can be obtained by minimizing a free energy functional", p.7). "first" is implied, not stated; OK with "introduced".
- Draft L273 "Other studies have bounded α ∈ [10⁻³, 10⁻¹]" (no citation): Kaempfer & Plapp Table I states this range (p.4), sourced to Libbrecht (2006, unpublished departmental). Cite K&P (and Libbrecht 2005/2013 primary data) for the range.
- Draft L316 Eq. 3 (Pvs): K&P Eq. 2 confirms product form K5 ln T; draft "+ k5 + ln(T+T0)" is a typo (contradicted as written).
- Draft L339 "Our model reproduces, under the conditions outlined by K&P, ... Gibbs–Thomson": supported in form (Eq. 11; conditions Eqs. 42–46).
- Draft L494 "quantitatively validate their model in 1D and 2D by comparing them to a pseudoanalytical solution and bubble migration experiments": supported (p.11–13), but the 2D agreement is order-of-magnitude (p.13) — "quantitatively" should be "to order of magnitude in 2D". "a tactic that has since become popular" — unsupported by this paper.
- Draft L494 "governed by the variations in a double-well free energy potential": supported (p.7).
- Outline L710–718 (thin-interface limit): K&P use the sharp-interface limit, not a thin-interface limit — "the intrinsic properties of the ice-vapor interface ... justify the use of the simpler sharp-interface limit" (p.3). If the outline says "under the conditions outlined by K&P", those conditions are the sharp-interface constraints (Eqs. 42–46), not thin-interface.
- Outline L437 (Kaempfer & Plapp use α ~ Libbrecht): supported (Table I ref [31]).

## Candidate additions
- Methods (time scaling): cite K&P Sect. II D for scaling diffusion and source together, if the paper uses τ_sub-based acceleration (consistent with project memory).
- Methods Table 1: K&P Table I is the in-bib source for d0, γ, κ, Kj, Dv exponent.

# fourteauMacroscopicWaterVapor2021 — Fourteau, Domine & Hagenmuller, TC 15, 389–406 (2021)
DOI (PDF p.1): 10.5194/tc-15-389-2021. Read in full (18 pp).

## Summary
Theory + steady-state 3D FEM (ElmerFEM) of vapor diffusion in snow under a macroscopic T gradient (5–200 K/m, mean 258 K), on an idealized 3-sphere structure and 6 µCT snow samples (125–380 kg/m3). Surface kinetics via Hertz–Knudsen with sticking coefficient α from 0 to 1, infinitely fast, and nonlinear α = exp(−σ0/σ), σ0 = 0.01. Result: macroscopic effective vapor diffusivity is always < D0 (free air); faster kinetics raise it but loss of pore space wins.

## Findings relevant to the manuscript
- Effective vapor diffusivity always below free air: "the macroscopic vapor flux in snow is less than the flux in free air under the same water vapor gradient" (p.10, printed 398).
- Transition between slow and fast kinetics regimes at α ≈ 1e-3: "The transition between the fast- and slow-surface-kinetics regimes occurs for values of α around 10−3" (p.10/398). Relevant to the paper's α_c = 1e-3 choice (our aggregates sit at the regime boundary, consistent with outline's Λ ≈ 2.8).
- α not constant: "In general α is not a constant and depends on the local vapor saturation as well as the crystallographic properties of the underlying ice crystal (Saito, 1996; Libbrecht and Rickerby, 2013)" (p.8/396).
- Uses α = exp(−σ0/σ) with σ0 = 0.01 "determined for the attachment of vapor to the basal and prismatic facets of ice crystals (Libbrecht and Rickerby, 2013) and might not properly apply for the entirety of ice surfaces in snowpacks" (p.8/396). Supports outline's \tA bullet giving α_c = exp(−σ0/σ).
- Separation of scales / REV framing cites Auriault 1991 and Auriault et al. 2010 (the book; note the book is listed as 2010 here, Wiley Hoboken) (p.3/391; refs p.17).
- k_ice = 2.34, k_air = 0.024 W/m/K (Riche & Schneebeli 2013) (p.8/396) — the manuscript uses 2.29 and 0.02; not an error, different sources.
- Hansen & Foslien effective conductivity including latent heat discussed (Appendix C) — this paper is about vapor flux, not about apparent conductivity. It does NOT itself quantify how much vapor transport adds to k_eff.
- Does not treat curvature-driven (isothermal) transport: "in a homogeneous snow layer curvature differences are distributed isotropically ... and thus do not result in a net movement of water vapor" (p.5/393).
- Detailed snow models (Lehning 2002; Vionnet 2012) "assume that all mass transport results from diffusion" (p.13/401).

## Claims supported / contradicted / qualified
- Outline L565–568 ("Fourteau et al. (2021) show that vapor transport under a gradient adds to the apparent conductivity"): this paper is NOT the source for that; it is about effective vapor diffusivity. Outline L1357 correctly says the apparent-conductivity point is fourteauImpactWaterVapor2021 — check that one.
- Outline L444–448 (\tA α_c = exp(−σ0/σ)): supported as a form used in snow modelling, citing Libbrecht & Rickerby 2013 (p.8/396).
- Supports "α_c uncertain / not constant" framing (Introduction paragraph on attachment kinetics; Section 5.1 effective α_c).

## Candidate additions
- Could be cited for "the transition between attachment-limited and diffusion-limited behaviour in snow lies near α ≈ 10^-3" (p.10/398) — directly relevant to Section 5.1's Λ argument; but their quantity is macroscopic vapor diffusivity under a gradient, so cite only as context for the regime boundary, not for our isothermal case. Borderline; prefer not unless the regime statement is made.

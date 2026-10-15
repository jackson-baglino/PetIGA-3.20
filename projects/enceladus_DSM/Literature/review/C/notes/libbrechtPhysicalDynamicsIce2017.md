# libbrechtPhysicalDynamicsIce2017 — Libbrecht, Annu. Rev. Mater. Res. 47, 271–295 (2017)
DOI (PDF p.1): 10.1146/annurev-matsci-070616-124135. Read in full (27 pp incl. journal contents pages).

## Summary
Review of ice crystal growth from liquid water and from vapor near the triple point: morphology diagram, equilibrium energies, Gibbs–Thomson, attachment kinetics by layer nucleation (α = A exp(−σ0/σ)), edge-sharpening instability, air-dependent kinetics, spherical/cylindrical/Ivantsov diffusion solutions, numerical methods. No snow metamorphism, no sintering.

## Findings relevant to the manuscript
- Hertz–Knudsen / Wilson–Frenkel: "vn = α vkin σsurf, where vkin = (csat/cice)(kTsurf/2π mmol)^1/2 ... The dimensionless factor α, called the attachment coefficient ... For an infinite planar surface, we must have 0 ≤ α ≤ 1 on average" (p.5, printed 275). Supports draft L273 form.
- Layer nucleation: "On faceted surfaces, the attachment kinetics are often limited by the nucleation of new molecular terraces" (p.8/278); vapor case "vn ≈ α vkin σsurf ≈ A vkin σsurf e^{−σ0/σsurf}" (p.9/279). Supports outline L439–448 rewording (nucleation is of layers, not molecules) and the α_c = exp(−σ0/σ) form (with A ≈ 1).
- Where the measurements come from: "The preferred experimental method is to observe the growth of small ice prisms on a substrate in near-vacuum conditions ... A careful elimination of systematic effects finally allowed suitable measurements of αbasal and αprism over a broad range of conditions (65)" — ref 65 = Libbrecht & Rickerby 2013 (p.10–11/280–281). So the review reports, it does not "lead an experiment"; and both basal and prism were measured. Supports outline L426–433.
- Rough surfaces: "Nonfaceted surfaces (also called rough surfaces) generally have quite fast kinetics, suggesting αrough ≈ 1" (p.10/280). Relevant: rounded sintering grains are not faceted, so α near the upper end is plausible (bears on Molaro α_c = 0.1 choice and on α_c's uncertainty).
- In air: "Observations of the growth of thin ice needles at T = −5 °C indicate that αprism ≈ 0.002 in air, much smaller than αprism values measured at low pressures" (p.13/283).
- Surface energy: "γsv ≈ γsl + γlv ≈ 0.109 J/m2" (Antonow's relation, p.7/277). This is the source of the manuscript's γ_ia = 0.109 J/m² (Table 1) — a good citation for Table 1.
- Capillary length "dsv = γsv/csolid kT ≈ 1 nm" (p.8/278) — matches Table 2's d0 ≈ 1.0 nm.
- Shape relaxation time: "τequil,vapor ≈ R²/2dsv α vkin for ice in vapor" (p.8/278). With β = 1/(α vkin), this is R²β/(2 d0): the same form as the outline's τ_R̄ = R̄²β/d0 (L723–727). Closely matching citation for that time scale.
- Series attachment/diffusion (spherical growth, Eq. 2, p.15/285): αdiff = X0/R with X0 = csat Dv/(cice vkin); "X0 ≈ 0.15 µm in air at −15 °C"; useful α information only if α ≪ αdiff. Note X0/α = (Dv/α) sqrt(2πm/kT) = the outline's L* (L1240), so Λ = L*/R̄ = αdiff/α. This is the published source of the Λ group used in Sections 2.2/5.1/5.2.
- "Attachment kinetics is typically the dominant factor determining growth behavior when the supersaturation is low, or the crystal size is small, or there is little background gas pressure" (p.5/275).
- Gibbs–Thomson term usually negligible vs kinetic term in growth experiments "except when σ∞ is quite small" (p.16/286) — in sintering at saturation σ∞ is tiny, so curvature drives everything (consistent with our setting).
- Latent heat of deposition Lsv ≈ 2.8e6 J/kg (p.7/277) — Table 1 has 2.83e6.
- α depends on temperature through σ0(T) and A(T), and on air pressure (pp.11–14).

## Claims supported / contradicted / qualified
- Draft L273 "Libbrecht (2017) led an experiment in near vacuum ... for crystal growth on the basal plane": wrong-citation / wrong description. It is a review; the measurements are Libbrecht & Rickerby 2013 on basal and prism facets. Supports outline L426–433.
- Draft L276 lists Libbrecht 2017 among DSM works: wrong — not about metamorphism. Supports outline L492–497.
- Draft L276 "Many experiments attempt to parameterize kinetic growth rates in relation to temperature and vapor saturation (Colbeck 1983; Libbrecht 2017)": Libbrecht 2017 reviews exactly this (α(σ,T) parameterization, Fig. 6) — acceptable as a review citation; the measurement paper is Libbrecht & Rickerby 2013.
- Draft L273 "Other studies have bounded α_c ∈ [1e-3, 1e-1]": not stated here. Libbrecht gives α ≤ 1, α_rough ≈ 1, α_prism ≈ 0.002 in air at −5 °C. Does not support the bound.
- Outline L738 Cite box lists Libbrecht 2017 for "range of α_c": weak — it gives the physics and examples, not a range. Prefer Braun 2024 / Libbrecht & Rickerby 2013.
- Draft L273 "Vapor molecules do not always nucleate on icy surfaces due to their hexagonal crystallographic structure": contradicted in wording; Libbrecht: layer (terrace) nucleation on faceted surfaces (p.8/278).

## Candidate additions
- Table 1: cite Libbrecht 2017 p.7/277 for γ_ia = 0.109 J/m² (Antonow's relation) — currently no source given.
- Section 2.2 time scale: cite Libbrecht 2017 p.8/278 for τ ≈ R²/(2 d_sv α v_kin) as the shape-relaxation time (form of τ_R̄).
- Section 2.2/5.2 Λ: cite Libbrecht 2017 Eq. 2 (p.15/285) for the attachment/diffusion series form and X0.

# libbrechtPhysicalDynamicsIce2017 — Libbrecht (2017), Annu. Rev. Mater. Res. 47, 271–295

**Summary.** Review of ice crystal growth from liquid water and from vapor near the triple point: morphology diagram, equilibrium energies, premelting, Gibbs–Thomson, attachment kinetics (layer nucleation, measured A(T), σ0(T) on basal and prism facets), edge-sharpening instability, air-dependent kinetics, diffusion-limited growth (spherical, cylindrical, Ivantsov), numerical methods. Not about snow metamorphism or sintering.

**Findings relevant to the manuscript**
- Kinetic velocity, PDF p.5 (275): "vn = α vkin σsurf, where vkin = (csat/cice)(kTsurf/2πmmol)^{1/2} ... The dimensionless factor α, called the attachment coefficient"; "0 ≤ α ≤ 1 on average". (Libbrecht calls this the Wilson–Frenkel form.) ⇒ β_sub = 1/(α vkin) = (ρice/ρvs)√(2πm/kT)/α — matches the manuscript's Eq. for β_sub.
- Surface energy, p.7 (277): "surface wetting measurements suggest that Antonow’s relation is a reasonable approximation near the triple point: γsv ≈ γsl + γlv ≈ 0.109 J/m²" — this is the value in manuscript Table 1 (γ_ia = 0.109 J/m²).
- Capillary length, p.8 (278): "dsv = γsv/csolid kT ≈ 1 nm" — matches Table 2 d0 ≈ 1.0 nm.
- Relaxation time, p.8: "The time needed to relax to the ice equilibrium shape is τequil,vapor ≈ R²/2dsvαvkin for ice in vapor" ⇒ ≈ R²β/(2d0) — the same form as the outline's τ_R̄ = R̄²β/d0.
- Layer nucleation, p.9 (279): vn ≈ A_sv vkin σ e^{−σ0/σsurf}; p.11 (281) Fig. 5: "αbasal(σsurf) = Abasal exp(−σ0,basal/σsurf)"; measurements are ref (65) = Libbrecht & Rickerby (2013): "The preferred experimental method is to observe the growth of small ice prisms on a substrate in near-vacuum conditions" (p.10, 280).
- Rough surfaces: "Nonfaceted surfaces (also called rough surfaces) generally have quite fast kinetics, suggesting αrough ≈ 1" (p.10).
- Air: "Observations of the growth of thin ice needles at T = −5 °C indicate that αprism ≈ 0.002 in air, much smaller than αprism values measured at low pressures" (p.13, 283).
- When attachment dominates, p.5: "Attachment kinetics is typically the dominant factor determining growth behavior when the supersaturation is low, or the crystal size is small, or there is little background gas pressure".
- Diffusion vs kinetics, p.15 (285): spherical solution with "αdiff = X0/R with X0 = csatDv/(cicevkin)"; "it is quite difficult to obtain useful information about α unless α ≪ αdiff ... (for example, X0 ≈ 0.15 µm in air at −15 °C)". In air the Gibbs–Thomson (dsv) term "is usually negligible compared with the kinetic term, except when σ∞ is quite small" (p.16).
- Gibbs–Thomson in vapor: ceq ≈ csat(1 + dsvκ ...) (p.8).

**Claims supported / contradicted / qualified**
- Draft L273 "Libbrecht (2017) led an experiment in near vacuum to quantify ... σ0 ... for crystal growth on the basal plane": **wrong-citation** — the 2017 paper is a review; the measurements are Libbrecht & Rickerby (2013) (ref 65), and cover both basal and prism facets. Outline L426–433 correct.
- Draft L273 HK law vn = α_c v_kin σ_surf: **supported** (p.5), Libbrecht's notation.
- Draft L273 "Vapor molecules do not always nucleate on icy surfaces due to their hexagonal crystallographic structure": **wrong as phrased** — Libbrecht: on faceted surfaces growth is limited by nucleation of new molecular terraces (p.8–9); outline L439–443 rewording supported.
- Draft L273 "condensation coefficient has a non-trivial temperature dependence through σ0(T)": **supported** (Fig. 6, p.12); also depends on σ (α = A e^{−σ0/σ}) — outline L444–448 supported.
- Draft L276 Libbrecht 2017 in DSM list and as experiment parameterizing kinetic growth with temperature and saturation: **partially** — it is a crystal-growth review, not DSM (outline L492–497 correct). It does parameterize α(T, σ) for faceted growth from measurements.
- Draft L273 "α_c ∈ [10⁻³, 10⁻¹]" (no citation): Libbrecht 2017 does not state this range; it gives α ≤ 1, rough α ≈ 1, and α_prism ≈ 0.002 in air at −5 °C (p.13).
- Manuscript Table 1 γ_ia = 0.109 J/m²: **source found** — Libbrecht 2017 p.7 (Antonow's relation). Cite it (the draft gives no source; Molaro uses 0.06).
- Outline §2.2 "Two time scales of the same form ℓ²β/d0 ... τ_R̄ the time for a grain to change shape": **supported in form** by Libbrecht p.8 τequil ≈ R²/(2 dsv α vkin) (factor 2).
- Outline L459–463 / L1194–1197 (Λ = L*/R̄; attachment and diffusion comparable for α_c = 10⁻³, 50 µm grains): Libbrecht's α vs α_diff = X0/R criterion (p.15) is the published form of this comparison (X0 ≈ 0.15 µm at −15 °C ⇒ α_diff ≈ 3×10⁻³ for R = 50 µm, comparable to α_c = 10⁻³; consistent with Λ ≈ 2.8).

**Candidate additions.** (1) Cite Libbrecht (2017, p.7) for γ_ia = 0.109 J/m² in Table 1. (2) Cite Libbrecht (2017, p.8) for the shape-relaxation time R²/(2 d0 α vkin) when introducing τ_R̄ in §2.2/§5.2. (3) Cite Libbrecht (2017, p.15) α vs α_diff = X0/R for the attachment/diffusion comparison (Λ). (4) p.13: α_prism ≈ 0.002 in air at −5 °C supports treating α_c as an effective, environment-dependent parameter.

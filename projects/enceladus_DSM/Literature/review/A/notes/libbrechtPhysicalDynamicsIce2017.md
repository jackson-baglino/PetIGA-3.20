# libbrechtPhysicalDynamicsIce2017 — Libbrecht 2017, Annu. Rev. Mater. Res. 47:271–295 (doi 10.1146/annurev-matsci-070616-124135)

Whole paper read (pp. 271–295 + contents pages).

**Summary.** Review of ice growth from liquid and from vapour: Wilson–Frenkel form v_n = α v_kin σ_surf, equilibrium energies (latent heats, surface energies, step energies), premelting, Gibbs–Thomson, polynuclear attachment kinetics α = A exp(−σ0/σ_surf), review of the Libbrecht & Rickerby (2013) measurements (its ref. 65), edge-sharpening instability, air-dependent kinetics, spherical/cylindrical/Ivantsov diffusion-limited solutions, numerical methods (phase field, cellular automata). It is a review — it reports no new experiment.

**Findings relevant to the manuscript**
- p.275 (Sect. 2.2): v_n = α v_kin σ_surf with v_kin = (c_sat/c_ice)(kT/2πm)^1/2; "For an infinite planar surface, we must have 0 ≤ α ≤ 1 on average".
- p.277 (Sect. 3.1): "Antonow's relation is a reasonable approximation near the triple point: γsv ≈ γsl + γlv ≈ 0.109 J/m2"; "The ice surface energies are likely isotropic to a few percent". L_sv ≈ 2.8×10^6 J/kg.
- p.278: d_sv = γ_sv/(c_solid kT) ≈ 1 nm (capillary length in vapour). Relaxation to equilibrium shape τ ≈ R²/(2 d_sv α v_kin).
- p.280 (Sect. 4.3): "The preferred experimental method is to observe the growth of small ice prisms on a substrate in near-vacuum conditions ... Faceted basal and prism surfaces ... low attachment coefficients (α_basal, α_prism ≪ 1) ... Nonfaceted surfaces (also called rough surfaces) generally have quite fast kinetics, suggesting α_rough ≈ 1."
- p.281: "A careful elimination of systematic effects finally allowed suitable measurements of α_basal and α_prism over a broad range of conditions (65)" — ref. 65 = Libbrecht & Rickerby 2013.
- p.283: α_prism ≈ 0.002 in air at −5 °C (needles).
- p.288: phase-field models apply to weakly anisotropic systems; "Cellular automata methods ... have been the most successful for modeling ice growth from water vapor".
- No numeric statement of a range "10^-3 to 10^-1".

**Claims it supports / contradicts / qualifies**
- Draft L273 "Libbrecht (2017) led an experiment in near vacuum ... for crystal growth on the basal plane": **wrong-citation** — the 2017 paper is a review; the experiment (substrate, reduced pressure, basal AND prism) is Libbrecht & Rickerby 2013 (its ref. 65, p.281). Outline's correction (cite L&R 2013) is right.
- Draft L276 lists Libbrecht 2017 among works parameterising kinetics in T and supersaturation: **supported** as a review (Fig. 6 A(T), σ0(T)), though primary data are L&R 2013.
- Table 1 γ = 0.109 J m^-2: **supported** by Libbrecht 2017 p.277 (Antonow estimate) — this is the natural citation for the surface energy value (K&P 2009 also). L_sv 2.8×10^6 J/kg consistent.
- "α_rough ≈ 1" (p.280) is relevant to the choice of constant α_c ≤ 0.1: Libbrecht expects fast kinetics on rough (non-faceted) surfaces; the sintering grains here are rounded. So α_c = 10^-3–10^-1 is low relative to Libbrecht's rough-surface expectation; the manuscript should present α_c as an effective calibration parameter (as Braun 2024 / Bouvet 2022 do) and note that a larger α_c moves toward diffusion-limited regime. Qualifies α_c discussion (partial).
- Isotropic surface energy assumption in the model: **supported** p.277 ("likely isotropic to a few percent").

**Candidate citations / claims**
- Section 2 (Table 1 note): "γ = 0.109 J m^-2 from Antonow's relation (Libbrecht, 2017)".
- Section 2.2: "Rough, non-faceted ice surfaces are expected to have α close to 1 (Libbrecht, 2017), so the α_c values used here are effective coefficients" — only if the authors agree to frame it this way.

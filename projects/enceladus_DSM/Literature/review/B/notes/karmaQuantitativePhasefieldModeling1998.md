# karmaQuantitativePhasefieldModeling1998 — Karma & Rappel (1998), Phys. Rev. E 57(4), 4323–4349

**Summary.** The thin-interface limit for pure-melt solidification: an interface thickness W small compared with the pattern scale but not with d0. It gives capillary length d0 = a1·W/λ and kinetic coefficient β = a1(τ/(λW) − a2·W/D), so interface kinetics can be made to vanish. Includes calibration of grid anisotropy and converged 2-D and 3-D dendrite simulations. Not cited in the manuscript.

**Findings relevant**
- p.7 (Table I): thin-interface constants for the ±1 well with the standard tanh profile: a1 = 0.8839 and a2 = 0.6267 ("IVF") or 0.3981 ("VF"), depending on the interpolation function h(ψ).
- Conclusion: the thin-interface limit gives "accurate quantitative solutions" with "a much less stringent restriction" on W/d0.

**Claims supported / qualified**
- Tables.tex footnote c (a1 = 5, a2 = 0.1581) and outline L711–718 (Δβ = a1·a2·ε(1/D_T + 1/D_v)·ρ_ice/ρ_vs): the K&R/K&P constants are 0.8839 and 0.6267 for a W-scaled ±1 well. The manuscript's 5 and 0.1581 come from Moure & Fu's Supporting Information (cited so in src/material_properties.c: "M&F SI footnote 1") for their 0–1 well with profile 0.5(1 + tanh(x/2ε)). The SI is not in the library, so the values are **not-checkable here**. The source code states the two sets are equivalent once rescaled (a1·a2 ≈ 0.79 versus 0.78). Fix: the footnote should cite Moure & Fu (2024, SI) for 5 and 0.1581, and note that they correspond to Karma & Rappel's 0.8839 and 0.6267 under ε = W/√2 and the 0–1 well.
- K&P 2009 Eq. 41 (thin-interface Δβ) builds on K&R. Citing K&R 1998 for the thin-interface correction is closely matching.

**Candidate additions.** Cite Karma & Rappel 1998 in the "thin-interface limit" paragraph (outline L710) as the origin of the a1/a2 correction.

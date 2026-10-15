# flatauPolynomialFitsSaturation1992 — Flatau, Walko & Cotton (1992), J. Appl. Meteorol. 31, 1507–1513

**Page convention.** PDF p.1 is the JSTOR cover; PDF p.N = journal p.(1505+N), so PDF p.2 = 1507.

**Summary.** Sixth- and eighth-order polynomial fits in (T − T0) to Wexler's (1976, 1977) and Hyland–Wexler's (1983) saturation-vapour-pressure expressions over water and ice, plus their temperature derivatives. Wexler's reference expressions and coefficients are restated in §2 and Table 1. OCR has digit errors in the tables, so the numbers below were checked on rendered pages 2–3.

**Findings relevant to the manuscript**
- p.2 (1507) Eq. 2.3, rendered: ln(e_i,sat) = [k0 + {k1 + k5 ln(T) + [k2 + (k3 + k4·T)T]T}T]/T, i.e. k0/T + k1 + k2·T + k3·T² + k4·T³ + **k5·ln T**, with T in K and e in Pa.
- p.3 (1508) Table 1 (Wexler, ice), rendered: k0 = −0.58653696×10⁴, k1 = 0.2224103300×10², k2 = 0.13749042×10⁻¹, k3 = −0.34031775×10⁻⁴, k4 = 0.26967687×10⁻⁷, k5 = 0.6918651×10⁰.
- p.3 (1508): "where T0 is a conversion between degrees Celsius and kelvins. It is T0 = 273.15 for Wexler and Hyland–Wexler". This T0 belongs to Flatau's polynomial (3.1), not to Eq. 2.3, but it supports the draft's T0 = 273.15 K.
- Hyland–Wexler ice (Eq. 2.6, Table 2) is a different form, with m0…m6 and validity 173.16–273.16 K. ASHRAE's handbook formula follows Hyland–Wexler (from general knowledge, not in this PDF).

**Claims supported / qualified**
- Draft L314–316 Eq. 3 and its citation of Flatau & Wexler: **partially**. The functional form and coefficients are in Flatau Eq. 2.3 and Table 1, but the draft's "+ k5 + ln(T + T0)" is a typo for "k5 ln(T + T0)". The sum Σ_{j=0}^{4} k_j (T + T0)^{j−1} matches Eq. 2.3 when T is in °C.
- Draft L316 ("values for k_j may be found in Moure"): better to cite Flatau Table 1 (or Wexler 1977) directly; they hold the primary numbers. The solver (src/material_properties.c RhoVS_I) uses Flatau Table 1's ice coefficients rounded to 4 significant figures (K0 = −0.5865e4 … K5 = 0.6918). Recomputing gives a difference of 0.04–0.08 % in P_vs between −5 and −80 °C, which is negligible.
- Tables.tex L72 "ASHRAE polynomial": **wrong-citation / inconsistent**. The code and the draft text use the Wexler (1977) form as restated by Flatau (Eq. 2.3, Table 1), not ASHRAE's Hyland–Wexler form. The note should read "Wexler (1977), coefficients as in Flatau et al. (1992), Table 1".

**Candidate additions.** Quote Table 1 as the source of the k_j and state that T is in K.

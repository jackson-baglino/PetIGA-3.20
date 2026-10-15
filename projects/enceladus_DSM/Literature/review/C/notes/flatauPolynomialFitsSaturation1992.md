# flatauPolynomialFitsSaturation1992 — Flatau, Walko & Cotton, J. Appl. Meteor. 31, 1507–1513 (1992)
No DOI in the PDF (JSTOR stable URL 26186606). Read in full (8 pp incl. JSTOR cover); Eq. 2.3 checked on the rendered PDF p.2 (printed 1507).

## Summary
Sixth- and eighth-order polynomial fits (relative- and absolute-error norms) to Wexler (1976, 1977) and Hyland–Wexler (1983) saturation-vapor-pressure expressions over water and ice, and their temperature derivatives; speed comparisons on a CRAY-YMP.

## Findings relevant to the manuscript
- Wexler (1977) SVP over ice (Pa), Eq. 2.3 (p.2/1507): ln(e_i,sat) = [k0 + {k1 + k5 ln(T) + [k2 + (k3 + k4 T)T]T}T]/T, i.e. ln e = Σ_{j=0}^{4} k_j T^{j−1} + k5 ln T, T in K. Coefficients (Table 1, p.3/1508): k0 = −0.58653696e4, k1 = 0.2224103300e2, k2 = 0.13749042e−1, k3 = −0.34031775e−4, k4 = 0.26967687e−7, k5 = 0.6918651e0.
- Hyland–Wexler (1983) over ice (Eq. 2.6): ln e = Σ_{j=0}^{5} m_j T^{j−1} + m6 ln T, valid 173.16–273.16 K; Hyland & Wexler 1983 was published in ASHRAE Transactions (reference list p.7).
- Wexler's formulations "are based on more recent experimental data than that of Goff–Gratch" (p.2).

## Claims supported / contradicted / qualified
- Draft Eq. (3) (L314): P_vs = exp(Σ_{j=0}^{4} k_j (T+T0)^{j−1} + k5 + ln(T+T0)) — WRONG as printed: the source form is "+ k5 ln(T + T0)" (product, not sum). Confirms outline L675–676.
- Draft L316 cites Flatau 1992 and Wexler 1977 for Eq. (3): the five-plus-one-coefficient form matches Wexler 1977 (Flatau Eq. 2.3). Supported once the typo is fixed.
- Table 1 note b says "ASHRAE polynomial": inconsistent with the text. The ASHRAE paper is Hyland & Wexler 1983, a different (m0–m6) formula; Flatau's polynomial fits are a third option. The note should name the formula actually used (check the solver/Moure & Fu coefficients); supports outline L677–678.

## Candidate additions
- None.

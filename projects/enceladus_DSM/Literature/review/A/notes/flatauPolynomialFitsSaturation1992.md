# flatauPolynomialFitsSaturation1992 — Flatau, Walko & Cotton 1992, J. Appl. Meteor. 31(12):1507–1513 (no DOI in PDF; JSTOR 26186606)

PDF p.1 = JSTOR cover; journal p.1507 = PDF p.2.

**Summary.** Sixth- and eighth-order polynomial fits (relative- and absolute-error norms) to the Wexler (1976, 1977) and Hyland–Wexler (1983, ASHRAE Trans.) saturation-vapor-pressure expressions over water and ice and their temperature derivatives; error analysis and CRAY timing.

**Findings relevant to the manuscript**
- **Wexler (1977) ice formula** (Flatau Eq. 2.3, p.2): ln(e_i,sat) = [k0 + (k1 + k5 ln T + [k2 + (k3 + k4 T) T] T) T]/T, i.e. ln e = Σ_{j=0}^{4} k_j T^{j−1} + k5 ln T (Pa, T in K). Coefficients (Table 1, p.3): k0 = −0.58653696×10^4, k1 = 0.2224103300×10^2, k2 = 0.13749042×10^-1, k3 = −0.34031775×10^-4, k4 = 0.26967687×10^-7, k5 = 0.6918651×10^0.
- Hyland–Wexler (1983) ice (Eq. 2.6): ln e = Σ_{i=0}^{5} m_i T^{i−1} + m6 ln T, valid 173.16–273.16 K (Table 2) — published in ASHRAE Transactions.
- Polynomial fits e_sat = Σ a_i (T − T0)^{i−1}, T0 = 273.15 for Wexler/Hyland–Wexler (p.3); ice sixth-order fit valid −50 to 0 °C (Table 3), eighth-order −90 to 0 °C (Table 4); Hyland–Wexler eighth-order ice −75 to 0 °C (Table 5). Units of the polynomial fits: mb (a1 ≈ 6.11).

**Claims it supports / contradicts / qualifies**
- **Draft Eq. (3) (L314) "P_vs(T) = exp(Σ_{j=0}^4 k_j (T+T0)^{j−1} + k5 + ln(T+T0))": typo confirmed.** The Wexler (1977) form, as written by Flatau (Eq. 2.3), has **k5·ln(T)**, not k5 + ln(T). Outline L675–676 suspicion is correct.
- Draft L316 cites Flatau 1992 and Wexler 1977 for Eq. (3): supported for the functional form (six coefficients k0–k5 = Wexler 1977).
- **Table 1 note "(ASHRAE polynomial)"**: ASHRAE refers to Hyland & Wexler (1983, ASHRAE Trans.), whose ice formula has **seven** coefficients (m0–m6), or to Flatau's polynomial fits — neither is the six-coefficient Wexler (1977) form written in Eq. (3). The text and Table 1 must name one formulation; check which the solver uses (Moure & Fu 2024 give the k_j — see that note). Outline L677–678 correct to flag it.
- Units: Wexler/Hyland–Wexler give Pa; Flatau's polynomial fits give mb — if a Flatau polynomial is used, P_atm in Eq. (2) must be in the same unit.

**Candidate citations / claims**
- None new. Fix Eq. (3) and the Table 1 note.

**Verified on rendered page (mutool, PDF p.2):** Eq. 2.3 reads ln(e_i,sat) = [k0 + {k1 + k5 ln(T) + [k2 + (k3 + k4 T) T] T} T]/T — i.e. k5·ln T. Typo verdict stands.

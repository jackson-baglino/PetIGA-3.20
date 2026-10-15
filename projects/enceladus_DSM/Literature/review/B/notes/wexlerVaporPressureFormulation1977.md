# wexlerVaporPressureFormulation1977 — Wexler (1977), J. Res. NBS 81A(1), 5–20

**Summary.** Thermodynamic (Clausius–Clapeyron integration with virial EOS, ice expansion, Giauque–Stout heat capacity) formulation of the saturation vapor pressure over hexagonal ice from the triple point to −100 °C (IPTS-68). Anchored to the triple-point pressure 611.657 Pa. Gives a 6-term closed form (eq 54), a simpler empirical 3-term alternative (eq 63) and a table at 0.1 K steps (Table 6).

**Findings relevant to the manuscript**
- Eq. (54), PDF p.7 (journal p.11): "ln p = Σ_{j=0}^{4} K_j T^{j−1} + K_5 ln T" (T absolute, kelvin). The last term is K_5·ln T, a product — not "k_5 + ln T".
- p.1 abstract: uncertainty "varies from 0.016 percent at the triple point to 0.50 percent at −100 °C".
- Table 3 / Table 6 values (p.7, p.13): p = 103.276 Pa at −20 °C, 38.0239 Pa at −30 °C, 12.8486 Pa at −40 °C, 611.154 Pa at 0 °C.
- p.12 (journal p.16): eq (54) "is the preferred equation because it has a rational thermodynamic basis"; eq (63) is empirical.
- Latent heat of sublimation ~2830 J/g (p.8: "l (~2830 J/g)") — consistent with Table 1 L_sub = 2.83e6 J/kg.
- Wexler cites the ASHRAE Guide only as a user of Goff–Gratch (ref [75], p.16); the formulation itself is not "ASHRAE".

**Claims supported / contradicted**
- Draft L314 Eq. (3): form `Σ k_j (T+T0)^{j−1} + k_5 + ln(T+T0)` — **contradicted in form**: Wexler eq (54) has `k_5 ln(T+T0)`. Outline (L675) already flags this; confirmed.
- Draft L316 cites Flatau & Wexler for Eq. (3): Wexler is the right source for this functional form (Flatau gives polynomial fits — check there). Table 1 footnote b says "ASHRAE polynomial": inconsistent with the text; the ASHRAE Handbook form for ice is the Hyland–Wexler form, but no local source states that, so the manuscript should cite Wexler (1977) for the form, and say which coefficients are used.
- Table 1: L_sub 2.83e6 J/kg consistent with p.8.

**Candidate additions.** None beyond correcting Eq. (3) and making text/Table-1 note agree (cite Wexler 1977 for eq. 54 form; drop "ASHRAE" or cite the ASHRAE source properly).

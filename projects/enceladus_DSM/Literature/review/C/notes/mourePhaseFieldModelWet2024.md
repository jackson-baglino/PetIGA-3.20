# mourePhaseFieldModelWet2024 — Moure & Fu, Cryst. Growth Des. 24, 7808–7821 (2024)
Read in full (14 pp; SI not in library). ACS text layer garbles equations (checked forms against draft).

## Summary
Non-variational three-phase (ice, water, air+vapor) phase-field model for wet snow metamorphism (sublimation, evaporation, melting, nucleation terms), with temperature and vapor-density equations; 1D and 2D tests (0.2 mm domain, 400×400). Reduces to the two-phase Kaempfer & Plapp dry model when φw = 0.

## Findings
- Title/scope: "we focus on wet snow metamorphism" (abstract); "Under thermal conditions when the liquid water phase is absent (φw = 0), our model reduces to the two-phase model for dry snow metamorphism first proposed in Kaempfer and Plapp" (p. ~7).
- Pvs^I(T) = exp(Σ_j k_j (T+T0)^{j−1} + k5 ln(T+T0)) (Eq. 8, refs 56, 57); text prints "T0 = 273.15 °C" (should be K). k_j in SI Table S2 (not in library).
- Dv(T) = Dv0 ((T+T0)/T0)^1.81 (Eq. 15) — matches draft.
- Time scaling: "multiplying the right-hand side of the T and ρv equations ... by a time-scale factor (ξT and ξv)" following Kaempfer & Plapp; ξv = 10⁻³ (Sect. 3.1) — i.e., diffusion and source scaled together.
- Gibbs–Thomson recovered only "under certain parameter regimes"; the K&P condition needs ε < 10⁻⁸ m in their wet case, they used ε < 10⁻⁶ m, so "less accurate in capturing interface kinetics along the air−ice ... interfaces" (Sect. 3.2).
- 2D dry example: four grains at −5 °C "quickly sinter (before t1 = 40 s)", coarsening slows, bubble trapped (Fig. 3).

## Claims
- Draft L282/L290 "updated phase-field model introduced by Moure & Fu ... largely based on": supported, but it is a wet-snow model whose two-phase reduction is used (outline L498 correct).
- Draft L316 "values for k_j may be found in Moure & Fu": partially — they are in the SI (Table S2), not in the article; and the article's own K&P Table I gives them (in-bib). Cite K&P Table I or Wexler.
- Draft Eq. 3 typo ("+ k5 + ln"): contradicted by Moure Eq. 8 (product).
- Draft L339 "Our model reproduces, under the conditions outlined by K&P, ... Gibbs–Thomson": qualified — Moure state the condition is met only for small ε; the paper should state its ε against the K&P bound (W ≪ 5×10⁻⁷–5×10⁻⁵ m, depending on β).

## Candidate additions
- Methods: cite Moure Sect. 3.1 / K&P Sect. II D for scaling diffusion and source together.

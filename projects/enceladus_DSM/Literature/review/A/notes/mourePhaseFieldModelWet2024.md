# mourePhaseFieldModelWet2024 — Moure & Fu 2024, Cryst. Growth Des. 24:7808–7821 (doi 10.1021/acs.cgd.4c00539)

Whole paper read (pp. 7808–7821). Eqs. 7–9 verified on rendered p.7810.

**Summary.** Nonvariational three-phase (ice/water/air) phase-field model with temperature and vapor density, Allen–Cahn kinetics plus α-terms for solidification, sublimation, evaporation; piecewise-constant mobility; K&P time scaling (ξ_v = 10^-3, ξ_T = 1 or 10^-2); FEM (bilinear), generalized-α time integration, Newton. 1D freezing/nucleation, 2D dry four-grain example (−5 °C, 0.2×0.2 mm², 90 min), 2D wet snow (vapor, temperature, LWC). No experimental comparison.

**Findings relevant to the manuscript**
- p.7810 Eq. 7: ρ_vs^J(T) = ρ_a 0.62 P_vs^J/(P_a − P_vs^J). Eq. 8 (rendered): P_vs^I(T) = exp(Σ_{j=0}^{4} k_j (T+T0)^{j−1} + k_5 ln(T+T0)), "where T0 = 273.15 °C" (T in °C), citing Flatau et al. 1992 and Wexler 1977 (refs 56, 57); k_j "listed in Table S2 in the Supporting Information".
- p.7811 Eq. 15: D_v(T) = D_v0 ((T+T0)/T0)^1.81, citing Massman 1998 (ref 58).
- p.7812 Sect. 3.1: time scaling "we leverage the procedure explained in Kaempfer and Plapp ... multiplying the right-hand side of the T and ρv equations ... by a time-scale factor (ξT and ξv)"; ξ_v = 10^-3.
- p.7813 Sect. 3.2: K&P condition ε < Dβ′ "then no antitrapping current is needed"; air–ice requires ε < 10^-7 m (garbled; text says most restrictive ε < 10^-8 m); "we consider a value ε < 10−6 m ... our simulations are less accurate in capturing interface kinetics along the air−ice and air−water interfaces".
- p.7815 Sect. 4.2: "Under thermal conditions when the liquid water phase is absent (φw = 0), our model reduces to the two-phase model for dry snow metamorphism first proposed in Kaempfer and Plapp" (Eqs. 30–32).
- p.7818: "A key characteristic of dry snow metamorphism is the reduction in the microstructure surface area while maintaining an almost constant ice volume".
- p.7809: "Currently, the model does not account for granular compaction".

**Claims it supports / contradicts / qualifies**
- Draft Eq. 3 (L314) "k_5 + ln(T+T_0)": **typo** — Moure Eq. 8 has k_5 ln(T+T0) (product). Also draft L316 "T_0 = 273.15 K" is right (Moure writes "°C", a slip on their side) and T is in °C in this form.
- Draft L316 "The values for k_j may be found in Moure" — they are only in Moure's Supporting Information (Table S2), which is not in the library; the original source is Wexler (1977) (Flatau 1992 gives polynomial fits instead). Cite Wexler 1977 for the k_j (and list them in an appendix or Table 1 note).
- Draft Eq. 9 (L332) D_v exponent 1.81: matches Moure Eq. 15; Moure's source is Massman (1998), which is not cited in the draft. Suggest citing Massman 1998 (not in library; cannot verify value — flag as not-checkable) or "as in Moure & Fu (2024)".
- Draft L282 "updated phase-field model introduced by Moure": **partially** — in the dry case Moure's model is K&P 2009 (p.7815). Say "the dry-snow phase-field model of Kaempfer & Plapp (2009), in the implementation of Moure & Fu (2024)". The \xi time scaling is K&P's procedure (Moure Sect. 3.1).
- Draft L276 (Moure in list of DSM efforts): fine, but Moure is a wet-snow model; its dry case is a four-grain demonstration.
- Constant porosity in the model (fixed ice volume under isothermal coarsening): consistent with Moure p.7818 (a property of the vapor-only model, no compaction p.7809) — the manuscript should present it as a model property, not as a property of real snow (cf. Kaempfer & Schneebeli 2007, Löwe 2011 densification).
- Interface-width compromise: Moure accept ε larger than the K&P bound for air–ice with reduced kinetic accuracy — relevant to the manuscript's thin-interface discussion (eps note in project CLAUDE.md).

**Candidate citations / claims**
- Sec 2.1: "With no liquid present the model of Moure & Fu (2024) reduces to that of Kaempfer & Plapp (2009)".
- Table 1 note: "D_v ∝ T^1.81 (Massman, 1998, as used by Moure & Fu, 2024)".

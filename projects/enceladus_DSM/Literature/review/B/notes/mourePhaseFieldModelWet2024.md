# mourePhaseFieldModelWet2024 — Moure & Fu (2024), Cryst. Growth Des. 24, 7808–7821

**Summary.** Non-variational three-phase (ice–water–air) phase-field model with temperature and vapor density: Allen–Cahn kinetics with a triple-well potential plus source terms for solidification, sublimation and evaporation; phase-dependent mobility (piecewise constant on the ternary diagram); time-scaling of heat/vapor equations (ξ_T, ξ_v) after Kaempfer & Plapp; 1-D and 2-D FE simulations (bilinear elements, generalized-α, Newton, adaptive Δt) of freezing, nucleation, dry (4 grains, −5 °C, 0.2 mm box, 90 min) and wet metamorphism.

**Findings relevant to the manuscript**
- Title/scope: "A Phase-Field Model for Wet Snow Metamorphism"; abstract p.1: "Our three-phase model recovers the corresponding two-phase transition model when one phase is not present".
- Two-phase dry model, PDF p.8 (journal 7815), Eqs. (30)–(32): ∂φi/∂t = −3M0[(1/ε)φi(1−φi)(1−2φi) − ε∇²φi] + α_sub φi²φa²(ρv−ρvs^I)/ρi; (1/ξ_T)ρ(φ)c_p ∂T/∂t = ∇·[K∇T] − ρ(φ)L_sub ∂φa/∂t; (1/ξ_v)∂(φaρv)/∂t = ∇·[φa Dv∇ρv] + ρ_SE ∂φa/∂t. "first proposed in Kaempfer and Plapp".
- Saturation vapor pressure, p.3 (7810), Eq. (8): P_vs^I = exp(Σ_{j=0}^{4} k_j (T+T0)^{j−1} + k5 ln(T+T0)), "the values of the parameters k_j ... are listed in Table S2 in the Supporting Information"; cites Flatau (ref 56) and Wexler (ref 57). Eq. (7): ρ_vs = ρ_a·0.62·P_vs/(P_a − P_vs). Note Moure writes "T0 = 273.15 °C" (unit typo).
- Diffusivity, p.4 (7811) Eq. (15): Dv = Dv0((T+T0)/T0)^1.81 citing Massman (1998, ref 58).
- ε, p.3: "The parameter ε represents the width of the phase-field interface" (Moure's own wording).
- Temporal scaling, p.5 (7812): "we take ξv = 10−3 and ξT = 1 or ξT = 10−2".
- Sharp-interface limit, p.6 (7813): "if we impose ε < Dβ′ ... then no antitrapping current is needed to correctly recover the sharp-interface limit" (from Kaempfer & Plapp); Gibbs–Thomson for air–ice (Eq. 20) recovered "under certain parameter regimes" (p.5).
- Scope, p.2: "the model does not account for granular compaction, melt film imbibition, or melt percolation".
- Dry 2-D example: "ice grains quickly sinter (before t1 = 40 s) and the coarsening process gradually slows down" (p.8); "A key characteristic of dry snow metamorphism is the reduction in the microstructure surface area while maintaining an almost constant ice volume" (p.11).

**Claims supported / qualified**
- Draft L282/L290/L494 (model from Moure & Fu; "two-phase equivalent"): **supported** (Eqs. 30–32). Outline L498–500 (say it is a wet-snow model, two-phase form used): **supported** (title, abstract).
- Draft L300–308 Eq. (1): matches Moure Eqs. 30–32 except the draft omits 1/ξ_T and 1/ξ_v (outline L663–664 correct to flag ξ_v) — **partially**.
- Draft Eq. (1b) latent term ρ(φ)L_sub: identical to Moure Eq. (31). The outline's correction (L657–659: should be ρ_ice "as the solver does") is a change *from* the published model; the paper should say so (e.g., "we use ρ_ice in the latent term, unlike Moure & Fu (2024) Eq. 31"), not present it as a typo.
- Draft L308 "ε represents the interface width": this is Moure's wording (p.3) but contradicts the project convention (ε is the interface parameter; visible band 9.2 ε) — outline correction stands.
- Draft L314 Eq. (3) with "+k5 + ln(T+T0)": **contradicted** — Moure Eq. (8) (and Wexler eq. 54) have k5·ln(T+T0).
- Draft L316 "The values for k_j may be found in Moure & Fu": **partially** — they are in Moure's Supporting Information Table S2 (not in the local library); say "Supporting Information of".
- Draft L332 Eq. (5) Dv(T) with exponent 1.81: matches Moure Eq. (15); Moure attributes it to Massman (1998). Add that citation (not in library; flag).
- Draft L297 "temperature (units: °C)": Moure also uses °C with T+T0 (p.2) — outline's unit-consistency point stands.
- Draft L339 "under the conditions outlined by Kaempfer & Plapp": Moure p.6 states the condition (ε < Dβ′) — the manuscript should state it explicitly.

**Candidate additions.** (1) Cite Massman (1998) for Dv(T) (as Moure does). (2) State the sharp-interface condition ε < D·β′ from Kaempfer & Plapp/Moure p.6 in §2.2. (3) Note departures from Moure & Fu explicitly (ρ_ice latent term, thin-interface Δβ, no ξ_T scaling).

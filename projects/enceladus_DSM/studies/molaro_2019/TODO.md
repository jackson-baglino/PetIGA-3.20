# Grain-pair sintering validation — to-do list

Covers the Molaro (2019) comparison, the historical datasets and the Demmenie
(2025) exponent test. Tick items as they close, with the date and where the
result lives, and commit. Jackson submits every HPC run; the assistant only
prepares the commands.

## Framing for the manuscript (decided 2026-10-02)

- **Molaro runs: α_c = 0.1 is a tuning choice.** It is the upper end of the
  literature band [1e-3, 1e-1], taken to best reproduce the data. Say so
  explicitly ("parameter tuning"). The final runs are
  `GrainPairSintering/batch_2026-09-08__17.20.46_molaro_T-20_round2`
  (second attempt) and `..._T-5_h0.99674_2h_a1e-1_dirichlet` (first attempt).
- **Demmenie run: α_c = 1e-3 has a theoretical reason.** It is the
  attachment-limited regime, where Kuczynski's t^(1/3) is derived.
- **The remaining Molaro rate deficit is surface diffusion.** That mechanism
  is in SA81 and Molaro's model, and absent from ours. Molaro's own weak
  point is the vacuum Hertz–Knudsen vapour term (their Appendix A).
- **Thomas et al. (1994): exponents only, no run**
  (`studies/historical_sintering/README.md`). Their a ≈ 0.33–0.40 implies at
  least saturated conditions. Do not claim super-saturation without more
  evidence.

## Now

- [x] **Back-of-the-envelope neck ODE** (2026-10-02),
  `studies/molaro_2019/envelope_ode/`.
  - No literature α_c reproduces the Molaro rate (α_c = 0.1 gives 11.6 % rms,
    even fully saturated).
  - The data's a = 0.20 is the saturated, diffusion-limited value.
  - Predicted saturated a: ≈ 0.20 at α_c = 0.1, ≈ 0.30 at 1e-3.
- [ ] **Per-face boundary-condition option** in `src/enceladus_main.c`: keep
  the x = 0 contact plane natural Neumann (mirror) while the outer walls take
  Dirichlet ρ_v / T. Build locally and add a unit-style check.
- [ ] **Quarter-domain geometry and opts** for the saturated
  equal-grain Demmenie run:
  - α_c = 1e-3, one grain, mirror at x = 0, axis at y = 0.
  - Pre-necked start (r₀ above the resolution floor).
  - Walls at ρ_v = ρ_vs(T)·(1 + 2d₀/R).
  - Output every 1–2 min equivalent.
  - The ODE predicts ~96 h simulated to reach x/R = 0.35.
  - **Open: the temperature and grain radius.** Molaro-like (−20 °C,
    R = 84.4 µm) or Demmenie-like (−3 °C, R = 500 µm)?
- [ ] **Check `neck_width.py` on a mirror-domain snapshot** (the neck sits on the boundary).
- [ ] **Submit** (Jackson). On hold while HPC priority is penalised.
  Budget ≲ $15.
- [ ] **Compare the run's a against Demmenie** (0.26–0.33) and the ODE
  prediction (0.29–0.30).

## Deferred / dropped

- Undersaturated control arm for the quarter domain (h = 0.99715): dropped
  2026-10-02, to save HPC priority. Equal grains are not expected to shift a
  materially.
- Saturated α_c = 0.1 quarter-domain arm: superseded by the α_c = 1e-3 run.
  The ODE predicts a ≈ 0.20 there.
- Thomas (1994) and Kingery (1960) replication runs: dropped
  (cost / mechanism).

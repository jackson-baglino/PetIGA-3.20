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
  - The ODE and the Fig. 2 run differ on how much undersaturation reaches the
    neck (ODE best fit ≈ 0; run 1.6e-4). Both agree the ambient was
    undersaturated (shrinkage gives s_∞ ≈ 3.6e-3). Explain both in the
    manuscript rather than reconcile them (decided 2026-10-02).
- [x] **Per-face boundary-condition option** (2026-10-02): `-bc_mirror x0`
  (any of x0,x1,y0,y1,z0,z1). Named faces stay natural Neumann (mirror) while
  `-flag_BC_rhovfix` / `-flag_BC_Tfix` pin the others. Unknown face names abort.
  The BC banner labels them "(mirror)". Builds clean. **Not yet exercised by
  a run**: check the BC table in `outp.txt` of the first job.
- [x] **Quarter-domain geometry and opts** (2026-10-02):
  - `inputs/geometry/molaro/molaro_2D_L220x188um_eps0.12um_axisym_T-20eq_R84um_r14um_mirror.opts`
  - `inputs/experiment/molaro/molaro_T-20_hsat_100h_a1e-3_dirichlet_mirror.opts`
  - Batch file: `batches/demmenie_mirror_T-20.txt`.
  - Settings: −20 °C, R = 84.4 µm (the Molaro pair's effective radius),
    α_c = 1e-3, walls h = 1 + 2d₀/R = 1.0000241, pre-necked r₀ = 14 µm,
    ε = 1.18e-7 m (17.9 M DoF), t_final = 100 h, hourly output.
  - dtmax = 5·τ_sub = 546 s keeps it at ~$9 (2·τ_sub would be ~$20; the
    ladder puts 5·τ_sub within 1 %).
  - Temperature and radius chosen Molaro-like rather than Demmenie-like
    (−3 °C, 500 µm): the kinetic-limit exponent does not depend on either,
    and this keeps the ODE and Fig. 2 run directly comparable.
- [x] **`neck_width.py --mirror-x0`**: reads the neck from the x = 0 column.
  `run_batch_measure.sh` passes it automatically when the opts carry
  `-bc_mirror x0`. `grain_shrinkage.py` still assumes two grains, so ignore
  its output for this run.
- [x] **Submit the saturated Demmenie mirror run** (Jackson). DONE 2026-10-07: ran 100 h (job 4066835).
  Released
  2026-10-06 (user), queued with the last k_eff stages. Re-verified that day:
  h = 1 + 2 d0/R = 1.0000241 (d0 = 1.0152e-9 m), set by `-humidity` as BOTH
  the initial pore vapour and the Dirichlet walls; 17.9 M DoF, 90 ranks;
  tau_sub = 108.9 s. Add `-- --time=0-12:00:00` to the command below.
  - Single job, ~$9, budget ≲ $15. Everything is committed and ready.
  - Command:
    `./scripts/HPC/submit_batch.sh --tag demmenie_mirror --tests-file studies/molaro_2019/batches/demmenie_mirror_T-20.txt --out-root /resnick/groups/rubyfu/jbaglino/simulation_outputs`
  - First check: the BC table in `outp.txt` shows x = 0 as "(mirror)".
- [x] **Compare the run's a against Demmenie** (0.26–0.33) and the ODE
  prediction (0.29–0.30). 2026-10-07: a = 0.19–0.23 by every fit form
  (free-t0 0.19–0.20 whole record, 0.22–0.23 above 40–50 um; Kuczynski
  m = 5.1–5.3; local slope 0.14 at 10 h rising to 0.21 at 100 h). Grain
  radius +0.06 %, so saturation held. BELOW both targets.
- [x] **Gibbs–Thomson prediction for this run (2026-10-09).**
  `demmenie/predict_neck_growth.py` → `demmenie/prediction/neck_prediction.png`.
  The neck ODE (σ = d0κ + βv_n, c = 0.275 from the envelope study, not
  refitted) gives 60.6 µm at 100 h against the run's 63.2 µm, and a = 0.29
  against the run's 0.21. The run is strongly attachment-limited (diffusion /
  attachment resistance 0.004–0.02), so the "mixed regime" idea is OUT. The
  difference is early: the run's neck velocity is ~2.5× the prediction at
  w = 35 µm and converges to it above ~55 µm. At 35 µm the fillet radius is
  only ~2 interface bands (6 at 63 µm): a finite-ε effect is the leading
  candidate. Thin-interface term: OFF in this run (β_realised/β_requested =
  1.000). At Demmenie's own conditions (−3 °C, mm spheres, 2.5 h) the vapor
  route alone gives a = 0.32–0.33 for α_c ≤ 0.01, 0.25–0.26 at 0.1, 0.19–0.20
  at 1: their measured range is reachable by vapor transport alone.
- [ ] **The −20 °C Molaro run (2026-09-08) predates the thin-interface fix
  (5ecdd236, 2026-09-09):** its τ_sub is 1.3403 s against 1.0891 s kinetic-only,
  so β_realised = 1.23 × requested (effective α_c ≈ 0.081). The −5 °C run
  (2026-09-29) is at 1.000. Manuscript Fig. 3 therefore mixes the two. Rerun
  −20 °C with the current solver, or state it; Table S3's τ_sub = 1.34 s is the
  inflated value.
- [ ] **ε-refinement of the saturated pair** (finer ε, same set-up): does the
  early excess shrink and the exponent rise toward 0.29? Also the one check
  the gt_deficit note leaves open.
- [ ] **Explain the low exponent.** Candidates to test, none checked yet:
  the pre-necked start (r0 = 14 um, no physical t = 0); the slope still
  rising at 100 h (not yet asymptotic); neck/grain ratio already 0.17–0.37
  (outside the small-neck limit the 1/3 law assumes); how the ODE prediction
  was fitted (same form and window?).

## Deferred / dropped

- Undersaturated control arm for the quarter domain (h = 0.99715): dropped
  2026-10-02, to save HPC priority. Equal grains are not expected to shift a
  materially.
- Saturated α_c = 0.1 quarter-domain arm: superseded by the α_c = 1e-3 run.
  The ODE predicts a ≈ 0.20 there.
- Thomas (1994) and Kingery (1960) replication runs: dropped
  (cost / mechanism).

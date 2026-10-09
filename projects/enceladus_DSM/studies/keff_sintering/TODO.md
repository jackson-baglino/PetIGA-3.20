# k_eff sintering campaign — to-do list

The working list for the manuscript runs. **Read this first in a new session.**
The history for writing -- every finished run, where its data is, what it
showed -- is [`RECORD.md`](RECORD.md); update it whenever a task finishes.
Tick items as they close, add the date and where the result lives, and commit.
Background and reasoning: `CAMPAIGN.md` (plan and results by stage),
`studies/rve_anisotropy/README.md` (domain size and packing bias),
`.claude/ACTIVITY_LOG.md` (what happened, session by session).

**Where the data lives (since 2026-10-01).** Every manuscript run is in ONE
folder on the cluster:
`.../simulation_outputs/enceladus_DSM/keff_sintering_campaign/<geom>__<exp>/`.
Each submission leaves a record in `stages/<stage>__<timestamp>/` (manifest,
stage file, job ids, inputs/src snapshot). An existing run folder is never
resubmitted over. Download a stage with
`./scripts/HPC/fetch_stage.sh <stage file>` (tables; `--full` or
`--full-run <substr>` for snapshots). It mirrors into
`~/SimulationResults/HPC_results/enceladus_DSM/keff_sintering_campaign/`.
Tests (scaling etc.) go one folder per test:
`.../enceladus_DSM/scaling_<date>[_suffix]/`.

**Standing rules**
- Every manuscript run goes through `scripts/HPC/submit_keff_production.sh`.
  That script owns every option; stage files list runs only. Never hand-type
  `--extra-opts` for these runs.
- Submit in stages. A stage's check must pass before the next is submitted.
- Jackson launches every run; the assistant only prepares commands.
- Push before submitting. The script refuses to run from an uncommitted or
  unpushed checkout.

---

## Now: batch 3, the porosity × temperature matrix (125 runs)

- **Packings:** `inputs/packings/keff_LR40/`, φ {0.275, 0.325, 0.375, 0.425,
  0.475} × 5, seeds 1601–2005, rebuilt 2026-10-01 with the seam and
  percolation gates only (option B, below). The gated build that 3a ran on is
  `inputs/packings/keff_LR40_gated/` (seeds 301–905).
- **Temperatures:** T {−5, −10, −20, −30, −40} °C.
- **Index:** `batch3_phi_T.txt` (never submitted whole).

- [x] **3a shakedown** — ran 2026-09-28 (commit f4179a9, 401 ranks), analysed
  2026-09-29. Results: `~/SimulationResults/HPC_results/enceladus_DSM/GrainPackingSintering/keff_b3a_shakedown/`
  (`compare/` holds the overlays). Cost **$195** (estimated $60–100).
  - [x] all 7 reach 30 d; health clean (mass drift ≤ 5e-7, 0 KSP failures,
    symmetry ≤ 1e-11, 0 bound trips)
  - [x] `k_eff.csv` every 5 steps; every step at −40 °C (108 samples)
  - [x] opening frame at t = 1.26 s; 45 snapshots, 8.1 GB/run
  - [x] **s/step at 401 ranks is no faster** (15–28 s vs ~20 at 241), and the
    k_eff solve is SLOWER: 0.85 s/iteration vs 0.49 at 241 ranks. So 60k
    DoF/core costs 1.66× for nothing.
    → **DECISION PENDING: revert `scripts/lib/alloc.sh` to 100k**
  - [x] −5 °C used 16.6 h of 24 (≈11–13 h expected at 241 ranks)
  - [x] φ 0.475 stable. k_iso falls monotonically with φ (0.86 → 0.31 W/m/K at
    30 d). Rise 1→30 d: +25.9 / +25.0 / +24.7 / +22.6 / +19.7 %.
  - [x] T collapse on the 0.325 packing: k_iso at matched SSA agrees to 0.06%;
    the speed-up equals the τ_sub ratio to ±1.5% (−5: 3.83×; −40: 0.128×)
  - Finding: at φ 0.475 (ice percolates in y only) k_xx rises +6% against
    +26% for k_yy; k_xx/k_yy goes 0.50 → 0.37.
- [x] **k_eff cadence → SSA trigger** (2026-09-29).
  - The 3a curves had visible corners at every 5 steps (0.49 d = 5.4 τ_sub
    apart). The new setting samples every 0.1% drop in SSA after 11 τ_sub,
    every 5 steps before that, with a 20 τ_sub backstop.
  - The same option at every temperature; the −40 rule is gone.
  - Max error 0.02% of the plotted range, 60× below the smallest real kink
    (the SSA ≈ 15,300 event on seed 301, ≥ 1.3%).
  - Samples per run ~308 / 179 / 46 at −5 / −20 / −40.
  - Code: `src/keff_sample.c` (KeffDue). `studies/keff_sintering/predict_cadence.py`
    mirrors it and reproduced all 7 3a schedules on the old path; run
    `--check` on the rerun.
- [ ] **Rerun 3a** with the final options, on the option-B packings
  (seeds 1601/1701/1801/1901/2001). First submission (2026-10-01, commit
  d7452d5, 1.09 tau_sub) CANCELLED while still pending: dtmax went to
  2 tau_sub on 2026-10-02, and jobs read inputs/ at start, so they would have
  run at 2 tau_sub under a manifest saying otherwise. Resubmit at the 2 tau_sub
  commit. Now 8 runs: + gated seed301 at -20 C, the paired dtmax check against
  the first shakedown's 1.09 run (pass: k_iso and SSA within 0.5% after
  11 tau_sub). Expected sampling cost of 2 tau_sub (batch-2 data at stride 2):
  max k(SSA) interpolation error 0.03-0.18%, still 7x under the smallest
  real kink.
  - 2026-10-03: 3 of 8 back and clean (see RECORD.md). Still to read: the
    0.375/0.425/0.475 runs at -20, 1701 at -5 (wall time), and the gated-301
    dtmax pair.
  - Changed for 3b: -keff_freq 5 -> 1 (first-day corners at 2 tau_sub);
    manifest run-list bug fixed.
  - 2026-10-03: user's call to submit 3b and batch_rve WITHOUT waiting for the
    rest of 3a (fair-share queue waits are long; the 3 runs read were clean).
    Risk accepted: if the gated-301 dtmax pair fails, every 2 tau_sub run
    (3a rerun, 3b, rve; ~$170) is redone at 1.09. src changed since 3a
    (458acdb, -bc_mirror): default off and periodic runs have no Dirichlet
    faces, so results are unaffected. The first 3a is kept as
  the shakedown record, not pooled. Check it as 3a was, plus
  `predict_cadence.py --check` on every run.
- [ ] **Before 3b** (cost, not correctness), in this order:
  - [x] **DoF/core scaling test**, ran 2026-09-29 (`scaling/README.md`).
    - The fewest ranks is cheapest: 198k DoF/core (121 ranks) is 3–5×
      cheaper than 60k.
    - k_eff is 15–30× slower at ≥ 8 nodes, from communication; the suspect
      is the mixed node-type constraint.
    - [ ] follow-up (~$5): bracketed `--constraint` at 241 ranks; 300k/400k
      targets.
    - [x] `scripts/lib/alloc.sh` set to 200k DoF/core (121 ranks) on
      2026-09-30.
    - [x] Memory (2026-09-30): the flat 1G OOM-killed the 61-rank jobs, and
      the 121-rank runs had peaked at 860 MB (84%). `scripts/lib/alloc.sh`
      `mem_per_cpu` now sizes `--mem-per-cpu` from measured MaxRSS:
      peak = 0.20 GB + 2.35 GB/MDoF-per-rank + 8 B × total DoF (rank 0),
      × 1.5. That gives 2G for production, 3–4G for L/R 80. Not billed.
      Note: the solver's "memory after setup" guard reads ~22% where the
      real peak is 84%, so don't trust it for sizing.
    - [x] test 2 (300k/400k), 2026-10-01: cost flat to ±10% from 61 to 161
      ranks. Kept 200k (121 ranks); 400k is the queue fallback.
      `scaling/README.md`.
    - [~] test 1 (bracketed constraint): cancelled after a day pending on
      Priority; moot at ≤ 6 nodes. Also check the login node's CPU
      (`lscpu`): the build is `-march=native`, so a binary built on an Ice
      Lake login node can SIGILL on Skylake/Cascade Lake nodes. That is the
      probable cause of the old "job starts, nothing runs" hangs.
  - [ ] k_eff solver benchmark: `solver_benchmark/` (~$10, replay only).
    Now LOW value: at ≤ 6 nodes a sample is ~5 s, ~15 min per −20 °C run.
    Optional. Adopt a faster setting only if max |Δk|/k_iso < 1e-5, and record
    it in the production script before 3b.
  - Projected cost of all 125 runs at 200k DoF/core with the SSA cadence:
    roughly $0.9k. That's per-temperature $/run from `scaling/README.md`
    (−30/−40 °C are cheaper still) × 25 runs each.
- [x] **3b** done 2026-10-04 (RECORD.md). The fabric-anisotropy check FAILED
  as written: k_yy overtakes k_xx above φ ≈ 0.35 (open question below).
- [x] **Convergence study** done 2026-10-04: rise size-independent (RECORD.md).
  - Check: seed scatter of the k_iso rise per φ (~3–4% expected at 0.325);
    the φ trend is larger than the scatter.
  - Check: k_xx/k_yy per φ follows the contact fabric (1.03 at 0.275 →
    ~1.1 at 0.425).
- **Reordered 2026-10-03 (user):** 3c–3e carry seeds 1–3 of each φ; seeds 4–5
  at −5/−10/−30/−40 °C moved to `batch3f_seeds45.txt` (40 runs), last and
  optional. Why: T is a per-packing time rescaling, so the −20 °C column
  (5 seeds) gives the seed-mean k(SSA); 3 packings per φ test the collapse
  and supply the time axis. Run 3f only if a φ fails the collapse or a
  reviewer asks for 5 everywhere.
- [x] **3c, 3d, 3e** done 2026-10-06 (RECORD.md): collapse holds at every φ.
- [x] **3r redo** `batch3r_redo.txt` (9 runs) and **3f** `batch3f_seeds45.txt`
  (40 runs): done 2026-10-07, all 49 clean; matrix complete at 125 (RECORD.md).
  Queued 2026-10-06. Before submitting 3r, move the old folders
  aside with `scripts/HPC/supersede_runs.sh` (header of the stage file).
  Watch φ 0.425 seed1902 at −40 °C: it hung twice.
  - [x] Old local copies of the nine redone runs moved to
    `LOC/keff_sintering_campaign/superseded_before_3r/` before the fetch
    (~30 GB; delete once nobody needs the old attempts).
  - [ ] Manuscript figures rebuilt on the full set are in
    `LOC/Figures2/` — the iCloud `Manuscript/Figures` folder was not writable
    from the session (macOS permission), so it still holds the 2026-10-06
    versions and the pre-swap numbers for Figs. 4 and 5.
- [ ] **After each stage:**
  - [ ] download the tables (see the `rsync` recipe in `ACTIVITY_LOG.md`, 2026-09-26)
  - [ ] run `postprocess/plot_keff.py` per run
  - [ ] run `postprocess/compare_keff.py` per batch
  - [ ] run `health_check.py`
  - [ ] record the measured cost here

## Supplement: k_eff domain-size convergence (planned 2026-09-30)

Does k_eff tend to one curve as the domain grows, and is L/R 40 close enough?

- **Design:** φ 0.325, −20 °C. L/R {20, 30, 40, 56, 80} with {8, 6, 5, 4, 3}
  seeds, all UNGATED (seam and percolation gates only). The L/R 40 point is
  the production 0.325 set (1701–1705, from the 3a rerun and 3b), now the
  same recipe; rve seeds 1501–1505 are built but not run. Plus the 5 gated
  pre-B packings (301–305) in `batch_rve.txt`, to measure the gate effect.
  Statistical-RVE test: the seed-mean curve stops moving with L and the seed
  scatter shrinks ~1/L. (A bigger periodic cell is a new realization, so
  there is no "same packing, larger".)
- [x] Packings (2026-09-30): `make_rve_packings.sh` ->
  `inputs/packings/rve_phi0.325/`, 26 built.
  - The first build switched void gates at L/R 40 and z_band jumped at the
    switch, so it was rebuilt with the homogeneity gates off at every size.
- [x] Opts at −20 °C + stage file `batch_rve.txt`, which passes the
  production script's checks.
- [x] **Allocation for L/R 56/80** (2026-10-01): `MAX_NODES_PER_JOB=6` in
  `scripts/lib/alloc.sh`, applied by `submit_batch.sh`. They run on 192
  ranks (250k / 500k DoF/core) instead of 9 / 16 nodes.
- [ ] Submit `batch_rve.txt` (26 runs, ~$155). Independent of batch 3; any
  time after the 3a rerun passes.
- [ ] Analysis: `rve_convergence/analyze_rve.py <rve stage> <3a> <3b>`. L/R 40
  is close enough if its gap to the L/R 80 mean is inside the gap's standard
  error, at 11, 100 and 330 τ_sub. Also report the gated-vs-ungated L/R 40
  k_eff gap (methods: why B).

- [x] **DECISION: production gates → option B** (user, 2026-10-01).
  - The homogeneity gates (void, density CV, half-domain asymmetry) test an
    extreme over the domain; at L/R 40 they filtered ordinary realizations.
    Gated z_band 3.32 ± 0.07 vs 3.48 ± 0.08 ungated at φ 0.325.
  - `keff_LR40` rebuilt with the seam and percolation gates only (all 25 on
    their base seed); z_band at 0.325 is now 3.46. Seed-mean z_band rose at
    0.275/0.325/0.375 by 0.14–0.20, unchanged at 0.425/0.475. φ 0.475
    percolates in both axes in 4 of 5 (1 of 5 gated).
  - Stage files remapped (301→1701, 601→1601, 701→1801, 801→1901,
    901→2001, …).

## Manuscript argument (agreed with Jackson, 2026-10-06)

Framed after Molaro et al. (2019, JGR Planets) and Choukroun et al. (2020,
GRL); both PDFs are in `Literature/`.

**Lead claims**
1. **k is a function of microstructural state, and temperature only sets the
   clock.** One sintering age θ = t/τ_sub(T) collapses every temperature
   (k at matched SSA within 0.27%, 59 paired runs). The novelty is the
   non-dimensionalization: it extrapolates to temperatures we did not run,
   as Molaro's later figures do for neck growth.
2. **After the transient, k grows as a weak power of age**, k ∝ θ^0.066
   (≈ +17% per decade), nearly the same for every well-connected packing.
   Weaker than claim 1: the porosity collapse is looser, and it is lost at
   φ ≥ 0.425, attributed to poor connectedness, which 2D amplifies.
3. **Level before drift**: porosity and packing set k; sintering is a slow
   correction on top. Claimed for φ ≤ 0.375 only.

**Minor (1–2 sentences each)**
- Anisotropy is constant in time: expected for non-densifying sintering, so
  it is a consistency check on the homogenization, not a finding.
- Small sensors: packing scatter ±9% in absolute k at 2 mm; the change with
  age is size-independent.

**Do not claim**
- Anything from φ > 0.40 beyond "connectivity is being lost"; those cases
  need 3D and larger domains.
- Absolute conductivities (2D, air-filled pores).
- "Rise from day 1" as width-independent.

**Enceladus framing**
- Target the near subsurface (~1 m): the overlying ice seals the pore vapour,
  which is what our closed periodic cell models, and gradients are weak
  there (little insolation at the south pole).
- Attachment-limited kinetics is what lets the result leave the laboratory:
  with Knudsen transport in a vacuum pore of size d, L* = D_K·β_HK =
  4d/(3α_c), so the regime holds for any α_c ≪ 1, at every temperature.

## Analysis tasks (2026-10-06)

- [ ] **Manuscript figures** — plan in `MANUSCRIPT_PLAN.md`. Build now,
  rebuild when 3f/3r land (Jackson, 2026-10-06).
  - [x] Fig. 6 merged (snapshots + collapse): `figures/fig_keff_collapse.py`,
    copied to `Figure6__KeffEvolution/Figure6_keff_collapse.*`.
  - [x] Figures renumbered 2026-10-06 (2 = homogenization method, 3 = Molaro,
    4 = k_eff evolution, 5 = gallery, 6 = state law, 7 = timescales); folders
    renamed to match.
  - [x] Fig. 2 methods figure complete (corrector replay of the master run,
    4 steps, `<master run>/corrector/`).
  - [x] One tables file `tables/tables.tex` (material, k_eff fixed, k_eff
    temperature, grain pair), one command per table.
  - [x] Fig. 1 = (a) mechanisms + (b) `aggregate_strip` (decided 2026-10-06; no setting panel).
  - [x] Fig. 3 variant with grain shrinkage as (c)(d): `Figure3_molaro_validation`.
  - [ ] Fig. 5 gallery, Fig. 7 state law, Fig. 8 timescale map: samples in
    `LOC/keff_sintering_campaign/compare/figure_samples/`
    (`figures/sample_figures.py`). Jackson to give figure numbers/folders.
  - [ ] Fig. 7 message agreed? (k is a power law of SSA; porosity sets the
    level). Fig. 8 waits on the surface-diffusion literature check.
  - [ ] Methods figure (§2.3), `figures/fig_Figure2_homogenization_method.py`:
    panels (a) cell and (b) neck zoom drawn; (c) corrector and (d) heat flux
    wait on ONE local replay of the master run. `-keff_write_corrector` was
    declared but never implemented; implemented 2026-10-06 (writes
    `igakeff.dat` + `t_vec_<step>_<m>.dat` beside the CSV). Command in the
    script header. Jackson runs it.
  - [x] Figure folders and LaTeX parameter tables created in the manuscript
    folder (list in `MANUSCRIPT_PLAN.md`). Figure numbers 4–6 are provisional;
    the methods figure's number is open.
  - [ ] Fig. 7: choose `Figure7_state_law` or `Figure7_state_law_alt` (adds k_eff(t) by
    porosity, pairing with the gallery).
- [ ] **Rewrite §3.2, §4, §5, key points and abstract** (plan in
  `MANUSCRIPT_PLAN.md`).

- [ ] **Surface diffusion at ~180 K: literature review before any claim.**
  The two framing papers point the other way from "negligible when the
  quasi-liquid layer is gone": Choukroun measured Q = 24.3 ± 3.3 kJ/mol over
  193–243 K and attributes it to surface self-diffusion (~23 kJ/mol, Nasello
  2007), explicitly not vapour (~51); Molaro has surface diffusion leading
  while necks are small. Our clock has Q ≈ 50 kJ/mol (corrected 2026-10-08 from 48). A lower Q wins at LOW
  temperature, so the vapour route should lead when warm and lose when cold.
  Read Nasello 2007 and the snow-sintering sources Molaro cites (Maeno &
  Ebinuma 1983; Löwe 2011; Vetter 2010) and decide how to state it.
- [ ] **THIN-INTERFACE CORRECTION: the manuscript runs are not consistent, and
  the project default disagrees with the user's intent (found 2026-10-09).**
  Audit: `thin_interface_audit.py` → `thin_interface_audit.txt`.
  - 125 aggregate runs: **OFF** in all (banner `-thin_iface_corr 0`). ON would
    add 1.9 % to τ_sub at every temperature, so the collapse, the state law and
    every dimensionless result are unaffected; only the clock shifts by 1.9 %.
  - Molaro pair −20 °C (2026-09-08): **ON** (pre-flag code), β realised 1.23 ×
    requested, effective α_c 0.081.
  - Molaro pair −5 °C (2026-09-29): **OFF**; ON would add 21.5 %.
  - Saturated pair (2026-10-06): **OFF**; ON would add 0.2 %, so the exponent
    result there does not depend on this.
  - History: enceladus removed the terms 2026-09-09 (5ecdd236, default 0);
    lunar switched its default back to 1 on 2026-09-13 (f7cdfbe5) because the
    OFF argument does not bound the spurious surface diffusion of a one-sided
    model without an anti-trapping current. enceladus was never brought in
    line. The user wants the correction ON (stability with the contact-angle
    model; precedent in Kaempfer & Plapp and Moure & Fu).
  - **DECIDED 2026-10-09 (user): the model is with the correction ON.** The 125
    aggregate runs are kept and relabelled: off at α_c = 1e-3 is the same
    calculation as on at **α_c = 1.020e-3** (exactly: the term enters only
    through τ_sub; 1.0199–1.0200e-3 across the five temperatures). τ_sub and
    every result are unchanged. `tables/tables.tex` updated (α_c, β_sub, new Δβ
    row and formula). The "1.019" quoted first was an estimate with an inferred
    thermal diffusivity; 1.020 uses the solver's own (D_T = 7.78e-6 m²/s).
  - Saturated pair: off at 1e-3 = on at 1.002e-3; no rerun.
  - Molaro −5 °C pair: off at 0.1 = on at 0.129, so it does not match the
    −20 °C pair (on at 0.1). Rerun with `-thin_iface_corr 1` (~12 h, ~$70) or
    keep and state it: the neck ODE puts the change in neck growth at −2 to
    −4 %. USER TO DECIDE. Table S3 carries a [REVISIT] note until then.
  - [x] 2026-10-09 (user): the `-thin_iface_corr` switch is REMOVED from the
    enceladus solver; the terms are always included. lunar still has the
    switch (default on) and was not touched.
  - **Consequence for reproducing the campaign:** the experiment files
    (`inputs/experiment/snow/snow_T*_h1.00_30d.opts`) still carry β_sub0 for
    α_c = 1e-3. With the terms now always on, resubmitting a stage file would
    run α_c = 1e-3 WITH the term, which is not the campaign (that is 1.020e-3
    with it). The β_sub0 values and `generate_study_opts.py --alpha-c` need
    changing before any campaign run is repeated. NOT changed; user to decide.
  - [ ] Rerun the Molaro −5 °C pair with the term (user, 2026-10-09).
  - [x] α_c = 1.020e-3 written into the tables, outline, plan and run table.
- [ ] **Decide the manuscript's clock: τ_sub (ε-dependent) or τ_R = R̄²β_sub/d₀.**
  Raised by the user 2026-10-08. ε = R̄/50 is a discretization choice and
  τ_sub = ε²β_sub/d₀ is the diffuse interface's relaxation time. The eps ×2
  runs show the physical late-time rate does NOT follow τ_sub: τ_sub ×4.00,
  yet k_eff rises 8.1/9.4/8.8 % (ε 1 µm) vs 8.0/9.6/8.7 % (2 µm) from day 8
  to 30. So θ = t/τ_sub is a good clock across temperature at fixed ε, but its
  numbers (θ_r = 30, θ = 331, "11 τ_sub") are tied to ε = 1 µm; the same state
  is θ = 83 at ε = 2 µm. The ε-free age is t/τ_R = θ/2500, from the
  sharp-interface law β v_n = u − d₀κ. Options: restate figures in t/τ_R, or
  keep θ and say once that τ_sub = τ_R/2500 here. Derivation:
  `timescale_map/timescale_map_derivation.pdf`, Section 3.
- [x] **Figure 8 (timescale map): DROPPED 2026-10-08 (user).** Replaced by a few
  sentences in the discussion; no Enceladus timescales in years. The paper is
  reframed around laboratory-observable sintering (MANUSCRIPT_PLAN.md,
  "Framing"). Notes kept for the record:
  Derivation and sources: `timescale_map/timescale_map_derivation.pdf`
  (numbers from `timescale_map/check_timescale_map.py`). Open items found
  while writing it up (2026-10-08):
  - [x] annotations checked against Literature/ (2026-10-08): plume grains
    0.1–5 µm and the Choukroun marker confirmed; surface corrected to 50–80 K;
    "fractures 175–185 K" replaced by a Tiger Stripes tick at 180 K. Primary
    sources to cite: Kempf 2010, Southworth 2019, Howett 2010, Spencer &
    Nimmo 2013 (all via Choukroun 2020);
  - [x] the map now uses the model's own ρ_vs (τ_sub ∝ √T/ρ_vs), so it equals
    the solver's τ_sub at the simulated temperatures; no pressure formula.
    The clock's activation energy is 50 kJ/mol (the "48" in older notes came
    from an ideal-gas ρ_vs);
  - Molaro's Table 6 gives τ ∝ R² above ~10 µm, independent support for the
    grain-size axis; their timescales at 180 K are ~100× shorter than ours;
  - R² scaling is an argument, never run at a second grain size, and needs
    R ≪ L* = D_v β_HK (139 µm at 253 K; R = 50 µm is already 0.36 L*);
  - the Choukroun marker compares different end states (their 10 MPa vs our
    θ = 331); 13.6 yr vs 15 yr is not agreement.
- [ ] **Compare our clock with Choukroun's and Molaro's numbers.** At 180 K
  and R = 6 µm our vapour route reaches the 30-day state in ~14 yr;
  Choukroun gets "very consolidated" (10 MPa) in ~15 yr. At 80 K they
  diverge (ours: never; theirs: ≥ 100 Myr to 1 MPa). Reproduce Molaro's
  timescale-vs-T-and-grain-size figure with our τ_sub and overlay.
- [ ] **F(θ, connectivity)**: first pass in `master_curve/` (README). Next:
  a connectivity measure that carries to 3D; refit with 3f and 3r; SSA form.
- [ ] **Vacuum-pore k_eff** (k_void → 0): a second k_eff solve on the
  snapshots we already have. **PARKED (Jackson, 2026-10-06): on the list,
  not started.** No new sintering runs; local.
  - Why: on Enceladus the pores are empty, so all heat crosses the necks and
    the sensitivity to sintering should be far larger than our +17% per decade.
  - Feasible on this Mac (64 GB, 12 cores): the corrector is 8 M unknowns, a
    few GB. The k_eff solve is cheapest on FEW cores (530 core-s at 61 ranks,
    630 at 121, 15–30× slower across ≥ 8 nodes), so ~1 min per snapshot with
    air; longer as k_void falls (contrast). Split the list across the two
    Macs; do not link them by MPI.
  - Exactly zero is singular: sweep k_void = 2e-2, 2e-3, 2e-4 W/m/K on one
    snapshot and see whether k_eff plateaus.
  - Unresolved necks carry all the heat in vacuum, so expect a strong
    interface-width dependence: replay the eps 1 µm and 2 µm pairs FIRST. If
    they disagree badly, the vacuum numbers are not publishable from these
    meshes.
  - First timing test (untried locally; four snapshots of seed 1701, −40 °C):
    ```
    R=~/SimulationResults/HPC_results/enceladus_DSM/keff_sintering_campaign/packing_2D_phi0.325_Rave50um_LR40_seed1701_L2mm_eps1000nm_perxy_T-40__snow_T-40_h1.00_30d
    ./scripts/Studio/run_enceladus.sh packing_2D_phi0.325_Rave50um_LR40_seed1701_L2mm_eps1000nm_perxy_T-40 \
        snow_T-40_h1.00_30d vac2e-3 -- -keff 1 -keff_replay "$R" -thcond_air 2e-3 \
        -keff_csv "$R/k_eff_void2e-3.csv" -keff_ksp_type cg -keff_pc_type gamg
    ```
  - Then: script the sweep; refit F on the vacuum values.
- [ ] **Is the rise really highest at φ 0.325?** In the day-1 → day-30 table
  0.325 ≥ 0.275 in 4 of 5 temperatures, but those four columns are the SAME
  three packings per porosity (T only rescales time), so it is one comparison
  of 3 vs 3, with gaps (0.6–1.7 points) inside the seed scatter. At −20 °C
  (5 packings) it is a tie, 27.4 vs 27.2%. Recheck with 3f. If it holds,
  count closed pores in the snapshots (Jackson's idea: pores sealing early
  limit sintering in dense packings). Against it so far: SSA falls FASTEST
  at φ 0.275 (exponent −0.085 vs −0.080 at 0.475). A φ 0.225 point (Optional)
  would show whether the rise turns over.
- [ ] Anisotropy "reversal" = k_xx/k_yy crossing 1 with porosity (1.02 at
  φ 0.275, 0.64 at 0.475) against a contact fabric that leans the other way.
  Confined to φ > 0.40, so out of the paper's claims; keep as a note.

## Optional — ready if reviewers ask, or if the write-up needs support

Not staged (Jackson, 2026-10-06). Each is cheap and independent.

- [ ] **α_c generality.** τ_sub ∝ 1/α_c, so α_c should only rescale the
  clock while R ≪ L* = D_v·β_HK ∝ 1/α_c (139 µm at 1e-3; 14 µm at 1e-2).
  Seed 1701 (or 1701–1703) at −20 °C with α_c = 1e-4 (deeper in the regime;
  ~$2 each) and 1e-2 (R > L*: where the collapse should break; ~$11 each).
  Needs new experiment files (`--alpha-c`) and a stage file.
- [ ] **eps = 0.5 µm** on seed 1701 (~$40): is the production early-time
  curve converged?
- [ ] **A porosity below 0.275** (e.g. 0.225 × 3 seeds, −20 °C): does the
  relative rise turn over at low porosity? At present 0.275 and 0.325 tie.

## Open questions and follow-ups

- [ ] **Hangs at −30/−40 °C** — diagnosed 2026-10-06: a cluster-side fault,
  not the packing. All four hung jobs ran 4–6× slower than healthy ones
  BEFORE stopping (k_eff solve 11–15 s vs a median 2.7 s over 55 healthy cold
  jobs, same iteration counts), on hpc-19-22, hpc-20-35, hpc-19-28 and
  hpc-34-38. Seed 1902 hung twice on two different slow node sets. Same work,
  slower per iteration, then a stall: the nodes or their interconnect.
  Mitigation in place: the 30-min stall watchdog. φ 0.425 seed 1902 at −40 °C
  is rerun in 3r. If hangs recur: requeue on stall, or the bracketed
  `--constraint` (one node type per job).
- [ ] **eps = 0.5 µm check** (~$40, one run on seed 1701): is the production
  early-time curve converged? eps ×2 changed the rise from day 1 by 9 points
  and nothing after day 8 (`eps_sensitivity/README.md`). Decide before
  writing the "rise" numbers.
- [ ] **Baseline for the quoted rise.** 11 τ_sub (1 d) sits inside the
  width-dependent transient. Consider a later baseline, or state it.
- [ ] Early-time corners of 1–2% remain on every-step runs (the ×1.33 ramp
  gives 8 steps per decade of time). `-factor 1.1` would cut them ~10× for
  ~35% more steps; not changed mid-campaign.

- [ ] **Anisotropy reversal at high φ** (2026-10-04): seed-mean k_xx/k_yy
  1.02 → 0.64 from φ 0.275 to 0.475, against the contact fabric; not the
  chord lengths either. Test a directional backbone measure (spanning-cluster
  mass or max-flow along x vs y). `studies/rve_anisotropy/README.md`.
- [ ] **Final-step rollback ends the run early** (up to 0.25 d, loses the
  30 d snapshot): a CFL/bounds rejection of the step that crosses t_final is
  deferred to a pre-step that never runs. Fix after the campaign (loop
  TSSolve while a rollback is pending); data unaffected.

- [ ] **k_eff solve cost.** 50–84 s per sample in batch 2, against 4 s on
  rev64 for the same solve at similar iteration counts, with 10× swings
  within a run. Looks environmental; unexplained.
  - Proposed: a small replay benchmark (tolerance 1e-7, GAMG threshold,
    hypre, pc_freeze) with `-log_view`. A few dollars, no new simulation.
  - Awaiting a go-ahead.
- [ ] **rev64 seeds 2–4** (2026-09-21, jobs 3264515–3264517). Are they on
  scratch? If so, fetch `k_eff.csv` + `SSA_evo.dat`: they test the
  sintered-state RVE directly.
- [ ] **Seam relaxation reach.** It moves 86 grains by more than 0.1 R (up
  to 1.2 R), reaching ~30% of L from the seam; seen in the deposition movie.
  - The generator's comment says only grains near the seam move. Correct the
    comment; decide whether it matters.
- [ ] **φ 0.475 connectivity at the band.** The generator's percolation test
  uses the sharp geometry, but the solver joins gaps smaller than 9.2·eps.
  Measure connectivity at the band before writing about 0.475.
- [ ] Figures normalized by k_eff(t = 0), e.g. the AGU snapshot set, are a
  post-processing choice (user, 2026-09-29), not a problem.
  - k(0) is the unrelaxed analytic IC: the first ~3.7 τ_sub are relaxation,
    worth ~+35–50% in k.
  - Say in the caption that the normalization includes that relaxation.
  - Before comparing porosities that way, check whether k(11 τ)/k(0) differs
    with φ.
- [ ] HPC support: job 3606286 (3a, −5 °C) stuck in COMPLETING on hpc-34-37
  since ~06:30 2026-09-29; billing should have stopped at job end (check `sacct`).
- [ ] **Kink near day 7.5** in seed 1, −20 °C (batch 2): the time step is
  flat there, so it's likely a topology event. Look at the snapshots.
- [ ] **Anisotropy direction vs real snow.** Our packings give k_xx > k_yy.
  - Check against measured snow anisotropy before claiming anything.
  - Cite only closely matching work, with DOIs verified via Crossref.
- [ ] Correct `studies/packing_design/README.md`: "vertical load paths"
  predicted k_yy > k_xx; the data say the opposite.
- [ ] `effective_thermal_cond/docs/tensor_conductivity_law.tex`: cite Nicoli,
  Plapp & Henry 2011 (prior art for the tensor law).
- Sequential metamorphism-then-k_eff (a separate replay job): considered
  2026-09-30 and not worth it at ≤ 6 nodes. k_eff is then ~10% of the cost,
  a split saves <$0.10/run, and it adds 30–55 GB of snapshots plus a second
  queue wait. Revisit only if k_eff becomes the dominant cost again.

## Later (CAMPAIGN.md stages 6–7)

- [x] **DECIDED (user, 2026-10-04): eps ×2 only**, staged as
  `batch_eps2.txt` (3 runs, < $2). All other arms dropped.
- [ ] Submit and read `batch_eps2.txt`.
- [~] (superseded) which sensitivity arms, if any.
  Criterion (user): run only what gives the manuscript data. Recommendation:
  keep eps ×2 (3 runs, < $2 total: the one reviewers will ask for, since the
  interface width is the model's main numerical choice and the tensor-law
  ladders tested static k only, not the evolution); α_c 1e-2 only if the
  paper discusses where the T collapse breaks; skip R_ave ×0.5 (pure time
  rescaling, studies/grain_size_scaling), R_ave ×2 and σ_ln 0.2.
- [ ] Sensitivity arms from the stage-5 centre point, 3 seeds each:
  - α_c {1e-4, 1e-3, 1e-2}: may change the kinetic regime, and with it the
    T collapse;
  - R_ave ×0.5 / ×2 at fixed L/R;
  - σ_ln 0.2;
  - eps ×2.
- [ ] Analysis:
  - partial correlation of k_eff with SSA at fixed porosity;
  - anisotropy from seed-mean k_xx/k_yy (never per seed);
  - per-snapshot metrics (Euler characteristic, chord lengths, D_eff);
  - lever ranking with seed error bars.

## Manuscript obligations (things the text must say)

- [ ] **Necks are under-resolved**, by design: floor r/R = √(12·eps/R) ≈ 0.49.
  No neck-growth exponent is claimed.
- [ ] **In this model, temperature is a pure time rescaling by τ_sub**
  (attachment-limited, α_c = 1e-3). State it as a model property, not a
  result about snow.
- [ ] **Anisotropy comes from the drop-and-roll deposition rule.** Report it
  as a seed mean; k_xy is a zero-mean fluctuation that shrinks as 1/L.
- [ ] **The y-seam gate** (≥ 0.76 of interior contact density) and its
  residual bound: k_yy at most ~1–1.5% low at L/R 40.
- [ ] **φ 0.475 is past 2D solid percolation**: it marks the limit of a 2D
  snow analogue.
- [ ] **The baseline is t = 1 d**, not t = 0 (IC relaxation). Compare
  temperatures at matched t/τ_sub or matched SSA.
- [ ] **Pore connectivity:** fragmented at t = 0, D_eff ≤ 0.0014 D_v.
  Long-range ripening is suppressed in 2D.

## Done (most recent first)

- 2026-09-30 — One normalization for every k_eff figure: `plot_keff.py`,
  `compare_keff.py` and `plot_keff_snapshots.py` open on the first k_eff
  sample with t >= 1 s, call it t = 0, and divide by its measured values
  (k_eff,0, SSA_0). The 11 τ_sub interpolated baseline ("b") and the greyed
  relaxation are gone; seed means cover only the span every seed does.
  Opening at 1.7 s vs the snapshots' 126 s differs by 0.12 % in k_iso.
  Colours: single runs cmocean `deep` (k_xx dashed, k_yy dotted, k_iso
  solid), temperature `thermal`, porosity amp-to-black.
  `coefficient_fix/compare_laws.py` keeps its own `--baseline-days` (a study
  record, not a figure script).

- 2026-09-29 — 3a analysed. Fixed on the way:
  - `health_check.py` keyed runs by seed, so it skipped the same packing at
    other temperatures;
  - `compare_keff.py` and `plot_keff.py` baselines are now 11 τ_sub, and
    interpolated. With −40 °C present, "1 d" was 1.4 τ_sub, inside the
    relaxation;
  - `run_batch_postprocess.sh` takes POSTPROCESS_DIR/PYTHON overrides.

- 2026-09-28 — Production submit script `scripts/HPC/submit_keff_production.sh`;
  this list.
- 2026-09-28 — Batch 3 staged into 3a–3e; 125 opts generated.
- 2026-09-28 — Deposition movie + storyboard
  (`postprocess/make_deposition_movie.py`); outputs in
  `~/SimulationResults/HPC_results/enceladus_DSM/GrainPackingSintering/presentation/`.
- 2026-09-28 — Packings: φ 0.275–0.475 by 0.05 × 5, unique seeds, y-seam
  gate, porosity-salted RNG.
- 2026-09-28 — RVE and anisotropy check from existing data
  (`studies/rve_anisotropy/`): L/R 40 is adequate; the seam is the only
  generator artifact.
- 2026-09-27 — Batch 2 (warm end, 9 runs) analysed: the dt check passes, T is
  a time rescaling, it cost $342; k_eff every 5 steps is enough.
- 2026-09-26 — k_eff plots (absolute, normalized, off-diagonal), compare
  plots; `k_eff.csv` naming fix; t ≥ 1 s opening frame.

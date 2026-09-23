# Campaign: how sintering changes k_eff

The executable plan, with a **stop rule at every stage**. `PLAN.md` holds the
study design and `ROADMAP.md` the earlier staging; this file supersedes both
wherever the tensor conductivity law changed the answer, and says so each time.

**The deliverable is comparative** — how much sintering changes `k_eff`, and
which levers matter most — not absolute `k_eff` in W/m/K.

**Scope is steady-state.** The cell problem, the surface excesses and the
tensor law concern the corrector problem on a frozen `φ`. Nothing here may be
used to justify an interpolation change in the transient metamorphism model
(`ρc`, latent-heat placement). If something there looks analogous, it is a
separate question.

## How to read a checkpoint

Each stage names **one number**, the threshold it must clear, and what happens
if it does not. Costs are at the tier-1 Resnick rate, $0.012/core-hour
(`scripts/HPC/hpc_cost.sh`); cost ≈ ranks × wall-hours × rate. The pilot's
measured **~24 s/step at 241 ranks** is the basis for every wall-clock estimate
below, and stage 2 replaces it with a measurement at other temperatures.

| stage | what | cost | running | figure |
|---|---|---|---|---|
| 0 | tensor law + both ladders | done | — | ✔ committed |
| 1 | pilot replay under tensor | ~$19 | ~$19 | ✔ |
| 2 | batch 1 — cold end of the T axis | ~$126 | ~$145 | ✔ |
| 3 | −10 °C | ~$321 | ~$466 | ✔ **SLIDES** |
| 4 | production packings | local | ~$466 | — |
| 5 | main matrix | measured at stage 2 | — | ✔ **SLIDES** |
| 6 | sensitivity arms | ~12 runs | — | — |
| 7 | analysis | — | — | ✔ **SLIDES** |

---

## Stage 0 — the conductivity law · DONE

`-keff_interp` defaults to `tensor` since 2026-09-23: arithmetic along the
interface, harmonic across, `n = ∇φ/|∇φ|`. Both branches satisfy `K(0)=K_air`
and `K(1)=K_ice`, so the law changes nothing outside the band and cancels the
band's first-order normal-flux bias.

| ladder | arith | tensor |
|---|---|---|
| planar slab, `studies/keff_sharp_limit/verification/` | `k_⊥` +59.4% → +3.8% | flat at the sharp value to **4e−7** |
| disk array (Rayleigh), `studies/keff_sharp_limit/disk/` | +18.19% → +1.93%, matching the no-free-constant dipole prediction to −1.9% | +1.212% → **+0.006%**, order **≈ 2.5** |

Writeup: `effective_thermal_cond/docs/tensor_conductivity_law.tex`.

**Every run before 2026-09-23 used `arith`**, the 2026-09-16 pilot included.

---

## Stage 1 — replay the pilot · GATES EVERYTHING

Neither ladder has a **contact**, and the contact is the point: the packing's
grains are seated in exact tangency, so at `t = 0` the path between two grains
is a single point and the arithmetic band bridges it with conducting material.
That inflates `k_eff(0)` and suppresses the relative rise — the headline. The
disk ladder's ~+8% at the pilot's `eps/R_ave = 0.02` is a **lower bound** on
the packing, not an estimate.

```bash
./scripts/HPC/submit_keff_replay.sh --dry-run --laws "tensor sharp" \
    --roots $SCRATCH/enceladus_DSM/batch_2026-09-16__13.16.11_pilot_keff \
            $SCRATCH/enceladus_DSM/packing_2D_pilot_phi0.325_Rave50um_LR40_seed*_L2mm_eps1000nm_perxy_T-20
venv_enceladus/bin/python studies/keff_sintering/coefficient_fix/compare_laws.py <batch>
```

Give both roots — all four seeds were resumed, and a resume leg lands in the
single-run tree, not the batch parent. ~$35 for 2 laws × 4 seeds.

**Read:** `tensor` vs `sharp` at `t_end`, per seed.

- **agree within ~5%** → adopt the tensor rise as the campaign headline and
  rescale stage 2's threshold by the same factor. Continue.
- **disagree by more than ~5%** → **stop.** Two unrelated routes to removing
  one bias must agree; a gap is a result in its own right and no production
  run is worth submitting until it is understood.

Also read, and record in this file when known:

- **the rise**, `arith +19.6% (sd 2.1) → tensor ?` — the correction factor.
- **`k_eff(0)`.** `ROADMAP.md:260` bans `k_eff(t=0)` as a baseline *because*
  arith pre-welded tangent grains. If tensor `k(0)` drops materially, that ban
  is obsolete; state here which baseline the campaign uses.
- **`ksp_its`**, against the disk ladder's ~2.5× arith. A larger rise on the
  packing is a real cost signal for every stage below.

> **Withdrawn here:** `PLAN.md:71-79`'s `R_feat/R_ave = 1/50` recommendation.
> Its case was a ~33%-vs-~15% `k_eff` eps bias, which was an *arithmetic-law*
> property and is gone. See "Necks are not resolved" below.
>
> **Also withdrawn:** `ROADMAP.md` Phase 2 / `PLAN.md` Stage 2's eps-correction
> runs were sized to measure the `O(eps)` `k_eff` bias the tensor law nulls.
> What survives is whether the *evolution* differs at smaller eps — which a
> replay cannot answer and which needs one real run, not a ladder. Deferred to
> stage 6.

---

## Necks are not resolved, by design

**The only length-scale ratio this campaign imposes is the RVE one,
`L/R_ave = 40`**, which is measured: seed scatter in `k_eff` falls from 12–13%
at `L/R_ave` 10–20 to 4.3% at 40. It must hold at **`t_final`**, not `t = 0` —
coarsening lowers it during the run.

`R_feat/R_ave` is **not** a resolution requirement. It is the cost knob: it
sets `eps = safety·R_feat`, hence `Nx`, hence the bill, and nothing else
depends on it. The campaign uses **1/25** (`eps` = 1 µm, `Nx` = 2829, 24 M DOF,
241 ranks) because that is what is affordable across ~80 runs.

At that mesh a neck is only represented above `r/R = √(12·eps/R) = 0.49`, so
**most of every trajectory's neck is below the floor.** That is accepted, not
fixed: resolving necks at this domain size is not something these simulations
can be asked for, and buying it would cost 4× the campaign for a quantity the
manuscript does not report. `k_eff`, SSA and porosity are the reported
quantities and none of them needs a resolved neck.

**Manuscript obligation.** State the floor as a measured number and its
consequence — necks are under-resolved, so no neck-growth exponent is claimed
from these runs, and the sintering *rate* inherits an error of unquantified
sign at early times. Do not present it as a caveat on `k_eff` itself; the
conductivity is a property of the field as it stands, whatever produced it.

Do not reopen this. If a later question genuinely needs neck radius, it needs
a different, smaller-domain experiment — not a finer mesh on an 80-run matrix.

---

## Stage 2 — batch 1, the cold end · ~$126

The axis the 2025-09 sweep destroyed, rebuilt with `alpha_c = 1e-3` constant
and `beta_sub0` therefore spanning 8.4× across the axis, one `eps` and one
`-eps_valid_temp` per temperature.

```bash
./scripts/HPC/submit_batch.sh --tag keff_T_batch1 \
    --tests-file studies/keff_sintering/batch1_T.txt \
    --extra-opts "-keff 1 -keff_step0 1 -keff_t_interv 5.0e4 \
                  -keff_ksp_type cg -keff_pc_type gamg"
```

−40 °C and −30 °C, seeds 1 and 2 — **two temperatures × two seeds, not three
× one**. The pilot already fixes the scatter at −20 °C; what is missing is
scatter at the cold end, which is what makes a temperature difference a result
rather than a realization.

| T | ρ_vs vs −20 | est. steps | est. wall @241 | legs | est. cost |
|---|---|---|---|---|---|
| −40 | 0.124× | ~825 | ~5.5 h | 1 | ~$16 ea |
| −30 | 0.368× | ~2 450 | ~16 h | 1 | ~$47 ea |

Note the mechanism: the pilot's realised `dt` is **median 200 s against
`dtmax` 8526 s, only 2.4% of steps at the cap**, so the step is `-dtCFL`-limited,
not `dtmax`-limited. Cold is expected to be cheaper because `dtCFL` tracks
interface velocity, which scales with ρ_vs — an *expectation* until this batch
measures it. **Record the measured s/step here; stages 3 and 5 are priced off
it, not off the table above.**

**Read:** the `k_eff` rise at −40 and −30 against the −20 ensemble, all on the
tensor scale.

- **cold-end rise separates from −20 by more than the seed scatter** →
  temperature is a real lever. Continue to stage 3.
- **it does not** → 30 days is too short at −40 °C for the effect to develop
  (ρ_vs is 8× lower, so the same wall-clock buys 8× less sintering). Do **not**
  buy the matrix. Either extend `t_final` at the cold end, or reformulate the
  axis at matched sintering state, and re-read this checkpoint.

---

## Stage 3 — −10 °C · ~$321 · **SLIDES 1**

Deferred out of batch 1 because at 2.5× the −20 °C vapour density it is
~16 700 steps, ~111 h at 241 ranks, and **five 24 h restart legs** — more than
the rest of the axis combined. Use `scripts/HPC/resume_batch.sh` between legs;
`merge_restart_legs.py` → `thin_snapshots.py` → `health_check.py` after, in
that order (thinning before merging frees nothing — the snapshots are
hardlinks).

Only worth submitting if stage 2 cleared. Completes the four-point axis.

> **SLIDES 1 — "the temperature lever is real and we can resolve it"**
> `k_eff` rise vs T with seed error bars; the audit table showing why the old
> sweep could not have seen it; the conductivity-law correction from stage 1.
> The ask: sign-off on the condition matrix before stage 5.

---

## Stage 4 — production packings · local, minutes

φ ∈ {0.25, 0.30, 0.35, 0.40} × 5 seeds, `L/R_ave` ≥ 40 at **`t_final`**, not
`t = 0` — coarsening lowers it during the run. Use `--max-void-per-L`, never
the mean-radii gate, when domain size varies. Porosity ceiling is **0.40**: the
solid stops percolating between 0.40 and 0.45, above which `k_eff` is not the
conductivity of a connected medium. Gate: solid percolates both axes on every
packing.

---

## Stage 5 — the main matrix · **SLIDES 2**

φ {0.25, 0.30, 0.35, 0.40} × T {−40, −30, −20, −10} × 5 seeds, minus what
stages 2–3 already cover. **Price it from stage 2's measured s/step**, not from
this file's estimates.

**Read:** between-condition difference against within-condition SEM. If the
matrix is unaffordable, cut temperatures or porosities — **never seeds**. The
seed count is what makes a trend a result rather than an anecdote, and the SEM
is already the limit on what can be claimed.

---

## Stage 6 — sensitivity arms · ~12 runs

One at a time from the stage 5 centre point, 3 seeds each:

- **`alpha_c`** {1e-4, 1e-3, 1e-2} — the one genuine free parameter, and it
  moves `L*` across the grain scale, so it may change the *mechanism*.
- **`R_ave`** {×0.5, ×1, ×2} **at fixed `L/R_ave`** — the packings are then
  geometrically identical and only the physics changes.
- **`σ_ln`** {0.2, 0.5} — how much the suppressed-ripening 2D artifact matters.
- **`eps` ×2 on one condition** — what survives of the withdrawn stage-1
  eps work: whether the *evolution*, not `k_eff`, differs at smaller eps.

---

## Stage 7 — analysis · **SLIDES 3**

1. **Make the SSA claim falsifiable.** Regress `k_eff` on porosity alone, then
   porosity + SSA, and report the **partial** correlation — pooled *and within
   each condition*. A relation that appears only when porosity varies is weaker
   evidence than one that holds at fixed porosity.
2. **Anisotropy from the tensor eigenvalues**, not `k_yy/k_xx`: `|k01|/k00`
   reaches 0.081, so the principal axes are not the grid axes. It is flat in
   time and seed-dependent — inherited from the packing, not made by sintering.
3. **Per-snapshot metrics**: Euler characteristic, chord-length distributions
   (what Calonne/Torquato models are built on, so the sharpest test of 1),
   `D_eff`. **Not neck radius** — it is below its resolution floor for most of
   every trajectory at this mesh; see "Necks are not resolved".
4. **Lever ranking** with seed-scatter error bars — the actual deliverable.

The strongest result already in hand: `k_eff` rises while the ice fraction does
not move at all (`phi_bar` constant to seven figures in every pilot run). Mass
is conserved in a closed box, so density is constant *by construction* and a
density-only parameterisation predicts a flat line. It is not that density
fails to explain the rise — density **cannot** be the explanation.

---

## Standing rules

- Runs are launched by Jackson, never by the assistant.
- 2+ HPC jobs go through `submit_batch.sh`; chained single submits race in
  `obj/` and give "Stale file handle".
- Push before submitting, so the cluster checkout matches.
- `merge` → `thin` → `health_check` on every batch, in that order.
- Production runs are 2D. The tensor law is verified in 2D; 3D is untested.

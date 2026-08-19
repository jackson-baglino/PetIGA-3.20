# Molaro 2019 — run ladder

**One simulation at a time.** Each rung answers one question, and its answer
changes what the next rung should be. A 5-arm humidity grid submitted blind
spends ~5x the core-hours to learn what two adaptive runs learn, because
shrinkage is very nearly linear in the wall undersaturation — so a measured
result gives an exact Newton step to the next humidity rather than a bracket.

Costs below assume ~220 core-h per 6 h dom3 arm; check the first run before
trusting them.

---

## Rung 0 — IC check (10 steps, ~1 core-h)

```bash
./scripts/HPC/submit_enceladus.sh \
    molaro_2D_L450x225um_eps0.26um_axisym_T-20pair_tangent_dom2 \
    molaro_T-20_h1.0000_2h_a1e-1_dirichlet icheck --time=0-01:00:00 \
    -- -ts_max_steps 10 -outp 1 -t_out_log 0
```

`-t_out_log 0` is load-bearing — the experiment file sets `-t_out_log 80`,
which would otherwise swallow most of the 10 steps.

On the cluster: `python postprocess/plot_fields.py --dir <run> --steps 0 1 10`,
then rsync those three `.vts` plus `pf.pvd`.

**Check:** the banner prints `-ic_grain_union 1` (the *default* is the additive
form, which sums two tanh tails into a spurious bridge at exact tangency);
φ = 0.5 is tangent at exactly one point; the intermediate-φ band at the cusp is
~4.4ε ≈ 1.14 µm axially (correct diffuse image of a cusp, not a defect); by
step 10 the O(ε) bridge has filled without spurious structure spreading along
the contact plane.

**Gate:** if the IC is wrong, nothing below is worth running.

---

## Rung 1 — pre-necked start at r = 14 µm (dom2 first)

**The CFL ladder is done and it split into two findings.**

*Fixed:* the φ > 1 overshoot is a timestep artifact. Peak excursion vs
`-dtCFL_dphimax`, four dom2 arms through the merge:

| `dphimax` | peak overshoot | ρ_v range / ρ_vs |
|---|---|---|
| 0.20 | +8.0e-03 | [−6.69, +8.17] |
| 0.05 | **+6.7e-16** | [+0.99, +1.00] |
| 0.02 | +2.2e-16 | [+0.99, +1.00] |
| 0.01 | +2.2e-16 | [+0.99, +1.00] |

`dphimax = 0.05` is now the default. 0.2 let the front cross 1.13 elements per
step; 0.05 is 0.28.

*Not fixed by dt:* a **trapped void**. Present in all four arms, forming at
t ≈ 150–165 s at r ≈ 5 µm, sealed over by ice, and persisting and deepening
(0.69 → 0.55 over 130 s). A 20× tightening moved the trapped minimum only
0.498 → 0.551 — it is not a temporal artifact.

The cause is spatial. At r = 5 µm the fillet radius is
`ρ = r²/(2(R−r)) = 0.153 µm = 0.59 ε` — **thinner than the interface itself**.
The whole tangent-start traverse (r = 0 → 16.4 µm) sits below the model's own
resolution floor. Fixing that by refinement needs ε = 2.5e-8, ~1044× the cost.

**So: start above the floor instead.** ε = 1.18e-7 and a pre-necked r = 14 µm,
slightly under Molaro's first measured 16.4 µm so the comparison window is still
derived rather than assumed. Two criteria independently pick that ε:

```
delta_ice < 1                     -> eps <= 1.1802e-07   <- binds
start above the resolution floor  -> eps <= 1.8828e-07
```

At ε = 1.18e-7: floor = 11.08 µm (start is 1.26× above, their first neck 1.48×),
`δ_ice = 1.000`, `δ_mean = 0.158`, `δ_vap = 0.067`, `τ_sub = 1.33 s`.

| arm | Lz × Lr (µm) | Nz × Nr | Melem | M DoF | ranks @100k |
|---|---|---|---|---|---|
| `..._r14um_dom2` | 450 × 225 | 5394 × 2697 | 14.55 | 43.7 | 437 |
| `..._r14um_dom3` | 680 × 340 | 8150 × 4075 | 33.21 | 99.7 | 998 |

```bash
./scripts/HPC/submit_enceladus.sh \
    molaro_2D_L450x225um_eps0.12um_axisym_T-20pair_r14um_dom2 \
    molaro_T-20_h1.0000_8h_a1e-1_dirichlet r14dom2
```

**Check, in order:** (1) no trapped void — `min φ` inside the neck stays at 1
and no bracketed dip appears; (2) `*** PHASE GUARD TRIPPED ***` never prints;
(3) the neck passes 16.4 µm so the comparison window is reached; (4) grab the
SLURM `.o` for wall time — it is now copied into the run folder automatically,
and the ~10× cost estimate is still unvalidated.

dom3 only after dom2 is clean — 998 ranks is 32 nodes.

## Rung 2 — domain convergence (~2 runs, ~600 core-h)

Repeat rung 1 at `dom3`, then at `dom4` **only if dom3 disagrees with dom2 by
more than the ~9 % the (a/L)³ estimates predict between them**.

| arm | L/a_eff | Lz × Lr (µm) | Melem | ranks @100k | predicted neck err |
|---|---|---|---|---|---|
| dom2 | 2.01 | 450 × 225 | 3.03 | 91 | 12.5 % |
| dom3 | 3.03 | 680 × 340 | 6.92 | 208 | 3.7 % |
| dom4 | 4.01 | 900 × 450 | 12.12 | 364 | 1.6 % |

**Decide:** the production domain is the smallest arm whose neck curve is
converged. Everything below runs on it.

---

## Rung 3 — humidity, by iteration (~2 runs, ~440 core-h)

Start at the derived estimate:

```bash
./scripts/HPC/submit_enceladus.sh <production-geom> \
    molaro_T-20_h0.9973_6h_a1e-1_dirichlet h9973
```

Then correct. Grain recession is linear in the wall undersaturation, so the
measured `dR_large_pct` gives an exact Newton step:

```
(1 - h)_next  =  (1 - h)_now  *  (-2.93) / dR_large_pct_measured
```

**Target: −2.93 %** (the large grain's least-squares slope over 78 min, which is
the −3 % the paper's caption quotes). Do **not** fit the small grain — its total
change is inside its own ±2 % error bar.

Pre-built rungs, if the correction lands near one: `h0.9989` (0.4×), `h0.9946`
(2.0×), `h0.9913` (3.2×). Otherwise copy the nearest and edit `-humidity`.

**Also check:** the `h1.0000` control from rung 1/2 should have essentially the
same *neck* curve as the fitted arm. Neck growth is driven by an internal
curvature difference and shrinkage by the external wall value, so the two
should be nearly decoupled. If they are not, the domain is still too small and
rung 2 concluded early.

---

## What "done" looks like

From `summary.csv`, on the production domain at the fitted humidity:

| quantity | target | source |
|---|---|---|
| `neck_w_at_78min_um` | 64.78 | their last measured neck |
| `t_star_s` | — | the clock shift; report it, do not fit it |
| `dR_large_pct` | −2.93 % | large-grain fit; caption says −3 % |
| `dR_small_pct` | report only | inside their own error bar |

The model is expected to reach **~50 % of the observed neck growth** and no
more, at any α_c — see `alpha_c_estimate.csv` and
`docs/molaro_validation_synthesis.md` §4. Landing there is the success
criterion; landing at 100 % would mean something is wrong.

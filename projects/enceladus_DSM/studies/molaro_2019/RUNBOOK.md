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

## Rung 1 — dom2, saturated (~2 h wall, ~70 core-h)

```bash
./scripts/HPC/submit_enceladus.sh \
    molaro_2D_L450x225um_eps0.26um_axisym_T-20pair_tangent_dom2 \
    molaro_T-20_h1.0000_6h_a1e-1_dirichlet dom2
```

The cheapest arm that produces a real neck curve. Two things to read:

1. **`plots/timestep.png` — is `dtCFL` binding?** `-dtCFL` is now *enforced*:
   a step whose measured `||dphi||_inf` exceeds `-dtCFL_dphimax` is rolled back
   and retried, so `max |dphi| <= 0.2` holds on every step the run keeps. That
   is why `-dtmax` is loose (1e4). Confirm the limiter is what binds:
   - **`grep -c 'Interface-CFL cap' <run>/outp.txt` is large** → the limiter is
     governing dt. Expected.
   - **dt pinned at 1e4 with no CFL lines** → the limiter is not firing; put
     `-dtmax 2.0e2` back before spending anything else.
   - **`grep -c 'Interface-CFL violated'` is more than a few percent of steps**
     → rollbacks are thrashing and each one is a wasted step. Lower
     `CFL_FORWARD_MARGIN` in `src/monitoring.c` (currently 0.9).
   - Watch the SNES counts: climbing toward `-NRmax 15`, or visible rejections,
     means dt is too large whatever the CFL limiter says.

   This single plot sets the cost of everything after it.

2. **Does the neck reach 64.78 µm inside 6 h?** Predicted t\* ≈ 0.42 h and
   t_end ≈ 3.16 h. If it overshoots, raise `-t_final` on the next rung; if it
   finishes in 1 h, lower it.

```bash
bash postprocess/run_batch_measure.sh <run_dir>      # works on a single run too
```

---

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

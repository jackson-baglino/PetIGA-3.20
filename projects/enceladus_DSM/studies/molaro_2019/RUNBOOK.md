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

## Rung 1 — CFL ladder through the neck merge (~4 short runs, ~40 core-h)

job1062679 (dom3, `-dtCFL_dphimax 0.2`) died at t = 457 s. Cause: the neck
**merged** at t = 126 → 142 s (steps 45 → 46), φ on the axis at the neck plane
went 0.9659 → 1.00200, and a ring of φ > 1 appeared at r = 1–7 µm and never
relaxed — peaking at **1.0143**. Once an accepted state passed the guard band
every residual evaluation domain-errored and no dt could recover, so it burned
`-max_rej 50` and exited `DIVERGED_STEP_REJECTED`.

φ outside [0,1] by more than ~0.01 means the numerics are wrong, so the fix is
to make the excursion not happen — not to tolerate it. The prime suspect is
front motion per step:

| `dphimax` | front move/step | elements | in units of W | dt (from the 42.2 s observed at 0.2) |
|---|---|---|---|---|
| 0.20 | 0.2068 µm | **1.131** | 0.566 | 42.2 s |
| 0.05 | 0.0517 µm | 0.283 | 0.141 | 10.6 s |
| 0.02 | 0.0207 µm | 0.113 | 0.057 | 4.2 s |
| 0.01 | 0.0103 µm | 0.057 | 0.028 | 2.1 s |

A front crossing **more than one element per step** on a C¹ quadratic basis is
where Gibbs overshoot starts, and 0.2 permits exactly that. Default is now 0.05;
this ladder measures what is actually required.

Run on **dom2** (cheap) with the merge-window experiment — `t_final = 300 s`
covers the merge with margin, so there is no reason to pay for 6 h to study it:

```bash
for d in 0.20 0.05 0.02 0.01; do
  ./scripts/HPC/submit_enceladus.sh \
      molaro_2D_L450x225um_eps0.26um_axisym_T-20pair_tangent_dom2 \
      molaro_T-20_h1.0000_merge_a1e-1_dirichlet cfl${d/./} \
      -- -dtCFL_dphimax $d
done
```

**Measure:** peak `max|φ − clamp(φ,0,1)|` through the merge. Healthy is near
machine precision; job1062679 reached 1.4e-2. The run now reports it itself —
any guard trip prints `*** PHASE GUARD TRIPPED ***` with the worst φ, so a
failure names its own cause instead of leaving a bare PETSc trace.

**If tightening `dphimax` alone does not drive the excursion down**, the next
suspect is spatial, not temporal: `h = ε/√2` is `W/2`, at the loose end of the
Karma–Rappel `dx/W = 0.4–0.8` band. Re-run the best `dphimax` at `h = ε/2`
(`Nx`, `Ny` ×1.41) before concluding anything about the model.

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

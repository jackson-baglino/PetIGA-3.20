# Contact-angle sweep — batch 2026-09-12 08:01

11 runs, all reached `t_final = 5.616e6 s` (65 days) and reported
`Solution completed`.

## Verdict

**The angle sweep passes decisively. The ε-convergence point at 1/100 is
invalid (solver stall), and the sessile runs are not converged — neither is a
physics failure, but neither is usable as it stands.**

## 1. Angle sweep — ✅ PASS

Channel, ε = 1.50 µm (ε/R = 1/50). θ_inf from a free three-parameter fit
`θ(t) = θ_inf + A·exp(−t/τ)`; the asymptote is fitted, never assumed.

| θ_Young | θ(t_end) | θ_inf | error | τ [d] | spread over 4 estimates |
|---|---|---|---|---|---|
| 30° | 32.685 | 30.079 ± 0.015 | **+0.079°** | 18.0 | 0.0000 |
| 60° | 60.970 | 60.035 ± 0.003 | **+0.035°** | 14.6 | 0.0000 |
| 90° | 90.500 | 90.007 ± 0.001 | **+0.007°** | 14.0 | 0.0000 |
| 120° | 120.223 | 120.001 ± 0.000 | **+0.001°** | 15.5 | 0.0000 |
| 150° | 149.128 | 150.020 ± 0.001 | **+0.020°** | 20.2 | 0.0000 |

**Max error 0.079°, RMS 0.040°**, against a 2° acceptance criterion — better
than the criterion by a factor of 25. The four independent estimates (two
menisci × two walls) agree to below the printed precision.

The extrapolation is trustworthy here, which matters because none of these runs
reached equilibrium directly. θ_inf is stable across fit windows:

| case | >30% | >40% | >50% | >60% | drift |
|---|---|---|---|---|---|
| θ=30 | 30.28 | 30.08 | 29.99 | 29.96 | 0.32 |
| θ=60 | 60.07 | 60.04 | 60.02 | 60.01 | 0.06 |
| θ=90 | 90.01 | 90.01 | 90.01 | 90.00 | 0.01 |
| θ=120 | 120.00 | 120.00 | 120.00 | 120.00 | 0.01 |
| θ=150 | 150.03 | 150.02 | 150.01 | 150.01 | 0.02 |

Errors are all *positive* and largest at the extremes (30°, 150°), which is the
signature of a small residual O(ε/R) diffuse-interface bias rather than noise.

## 2. ε convergence — ⚠️ one point valid, one invalid

| ε/R | error | status |
|---|---|---|
| 1/25 | +0.067° | ok |
| 1/50 | +0.035° | ok |
| 1/100 | +45.239° | **INVALID — solver stalled** |

1/25 → 1/50 halves the error, consistent with an O(ε/R) bias. The 1/100 point
is not a measurement at all:

- **Only 9 distinct snapshots out of 61.** The last 52 `sol_*.dat` files are
  byte-identical; θ is bit-identical at 105.239° from day 10 to day 65.
- A phase-field bounds violation at steps 383–385 cut dt from 2.0e3 → 5.9e-1,
  followed by 23 `DIVERGED_LINE_SEARCH iterations 0` — the line search failing
  immediately, i.e. the state cannot move at all.
- The clock kept advancing to 65 days at dtmax the whole time, so the run
  *reported success*. Nothing in the solver's exit status flags this.

Note `tot_ice` in `SSA_evo.dat` is constant to all printed digits in **every**
run including healthy ones — the vapour holds ~1e-6 of the total water, so ice
volume cannot change within `%e` precision. It is not a freeze diagnostic. The
snapshot-fingerprint count is.

### Cause, confirmed — and it was the wall term, not dtmax

The first diagnosis (dtmax too coarse) was wrong. A rerun at
`dtmax/tau_sub = 0.10` stalled again, at step 2653 of 25584.

**The real cause is a sign-flip instability in the wall term.**
`h'(phi) = 6*phi*(1-phi)` changes sign for `phi < 0`, so the wall residual
`-3*M*cos(theta)*h'/6*N` flips with it and drives phi *further* out of range.
And nothing opposes it: the bulk `f1` and `loc` are evaluated at the **clamped**
`phi_c`, so for `phi < 0` both are identically zero — there is no restoring
force outside [0,1] by construction. An unclamped wall term was the only thing
acting out there.

The evidence is a clean boundary layer, not a bulk effect:

| distance from wall | 0 | 1.4 ε | 2.1 ε | 4.2 ε | 14 ε | interior (rows 20–170) |
|---|---|---|---|---|---|---|
| φ_min | −0.0521 | −0.0294 | −0.0204 | −0.0069 | −4.7e-5 | **−4.7e-5** |

Control points below −1e-3 go from **0 at step 1779 to 4245 at step 2211** — an
exponential runaway from a −2e-6 seed, exactly as a positive feedback predicts.

Once phi crossed `-phase_lo`, the residual's domain guard zeroed the residual,
SNES read ‖F‖ = 0 as converged at iteration 0, and the solution froze for 23,370
steps while the clock ran to `t_final` and the run **reported success**.

**Fixed** by clamping phi to [0,1] before evaluating `h'` (the Jacobian clamps
identically and carries the chain-rule factor, so it stays exact; all 24 gates
still pass), plus a stall detector that aborts when ‖U‖ is bit-identical for 500
consecutive steps.

Ruled out: `-thin_iface_corr` was off, but enabling it moves `tau_sub` by only
2.9% at this ε. A sign-flip instability is not a rate problem.

### Mesh and dtmax, checked anyway

The mesh was **not** at fault. All four geometry files match `comp_eps.py`
exactly (eps, Nx and Ny), and `dx/eps = 0.707 = 1/√2` independently confirms the
mesh rule. ε/R = 1/100 was correctly sized.

The timestep was. `tau_sub` scales as **ε²** — measured across this family:

| ε [µm] | τ_sub [s] | dtmax | dtmax/τ_sub | outcome |
|---|---|---|---|---|
| 3.00 | 3.515e4 | 2.0e3 | 0.057 | clean |
| 1.50 | 8.788e3 | 2.0e3 | 0.228 | clean |
| 0.75 | 2.197e3 | 2.0e3 | **0.910** | stalled |

A single step spanned 91% of the entire interface relaxation time. The
structural cause was that `-dtmax` lived in the **experiment** file, which is
shared across all three resolutions, while ε lives in the **geometry** file —
and options are applied geometry-then-experiment, so the experiment's value
overrode everything.

**Fixed two ways.** `-dtmax` moved into the geometry files at `τ_sub/10`
(2.0e3 / 8.79e2 / 2.20e2 for ε = 3.00 / 1.50 / 0.75 µm), never looser than a
value already proven for that geometry; and the solver now prints
`dtmax/tau_sub` at startup and warns loudly above 0.25.

## 3. Sessile drop — ⚠️ not converged, do not interpret

| θ_Young | θ(t_end) | θ_inf (fitted) | still moving at t_end |
|---|---|---|---|
| 30° | 49.76 | 43.2 ± 1.8 | −0.127 °/d |
| 60° | 68.92 | 66.2 ± 0.9 | −0.139 °/d |
| 120° | 110.75 | 115.6 ± 2.0 | +0.180 °/d |
| 150° | 127.97 | 133.5 ± 2.3 | +0.241 °/d |

Every drop is moving **toward** its Young target and none has arrived. The
drift rates are 2–10× the channel's, and τ = 30–45 days against a 65-day run —
only 1.4–2.2 time constants, where the channel got 3.2–4.6.

The fitted θ_inf is **not** a measurement here: it drifts by 2.8–9.2° depending
on the fit window (channel: 0.01–0.32°), so the single-exponential model does
not yet hold. Reporting −16.5° of "error" at θ=150 would be reading noise.

A sessile drop has no second wall to pin it and must transport mass further to
spread, so it relaxes more slowly than a confined bridge. This needs
`t_final ≈ 1.7e7 s` (≈200 days, 4.5τ), not 5.6e6.

## 4. Mass conservation — ✅

Sealed-box total water is constant to printed precision in every channel run at
ε = 1.5 µm, and drifts by ≤1.4e-4 % in the sessile runs.

## What to rerun

```bash
# eps/R = 1/100 -- dtmax now comes from the geometry file, no override needed
./scripts/HPC/submit_lunar.sh channel_2D_H100um_eps0.75um relax_T-20_theta60 eps100_retry

# sessile, long enough to actually converge
#   raise -t_final to 1.7e7 in inputs/experiment/contactangle/ first
```

Check any rerun with the fingerprint test before trusting it:

```bash
ls <run>/sol_*.dat | while read f; do shasum -a1 "$f"; done | awk '{print $1}' | sort -u | wc -l
```

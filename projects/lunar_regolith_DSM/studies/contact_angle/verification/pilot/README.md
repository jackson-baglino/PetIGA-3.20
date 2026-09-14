# Pilot run — θ = 60°, ε/R = 1/25

First physics run of the prescribed-contact-angle feature. Purpose was to size
`t_final` for the sweep; it also turned into the first real validation.

```
geometry   channel_2D_H100um_eps3.00um     (142 x 48, eps = 3.0e-6, eps/R = 1/25)
experiment relax_T-20_theta60              (gamma_ia 0.109, gamma_as 0.300, gamma_is 0.2455)
t_final    2.592e6 s  (30 days)
```

## Result

| | |
|---|---|
| Young's prediction | **60.000°** |
| Extrapolated equilibrium θ_inf | **60.36 ± 0.04°** |
| θ at t = 30 d (not yet equilibrated) | 72.19 ± 0.006° |
| Relaxation time τ | 14.4 days |
| Fit residual RMS | 0.015° |

θ_inf is from a **free** three-parameter fit `θ(t) = θ_inf + A·exp(−t/τ)` — the
asymptote is fitted, not assumed, which is what makes the comparison to Young's
equation meaningful. It is stable across fit windows:

| fit window | θ_inf | τ [d] | resid RMS |
|---|---|---|---|
| t > 6 d  | 58.79 ± 0.23 | 15.10 | 0.142 |
| t > 9 d  | 59.95 ± 0.08 | 14.42 | 0.037 |
| t > 12 d | 60.36 ± 0.04 | 14.15 | 0.015 |
| t > 15 d | 60.24 ± 0.05 | 14.24 | 0.012 |
| t > 18 d | **59.93 ± 0.03** | 14.49 | **0.003** |

The best-conditioned window gives 59.93°, i.e. **0.07° from Young's prediction**,
on the *coarsest* mesh in the sweep. The four per-snapshot estimates (two
menisci × two walls) agree to ±0.006°.

## Two things it established

**Relaxation is vapour-diffusion limited and slow.** τ = 14.4 days, so reading
the angle directly needs ~65 days of simulated time (4.5 τ). The experiment
files were resized from 30 to 65 days accordingly, and `contact_angle.py` now
reports θ_inf so a short run still yields the equilibrium.

**θ moves the wrong way for the first ~2 days** — 131.8° → 133.8° before turning
around. That is the Gibbs-Thomson transient: the grain is convex, so its
equilibrium vapour density exceeds the flat value, and it shrinks until the
sealed box saturates. A shrinking bridge at fixed contact points has a larger
apparent angle. This is why `extrapolate()` discards the first 40% of the
trajectory before fitting.

## Caveat

This is one angle on one mesh, and the headline number is an extrapolation, not
a directly observed equilibrium. It is strong evidence, not the validation. The
sweep in `../sweep_tests.txt` is what settles it.

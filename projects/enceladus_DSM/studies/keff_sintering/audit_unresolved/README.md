# Are the 2025-09 `unresolved_results` runs publishable?

**Short answer: the porosity sweep is a usable preliminary result. The
temperature sweep is not a result at all — its headline finding is an artefact
of the parameter choice, and it should not be shown to anyone.**

Nine completed 28-day runs live in
`~/SimulationResults/HPC_results/enceladus_DSM/unresolved_results/`
(2 mm × 2 mm, 441 grains, `R_ave = 45 µm`, humidity 0.98, seed 21):

- `2mm_results2` — temperature sweep, `T = −20/−25/−30/−35/−40 °C`, `φ = 0.24`
- `2mm_results3` — porosity sweep, `φ = 0.24/0.26/0.28/0.30`, `T = −25 °C`

`k_eff` tables for all nine already exist, from the old standalone
`effective_thermal_cond` code, under
`~/SimulationResults/effective_thermal_cond/`.

Reproduce with:

```
venv_enceladus/bin/python studies/keff_sintering/audit_unresolved/audit.py
venv_enceladus/bin/python studies/keff_sintering/audit_unresolved/plot_audit.py
```

Output is committed here as `audit_output.txt` and `audit.png`.

## Why this had to be reconstructed rather than read off

The code that produced these runs printed **no kinetics banner** — no `alpha_c`,
no `tau_sub`, no ε bound, nothing. `outp.txt` records only `humidity`. Every
kinetic number below was recovered from the `src/` copy each run directory
stages, which is the only reason the audit was possible at all. The current
solver prints all of it; that change was worth making.

## 1. The temperature sweep is confounded with `alpha_c`

`dry_snow_metamorphism.c:49` hardcodes `beta_sub0 = 1.4e5`, identical at every
temperature. The solver then forms `beta_sub = beta_sub0 / (rho_i/rho_vs(T))`.
**A constant `beta_sub0` is therefore not a constant condensation coefficient.**
The implied `alpha_c` is whatever makes `comp_eps.py` return `beta0 = 1.4e5`:

| T [°C] | `rho_vs` | implied `alpha_c` | `alpha_c · rho_vs` | `eps/eps_max` | `Nx` required | too coarse |
|---|---|---|---|---|---|---|
| −20 | 8.478e-4 | 5.67e-2 | 4.809e-5 | 1.68 |  5 427 |  4.8× |
| −25 | 5.197e-4 | 9.35e-2 | 4.859e-5 | 2.74 |  8 857 |  7.8× |
| −30 | 3.122e-4 | 1.57e-1 | 4.910e-5 | 4.57 | 14 747 | 12.9× |
| −35 | 1.835e-4 | 2.70e-1 | 4.962e-5 | 7.77 | 25 085 | 22.0× |
| −40 | 1.055e-4 | 4.76e-1 | 5.015e-5 | 13.51 | 43 651 | 38.2× |

`alpha_c` rises **8.4×** exactly as `rho_vs` falls **8.0×**. Their product — which
is what sets the sublimation flux — is **constant to 4.2% across the whole 20 °C
sweep**. The sweep cannot resolve a temperature effect because the temperature
effect was cancelled in the parameter choice.

And that is precisely what the `k_eff` tables report: the 28-day rise spans
**+47.78% (−40 °C) to +49.46% (−20 °C)** — 1.7 points, or 3.4% of the mean, over
20 °C. That residual is the same size as the 4.2% drift left in the cancelled
product and as the growing ε violation, and it runs opposite in sign to the
first. It is numerical, not physical.

Separately, `alpha_c` of 0.057–0.48 is **57× to 475× above the top of the
literature band** (Libbrecht 2017, Braun 2024: 1e-4 … 1e-3), and 0.48 is
approaching the kinetic ceiling of unity.

This failure mode is already guarded in the current pipeline —
`generate_study_opts.py:13-15` explains it, and `-eps_valid_temp` makes the
solver **abort** if `-temp` differs from the temperature ε was computed at by
more than 1 °C. These runs predate that guard and do exactly what it forbids.

## 2. Mesh and time step

One ε and one mesh were used for all nine runs:

- `eps = 8.757e-7 m`, `Nx = Ny = 1142`, and **`p = 1, C = 0` — linear elements**
  (confirmed against `sol_00000.dat` = 31 354 776 + 8 B, which matches
  `3 · 1143² · 8` exactly and rules out `p=2, C=1`).
- `h/eps = 2.00` against the K&P mesh rule `h = eps/√2 = 0.707` → **2.83× too
  coarse**, and that rule was calibrated for the `p=2, C=1` basis the current
  code uses, so 2.83× understates the deficit.
- The 1%–99% band `9.2·eps` = 8.06 µm spans **4.6 elements**, against the project
  standard of 7.5–10.

`dtmax` was `0.5 · t_interv` — half the *output* interval, with no reference to
the interface physics. Realised median `dt = 1334 s` is **6.4× `tau_sub` at
−20 °C** and 1.5× at −40 °C, against the lunar-calibrated safe ratio of **1.09**.
This is the same class of error as the `dtmax` bug fixed in the current pilot,
and it too is worst at the warm end of the sweep.

## 3. What is sound

- **Mass.** Ice loss is 0.007%–0.062%, monotone in T, tracking `rho_vs`. That is
  the physical response to humidity 0.98, not a conservation leak.
- **Completion.** All nine reached exactly 28.00 days with 100 snapshots.
- **REV.** `L/R_ave = 44.4` at `t = 0`, above the measured floor of 40. (It falls
  during the run as grains coarsen; not checked at `t_final` here.)
- **The magnitude of the rise is corroborated.** The old standalone code gives
  ~+50%; the new in-line `k_eff` with a sharp coefficient gives +57% on a
  different packing. Two independent codes agreeing to that tolerance says the
  ~+50% is not a coding error in either — though both inherit the same
  microstructure bias.

## 4. Verdict

| | usable? | why |
|---|---|---|
| Temperature sweep | **No** | The swept variable is confounded with `alpha_c` (8.4×) and with `eps/eps_max` (1.7→13.5×). Its null result is forced by construction. |
| Porosity sweep | **Preliminary** | Single T, single ε, single `beta_sub0` — every bias is a *common offset*. The monotone `k_eff(φ)` trend (0.894 → 0.621 at t=0) survives. |
| Absolute `k_eff` | **No** | 4.6 elements across the band, plus the arithmetic-vs-sharp coefficient question (17.8% vs 57%). |
| Sintering *rate* | **No** | `alpha_c` is 57–475× too large and `dt/tau_sub` is 1.5–6.4× over. |

**Recommendation:** treat these as a pilot that established the pipeline and the
rough magnitude, cite nothing from them quantitatively, and regenerate the
temperature sweep through `generate_study_opts.py` so each temperature gets its
own ε, `beta_sub0`, mesh, and `dtmax`. The porosity sweep may be reused as a
qualitative "trend is monotone and steep" statement while the resolved version
runs.

## 5. Cost note

The reason these were cheap is the reason they are wrong: `Nx = 1142` with
linear elements and `dt` 6× over the CFL limit. The resolved −40 °C case wants
`Nx ≈ 43 651`, which is **1 460× the degrees of freedom** — not a run that can be
done at fixed grain size. The way out is the `L/R_ave` scaling already adopted
for the new campaign (work in `L/R_ave`, keep `R_feat/R_ave` fixed), not a finer
mesh at 2 mm.

# k_eff sharp-limit verification

Measures the finite-`eps` bias of the phase-field homogenization against a
closed form, on the periodic ice-slab (laminate) cell.

The theory being tested is in
[`projects/effective_thermal_cond/docs/calonne_to_phasefield_equivalence.tex`](../../../../effective_thermal_cond/docs/calonne_to_phasefield_equivalence.tex),
sections 6 and 7.

## What question this answers

The note shows that our phase-field cell problem and Calonne et al.'s
sharp-interface cell problem are the **same** boundary-value problem — not in a
limit, but identically. The reason is that the sharp coefficient and the
phase-field interpolation are one function evaluated at two arguments:

```
K_star = K_i·1_{Ω_i} + K_a·(1 − 1_{Ω_i}) = K_a + (K_i − K_a)·1_{Ω_i} = K(φ)|_{φ = 1_{Ω_i}}
```

because the indicator `1_{Ω_i}` only ever takes the values 0 and 1. That part is
algebra and needs no test.

What *does* need a test is the residual bias at the finite `eps` we actually
run, where `φ^eps ≠ 1_{Ω_i}`. The note derives it in closed form. **This study
measures it.**

## The prediction

To first order the diffuse band acts as a sharp interface carrying two Gibbs
surface excesses per unit area — one tangential (a conductance), one normal (a
resistance):

| | definition | value | consequence |
|---|---|---|---|
| `Σ_t` | `∫[K(φ^eps) − K_star] ds` | **exactly 0** | `k_∥` is `eps`-exact |
| `Σ_n` | `∫[1/K(φ^eps) − 1/K_star] ds` | `−eps·C`, `C > 0` | `k_⊥` biased **high**, first order |

`Σ_t` vanishes because the logistic profile is antisymmetric, `σ(−u) = 1 − σ(u)`:
the arithmetic law moves no ice across the interface, so the cell's ice fraction
is exact at every `eps`. `Σ_n` does not vanish because `1/K` is convex in `φ` —
the band short-circuits what should be a series resistance. The constant is

```
C = (1/K_a − 1/K_i)·ln(K_i/K_a) = 234.959 m K W⁻¹   at K_i/K_a = 2.29/0.02
```

For the laminate this closes completely:

```
k_∥(eps)    = K_a + (K_i − K_a)·φ̄                       exact at every eps
⟨1/K⟩(eps)  = [φ̄/K_i + (1−φ̄)/K_a] − (n_Γ·eps/L)·C
k_⊥(eps)    = 1/⟨1/K⟩(eps)
```

`n_Γ = 2`: under `-periodic 1` the `y = 0` seam is a real second interface.

| eps | ⟨1/K⟩ | k_⊥ | bias |
|---|---|---|---|
| L/50 | 15.8200 | 0.063211 | +59.4% |
| L/64 | 17.8759 | 0.055941 | +41.1% |
| L/128 | 21.5471 | 0.046410 | +17.0% |
| L/256 | 23.3827 | 0.042767 | +7.9% |
| L/512 | 24.3005 | 0.041151 | +3.8% |
| sharp | 25.2183 | 0.039654 | — |

This is not merely a first-order truncation. The only approximation is the
overlap of neighbouring interfaces' exponential tails, `O(exp(−L/4·eps))`;
against direct quadrature the formula is good to 5e−5 relative at `eps = L/50`
and to machine precision by `eps = L/128`.

## Plot it in resistivity, not conductivity

`⟨1/K⟩` is **linear in eps**; `k_⊥` is a hyperbola. Fit and plot `1/k_⊥` and the
test has two independently predicted numbers and no fitted ones:

```
intercept = φ̄/K_i + (1−φ̄)/K_a        slope = (n_Γ/L)·C = 2 × 234.959 / L
```

Plotted as `k_⊥` the same data is a curve whose agreement can only be eyeballed,
which is how a 10% slope error passes review.

## Running it

```bash
cd /Users/jacksonbaglino/PetIGA-3.20/projects/enceladus_DSM
./studies/keff_sharp_limit/verification/verify_keff_sharp_limit.sh
```

Options: `--rungs "64 128 256"` to shorten the ladder, `--dry-run` to print the
`eps`/`Ny` schedule without running anything.

Each rung is `-keff_only`: one sample from the initial condition, `dim` scalar
Poisson solves, exit. There is **no time integration**, so even the 2560²  rung
is minutes. This runs locally; it is not an HPC job.

Outputs, all written next to this README:

- `keff_sharp_limit.csv` — one row per rung: the solver's `k_eff.csv` row plus
  `eps`, `Ny`, and the run directory it came from.
- `keff_sharp_limit.png` — the two-panel figure.
- `keff_sharp_limit.log` — raw solver output for every rung.

## The gates

Run automatically at the end; `plot_keff_sharp_limit.py` exits non-zero if any
fails, so the driver doubles as a regression check.

**Gate 0 — analytic self-check.** The module must reproduce `k_∥ = 1.155000` and
`k_⊥ = 0.0396538`, the two constants stated independently in the geometry file
[`inputs/geometry/iceslab/iceslab_2D_L1mm_eps20um_keff.opts`](../../../inputs/geometry/iceslab/iceslab_2D_L1mm_eps20um_keff.opts).
Checked *before* any run: failing it means the module is wrong, not the solver.

**Gate 1 — `k_00` flat**, at `φ̄K_i + (1−φ̄)K_a`, to 2e−3. The tangential excess
is exactly zero, so a drift here is a **real defect** in the cell solver or the
IC, not interface bias.

**Gate 2 — `1/k_11` linear in `eps`**, intercept and slope as above. This is the
test of the theory.

**Gate 3 — `k_01`, `k_10` ≈ 0.** Off-diagonal isotropy check on the cell solver.

## Two ways to get a meaningless answer

**Resolution must scale with `eps`.** The driver sets `Ny` to hold `eps/dy`
fixed at 5.12, anchored on the geometry file's reference point (`eps = 2e−5` at
`Ny = 256`). At *fixed* resolution an `eps` ladder confounds the interface bias
being measured with ordinary mesh convergence. A **bad slope with a good
intercept** is the signature — the gate prints this hint when it sees that
pattern. Check `EPS_PER_ELEM` before doubting the theory.

**Evaluate predictions at the measured `φ̄`,** not the nominal 0.5. The discrete
field is a spline fit to nodal values and its mean lands near but not on 0.5
(~4e−4 relative at `Nx = Ny = 256`). Feeding 0.5 in predicts `k_00 = 1.155000`
where the correct discrete answer is `1.155454`, and the gate flags a 4e−4
"solver error" that is really the IC's quadrature error. The scripts read `φ̄`
from the CSV column for this reason.

## Files

| | |
|---|---|
| `keff_laminate_analytic.py` | the closed forms; single source of truth for gate and plot. Run it directly for the self-check and the predicted ladder. |
| `verify_keff_sharp_limit.sh` | driver: loops `eps`, calls `run_enceladus.sh`, collects the CSV |
| `plot_keff_sharp_limit.py` | gates + figure; exits non-zero on failure |

## The tensor law (`-keff_interp tensor`)

Section 5 of the note gives the fix. Neither scalar interpolation can null both
excesses — arithmetic is exact tangentially, harmonic is exact normally. The
tensorial law

```
K(φ) = K_arith(φ)·(I − n⊗n) + K_harm(φ)·(n⊗n),    n = ∇φ/|∇φ|
```

uses each where its excess vanishes, giving `O(eps²)`. This is legitimate
because the conductivity cell problem is pure post-processing on a frozen `φ`
field, decoupled from the phase evolution — the interpolation used for `k_eff`
need not be the one used in the evolution equations.

It is implemented in `KeffPointCond`
([`src/keff_cell.c`](../../../src/keff_cell.c)) and selected with
`-keff_interp tensor` (default `arith`; `sharp` evaluates `K(H(φ−½))` on the
same mesh as a cross-check). On this laminate the prediction becomes a **flat**
`1/k_11` at the sharp intercept, because `1/K_harm` is affine in `φ`:

```bash
./studies/keff_sharp_limit/verification/verify_keff_sharp_limit.sh --interp tensor
```

writes `keff_sharp_limit_tensor.{csv,log,png}` and gates the fitted slope at 5%
of the arithmetic one.

The laminate cannot exercise the tensor's off-diagonal terms (its normal is
exactly `ŷ`). The curved-interface test that does is
[`../disk/`](../disk/README.md).

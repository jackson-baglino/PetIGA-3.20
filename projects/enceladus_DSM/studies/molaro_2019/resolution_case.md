# Is the expensive mesh worth it? — archived runs vs the 2026-09 campaign

Written 2026-09-10 to answer three questions: what the mesh in the archived
`BestParams` runs actually was, why an under-resolved neck cannot be rescued by
interpolating to φ = 0.5, and whether the ~500× more expensive current runs
earn their cost. Companion audit of the archived parameters:
[`bestparams_archive.md`](bestparams_archive.md).

## The short answer

**Yes — but two of the three improvements are free, and the expensive one is
stricter than it has to be.**

The archived runs fail on three independent counts. Only one of them is the
cost of a fine mesh:

| # | defect | cost to fix |
|---|---|---|
| 1 | mesh 2.8× too coarse **for their own `eps`** | ~8× elements — cheap, and was never a choice |
| 2 | additive IC inflates the t = 0 neck by **1.71×** | **free** (`-ic_grain_union 1`) |
| 3 | 5 of Molaro's 9 points below the fillet-resolution floor | this is the expensive one |

Defects 1 and 2 are not trade-offs. They are a mesh that does not match its own
`eps`, and an initial condition that measures the wrong thing. Fixing them
changes the answer and costs essentially nothing.

## You do not choose the mesh — `eps` does

This is the part that makes the cost argument sound rather than a preference.
The mesh follows from `eps` by `dx = eps/√2`, and `eps` is bounded by
Kaempfer & Plapp's thin-interface conditions. So "is the fine mesh worth it" is
really **"is the constraint on `eps` worth respecting"** — and that constraint is
what makes the model's own derived coefficients (`tau_sub`, `mob_sub`,
`alph_sub`) mean what the asymptotic matching says they mean.

The governing number is the thin-interface expansion parameter
`δ = a₁a₂·eps/(D·β_HK)`. `tau_sub` compensates it to first order, so the
residual error is O(δ²) — but only while δ < 1. Above that the "correction" is
larger than the term it corrects and the expansion is not an expansion.

| configuration | `eps` | `δ_ice` | mesh |
|---|---|---|---|
| archived `BestParams` | 9.096e-07 | **4.38** | 134 × 214 |
| `comp_eps.py` default for that α_c | 5.205e-07 | 2.51 | 659 × 1056 |
| **2026-09 campaign** | 1.18e-07 | **1.000** | 5394 × 2697 |

Worth being straight about: `comp_eps.py`'s own default lands at δ_ice = 2.5 and
says the ice channel is "intractable… M&F violate it by 100× in print". So
δ_ice < 1 is *stricter than the literature norm*, and the current campaign chose
it deliberately. The archived δ_ice = 4.38 is not a near miss of a strict
standard — it is ~1.7× past even the tool's permissive default.

Separately, and independently of `eps`: their `dx = 1.809e-6` against a sizing
rule of `eps/√2 = 6.43e-7` — the mesh is **2.81× too coarse for the `eps` they
themselves chose**, giving 4.6 elements across the φ = 0.01–0.99 band where the
current runs have 13.2. That one is just a mistake, not a budget decision.

## Why interpolating to φ = 0.5 does not rescue an unresolved neck

The natural objection: the φ = 0.5 contour is well defined at any resolution, so
why not just locate it precisely and read off the neck?

**Because the error is in the field, not in the level-set extraction.**
Interpolation gives an exact answer to the wrong question.

We have just proved the first half of that directly. The sub-cell crossing was
fixed on 2026-09-09 to interpolate in `logit(φ)`, which is *exact* for this
model's equilibrium profile at any spacing — verified to 1.6e-16 µm against a
synthetic profile. That removed a real artefact, and it moved the answers by
**±0.06 µm**. Measurement precision was never what was wrong.

What interpolation cannot fix:

**1. The initial condition measures the wrong surface.** The archived IC is
additive (`ice += 0.5 − 0.5·tanh(...)` summed over grains), so in the overlap
region both grains contribute and φ crosses 0.5 well outside either sphere.
Reconstructing their exact configuration:

| separation | exact lens neck | additive IC reports | bias | union IC |
|---|---|---|---|---|
| 172.95 µm (t = 0) | 18.83 µm | **32.27 µm** | **+13.44 (1.71×)** | 18.83 (+0.00) |
| 165 µm | 54.43 µm | 60.41 µm | +5.98 (1.11×) | 54.43 (+0.00) |
| 155 µm | 77.76 µm | 82.06 µm | +4.30 (1.06×) | 77.76 (−0.00) |
| 145 µm | 94.37 µm | 97.94 µm | +3.57 (1.04×) | 94.37 (+0.00) |

(32.27 against the 33.05 the run actually reported — the small difference is the
measurement convention, so this reproduces their IC correctly.)

The bias is **large early and small late**, so it does not cancel in a
growth-rate or exponent fit: it flattens exactly the early part of the curve
that sets the exponent. And no interpolation scheme removes it, because φ
genuinely *is* 0.5 there. The union IC has **zero** bias at any `eps` by
construction, which is why the current geometries set `-ic_grain_union 1`. This
fix is free.

**2. The dynamics below the floor are wrong, not just the measurement.** The
fillet radius is `ρ = r²/(2(R−r))`. When ρ approaches `eps` the two interfaces
bounding the fillet overlap, the double well cannot sustain them as separate
interfaces, and the Gibbs–Thomson driving force — which *is* the curvature —
saturates at ~1/eps. The model cannot represent a fillet sharper than its own
interface, however precisely you locate φ = 0.5.

The floor `sqrt(12·eps·R)` is just `ρ ≥ 6·eps` rearranged. Across Molaro's nine
measured points:

| neck width | ρ | ρ/eps archived | ρ/eps current |
|---|---|---|---|
| 32.81 µm | 1.906 µm | **2.10** | 16.2 |
| 33.60 | 2.010 | **2.21** | 17.0 |
| 37.16 | 2.523 | **2.77** | 21.4 |
| 44.68 | 3.859 | **4.24** | 32.7 |
| 46.84 | 4.313 | **4.74** | 36.6 |
| 54.18 | 6.125 | 6.73 | 51.9 |
| … | | | |
| 64.78 | 9.605 | 10.56 | 81.4 |

**Five of the nine — the whole first 18 minutes — are below the floor in the
archived setup. All nine clear it by 16× or more in the current one.**

That window is not incidental. It is where the growth is steepest and where the
exponent is determined, and the exponent is the claim the Molaro comparison
rests on.

**Honest limit on this point:** we have *directly measured* the failure at
ρ/eps = 0.59 — the 2026-08 dom2 runs opened a trapped void at r ≈ 5 µm that
persisted and deepened (φ 0.69 → 0.55) and moved only 0.498 → 0.551 under a 20×
tightening of the timestep, so it is spatial and not temporal. The archived runs
start at ρ/eps = 2.10, which is below the 6·eps criterion but well above the
0.59 where we have hard evidence. We have not measured the failure at 2.1. The
criterion says it is unsafe; we have not demonstrated how unsafe.

## How the violations show up in the dynamics

Not "are the constraints violated" but "what do they *do*". Four mechanisms,
each traceable to a number.

### 1. Half the interface kinetics is the mesh, not the physics

`tau_sub = eps·lambda·( beta_HK/a1 + a2·eps/D_th + a2·eps/D_v )`. The first term
is the physical attachment kinetics. **The other two are proportional to `eps`** —
they are the thin-interface correction, i.e. the model's own discretisation
entering the interface mobility. Their share:

| configuration | physics | **numerical** | interface is |
|---|---|---|---|
| archived, −20 °C | 49.8 % | **50.2 %** | **2.01× slower** than physics alone |
| archived, −5 °C | 79.9 % | **20.1 %** | 1.25× slower |
| `comp_eps` default, −20 °C | 63.4 % | 36.6 % | 1.58× slower |
| **2026-09 campaign, −20 °C** | **81.3 %** | **18.7 %** | 1.23× slower |
| proposed ε = 1.883e-7 | 73.1 % | 26.9 % | 1.37× slower |

Two consequences, and the second is the damaging one:

- At the archived ε the interface moves at **half** the speed the physical
  kinetics specify — the mesh parameter is setting the rate.
- **It is temperature-dependent.** The spurious slowdown is 2.01× at −20 °C and
  1.25× at −5 °C, a factor **1.6× difference between the two cases those runs
  existed to compare.** The artefact does not cancel between temperatures; it
  tilts the −20/−5 comparison directly.

And it reframes the fitted mobility: `mob_sub` was already 2.01× too slow from
discretisation, so **roughly half of the ×5 is undoing the mesh**, not a
statement about physical mobility.

### 2. The ÷100 on the source is very nearly the reciprocal of the BC error

This is the sharpest one. At `humidity = 0.70`:

| | archived | 2026-09 campaign | ratio |
|---|---|---|---|
| wall undersaturation `1 − h` | 0.30 | 2.875e-03 | **104×** stronger |
| `alph_scale` on the source | 0.01 | 1.0 | **100×** weaker |
| product | | | **1.04** |

Two ~100× errors in opposite directions. The measured ice loss confirms the
cancellation: the archived run loses **−2.87 % in area (−1.45 % equivalent
radius) over 2 h** — a plausible-looking number produced by a wall 104× too
undersaturated feeding a source 100× too weak.

So the fitted `alph_sub` factor is not a claim about attachment kinetics. To the
precision of this comparison it is **the reciprocal of the humidity error**. That
is why it looked like a good fit, and it is why the value does not transfer to a
run with a defensible wall.

### 3. The neck is fed by curvature relaxation, not vapour

With the source 100× weaker and the mobility 5× stronger, the Allen–Cahn term
dominates neck filling — and AC curvature motion is deliberately *not* coupled to
vapour (`docs/model_description.md` §3.4), so that ice does not have to come from
anywhere. The 2026-09 arm-3 run reproduced exactly this signature in the current
solver: 4× the wall undersaturation and 2.5× the grain recession of the untuned
arm, yet a 75.5 µm neck against 46.1 — inverting the ordering vapour-only
transport requires (`three_options/README.md` §6).

### 4. The IC transient, and what is *not* wrong

The additive IC is not an equilibrium profile, so the run opens by relaxing
toward one. Measured: **+1.32 µm of neck in the first 60 s**, 3.1 % of the total
2 h growth compressed into 0.8 % of the run, before the sintering rate settles.
Small in absolute terms, but it lands exactly where Molaro's data is densest.

Two things worth clearing, because a sound argument should not carry weak claims:

- **`Lambd = 1.0` is inert.** Every Λ term in the free energy is multiplied by
  the sediment fraction, and there is no sediment (`No sed grains`). Wrong in
  principle, zero dynamical effect.
- **`grad_T = 1e-4 °C/m` is negligible.** Over `Ly = 3.884e-4 m` that is a
  **3.9e-8 °C** difference across the domain. The workaround cost nothing.
- **No grid pinning is visible.** The neck series is not quantised to the mesh
  (radius/dy fractional part spans 0.001–0.992), so the measurement interpolated
  sub-cell and there is no stick-slip signature to point at.

## What the archived runs were good for

They were a **parameter search**, and as a search they were legitimate and
efficient. A 500×-cheaper model that gets the *direction* of a parameter's
effect right is the correct tool for scanning `mob_sub` × `alph_sub` space, and
it is how the ×5 / ÷100 pair was found at all. Two groups were explored
(Group1 5/0.01, Group2 5.5/0.0075), which is what a search looks like.

Other things that do not need the production mesh: domain-size convergence, the
wall-humidity calibration (a far-field property, insensitive to the neck), cost
estimation, and IC smoke tests.

The `grad_T = 1e-4` is in the same category — a workaround for a solver that
could not then run at exactly zero gradient, not a physics claim. The current
solver runs `grad_temp0 0,0,0` without complaint, so that constraint is simply
gone.

## The cost, and the lever

| | archived | current |
|---|---|---|
| elements | 28,676 | 14,547,618 (**507×**) |
| DoF | 86,028 | 43,642,854 |
| ranks | 12 | 437 |
| wall clock, 2 h sim | minutes | 186–342 min |
| **core-hours per arm** | ~1 | **1355–2491** |

Cost scales as `1/eps³` in 2D (mesh² × timestep): `(9.096/1.18)³ = 458×`,
matching the measured 507× in elements.

**The lever, if the budget matters.** The current `eps = 1.18e-7` is set by
`δ_ice < 1`. The *neck floor* — the constraint that actually protects the
comparison window — would allow `eps = 1.883e-7`:

- `δ_ice` = 1.60 (still 2.7× better than the archived 4.38, and better than
  `comp_eps`'s own 2.51 default);
- Molaro's first point sits at ρ/eps = **10.1**, comfortably above 6;
- mesh 3380 × 1690 — **2.6× fewer elements, ~4× less total cost**.

So the honest minimum for a defensible comparison is ~4× cheaper than what we
are running. Going back toward the archived mesh is not on the menu: it fails
the floor over half the comparison window and its mesh does not even match its
own `eps`.

## Recommendation

1. **Keep the current resolution for anything that gets published.** The
   δ_ice < 1 choice is defensible and the runs are already done.
2. **If cost becomes binding, drop to `eps = 1.883e-7`** (~4× cheaper) — not
   further. That is the neck floor, and below it the comparison window stops
   being resolved.
3. **Do not re-use the archived numbers as a quantitative result.** Their
   reported growth 33.05 → 76.30 µm starts from a neck inflated 1.71× by the
   IC; corrected, the true starting neck was ~18.8 µm and the growth was much
   larger than reported. Cite them as a parameter search, which is what they
   were.
4. The re-run argument does not rest on the mesh alone. It rests on: a mesh
   inconsistent with its own `eps`, an IC measuring the wrong surface, a
   comparison window half of which is below the model's resolution floor, a wall
   100× more undersaturated than the experiment supports, an imposed temperature
   gradient that was a solver workaround, `Λ = 1` with no triple junction, and
   2D planar cylinders standing in for spheres. **Any one of those alone would
   justify re-running; the mesh is simply the one that costs money.**

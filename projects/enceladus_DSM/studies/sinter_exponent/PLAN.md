# Does our model reproduce the r ~ t^(1/3) sintering exponent?

## Why

Demmenie, Woutersen & Bonn, *Ice Sintering by Sublimation and Condensation*,
J. Phys. Chem. Lett. **16**(8) 2104–2109 (2025),
[doi:10.1021/acs.jpclett.5c00050](https://doi.org/10.1021/acs.jpclett.5c00050).

Two ~1 mm-diameter ice spheres, side-by-side horizontal contact, no applied
load, in an air box held at the **saturation vapour pressure of ice** at
**−3 °C**, imaged from above for **2.5 h**. Four runs give

    r/R0 = [C_m(T)·t]^(1/m)   fitted as  r ~ t^alpha
    alpha = 0.29 / 0.33 / 0.30 / 0.26   (±0.01)

i.e. **m = 3 → evaporation–condensation**, which is precisely this model's
transport route. They also show an out-of-equilibrium control (deliberately
under-saturated) giving alpha ≈ 1/7, and argue that the century of literature
values scattered between 1/3 and 1/7 is an artifact of humidity control.

This is a **shape** observable, so it is independent of the ~50 % rate deficit
already documented in `docs/molaro_validation_synthesis.md`. It is the
cleanest external test of the model we have found.

> **Notation trap.** In their Table 1, `n` is a *mechanism index* running 1–4;
> the Kuczynski exponent is their `m`. Their mapping is viscous flow 1,
> sublimation–condensation 3, bulk diffusion 5, surface diffusion 7 — and they
> argue viscous flow should be **m = 1**, not the textbook 2, because Frenkel's
> derivation neither conserves mass nor localises dissipation at the neck.
> Quote `m`, never `n`.

## Scope

Two questions, in order of how much they cost:

1. **Do our existing Molaro simulations fit the trend?** No compute. *Done —
   see `README.md`.*
2. **Does a from-scratch replication of Demmenie's geometry give 1/3?** One
   pilot + one production run.

---

## The two things that decide everything

### The exponent is a property of the protocol, not of the curve

The same committed model curve gives `a = 0.087` fitted over the solver's dense
output and `a = 0.078` at the nine times Molaro sampled, against `a = 0.121`
once `t0` is freed the way Demmenie free it. A quoted exponent without its fit
form and window attached means nothing.

`postprocess/fit_neck_growth.py` therefore always reports three forms side by
side, plus the local log-slope:

| form | model | when it is the right one |
|---|---|---|
| `d_free` | `r = C(t+t0)^a`, t0 free | Demmenie's own protocol — the only like-for-like comparison to their alpha |
| `d_fixed` | `r = C t^a` | only when the clock zero is exact — **not** our runs, which start at r0/R = 0.09; kept as a diagnostic of how much the protocol matters |
| `kucz` | `r^m − r0^m = Kt` | the only meaningful form when r0 > 0 |
| local slope | `d ln r / d ln t` | shows whether a power law exists at all |

A fourth tool matters for the Molaro comparison: `--anchor-neck` re-zeros each
series' clock at a shared neck radius. Two curves on the same trajectory but
started at different r0 differ *only* by a time offset, so anchoring at equal
neck size overlays them correctly — it removes the offset without fitting it
away, and it is the convention Molaro et al. themselves used in their Fig. 12.

### The benchmark is not a flat 1/3

`alpha = 1/3` follows from the small-neck approximation `rho = r²/(2R)`. With
the exact rolling-ball fillet (`R² + c² = (R+rho)²`, `r = c - rho`), rho is
larger — 5 % at r/R = 0.1, 20 % at 0.3 — so the driving force `d0/rho` is
smaller and **a perfect kinetic-limited model sags**:

| r/R | 0.087 | 0.10 | 0.15 | 0.20 | 0.25 | 0.30 | 0.35 |
|---|---|---|---|---|---|---|---|
| ideal local slope | 0.326 | 0.324 | 0.320 | 0.314 | 0.309 | 0.303 | 0.296 |

So the target over the production window is **0.30–0.33**, not 0.3333 — still
well inside Demmenie's own 0.26–0.33 spread. `fit_neck_growth.py::ideal_slope`
computes this and draws it on the local-slope figure; benchmarking against the
flat line alone would make a correct model look progressively wrong at large
necks.

### There is a hard resolution floor, and it is why the grains must be big

The neck fillet has radius `rho ~ r²/(2R)`, and a diffuse interface resolves it
only once `rho` exceeds the ~`6*eps` visible band:

    r/R  >=  sqrt(12*eps/R)

| geometry | eps | R | floor | measured range |
|---|---|---|---|---|
| Molaro (committed) | 0.35 µm | 85 µm | **0.22** | 0.19 – 0.38 |
| Molaro, eps 0.60 µm | 0.60 µm | 85 µm | 0.29 | *entirely below* |
| Molaro, eps 0.86 µm | 0.86 µm | 85 µm | 0.35 | *entirely below* |
| **Demmenie (production)** | 0.32 µm | 500 µm | **0.087** | 0.09 → 0.35 |

**No exponent can be extracted from the Molaro geometry at any eps we can
afford.** Six times the grain radius at comparable eps buys a floor 2.5× lower
and an actual fit window — and 0.087 happens to land near the r/R0 ≈ 0.08
vapour/surface-diffusion crossover of Maeno & Ebinuma (1983), so the numerical
floor and the validity limit of a vapour-only model nearly coincide.

---

### The initial condition: start ON the trajectory, just above the floor

A pre-existing neck is **not** a problem in itself — it only means the clock
starts partway along the curve, which `d_free` (t0 free) and equal-neck
anchoring both absorb exactly. What matters is whether the initial state lies
*on* the trajectory. Two errors compete, and **r0/R = 0.09** is where they
cross:

- **Too small (tangent contact).** Below the floor the fillet is thinner than
  the interface, so the *dynamics* are wrong, not just the measurement. From
  tangent the run spends its first ~287 s there — and, worse, log-spaced output
  then puts only **11 of 50** snapshots inside the usable window, because 6 of
  the 8 decades covered are sub-floor.
- **Too large.** The union IC intersects two spheres, so its neck is a
  zero-radius crease while the true trajectory carries a fillet of radius
  `r0²/(2R)`. Building it is an off-trajectory transient
  `t_fill ≈ beta·rho²/(2·d0) = beta·r0⁴/(8·d0·R²)`, growing as **r0⁴**.
  Calibrated: for the committed Molaro IC (r0/R = 0.19 at −20 °C) it predicts
  701 s, and the `d_free` fit independently returned t0 = 864 s.

At r0/R = 0.09 the transient is 40 s (0.2 % of an 18 ks run) and every snapshot
lands in the window. `solve_d_union(5e-4, 5e-4, 9.0e-5)` → d = 9.959418e-04.

## Stage 1 — Reanalysis (done, no compute)

`bash studies/sinter_exponent/verification/run_exponent_fits.sh`

Fits the three committed Molaro eps arms and both validation CSVs. Results and
interpretation in `README.md`.

## Stage 2 — Pilot (~3 % of production)

```bash
./scripts/HPC/submit_batch.sh --tag sinter_pilot \
    --tests-file studies/sinter_exponent/pilot_batch.txt
```

**2a. Demmenie, production mesh truncated to 10 min** — the real 78.3M-dof
problem, ~180 steps. Deliberately *not* a coarse-eps pilot: the
kinetic/thermal/vapour split of tau_sub is set by the comp_eps safety factor
alone (64.7 % at 0.5, 47.9 % at 1.0), so a coarser arm sits in a **different
physical regime** and predicts nothing about the production run. Truncating
time changes neither the discretisation nor the regime.

**2b. Molaro from its own resolution floor**, 8 h — 765k nodes, ~14 core-h.

**Gate:** stability, and confirm the t = 0 neck is r0/R = 0.09 and the early
rate matches the `u³` extrapolation that sized `t_final`.

## Stage 3 — Production

```bash
./scripts/HPC/submit_enceladus.sh \
    sinter_2D_L2196um_eps0.32um_axisym_D1mm_r0p09 \
    sinter_T-3_h1.00_5h_a9.0e-2 demmenie
```

78.3M dof, ~5450 steps at dtmax 3.3. `t_final = 1.8e4` (5 h) rather than their
2.5 h, because their duration is a property of their apparatus:

| t_final | 2.5 h | 3.3 h | 5.0 h | 6.7 h |
|---|---|---|---|---|
| r/R end | 0.279 | 0.306 | 0.349 | 0.384 |
| window (decades) | 1.47 | 1.59 | 1.77 | 1.89 |

Report the fit over the whole window **and** over `t ≤ 9000 s`, the part
directly comparable to their runs. One snapshot is ~627 MB, so run
`neck_width.py` on the cluster and bring back only the CSV. Push before
submitting.

**Regime control, not an eps ladder.** Because eps and regime are entangled, a
coarse/fine pair cannot separate discretisation from physics. Instead:

```bash
./scripts/HPC/submit_enceladus.sh \
    sinter_2D_L2196um_eps0.32um_axisym_D1mm_r0p09 \
    sinter_T-3_h1.00_5h_a9.0e-2_dirichlet demmenie_dirichlet
```

pins T and rho_v at the outer faces, so the latent heat released at the neck
has a sink and the vapour reservoir cannot deplete — the sealed box otherwise
self-warms and starves, both of which bend the slope *downward*. Compare the
**exponent**, not the rate: if they agree, the sealed arm is fine and this is a
one-line robustness statement; if they differ, the Dirichlet arm is the more
faithful representation of their cooled-stage apparatus and is the one to
quote.

A genuine discretisation check, if wanted, is **refinement at fixed eps**
(h = eps/2 instead of eps/√2, ~4× cost) — that isolates the mesh with zero
regime change.

## Stage 4 — Synthesis

`docs/sintering_exponent_synthesis.md`, in the style of
`docs/molaro_validation_synthesis.md`. Must state:

- our alpha for the replication under **their** protocol, with the fit window
  and resolution floor quoted, against their 0.26–0.33;
- that the Molaro geometry cannot resolve an exponent, and why;
- whether the floor-start Molaro run recovers the same alpha as the 1 mm run —
  Kuczynski scaling says m should not depend on R, and if it does in our model
  then the two studies cannot be compared, which is the more interesting result;
- the standing caveat: **this model has no surface-diffusion mobility.**
  Demmenie's under-saturated alpha ≈ 1/7 is the surface-diffusion signature, so
  a vapour-only model reproducing 1/3 at equilibrium and *failing* to reproduce
  1/7 out of equilibrium is the honest expected outcome. One `-humidity 0.98`
  arm tests it cheaply.

## Known limits

- **Demmenie publish no tabulated r(t)** — Figure 3 is log–log only. The
  benchmark is their reported alpha, not a point-by-point overlay. Digitising
  the figure from the open-access copy (PMC11874037) is possible and would drop
  into `inputs/validation/` with no code change, but is out of scope here.
- **alpha_c = 9.007e-2** at −3 °C is the repo's Arrhenius fit, not a fitted
  value. It should move the prefactor `C_m(T)`, not `m`. If the exponent ever
  needs defending, re-run at alpha_c = 3e-2 (eps 9.54e-7, 3262 × 890,
  ~20 core-h) and check alpha is unchanged while C moves.
- **`Patm` is hardcoded** to 101325 Pa (`material_properties.c:160`), which is
  correct for both Demmenie's air box and Molaro's sealed cryostage. Worth a
  line in the write-up, since Molaro's own vapour-transport term is derived for
  vacuum.

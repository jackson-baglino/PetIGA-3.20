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
| `d_fixed` | `r = C t^a` | only when the clock zero is exact, i.e. our tangent-contact runs |
| `kucz` | `r^m − r0^m = Kt` | the only meaningful form when r0 > 0 |
| local slope | `d ln r / d ln t` | shows whether a power law exists at all |

### There is a hard resolution floor, and it is why the grains must be big

The neck fillet has radius `rho ~ r²/(2R)`, and a diffuse interface resolves it
only once `rho` exceeds the ~`6*eps` visible band:

    r/R  >=  sqrt(12*eps/R)

| geometry | eps | R | floor | measured range |
|---|---|---|---|---|
| Molaro (committed) | 0.35 µm | 85 µm | **0.22** | 0.19 – 0.38 |
| Molaro, eps 0.60 µm | 0.60 µm | 85 µm | 0.29 | *entirely below* |
| Molaro, eps 0.86 µm | 0.86 µm | 85 µm | 0.35 | *entirely below* |
| Demmenie strict | 0.32 µm | 500 µm | **0.087** | to be measured |
| Demmenie coarse | 0.64 µm | 500 µm | 0.124 | to be measured |

**No exponent can be extracted from the Molaro geometry at any eps we can
afford.** Six times the grain radius at comparable eps buys a floor 2.5× lower
and an actual fit window — and 0.087 happens to land near the r/R0 ≈ 0.08
vapour/surface-diffusion crossover of Maeno & Ebinuma (1983), so the numerical
floor and the validity limit of a vapour-only model nearly coincide.

---

## Stage 1 — Reanalysis (done, no compute)

`bash studies/sinter_exponent/verification/run_exponent_fits.sh`

Fits the three committed Molaro eps arms and both validation CSVs. Results and
interpretation in `README.md`.

## Stage 2 — Pilot (~$5)

```bash
./scripts/HPC/submit_batch.sh --tag sinter_pilot \
    --tests-file studies/sinter_exponent/pilot_batch.txt
```

**2a. Demmenie coarse arm** — 19.6M dof, tau_sub 16.3 s, ~550 steps.
**2b. Molaro from tangent contact**, 8 h — 774k nodes, ~14 core-h.

**Gate.** The strict arm's window does not open until r/R0 = 0.087. If the
pilot finishes below ~0.15 there is nothing to fit and `t_final` must go up
before stage 3 is worth submitting — cost is linear in steps, so this is cheap
to fix and expensive to get wrong.

Also read off the pilot: where the neck actually starts moving, so stage 3 can
set `-t_out_log_t0` there instead of spending half its snapshot budget on
IC relaxation at t < 1 s.

## Stage 3 — Production (~$60–110)

```bash
./scripts/HPC/submit_enceladus.sh \
    sinter_2D_L2200um_eps0.32um_axisym_D1mm_tangent \
    sinter_T-3_h1.00_2.5h_a9.0e-2 demmenie_strict \
    -- -dtmax 3.3 -t_out_log 40
```

78.5M dof, ~2700 steps. `-dtmax 3.3` and `-t_out_log 40` override the coarse
arm's values in the shared experiment file (one snapshot is ~630 MB, so run
`neck_width.py` **on the cluster** and bring back only the CSV).

Push before submitting so the cluster can pull.

**eps check.** Coarse vs strict is a 2× eps ratio. If the fitted alpha agrees
over the shared window, the exponent is not a discretisation artifact. If it
does not, add the `--safety 0.7` arm (eps 4.4486e-07, 6994 × 1908, tau_sub
6.75 s, ~1100 core-h) and look at the trend rather than at two points.

## Stage 4 — Synthesis

`docs/sintering_exponent_synthesis.md`, in the style of
`docs/molaro_validation_synthesis.md`. Must state:

- our alpha for the replication under **their** protocol, with the fit window
  and resolution floor quoted, against their 0.26–0.33;
- that the Molaro geometry cannot resolve an exponent, and why;
- whether the tangent Molaro run recovers the same alpha as the 1 mm run —
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

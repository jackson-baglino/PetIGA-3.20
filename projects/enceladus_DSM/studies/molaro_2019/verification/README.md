# Molaro 2019 verification gates

Two gates live here. Both write their results into this directory and exit
non-zero on failure.

---

# 1. `grain_shrinkage.py` verification

`python studies/molaro_2019/verification/verify_grain_shrinkage.py`
→ `verification.csv`, non-zero exit on failure.

## What it tests and why

Grain shrinkage is the observable the chamber humidity is fitted against, and
the signal is small: Molaro's large grain loses **2.93 %** of its diameter over
78 min. A systematic bias of even 1 % in the measurement would swallow a third
of it. There is no analytic check available on a real run — the grains are
deforming — so this synthesises a snapshot whose answer is known exactly.

Two tangent spheres (R = 72.5 / 101 µm) on the symmetry axis, carrying the
model's own equilibrium profile `phi = 1/(1 + exp(-(R-d)/eps))` at
ε = 2.5852e-07 and 4 points per ε. Tangent spheres do not overlap, so each
side's volume is exactly `(4/3)πR³` and the reported sphere-equivalent radius
must come back as R.

That exercises three things that fail silently:

- the `2πy` axisymmetric measure (a planar integral is off by orders of
  magnitude; a dropped `2π` by 6.28×)
- the split at the neck plane (no double-counting, no dropped mass)
- the volume → equivalent-radius conversion

## Result

| case | large | small |
|---|---|---|
| default split (radical plane) | +0.002 % | +0.004 % |
| split shifted +5 µm | −0.059 % | +0.167 % |
| split shifted −5 µm | +0.045 % | −0.113 % |

Tolerance 0.5 %. The ±5 µm cases show the answer is insensitive to the exact
plane once it is near the contact — the neck disc is ~1e-5 of a grain's volume.

## What it caught

**The default split plane was the midpoint between grain centres.** For unequal
radii that is not the contact point: for this pair it sits 14 µm away, and it
moved **+1.3 %** onto the small grain and −0.5 % off the large one — half the
shrinkage signal, as a pure artifact. Fixed to the radical plane

```
x = c0 + (d^2 + R0^2 - R1^2) / (2 d),    d = |c1 - c0|
```

which reduces to the contact point `c0 + R0` for tangent spheres and stays
correct once the grains overlap.

It also caught a second, smaller bug on the way: splitting the trapezoid rule
at an array index drops the segment straddling the split. `integrate_sides()`
now inserts an interpolated node at the plane, so `V_low + V_high` reproduces
the whole-domain integral. A resolution sweep (2/4/8 points per ε) puts the
residual at +0.007 % on the total, matching the analytic diffuse-tail
correction `(4π³/3)·R·ε² / V = 1.25e-4` — i.e. the integrator is exact to the
profile's own smearing and nothing more.

---

# 2. `-mob_scale` / `-alph_scale` verification

`bash studies/molaro_2019/verification/verify_kinetics_scaling.sh`
→ `kinetics_scaling.csv`, `kinetics_scaling.log`, non-zero exit on failure.

## What it tests and why

Option 3 of the three-option study scales the mobility `M_0` and the
phase-change rate `alph_sub` independently (×5 and ×1/100). The obvious knobs
for that, `-mob_sub` and `-alph_sub`, do **not** work: they set absolute values
on the scalar path in `main()`, and every Molaro run uses
`-alpha_pointwise 1`, under which `SubKinetics()` rebuilds both at every
quadrature point and overwrites them — silently, exactly as it does with
`-beta_sub0`. `-mob_scale` / `-alph_scale` are multipliers applied on *both*
paths instead.

Three things have to hold, and only a run can show them:

1. **The default is a no-op.** Omitting the flags and passing `1.0` must give
   identical kinetics, or every result predating the change is invalidated.
2. **The factors bite on the pointwise path** — this is the entire reason the
   options exist, and the failure mode being fixed is a silent one.
3. **The Jacobian stays exact.** The scaling is linear, so each derivative has
   to carry its own factor. Missing one in `SubKinetics()`'s tail would leave
   the residual right and the Jacobian wrong — which costs Newton iterations
   and dt, not correctness, and so would never announce itself.

Run on a deliberately tiny 32 × 16 mesh (with ε coarsened to match), because
`-snes_test_jacobian` builds the finite-difference Jacobian one column per DoF.
The physics is meaningless at that resolution; the derivative algebra is not,
and that is all this gate claims to test.

## Result

| check | measured | verdict |
|---|---|---|
| default is a no-op | `mob_sub` and `alph_sub` identical to the baseline | PASS |
| `-mob_scale 5.0` | ratio 4.999822 | PASS |
| `-alph_scale 0.01` | ratio 0.010000 | PASS |
| `-snes_test_jacobian` | worst `‖J − Jfd‖_F/‖J‖_F` = 5.54e-09 over 4 Newton iterations | PASS |

The ratio tolerance is 1e-3, not machine epsilon: the numbers are parsed from
the startup banner, which prints `%.4e`, so a ratio of two 5-significant-figure
values is only good to ~1e-4. That is why 5× reads as 4.999822 — print
precision, not an error in the scaling.

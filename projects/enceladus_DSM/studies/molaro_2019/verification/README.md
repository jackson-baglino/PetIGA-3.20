# `grain_shrinkage.py` verification

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

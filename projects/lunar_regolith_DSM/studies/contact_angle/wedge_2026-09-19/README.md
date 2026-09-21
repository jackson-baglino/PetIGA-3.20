# Prescribed contact angle in a wedge — batch 2026-09-19

Four runs, all clean: 61/61 distinct snapshots each, no bounds violations, no
divergences. Tapered domain of the August `wedge_bc` campaign, with the wall
free-energy term enforcing 60° and 120°, under both boundary-condition types.

## Results

| BCs | θ_Young | θ_inf | error | τ [d] | ice change |
|---|---|---|---|---|---|
| sealed | 60° | 60.039 | **+0.039** | 32.8 | 0.000 % |
| sealed | 120° | 119.953 | **−0.047** | 34.0 | 0.000 % |
| reservoir | 60° | 58.936 | −1.064 | 34.6 | **+21.3 %** |
| reservoir | 120° | 119.846 | −0.154 | 32.5 | **−21.8 %** |

The two BC types do exactly what was predicted of them, and the split is clean:

**Sealed recovers Young's angle to within 0.05°** on a curvilinear mesh, and the
ice volume does not move — a sealed box holds ~1e-6 of its water as vapour, so
only the shape relaxes.

**The reservoir grows the wetting case by 21.3 % and shrinks the non-wetting
case by 21.8 %**, near-symmetric, at the cost of ~1° of angle error at 60°. That
is about 3× the ±7.5 % the flat channel produced under the same BCs, which is
the taper doing its job: in a flat channel θ fixes the meniscus curvature
outright, whereas here the channel width varies along x, so the same θ gives a
curvature that varies with position and the ice has a gradient to move along.

The reservoir's angle error is strongly asymmetric — −1.064° at 60° against
−0.154° at 120°. Both are negative, i.e. both read below Young, which is what a
moving interface should do (σ = d₀χ + βv_n rather than d₀χ alone). Why the
wetting case pays ~7× more is not explained here.

## What this validates

The wall term and the measurement both work on a **curvilinear** mesh, and they
are independent of each other — the solver imposes cos θ through a boundary
integral, the measurement reads the angle off the φ = 0.5 contour geometry.
Agreement to 0.05° between them is not something to get by accident.

Before the runs, `-test_wall_measure` confirmed that PetIGA integrates the
boundary form over the **true** sloped-wall length (6.184658e-04 m =
Lx·√(1+0.25²) per wall) rather than the projected 6.0e-04, to 1.0e-15. That is
a different code path from the flat channel, where PetIGA takes a detS = 1
shortcut, and it had never been exercised.

## Measurement changes

`contact_angle.py` now handles ruled patches generally:

- coordinates come from the full physical `X, Y` arrays rather than a separable
  `x1d, y1d`, so rows need not lie at constant y;
- walls are affine curves `y = y₀ + s·x` read from `-wall_bot_*`/`-wall_top_*`,
  defaulting to the flat channel;
- the 5ε exclusion uses **perpendicular** distance to the wall line — a vertical
  offset would under-cut the near wall and over-cut the far one on a taper;
- the contact point comes from intersecting the fitted conic with the wall
  *line*, a quadratic that reduces to the old form when s = 0;
- φ probes use an exact inverse map, available because the patch is ruled:
  x depends on u alone, and v = (y − y_bot)/(y_top − y_bot).

Verified backward-compatible: the flat-channel sealed run re-measures to
60.074 / 60.004, identical to before the change, and all 11 synthetic angles
still return exact.

## A withdrawn gate

A synthetic sloped-wall case was added to
`verify_contact_angle_measure.py` and then **withdrawn**. It built the bridge in
an unsheared frame and mapped it onto the tapered walls, but that map scales y
by an x-dependent factor and does not preserve angles — so the synthetic never
had θ at the walls, and the gate was measuring its own distortion. It reported
errors up to 25° against a measurement that is demonstrably correct.

A correct construction exists (a circle centred on the wedge axis at distance c
from the apex meets both walls at cos θ = c·sin α/ρ, with c = 0 giving the
apex-centred 90° arc) but is not attempted here. Until it is, the sloped-wall
path rests on the solver comparison above.

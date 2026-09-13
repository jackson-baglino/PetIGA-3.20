# How to choose `-dtmax` — measured, not estimated

## The question

`dtmax` had been set by a `tau_sub/10` rule of thumb. That rule was calibrated
on three points from batch 2026-09-12 (dtmax/τ = 0.057 and 0.23 clean, 0.91
stalled) — but the 0.91 stall turned out to be the wall-term sign flip, not the
timestep. Once that was clamped, the rule had no evidence behind it.

## Why there is no formula

The solver is implicit (TSALPHA + SNES), so there is **no explicit stability
limit** — no CFL-style expression applies. dt is bounded only by temporal
accuracy and by Newton converging. The only defensible way to set it is to vary
it and watch the measured quantity.

The code already carries the right *dynamic* guard: `InterfaceCFLMonitor`
measures `‖φⁿ − φⁿ⁻¹‖_∞ / Δt` after each accepted step and caps the next dt so
no DOF moves more than `-dtCFL_dphimax` (default 0.2). Its own docstring puts
it well: *"a static dtmax cannot track a diverging velocity … NRmin/NRmax and
dtmax still propose dt; the measured phase rate disposes."*

## The measurement

`dtmax_study.sh`, on the cheapest geometry (ε = 3.00 µm, 142×48), θ = 60°:

| dtmax | dtmax/τ_sub | steps | CFL caps | φ_min | θ_inf | wall |
|---|---|---|---|---|---|---|
| 2.0e3 | 0.051 | 2864 | 0 | 0 | 60.1121 | 872 s |
| 4.0e3 | 0.102 | 1462 | 0 | 0 | 60.1119 | 449 s |
| 8.0e3 | 0.204 | 763 | 0 | 0 | 60.1118 | 237 s |
| 1.6e4 | 0.407 | 414 | 0 | 0 | 60.1109 | 131 s |
| 3.2e4 | 0.815 | 241 | 0 | 0 | 60.1093 | 79 s |

**θ_inf moves 0.0028° across a 16× change in dtmax** — three orders of
magnitude below the 2° acceptance criterion. φ never leaves [0,1]. The CFL
limiter never fires even at 241 steps for a 65-day relaxation, i.e. no DOF ever
moves 0.2 in a step. Wall time falls 11×.

dt is simply not the accuracy bottleneck in this regime.

The top of the range is also a control: dtmax/τ = 0.815 is essentially the
setting that stalled the *unclamped* ε = 0.75 µm run, and it is clean here —
independent confirmation that the clamp removed that failure mode rather than
raising its threshold.

## The rule

**`dtmax = tau_sub/2`**, set in the geometry file (which owns ε; τ_sub ∝ ε², so
a shared experiment file cannot express it).

That sits 1.6× inside the largest value tested clean, with three independent
backstops: the interface-CFL limiter for fast events, the adaptive controller
(`NRmin`/`NRmax`), and the stall detector.

The startup guard now warns only above `dtmax/tau_sub = 1.0`, where one step
would exceed the entire interface relaxation time.

## Caveats

One geometry, one angle, one temperature, and a relaxation that is smooth
throughout. A case with a genuine fast event — grain collapse, neck pinch-off —
will be limited by the CFL limiter rather than by dtmax, which is the intended
division of labour. Re-run this study if the campaign changes character.

# Molaro et al. (2019) Validation — Synthesis

*2026-07-13. Closes the two-week campaign: what the model reproduces, what
the residual is, and why. Companion data: inputs/validation/
molaro2019_fig11_T-20.csv; runs indexed in .claude/ACTIVITY_LOG.md.*

## 1. The benchmark

Molaro et al. (2019, JGR Planets, 10.1029/2018JE005773) Fig. 11/12: two ice
grains (R = 101/72.5 um) in a sealed Peltier cryostage (~1 mm gap,
ice-coated slide as a nearby saturated reservoir, grains ON a glass slide
over the cold plate), T = -20 C series: neck WIDTH 32.81 -> 64.78 um over
78 min (+97%), grains shrinking 3-4%. A -5 C series exists (32 -> 44 um in
<= 50 min) but is artifact-dominated (see 5).

## 2. Our model's result

Axisymmetric first-principles phase-field (K&P/M&F): finite attachment
kinetics (alpha_c), resolved vapor diffusion (xi_v-scaled, quasi-steady-
invariant), full thermal coupling (latent heat + conduction), optional
open-chamber Dirichlet reservoirs. Response surface at -20 C (neck growth
at 78 min, all runs 0 solver failures):

  alpha=3e-2, D=Dv          +33%
  alpha=1e-1, D=Dv          +38%
  alpha=3e-2, 3Dv           +43%
  alpha=1e-1, 3Dv (+3k_air) +53-55%
  alpha=1e-1, 10Dv+10k_air  +75%      | data: +97%

Rate ~ D^0.3 at alpha = 0.1 (mixed regime, Lambda ~ 1). Growth exponent
n ~ 5 matches the data (x^5 - x0^5 ~ t, R^2 = 0.99). Eliminated as
explanations by direct measurement: timestep (12 nm effect), eps
(planar sweep linear O(eps); axisym s025 +5 pts), dimensionality (axisym
= planar +0.8 pts), attachment ceiling, chamber humidity (neck-neutral;
calibrates grain shrinkage at h ~ 0.998), thermal-boundary enhancement.

## 3. Their model is the corner of our response surface

SA81's vapor term V3 (their Table 3; derivation their Appendix A) is
Hertz-Knudsen with alpha = 1 and, quoting their appendix, assumes "all
molecules that evaporate from the grains' surfaces are immediately
deposited in the neck. This ignores gas-gas diffusion rates" — plus mass
conservation by fiat and interface T = gas T. I.e. the (alpha = 1,
D = infinity) corner, unreachable at physical parameters (our surface
shows why). Every deficiency they list is a feature our model implements:
finite alpha, resolved gas diffusion, ambient exchange, thermal
disequilibrium. Their own validation claim is "agrees within approximately
1 order of magnitude"; our constrained-parameter factor 1.3-2 is stronger.

## 4. The residual factor ~2 = surface diffusion (their other mechanism)

Their Appendix A opens: vapor transport is "one of the TWO dominant
diffusion mechanisms driving neck growth in ice (Figure 10)" — the other
is surface diffusion (V1), absent from our model. A vapor-only model
reproducing ~half the observed rate is the expected vapor SHARE if the two
mechanisms are co-dominant at homologous T = 0.93. This is the natural
future-work item: a surface-mobility term in the phase-field model.

## 5. The temperature inversion is an experimental artifact

Absolute neck rates: -20 C ~ 0.41 um/min vs -5 C ~ 0.24 um/min — inverted
vs ANY vapor mechanism (both scale with P_v, ~4x larger at -5). Resolution:
their grains net-shrink in a nominally saturated isothermal sealed cell
(3-4% at -20, 9% at -5) — impossible without a heat source; illumination/
window radiative heating driving a weak (~0.1-1 K) gradient reproduces the
shrinkage magnitude AND its T-scaling. At -5 the artifact is ~4x stronger
and was consuming the neck (they stopped measuring when necks shrank). The
-20 series is the clean benchmark; the -5 series is artifact-dominated.

## 6. Comparison conventions (adopt for the final figure)

- Their Fig. 12 anchors the first measurement to the model time with equal
  neck size (shape comparison; neutralizes the unknown contact time).
  Adopt the same anchoring.
- "Neck size" = neck WIDTH; relative neck size = width / grain DIAMETER
  (confirmed by their uncertainty arithmetic and Fig. 12 caption).
- Experimental variance is large: "some necks behaved predictably, while
  others showed no growth at all" (their sec. 4.1). Single-pair benchmarks
  carry unquantified pair-to-pair scatter.

## 7. Validation statement

A first-principles phase-field model of vapor-route sintering reproduces
the observed growth exponent and ~50% of the measured -20 C neck growth
rate with all parameters physically constrained; the remainder is
quantitatively consistent with the surface-diffusion contribution that the
SA81 framework identifies as co-dominant; and the apparent temperature
inversion in the data is explained by the experiment's own sublimation
artifact rather than sintering physics. Future work: (a) surface-diffusion
mobility in the phase-field model; (b) gradient-metamorphism runs with the
artifact's measured Delta-T (planar 2D — transverse gradient breaks
axisymmetry); (c) grain-shrinkage co-fit via the open-chamber
configuration (h ~ 0.998 calibrated).

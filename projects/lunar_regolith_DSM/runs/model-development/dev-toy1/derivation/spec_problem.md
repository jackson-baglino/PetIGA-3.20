# M1 Scope: toy 1D heat conduction with a source (DRY RUN ONLY)

This is a plumbing test of the model line. Per the brief
(`workflows/model-development/problems/toy_heat1d.md`, lines 3–4), it is not a result and must not be
cited. No paper is cited here, so `citations.json` is an empty list. Every number below comes from the
brief (file and line), from the proposed reference case in choice C1, or from a computation shown in the
text. The computations are repeated in `scratch/m1_scales.py`, and its output matches the values given here.

## Problem

A rod occupies 0 ≤ x ≤ L. It has constant conductivity k, density ρ and heat capacity c_p, and a uniform
volumetric heat source q. Both ends are held at T(0, t) = T(L, t) = T_0, and the rod starts at
T(x, 0) = T_0 (brief, lines 8–9).

The model must answer these questions:
1. What is T(x, t) for t > 0?
2. Does the computed steady state equal the analytic parabola T_0 + q x (L − x) / (2k) (brief, line 17)?
3. How long does the rod take to reach that steady state? This sets t_end for the tests (see `## Scales`).

## Mechanisms

| Mechanism | In/out | Scaling estimate and source |
|---|---|---|
| Heat conduction (Fourier's law, constant k) | **In** | Brief, lines 8 and 21. The conduction term k ∂²T/∂x² is about k ΔT/L². It is the only term that removes heat, so at steady state it balances q exactly. That balance gives ΔT ~ qL²/k (computation in `## Scales`). |
| Uniform volumetric source q | **In** | Brief, lines 8 and 21. Its size is q, which drives the peak steady temperature rise ΔT_max = qL²/(8k) (computed in `## Scales`). |
| Heat storage ρ c_p ∂T/∂t | **In** | It is part of the transient heat equation, and the brief asks for T(x, t) (line 9). It balances conduction on the time scale L²/κ. |
| Temperature dependence of k, ρ, c_p | Out | The brief fixes them as constant (line 8). The relative change in k over the rod is about β ΔT_max, where β = d ln k / dT. For the reference case ΔT_max/T_0 = 0.125 (computed in `## Scales`). So if k varied like a power of T, the error would be of order 10⁻¹ times that power. No material is specified, so no value of β can be sourced; the assumption is accepted by decision of the brief (line 21). |
| Radiation, convection or lateral losses from the rod surface | Out | Brief, line 21: "out: everything else". No surface coefficient is given, so no number can be sourced; a 1D rod with fixed ends and no lateral flux is the stated system (line 8). |
| Phase change, advection, thermal expansion | Out | Brief, line 21. No latent heat, velocity or expansion coefficient appears in the system (lines 8–9). |

## Scales

The scales are built only from the parameters in the brief (line 8). Numbers use the proposed unit
reference case of C1 (L = 1 m, k = 1 W m⁻¹ K⁻¹, ρ = 1 kg m⁻³, c_p = 1 J kg⁻¹ K⁻¹, q = 1 W m⁻³,
T_0 = 1 K). These values are **provisional until the user decides C1**.

- **Diffusivity:** κ = k/(ρ c_p) (brief, line 13). Reference case: κ = 1/(1·1) = 1 m² s⁻¹.
- **Length:** x_c = L. Reference case: 1 m.
- **Time:** t_c = L²/κ. Reference case: 1²/1 = 1 s.
  - The slowest transient mode, sin(πx/L), decays as exp(−π²κt/L²). Its e-folding time is
    t_slow = L²/(π²κ) = 1/π² ≈ 0.1013 s.
  - That mode falls to 1 % of its initial value after t_99 = t_slow ln 100 ≈ 0.1013 × 4.605 ≈ 0.4666 s.
  - So the steady-state test needs t_end of a few tenths of t_c or more.
- **Temperature:** the base value is T_0 (brief, line 9). The source-driven rise scale is
  ΔT_q = qL²/k, obtained by balancing k ΔT/L² against q. Reference case: 1·1²/1 = 1 K.
  - The peak of the steady parabola is at x = L/2. Its value is
    T(L/2) − T_0 = q (L/2)(L/2)/(2k) = qL²/(8k). Reference case: 0.125 K.
- **Vapor density:** not applicable. There is no vapor in this toy problem (brief, lines 8 and 21).

**Dimensionless form.** Let ξ = x/L, τ = t κ/L² and θ = (T − T_0)/ΔT_q. Then ρ c_p ∂T/∂t = k ∂²T/∂x² + q
becomes
∂θ/∂τ = ∂²θ/∂ξ² + 1, with θ(0, τ) = θ(1, τ) = 0 and θ(ξ, 0) = 0.
This form has no parameters. M3 will write the step-by-step derivation into the ledger.

**The one dimensionless group** (brief, line 22) is Q = qL²/(k T_0) = ΔT_q / T_0, the source-driven rise
relative to the base temperature. Reference case: Q = 1·1²/(1·1) = 1.
- Q appears only if T is scaled by T_0 (θ' = T/T_0 gives ∂θ'/∂τ = ∂²θ'/∂ξ² + Q).
- It does not change the shape of the solution. It measures how large the temperature excursion is, and
  so how far the constant-property assumption is stretched: ΔT_max/T_0 = Q/8 = 0.125.

**Dominant terms.**
- On the time scale τ ~ 1, storage, conduction and the source are all O(1). No term is negligible.
- At steady state (τ ≳ 0.47, from t_99 above), storage vanishes and conduction balances the source
  exactly.

## Constitutive laws

| Law | Formula | Source | Validity range | Holds for this problem? |
|---|---|---|---|---|
| Fourier's law of conduction | heat flux j = −k ∂T/∂x (W m⁻²), with k constant | Brief, line 8 (constant conductivity k). This is a standard law; no paper is cited, because the brief says none is needed and the toy must not be cited (lines 3–4, 23). | Continuum media where gradients are slow compared with mean free paths. Constant k needs a small relative temperature change (β ΔT_max ≪ 1). | Yes, by decision of the brief. The relative excursion is ΔT_max/T_0 = Q/8 = 0.125 for the reference case. |
| Caloric law | internal energy per volume e = ρ c_p (T − T_ref), with ρ and c_p constant | Brief, line 8 | No phase change and constant c_p over the temperature range | Yes, by decision of the brief. There is no phase change (line 21). |
| Source | uniform q (W m⁻³), constant in x and t | Brief, line 8 | As stated | Yes |

There are no pressure-dependent laws in this problem.

## Notation

Symbols from the brief (line 13) are locked here, pending user approval.

| Symbol | Meaning | Units | Range |
|---|---|---|---|
| x | position along the rod | m | [0, L] |
| t | time | s | [0, t_end] |
| L | rod length | m | > 0 |
| T(x, t) | temperature | K | > 0. With q ≥ 0, T ∈ [T_0, T_0 + qL²/(8k)] (maximum principle) |
| T_0 | boundary and initial temperature | K | > 0 |
| k | thermal conductivity | W m⁻¹ K⁻¹ | > 0 |
| ρ | density | kg m⁻³ | > 0 |
| c_p | specific heat capacity | J kg⁻¹ K⁻¹ | > 0 |
| κ | thermal diffusivity, k/(ρ c_p) | m² s⁻¹ | > 0 |
| q | volumetric heat source | W m⁻³ | real. Assumed ≥ 0 for the range of T above |
| j | conductive heat flux, −k ∂T/∂x | W m⁻² | real |
| ΔT_q | source temperature scale, qL²/k | K | ≥ 0 |
| ξ, τ, θ | dimensionless x/L, κt/L², (T − T_0)/ΔT_q | 1 | ξ ∈ [0, 1], τ ≥ 0, θ ∈ [0, 1/8] |
| Q | the one dimensionless group, qL²/(k T_0) | 1 | ≥ 0 |
| n | outward unit normal at the ends | 1 | n = −1 at x = 0, n = +1 at x = L |

**Sign convention.** The outward heat flux at an end is j·n. Both ends have Dirichlet conditions, so no
flux condition is imposed. The flux at the ends is an output of the solution. At steady state each end
carries qL/2 outward, which is a computation from the parabola: −k ∂T/∂x at x = L equals q L/2.

## Out of scope

- **Temperature-dependent properties.** The brief fixes them as constant (line 8). The error they would
  cause scales with β ΔT_max (see `## Mechanisms`).
- **Radiation, convection, lateral losses, phase change, advection.** The brief puts "everything else" out
  of scope (line 21).
- **Multi-dimensional effects.** The system is a 1D rod (line 8).
- **The C solver and PetIGA.** The implementation is a Python module, `toy/heat1d.py` (lines 14–16).
- **The ice/vapor/regolith model.** This is a plumbing dry run (lines 3–4). Nothing here feeds the PhD
  model.
- **Performance work (M8).** The brief says to skip it (line 26).

## Choices for the user

**C1: Reference parameter values.** The brief gives no numbers for L, k, ρ, c_p, q, T_0.
- **Options:**
  - (a) the unit case L = k = ρ = c_p = q = T_0 = 1 in SI units, which gives Q = 1;
  - (b) a physical material, for example a metal rod, with values that would then need a source;
  - (c) a set of user-specified values.
- **Recommendation:** (a). The equation in θ has no parameters, so the unit case tests the plumbing fully.
  It also keeps the analytic peak at a clean 0.125 K, and it needs no external source. All reference
  numbers in this file assume (a).

**C2: Sign of q.** Should negative q (a sink) be allowed?
- **Recommendation:** allow any real q in the code, but test with q > 0.
- **Reason:** the parabola test (line 17) holds for any sign. Only the stated range of T assumes q ≥ 0.

Already decided by the brief, and not reopened here: the notation (line 13), finite differences with
backward Euler in a NumPy module (lines 14–16), and the parabola comparison test (line 17).

# M3 Derivation: toy 1D heat conduction with a source (DRY RUN ONLY)

This is a plumbing test of the model line. It is not a result and must not be cited (problem.md, lines 3–4).
The model is `heat1d-fd-be` (decisions.json, G:Selection). This file presents `ledger.json` in list order.
Each displayed equation carries its step ID. No paper is cited, as the brief asks (problem.md, line 23). So
every step is one of three things: an assumption taken from the brief or the spec (A1–A8), the decided
choice C1, or a derived step that names its operation and its premises.

Symbols and units are listed in `ledger.json` → `symbols`. All of them come from the notation table that was
approved at G:Scope (spec_problem.md, lines 80–96), except these helpers:
- the control-segment end points a, b, the internal energy e, T_ref, and the outward fluxes F_0, F_L;
- for the uniqueness argument: T_1, T_2, w, E;
- for the minimum principle: τ_1, θ_min, (ξ_*, τ_*);
- for the steady state: θ′, θ_s, T_s, c_0, c_1;
- for the transient: v, ℓ, r, β_ℓ, g, τ_g.

The symbols ψ, φ, eps and eps_s are named only to say that they are absent.

**Where the numbers come from.** Every number in this file was printed by one of two scratch scripts, both
run at try 3 with the project venv:
- `derivation/scratch/m3_series.py` covers S50–S57;
- `derivation/scratch/m3_try3.py` covers S13, S21, S22, S31, S56, S57, S32 and the §8 test note.

None of the numbers was computed by hand. The scripts are my own checks, not evidence for anyone else.

## 0. What changed in this revision (try 3, answer to G:Derivation)

The main session's note for try 3 (decisions.json, 2026-11-16T18:28:21-08:00) asks for five changes. No step
ID was added, removed or renumbered, and no displayed equation changed.

| Item in the note | Fix | Steps |
|---|---|---|
| `required_premises`: S12 needs S5 | S12 now lists S5. Its operation says why: the localization lemma needs a continuous integrand, and q is continuous because it is constant (S5). | S12 |
| `required_premises`: S21 and S22 need S6 | S21 and S22 now list S6. T is differentiable once in t and twice in x by S6, so θ is too, through the affine map S20. | S21, S22 |
| `tref_constant_for_S13` | S1 now states that T_ref is constant in x and t, and S13 cites S1 for it. `m3_try3.py` confirms the term this removes: with T_ref(t), ∂_t e = ρc_p(∂_tT − dT_ref/dt). S1's `if_wrong` says what that term would mean. | S1, S13 |
| (my own audit, not in the note) | S28 reuses the chain rule "as in S21–S22". It now lists S21 and S22 as premises, as S34 already lists S9. The verifier rated this pattern informative, so this is for consistency only. | S28 |
| D-b (a) | S6 kept as broadened (T, ∂_xT, ∂_x²T and ∂_tT continuous up to the ends for t > 0). | S6 |
| D-c (a) | The upper bound θ ≤ 1/8 stays out of the ledger (§5). | none |
| derivation.md up to date | S48–S57 are presented in §5 and §8. S32 is worded as in the ledger. The out-of-date remarks about coverage and `S32_magnitudes` are removed. The hand arithmetic is replaced by script output. | S32, S48–S57 |

Try 2 fixed failures 1–6 of the earlier results. Those fixes are still in place: the `from` lists, S6
regularity at the ends and at t = 0, θ ≥ 0 (S48–S50), S32 restricted to 0 < τ ≤ 0.1 and τ ≥ 1, and T_1 and
T_2 in the symbols table. The ledger order is the derivation order, so the IDs are not monotone in two
places: S48–S50 come before S38, and S51–S57 and S32 come after S45.

## 1. Assumptions

- **S1 (A1).** The rod is 0 < x < L and t > 0. The coefficients k, ρ and c_p are constant and positive.
  The reference temperature T_ref of the caloric law is constant in x and t. The rod is 1D and exchanges
  no heat through its sides.
- **S2 (A2).** Energy balance on every segment [a, b] ⊂ [0, L]. The only channels are conduction through
  the ends of the segment and the source:

  $$\frac{d}{dt}\int_a^b e\,dx = j(a,t) - j(b,t) + \int_a^b q\,dx \tag{S2}$$
- **S3 (A3).** Caloric law:

  $$e = \rho c_p (T - T_{ref}) \tag{S3}$$
- **S4 (A4).** Fourier's law:

  $$j = -k\,\frac{\partial T}{\partial x} \tag{S4}$$
- **S5 (A5).** q is a real constant, uniform in x and t. Decision C2 allows either sign.

  $$q\in\mathbb{R}\ \text{constant in } x \text{ and } t \tag{S5}$$
- **S6 (A6).** Regularity, which defines the class of solutions considered. It is broadened as accepted at
  D-b (a):

  $$T\in C^{2,1}\big((0,L)\times(0,\infty)\big);\quad T\in C^{2,1}\big([0,L]\times(0,\infty)\big);\quad T\in C\big([0,L]\times[0,\infty)\big) \tag{S6}$$

  Each part has a user:
  - The interior class serves the localization step S12 and the chain rule in S21 and S22.
  - The class up to the ends for t > 0 serves S9–S11 (with a = 0 or b = L), S34, S35 and S39.
  - Continuity on the closed set, including t = 0 and the corners, serves S36 and S48. It is compatible with
    the data, because the boundary and initial values agree at the corners (both are T_0).

  Derivatives are not required to be continuous at t = 0. At the corners (0, 0) and (L, 0) the data are
  incompatible in their derivatives when q ≠ 0 (selection.md, item 1). There ∂_tT jumps from 0 at the end
  to q/(ρc_p) inside the rod. S54 shows that the class is not empty: the series solution belongs to it.
- **S7 (A7).** End and initial data:

  $$T(0,t) = T(L,t) = T_0\ (t>0),\qquad T(x,0) = T_0,\qquad T_0 > 0 \tag{S7}$$
- **S8 (A8).** There are no grains, no ice, no vapor and no diffuse interface. The symbols ψ, φ, eps and
  eps_s do not appear.

Each assumption's `if_wrong` entry in the ledger says what would change if it failed.

## 2. Strong form

The integral balance is converted into a pointwise equation in four moves:

$$\frac{d}{dt}\int_a^b e\,dx = \int_a^b \frac{\partial e}{\partial t}\,dx \tag{S9}$$

S9 is the Leibniz rule with fixed limits. ∂e/∂t = ρc_p ∂T/∂t is continuous on [a, b] × (0, ∞), including
a = 0 and b = L, by S3 and S6.

$$j(a,t)-j(b,t) = -\int_a^b \frac{\partial j}{\partial x}\,dx \tag{S10}$$

S10 is the fundamental theorem of calculus, the 1D divergence theorem. ∂j/∂x = −k ∂²T/∂x² is continuous up
to the ends, by S4 and S6.

$$\int_a^b \Big(\frac{\partial e}{\partial t} + \frac{\partial j}{\partial x} - q\Big)dx = 0\quad \forall\, 0\le a<b\le L \tag{S11}$$

S11 is S9 and S10 substituted into S2.

$$\frac{\partial e}{\partial t} + \frac{\partial j}{\partial x} = q\quad\text{in }(0,L)\times(0,\infty) \tag{S12}$$

S12 is the localization lemma: a continuous integrand that integrates to zero on every subinterval is zero.
The integrand of S11 is continuous for two reasons:
- ∂_t e and ∂_x j are continuous by S6;
- q is continuous because it is constant (S5).

Next, the constitutive laws are substituted:

$$\rho c_p\frac{\partial T}{\partial t} + \frac{\partial j}{\partial x} = q \tag{S13}$$

S13 is S3 substituted into S12. ρ and c_p are constant, and T_ref is constant in t, all by S1. So
∂_t e = ρc_p ∂_tT. If T_ref depended on t, the extra term −ρc_p dT_ref/dt would appear (`m3_try3.py`).

$$\frac{\partial j}{\partial x} = -k\frac{\partial^2 T}{\partial x^2} \tag{S14}$$

S14 differentiates S4. The product-rule term (∂_x k)(∂_x T) is zero because k is constant (S1).

$$\boxed{\rho c_p\frac{\partial T}{\partial t} = k\frac{\partial^2 T}{\partial x^2} + q}\quad\text{in }(0,L)\times(0,\infty) \tag{S15}$$

S15 together with the data S7 is the full strong form. Its units are W m⁻³ in every term:
- ρ c_p ∂T/∂t: kg m⁻³ · J kg⁻¹ K⁻¹ · K s⁻¹ = W m⁻³;
- k ∂²T/∂x²: W m⁻¹ K⁻¹ · K m⁻² = W m⁻³;
- q: W m⁻³.

$$\kappa \equiv \frac{k}{\rho c_p} > 0 \tag{S16}$$
$$\frac{\partial T}{\partial t} = \kappa\frac{\partial^2 T}{\partial x^2} + \frac{q}{\rho c_p} \tag{S17}$$

S16 is the definition of κ. S17 is S15 divided by ρ c_p > 0.

## 3. Reduction where the grains are absent

$$\text{S15}\big|_{\psi\equiv1} = \text{S15} \tag{S18}$$

By S8, no grain indicator ψ appears in S15, so the grain-free reduction is the identity, term by term:
- storage ρ c_p ∂T/∂t is unchanged;
- conduction k ∂²T/∂x² is unchanged;
- the source q is unchanged.

**Anchor comparison.** This is informative only and not a ledger step; see the last section of this file.
- The heat block of today's solver is `R[a][1]` in `Residual_A1` (src/assembly.c, lines 183–185). Its weak
  form is ∫ [ρ c_p N ∂_t T + ξ_T k ∇N·∇T − pc ξ_T ρ_ice L_sub ∂_t φ N] = 0.
- k, c_p and ρ are linear mixtures in φ (src/material_properties.c, lines 40, 77 and 113). So with φ uniform
  and constant in time, ξ_T = 1 and dim = 1, they are constants and the latent term vanishes.
- Integrating by parts then gives ρ c_p ∂_t T = k ∂²_x T. This is S15 with q = 0.

The toy therefore adds exactly one term to the anchor's heat block, the uniform source q. The toy does not
use the C solver (spec_problem.md, line 109), so this is a consistency remark, not a requirement. The
verifier's `reduction_vs_Residual_A1` check reached the same result independently.

## 4. Nondimensional form

The reference case comes from the decided choice C1:

$$L=1\,\mathrm{m},\ k=1,\ \rho=1,\ c_p=1,\ q=1\,\mathrm{W\,m^{-3}},\ T_0=1\,\mathrm{K}\ \text{(SI)} \tag{S19}$$

S19 is choice **C1**. The main session decided it as (a) at G:Scope, for plumbing only. Its later notes
confirm that C1 stands for the toy.

The scales are the approved notation (spec_problem.md, lines 93–94). They need q ≠ 0:

$$\Delta T_q \equiv \frac{qL^2}{k},\quad x_c=L,\quad t_c=\frac{L^2}{\kappa},\quad \xi=\frac{x}{L},\quad \tau=\frac{\kappa t}{L^2},\quad T = T_0+\Delta T_q\,\theta(\xi,\tau) \tag{S20}$$

The chain rule converts the derivatives. It applies because T is differentiable once in t and twice in x
(S6), and S20 is an affine change of variables, so θ inherits that smoothness:

$$\frac{\partial T}{\partial t} = \frac{\kappa\Delta T_q}{L^2}\frac{\partial\theta}{\partial\tau} \tag{S21}$$
$$\frac{\partial T}{\partial x} = \frac{\Delta T_q}{L}\frac{\partial\theta}{\partial\xi},\qquad \frac{\partial^2 T}{\partial x^2} = \frac{\Delta T_q}{L^2}\frac{\partial^2\theta}{\partial\xi^2} \tag{S22}$$

`m3_try3.py` applies the chain rule symbolically and finds residual 0 for S21 and for both parts of S22.

Substituting S21 and S22 into S17:

$$\frac{\kappa\Delta T_q}{L^2}\frac{\partial\theta}{\partial\tau} = \kappa\frac{\Delta T_q}{L^2}\frac{\partial^2\theta}{\partial\xi^2} + \frac{q}{\rho c_p} \tag{S23}$$

Dividing by κ ΔT_q / L², which is nonzero because q ≠ 0:

$$\frac{\partial\theta}{\partial\tau} = \frac{\partial^2\theta}{\partial\xi^2} + \frac{qL^2}{\rho c_p\,\kappa\,\Delta T_q} \tag{S24}$$

The coefficient simplifies in two steps. First ρ c_p κ = k (S16). Then qL²/(k · qL²/k) = 1. This gives

$$\boxed{\frac{\partial\theta}{\partial\tau} = \frac{\partial^2\theta}{\partial\xi^2} + 1}\quad\text{in }(0,1)\times(0,\infty) \tag{S25}$$
$$\theta(0,\tau)=\theta(1,\tau)=0,\qquad \theta(\xi,0)=0 \tag{S26}$$

S26 maps the data S7 through S20: x = 0 ↦ ξ = 0, x = L ↦ ξ = 1, and T = T_0 ↦ θ = 0.

**Scaling valid for every q, including q = 0.** This is spec_problem.md, lines 55–57:

$$T = T_0\,\theta'(\xi,\tau),\qquad Q\equiv\frac{qL^2}{kT_0} \tag{S27}$$
$$\frac{\kappa T_0}{L^2}\frac{\partial\theta'}{\partial\tau} = \kappa\frac{T_0}{L^2}\frac{\partial^2\theta'}{\partial\xi^2} + \frac{q}{\rho c_p} \tag{S28}$$
$$\frac{\partial\theta'}{\partial\tau} = \frac{\partial^2\theta'}{\partial\xi^2} + Q,\qquad \theta'(0,\tau)=\theta'(1,\tau)=1,\ \ \theta'(\xi,0)=1 \tag{S29}$$
$$Q = \frac{\Delta T_q}{T_0},\qquad \theta' = 1 + Q\,\theta\quad(q\ne0) \tag{S30}$$

**The dimensionless group and its reference value.** `m3_try3.py` substitutes S19 into S16, S20 and S27
and prints κ = 1, t_c = 1, ΔT_q = 1 and Q = 1:

$$\kappa=1\,\mathrm{m^2s^{-1}},\quad t_c=1\,\mathrm{s},\quad \Delta T_q=1\,\mathrm{K},\quad Q=1 \tag{S31}$$

**Dominant terms: forward pointer.** All coefficients of S25 are 1. Which terms dominate depends on τ, and
the proof needs the steady state and the series solution. So that statement, S32, is in §8, after its
premises. In short: for 0 < τ ≤ 0.1 no term is negligible, and for τ ≥ 1 the solution is near steady.

The solver-scaling concern in the role brief does not arise here. That concern is SI residual blocks of
1e-13 and smaller, and in θ form every coefficient is 1 (S25).

## 5. Uniqueness, minimum principle, and the case q = 0

The θ scaling breaks down at q = 0, so I first show that the solution is unique and that this case is
trivial.

$$w\equiv T_1-T_2:\quad \partial_t w = \kappa\,\partial_x^2 w,\quad w(0,t)=w(L,t)=0,\quad w(x,0)=0 \tag{S33}$$

S33 follows from linearity: subtracting the problems for two solutions in the class S6 cancels q.

$$\frac{dE}{dt} = \kappa\int_0^L w\,\partial_x^2 w\,dx,\qquad E=\tfrac12\int_0^L w^2dx \tag{S34}$$

S34 multiplies S33 by w, integrates over (0, L) and differentiates under the integral (as in S9). ∂_tw is
continuous on [0, L] × (0, ∞) by S6.

$$\frac{dE}{dt} = -\kappa\int_0^L(\partial_x w)^2dx\le0 \tag{S35}$$

S35 integrates by parts on [0, L]. This is valid because w, ∂_xw and ∂_x²w are continuous up to the ends
for t > 0 (S6). The boundary term [w ∂_x w]₀^L vanishes because w = 0 at both ends.

$$T_1\equiv T_2 \tag{S36}$$

S36 takes three steps:
1. **E(t) → 0 as t → 0⁺.** w is continuous on the compact set [0, L] × [0, 1] (S6), so it is uniformly
   continuous there. Hence max_x |w(x,t)| → max_x |w(x,0)| = 0, and E(t) ≤ (L/2) max_x w² → 0.
2. **E ≡ 0.** By S35 and the mean value theorem, E(t) ≤ E(s) for 0 < s < t. Letting s → 0⁺ gives E(t) ≤ 0.
   Since E ≥ 0, E ≡ 0.
3. **w ≡ 0.** A continuous w whose square integrates to zero vanishes identically.

$$q=0\ \Rightarrow\ T\equiv T_0 \tag{S37}$$

T ≡ T_0 solves S17 with S7 when q = 0, and it lies in the class S6. By S36 it is the only solution.

**Minimum principle.** Take q ≠ 0, so that θ is defined. Fix any τ_1 > 0.

$$\theta_{min}(\tau_1)\equiv\min_{[0,1]\times[0,\tau_1]}\theta=\theta(\xi_*,\tau_*) \tag{S48}$$

S48 is the extreme value theorem. θ = (T − T_0)/ΔT_q is continuous on the compact set [0, 1] × [0, τ_1],
because T is continuous on [0, L] × [0, ∞) (S6) and S20 is a continuous change of variables.

$$\theta_{min}(\tau_1)<0\ \Rightarrow\ 0<\xi_*<1,\ 0<\tau_*\le\tau_1,\ \partial_\xi^2\theta(\xi_*,\tau_*)\ge0,\ \partial_\tau\theta(\xi_*,\tau_*)\le0 \tag{S49}$$

S49 combines three operations:
- **Location.** θ = 0 on ξ = 0, ξ = 1 and τ = 0 (S26), so a negative minimum cannot be on those three
  sides.
- **Second-derivative test in ξ.** ξ_* is interior, and θ(·, τ_*) is C² near it (S6), so ∂_ξ²θ ≥ 0 there.
- **One-sided difference quotient in τ.** For 0 < h < τ_*,
  [θ(ξ_*, τ_*) − θ(ξ_*, τ_* − h)]/h ≤ 0, because θ(ξ_*, τ_*) is the minimum. ∂_τθ exists at (ξ_*, τ_*),
  since τ_* > 0 (S6). Letting h → 0⁺ gives ∂_τθ ≤ 0. This covers τ_* = τ_1 as well as τ_* < τ_1.

$$\theta(\xi,\tau)\ge0\quad\text{on }[0,1]\times[0,\infty) \tag{S50}$$

S50 is a proof by contradiction. If θ_min(τ_1) < 0, S49 gives ∂_τθ − ∂_ξ²θ ≤ 0 at (ξ_*, τ_*). But S25 says
that this combination equals 1 > 0 everywhere in the open set. So θ_min(τ_1) ≥ 0 for every τ_1, which is
S50. Because the source is strictly positive (+1), the usual ε-perturbation of the weak maximum principle is
not needed. The scratch series agrees: its minimum of θ over a grid is −1.8e-17, which is roundoff
(`m3_series.py`).

$$\theta \text{ is the same function for every } q\ne0;\ \ T-T_0=\Delta T_q\,\theta \text{ has the sign of } q \tag{S38}$$

- **First clause (uniqueness).** S25–S26 contain no parameter. Two solutions θ_1 and θ_2 map through S20 to
  two solutions of S17 with S7 in the class S6, and these are equal by S36.
- **Second clause (sign).** ΔT_q = qL²/k has the sign of q, because k > 0 (S1). And θ ≥ 0 (S50). So where
  θ > 0, T − T_0 has the sign of q.

**Left out by decision D-c (a): the upper bound θ ≤ 1/8** of the spec range (spec_problem.md, line 94). No
step uses it. It would follow from the same kind of argument applied to v = θ_s − θ (S51), with the
ε-perturbation, because v's equation has no source. The main session decided to leave it out of the ledger.

## 6. Boundary and wall conditions

$$F_0=-j(0,t)=k\,\partial_xT(0,t),\qquad F_L=j(L,t)=-k\,\partial_xT(L,t) \tag{S39}$$

Both ends are essential (Dirichlet) conditions, from S7, and no flux condition is imposed. The outward
fluxes F_0 and F_L are outputs of the solution. They are defined pointwise for t > 0 because ∂_xT is
continuous up to the ends (S6). The normal is n = −1 at x = 0 and n = +1 at x = L.

$$\text{wall (contact-angle) condition: none} \tag{S40}$$

The work order says to start from today's validated wall term, the `FACE_WALL` branch of `Residual_A1`
(src/assembly.c, lines 81–87). That term acts only on the phase-field row `Rb[a][0]`; the heat row `Rb[a][1]`
is 0 (line 100). The toy has no solid wall and no phase field (S8), so there is nothing to convert.

## 7. Steady state (the brief's comparison test)

$$\theta_s''=-1,\qquad\theta_s(0)=\theta_s(1)=0 \tag{S41}$$
$$\theta_s=-\tfrac{\xi^2}{2}+c_1\xi+c_0 \tag{S42}$$
$$\theta_s=\frac{\xi(1-\xi)}{2},\qquad \max\theta_s=\theta_s(\tfrac12)=\tfrac18 \tag{S43}$$

These come from:
- S41: setting ∂_τθ = 0 in S25;
- S42: integrating twice;
- S43: matching the end values, θ_s(0) = c_0 = 0 and θ_s(1) = −1/2 + c_1 = 0.

$$T_s(x)=T_0+\frac{q\,x(L-x)}{2k} \tag{S44}$$

S44 is the back-substitution of S20 into S43:
(qL²/k)(x/L)(1 − x/L)/2 = q x (L − x)/(2k).
This is exactly the brief's parabola (problem.md, line 17). For q = 0 it gives T_0, which agrees with S37.

$$F_0=F_L=\frac{qL}{2},\qquad F_0+F_L=qL=\int_0^Lq\,dx \tag{S45}$$

From S44, ∂_x T_s = q(L − 2x)/(2k). So F_0 = k · qL/(2k) and F_L = −k · (−qL/(2k)). The two ends together
carry away exactly the heat that the source puts in, which is a global energy check. This agrees with
spec_problem.md, line 100.

## 8. Approach to steady state and dominant terms

$$v\equiv\theta_s-\theta:\quad \partial_\tau v=\partial_\xi^2 v,\quad v(0,\tau)=v(1,\tau)=0,\quad v(\xi,0)=\frac{\xi(1-\xi)}{2} \tag{S51}$$

θ_s does not depend on τ, and θ_s″ = −1 (S41). So
∂_τv = −∂_τθ = −∂_ξ²θ − 1 = −∂_ξ²θ + θ_s″ = ∂_ξ²v.
The end values are 0 − 0 (S26, S41). The initial value is θ_s − 0 (S26, S43).

The sine coefficients of θ_s need two integrals, each by integration by parts:

$$\int_0^1\xi\sin(\ell\pi\xi)\,d\xi=-\frac{(-1)^\ell}{\ell\pi},\qquad \int_0^1\xi^2\sin(\ell\pi\xi)\,d\xi=-\frac{(-1)^\ell}{\ell\pi}+\frac{2((-1)^\ell-1)}{(\ell\pi)^3} \tag{S52}$$

The working, line by line:
- **First integral.** ∫ξ sin = [−ξ cos(ℓπξ)/(ℓπ)]₀¹ + (1/(ℓπ))∫cos(ℓπξ) = −(−1)^ℓ/(ℓπ) + 0.
- **Second integral.** ∫ξ² sin = [−ξ² cos(ℓπξ)/(ℓπ)]₀¹ + (2/(ℓπ))∫ξ cos(ℓπξ).
  - The intermediate integral is ∫₀¹ ξ cos(ℓπξ) dξ = [ξ sin(ℓπξ)/(ℓπ)]₀¹ − (1/(ℓπ))∫sin(ℓπξ)
    = ((−1)^ℓ − 1)/(ℓπ)².
  - Substituting it gives the second formula of S52.

`m3_series.py` evaluates both integrals symbolically for integer ℓ, and both agree with S52.

$$\beta_\ell\equiv2\int_0^1\theta_s\sin(\ell\pi\xi)\,d\xi=\frac{2(1-(-1)^\ell)}{(\ell\pi)^3}=\begin{cases}4/(\ell\pi)^3&\ell\text{ odd}\\0&\ell\text{ even}\end{cases} \tag{S53}$$

Since 2θ_s = ξ − ξ², β_ℓ is the first integral of S52 minus the second:
−(−1)^ℓ/(ℓπ) + (−1)^ℓ/(ℓπ) − 2((−1)^ℓ − 1)/(ℓπ)³ = 2(1 − (−1)^ℓ)/(ℓπ)³.
This agrees with selection.md D9 and with `m3_series.py`, which gives β_1 = 4/π³ and β_2 = 0.

$$v(\xi,\tau)=\sum_{\ell\text{ odd}}\frac{4}{(\ell\pi)^3}\sin(\ell\pi\xi)\exp(-\ell^2\pi^2\tau),\qquad \theta=\theta_s-v \tag{S54}$$

S54 is proved in four parts:
1. **Each term solves S51.** sin(ℓπξ) vanishes at both ends and satisfies ∂_ξ² sin = −ℓ²π² sin.
   `m3_series.py` confirms that each term satisfies the heat equation symbolically.
2. **Convergence (Weierstrass M-test).** Each term is bounded by β_ℓ, and Σβ_ℓ < ∞. So the series converges
   uniformly on [0, 1] × [0, ∞), and its sum is continuous there.
   - For τ ≥ τ_0 > 0, the differentiated series (∂_τ, ∂_ξ, ∂_ξ²) have terms bounded by a constant times
     ℓ⁻¹ exp(−ℓ²π²τ_0), which are summable.
   - So term-by-term differentiation is valid, and the sum is C^{2,1} on [0, 1] × (0, ∞).
3. **Initial value.** At τ = 0 the sum is the sine series of θ_s. The odd 2-periodic extension of θ_s is
   continuous (θ_s(0) = θ_s(1) = 0) and piecewise C¹, so its Fourier series converges to it uniformly.
4. **Identification.** θ_s minus the series solves S25–S26 and lies in the class S6, through S20. By the
   uniqueness clause of S38, it is θ.

Parts 2 and 3 also show that the class S6 contains a solution. So the hypothesis S6 excludes nothing that
the problem actually produces.

$$\frac{\partial\theta}{\partial\tau}(\xi,\tau)=\sum_{\ell\text{ odd}}\frac{4}{\ell\pi}\sin(\ell\pi\xi)\exp(-\ell^2\pi^2\tau)\quad(\tau>0) \tag{S55}$$

S55 differentiates S54 term by term, which part 2 of S54 justifies. Since ∂_τθ = −∂_τv, each coefficient is
ℓ²π² · 4/(ℓπ)³ = 4/(ℓπ). `m3_series.py` prints −∂_τ(term) = 4 e^{−π²ℓ²τ} sin(πℓξ)/(πℓ).

$$0<\tau\le0.1:\quad \max_\xi\partial_\tau\theta\ \ge\ \partial_\tau\theta(\tfrac12,\tau)\ \ge\ g(\tau)\ \ge\ 0.474 \tag{S56}$$

S56 has two parts:
- **Alternating series.** At ξ = 1/2, sin(ℓπ/2) = (−1)^r for ℓ = 2r + 1. So S55 is an alternating series
  whose terms a_r = 4/((2r+1)π) exp(−(2r+1)²π²τ) decrease strictly and tend to 0. By the Leibniz bound, the
  sum is at least a_0 − a_1 = g(τ) = (4/π)e^{−π²τ} − (4/(3π))e^{−9π²τ}.
- **Calculus on g.** g′(τ) = −4π e^{−π²τ} + 12π e^{−9π²τ}. It vanishes only at τ_g = ln 3/(8π²), because
  g′ = 0 means e^{8π²τ} = 3. It is positive before τ_g and negative after. So on (0, 0.1],
  g ≥ min(g(0⁺), g(0.1)).

The numbers, all from `m3_try3.py`:

| Quantity | Value |
|---|---|
| τ_g = ln 3 / (8π²) = 1.098612 / 78.956835 | 0.0139141 |
| g′ at τ_g/2, τ_g, 2τ_g, 0.1 | +8.59, 0 (to 2e-15), −6.37, −4.68 |
| g(0⁺) = 8/(3π) | 0.848826 |
| (4/π) e^{−π²·0.1} = 1.273240 × 0.372708 | 0.474546 |
| (4/(3π)) e^{−9π²·0.1} = 0.424413 × 1.38777e-4 | 5.88987e-5 |
| g(0.1) | 0.4744875 |

So min over (0, 0.1] of g is g(0.1) = 0.4744875 ≥ 0.474. `m3_series.py` gives the true value
∂_τθ(1/2, 0.1) = 0.4744875 (the two differ by 5e-12). The verifier's `m3_series.py`/`S56` check reports
the same minimum, 0.4744875.

$$\tau>0:\ \max_\xi|\partial_\tau\theta|\le\tfrac4\pi e^{-\pi^2\tau}\Big(1+\tfrac{e^{-8\pi^2\tau}}{3(1-e^{-8\pi^2\tau})}\Big),\ \ \max_\xi|\theta-\theta_s|\le\tfrac4{\pi^3}e^{-\pi^2\tau}\Big(1+\tfrac{e^{-8\pi^2\tau}}{27(1-e^{-8\pi^2\tau})}\Big);\quad \tau\ge1:\ 6.59\times10^{-5},\ 6.68\times10^{-6} \tag{S57}$$

S57 combines three operations:
- **The bound on |∂_τθ|.** Use the triangle inequality on S55 with |sin| ≤ 1, and factor out e^{−π²τ}.
  For ℓ = 2r + 1 ≥ 3, ℓ² − 1 = 4r(r+1) ≥ 8r and 1/ℓ ≤ 1/3. Then sum the geometric series
  Σ_{r≥1} e^{−8rπ²τ} = e^{−8π²τ}/(1 − e^{−8π²τ}).
- **The bound on |v| = |θ − θ_s|.** The same steps on S54, with 1/ℓ³ ≤ 1/27.
- **Values for τ ≥ 1.** Both bounds decrease in τ, so their values at τ = 1 hold for every τ ≥ 1.

The numbers at τ = 1, from `m3_try3.py`:

| Quantity | Value |
|---|---|
| e^{−π²} | 5.17232e-5 |
| (4/π) e^{−π²} | 6.58560e-5 ≤ 6.59e-5 |
| 4/π³ | 0.129006 |
| (4/π³) e^{−π²} | 6.67261e-6 ≤ 6.68e-6 |
| e^{−8π²}, the size of the bracket's departure from 1 | 5.12e-35 |

At τ = 1, `m3_series.py` finds the series maxima over ξ equal to these bounds to every printed digit: the
ℓ = 1 term dominates, and its maximum is at ξ = 1/2.

$$0<\tau\le0.1:\ \max_\xi\partial_\tau\theta\ge0.47,\ \max_\xi|\partial_\xi^2\theta|\ge1,\ \text{source}=1\ \Rightarrow\ \text{no term is negligible};\quad \tau\ge1:\ \max_\xi|\partial_\tau\theta|\le6.6\times10^{-5},\ \max_\xi|\theta-\theta_s|\le6.7\times10^{-6}\ \Rightarrow\ \text{near steady: conduction balances the source} \tag{S32}$$

S32 is the dominant-term statement, and its wording is the same as in the ledger. It is built from four
pieces:
- **Coefficients.** All coefficients of S25 are 1. In S29 they are 1, 1 and Q = 1 (S31).
- **Storage** comes from S56 for 0 < τ ≤ 0.1, and from S57 for τ ≥ 1.
- **Conduction.** θ(0, τ) = 0 for all τ (S26), and ∂_τθ is continuous up to ξ = 0 (S6). So ∂_τθ(0, τ) = 0,
  and S25 gives ∂_ξ²θ(0, τ) = −1 for every τ > 0. Hence max_ξ |∂_ξ²θ| ≥ 1.
- **Transition.** For 0.1 < τ < 1, storage decays like e^{−π²τ} (S57). S32 makes no claim there. For
  example, at τ = 0.3 the S57 storage bound is 0.0659 (`m3_try3.py`).

**Consequence for the tests (G:Derivation D-d).** The reference run ends at τ_end = 1. There the remaining
transient is at most 6.67e-6 in θ (S57). Relative to the peak 1/8, that is 5.34e-5 (`m3_try3.py`). Storage
is at most 6.59e-5 of the source. So a τ_end-versus-S44 tolerance tighter than about 1e-5 in θ would fail
for this reason alone, before any discretization error. The direct steady solve chosen at G:Selection
(V1 (c)) avoids this problem.

**Not shown:** an upper bound on storage for small τ. S32 claims only that storage is not negligible
(≥ 0.47). It does not claim that storage is ≤ O(1). The bound of S57 grows without limit as τ → 0⁺, so it
cannot supply one.

## 9. Thin-interface correction

$$\text{thin-interface correction: does not apply} \tag{S46}$$

The model has no diffuse interface and no width eps (S8). So there is no finite-width error in a kinetic
coefficient to correct.

## 10. Sharp-interface limit to first order in the widths

$$\lim_{eps,\,eps_s\to0}\text{S15}=\text{S15};\quad O(eps),\ O(eps_s)\ \text{corrections}\equiv0 \tag{S47}$$

No interface width enters S15 (S8). The model is its own sharp limit, and every first-order correction is
identically zero.

## Claims from papers to confirm at the gate

None. The ledger cites no paper, as the brief requires. Every source is problem.md, spec_problem.md,
selection.md, decisions.json, or src/ files read in this stage.

## What I could not show, or did not do

- **Existence is shown only through the series (S54).** It relies on three standard theorems that are named
  but not proved here: the Weierstrass M-test, term-by-term differentiation, and uniform convergence of the
  Fourier series of a continuous, piecewise C¹ periodic function.
- **The upper bound θ ≤ 1/8** of the spec range is not in the ledger, by decision D-c (a) (see §5).
- **The storage upper bound for small τ** is not shown (see §8).
- **The g′ sign pattern in S56** is shown analytically: there is one root, at e^{8π²τ} = 3. The scratch run
  only samples it at four points, so the samples support the claim but do not prove it.
- **C1 and C2.** The main session approved them for plumbing only, and decisions.json records that this was
  not a user decision. The G:Derivation notes say that they stand for the toy, and I treat them as decided.
- **Scratch runs are mine only.** `m3_series.py`, `m3_try3.py` and `m3_ledger_check.py` were run at try 3.
  Their output supports my own confidence in the numbers. It is not evidence for the gate. The verifier's
  results_M3.json (try 2) already agreed on 0.4744875, 6.5856e-5 and 6.6726e-6. It has not yet been rerun on
  this revision.

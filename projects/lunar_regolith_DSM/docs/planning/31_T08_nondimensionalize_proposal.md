# 31. T08: nondimensionalize the coupled system (proposal for approval)

Date: 2026-11-19. Status: **proposed, nothing built or run.** Step 1 of the Lab Report 4 plan
(`project/LOG.md`, 2026-11-17 to 2026-11-19 entry).

Tags: **[checked]** verified by a script or in a file, with the place given; **[inferred]** my reasoning.

## 1. What T08 asks and what it feeds

The agent gets the three coupled equations of the current model, written with general material functions,
and nondimensionalizes them. It chooses the scales, writes the system in dimensionless variables, collects
the parameters into independent dimensionless groups, and evaluates each group at the reference state
(−20 °C, the sintering pair). It keeps every term; it reports which terms are small but does not drop them.

What T08 feeds:

- **T09 (the thin-interface limit).** T09 expands in $\varepsilon/\ell$ and needs the groups that set the
  kinetic and capillary terms. T08 fixes those groups and their sizes first.
- **The held T05 test (doc 30 §1).** See §6: in dimensionless form, $\rho_{vs}/\rho_i$ becomes a group
  $\beta$, which the grader can sample at random $O(1)$ values. An error of relative size $\beta$ then has
  nowhere to hide.
- **The M&F time scaling in the code** (`xi_v = 1e-3`, `src/lunar_main.c:163-170`). That scaling is
  justified only if vapor diffusion is much faster than interface motion. T08 measures the ratio.

## 2. The given system

These are the code's continuous equations (`src/assembly.c:177-191`), with $\xi_T = \xi_v = 1$ and no
clamps or floors (`air_lim`). Three items are left general on purpose, so that T08 does not prejudge a
decision:

| Kept general | Why |
|---|---|
| $C(\phi, T)$ | Doc 30 §5: the mixture rule against the code's product $\rho(\phi)c_p(\phi)$ is your open decision. |
| $L(T)$ | T04 derived $L(T) = h_v - h_i$; the code uses a constant. |
| $\alpha(T)$, with shape $\phi^2(1-\phi)^2$ | T09 decides the coupling; the three T07 runs gave constant $\alpha$, $\alpha \propto M$ and $\alpha(T)$. |

$$R_\phi = \partial_t\phi + M\Big[\tfrac{6\gamma}{\varepsilon}\phi(1-\phi)(1-2\phi) - 6\gamma\varepsilon\nabla^2\phi\Big]
- \alpha(T)\,\phi^2(1-\phi)^2\big(\rho_v - \rho_{vs}(T)\big)$$

$$R_T = C(\phi,T)\,\partial_tT - \nabla\cdot\big(k(\phi,T)\nabla T\big) - \rho_iL(T)\,\partial_t\phi$$

$$R_\rho = (1-\phi)\,\partial_t\rho_v - \nabla\cdot\big(D_v(T)(1-\phi)\nabla\rho_v\big) + (\rho_i - \rho_v)\,\partial_t\phi$$

**[checked]** All three parse with the graders' parser (`tasks/_shared/blind.py`), and every term of each
has one unit: $\mathrm{s^{-1}}$, $\mathrm{W\,m^{-3}}$, $\mathrm{kg\,m^{-3}\,s^{-1}}$.
$R_\phi$ is D1's form; $M = M_{\text{code}}/(2\gamma)$ and $\alpha = \alpha_{\text{code}}/\rho_i$ map it
onto A1 (doc 30 §2–3).

## 3. Reference values

These are the values of the reference run (`2026-11-18__15.21.44_dtlad60_T-20_sealed_ac1e-3_dt2tau_cmpD1`,
its `outp.txt` and opts files, and the code defaults at `src/lunar_main.c:128-139`). A script recomputed
each derived value: `t08_sizes.py`, kept in the session scratchpad, not in the repository.

| Symbol | Value | Source |
|---|---|---|
| $M$ | $2.24\times10^{-10}$ m$^3$ J$^{-1}$ s$^{-1}$ | $M_{\text{code}}/(2\gamma)$; printed 2.2397e-10 [checked] |
| $\gamma$, $\varepsilon$ | 0.109 J m$^{-2}$, $8.584\times10^{-7}$ m | outp.txt |
| $R_g$ (grain radius) | $1.5\times10^{-5}$ m | geometry opts (`-ice_grain_R`) |
| $\rho_i$, $T_0$ | 919 kg m$^{-3}$, 253.15 K | defaults; `-temp -20` |
| $\alpha_0 = \alpha_{\text{code}}/\rho_i$ | 850 m$^3$ kg$^{-1}$ s$^{-1}$ | $7.812\times10^5/919$ |
| $C_i = \rho_ic_i$, $C_g = \rho_ac_a$ | $1.801\times10^6$, 1400 J m$^{-3}$ K$^{-1}$ | defaults |
| $k_i$, $k_g$ | 2.29, 0.02 W m$^{-1}$ K$^{-1}$ | defaults |
| $L_0$ | $2.83\times10^6$ J kg$^{-1}$ | default |
| $D_0 = D_v(T_0)$ | $1.898\times10^{-5}$ m$^2$ s$^{-1}$ | `VaporDiffus`, $(T/273.15)^{1.81}$ law |
| $\rho_{vs,0}$ | $8.487\times10^{-4}$ kg m$^{-3}$ | `RhoVS_I`; $\rho_i/\rho_{vs}$ printed 1.0829e6 [checked] |
| $b_0 = d\ln\rho_{vs}/dT$ | 0.0961 K$^{-1}$ | `RhoVS_I`, central difference |
| $\sigma_0$ | $10^{-5}$ | typical supersaturation of the runs [inferred from experiment names] |

Sizes I computed for orientation (the agent is not shown them) [checked by script]:

- $\beta = \rho_{vs,0}/\rho_i = 9.2\times10^{-7}$;
- with $\ell = R_g$, time to diffuse vapor across a grain: $1.2\times10^{-5}$ s;
- heat across a grain: $1.8\times10^{-4}$ s through ice, $1.6\times10^{-5}$ s through gas;
- interface relaxation, $\varepsilon/(6\gamma M) = 5860$ s, which equals the code's printed
  `tau_sub` = 5860.4 s;
- vapor-to-interface time ratio: $2\times10^{-9}$.

**Correction to doc 30 §1.** The blind spot of the T05 grader is $\rho_v/\rho_i$. Doc 30 put it at
$10^{-9}$, using 200 K values. At −20 °C it is $9.2\times10^{-7}$, just under the grader's $10^{-6}$
tolerance. The T05 grader still samples at 200 K sizes, so its test result does not change. The margin
at our working temperature, however, is a factor of about 1.1, not 1000.

## 4. The work order (draft `task.json`)

Two things are new compared with T02–T07.

**(a) The answer is a system.** `result_expr` is a tuple `(R_phi, R_T, R_u)`.

**(b) Every scale, group and dimensionless function needs a machine-readable definition.** The grader
substitutes these definitions back. This needs schema changes, which are decision 1 in §7:

- `result.schema.json`: an optional `definition` field (a SymPy string) on each entry of `symbols`.
- `task.schema.json`: a new deliverable kind, `sympy_system`, and an optional `value` field on work-order
  symbols. The reference values then live in one place, the work order, which the grader also reads.
  The symbol registry will need the same field later.
- One line in the definition, `beth-mathematician.md` (lines 3 and 21): "one SymPy expression, a tuple
  when the work order asks for a system".

I chose not to restrict the vapor map to an affine one: the agent may write $\rho_v$ in terms of
$u$ and $\theta$, for example $\rho_v = \rho_{vs}(T)(1+\sigma_0u)$. Length and time scales are its choice.

```json
{
  "task_id": "T08_nondimensionalize_system",
  "goal": "nondimensionalize",
  "statement": "Below are the three coupled equations of the current phase-field model of ice and water vapor in the pore space of lunar regolith: the phase field phi (ice fraction), the temperature T and the vapor density rho_v. You derived each of them in earlier tasks; here they are written with general material functions so that no choice of material law or coupling is built in. Nondimensionalize the coupled system. (1) Choose one length scale ell for all three coordinates, one time scale t_c for all three equations, and a map for each field: x = ell*X (and the same for y, z), t = t_c*s, T = T_r + DT*theta(X, Y, Z, s), and rho_v written in terms of a dimensionless vapor variable u(X, Y, Z, s) (and theta, if you wish). Define every scale in terms of the work order's parameters. (2) Write each material function as a reference value times a dimensionless function of phi and theta. (3) Substitute, apply the chain rule, and divide each equation by one factor so that every term is dimensionless. Keep every term: do not drop any term, however small; a term's size is reported, not acted on. (4) Collect the parameters into dimensionless groups. Each group must be independent of the others (none is a product of powers of the others). (5) Evaluate every group, and every time scale you identify, at the reference state given by the symbols' values, and say which processes are fast and which slow, and which terms are small. Return the three dimensionless residuals as a tuple (R_phi, R_T, R_u), in that order, each in 3D Cartesian coordinates X, Y, Z and time s.",
  "symbols": [
    {"name": "x", "latex": "x", "meaning": "first position coordinate", "units": "m"},
    {"name": "y", "latex": "y", "meaning": "second position coordinate", "units": "m"},
    {"name": "z", "latex": "z", "meaning": "third position coordinate", "units": "m"},
    {"name": "t", "latex": "t", "meaning": "time", "units": "s"},
    {"name": "phi(x, y, z, t)", "latex": "\\phi", "meaning": "local volume fraction of ice (phase field), between 0 and 1", "units": "1"},
    {"name": "T(x, y, z, t)", "latex": "T", "meaning": "temperature", "units": "K"},
    {"name": "rho_v(x, y, z, t)", "latex": "\\rho_v", "meaning": "vapor density: mass of water vapor per unit volume of gas", "units": "kg m^-3"},
    {"name": "M", "latex": "M", "meaning": "interface mobility, constant", "units": "m^3 J^-1 s^-1", "value": 2.24e-10},
    {"name": "gamma", "latex": "\\gamma", "meaning": "energy per unit area of the ice-gas interface, constant", "units": "J m^-2", "value": 0.109},
    {"name": "eps", "latex": "\\varepsilon", "meaning": "the length that sets the width of the diffuse interface, constant", "units": "m", "value": 8.584e-7},
    {"name": "R_g", "latex": "R_g", "meaning": "radius of an ice grain: the size of the ice structure", "units": "m", "value": 1.5e-5},
    {"name": "rho_i", "latex": "\\rho_i", "meaning": "density of ice, constant", "units": "kg m^-3", "value": 919},
    {"name": "T0", "latex": "T_0", "meaning": "reference temperature", "units": "K", "value": 253.15},
    {"name": "alpha(T)", "latex": "\\alpha(T)", "meaning": "coefficient of the growth term: rate of change of phi per unit vapor excess rho_v - rho_vs(T), per unit volume of vapor", "units": "m^3 kg^-1 s^-1"},
    {"name": "C(phi, T)", "latex": "C(\\phi, T)", "meaning": "volumetric heat capacity of the ice-gas mixture", "units": "J m^-3 K^-1"},
    {"name": "k(phi, T)", "latex": "k(\\phi, T)", "meaning": "thermal conductivity of the ice-gas mixture", "units": "W m^-1 K^-1"},
    {"name": "L(T)", "latex": "L(T)", "meaning": "latent heat of sublimation per unit mass", "units": "J kg^-1"},
    {"name": "D_v(T)", "latex": "D_v(T)", "meaning": "diffusivity of water vapor in the pore gas", "units": "m^2 s^-1"},
    {"name": "rho_vs(T)", "latex": "\\rho_{vs}(T)", "meaning": "saturation vapor density over a flat ice surface at temperature T", "units": "kg m^-3"},
    {"name": "alpha_0", "latex": "\\alpha_0", "meaning": "alpha(T0)", "units": "m^3 kg^-1 s^-1", "value": 850.0},
    {"name": "C_i", "latex": "C_i", "meaning": "C(1, T0), pure ice", "units": "J m^-3 K^-1", "value": 1.801e6},
    {"name": "C_g", "latex": "C_g", "meaning": "C(0, T0), pure gas", "units": "J m^-3 K^-1", "value": 1400},
    {"name": "k_i", "latex": "k_i", "meaning": "k(1, T0), pure ice", "units": "W m^-1 K^-1", "value": 2.29},
    {"name": "k_g", "latex": "k_g", "meaning": "k(0, T0), pure gas", "units": "W m^-1 K^-1", "value": 0.02},
    {"name": "L_0", "latex": "L_0", "meaning": "L(T0)", "units": "J kg^-1", "value": 2.83e6},
    {"name": "D_0", "latex": "D_0", "meaning": "D_v(T0)", "units": "m^2 s^-1", "value": 1.898e-5},
    {"name": "rho_vs0", "latex": "\\rho_{vs,0}", "meaning": "rho_vs(T0)", "units": "kg m^-3", "value": 8.487e-4},
    {"name": "b_0", "latex": "b_0", "meaning": "d ln rho_vs / dT at T0 (Clausius-Clapeyron slope)", "units": "K^-1", "value": 0.0961},
    {"name": "sigma_0", "latex": "\\sigma_0", "meaning": "typical supersaturation (rho_v - rho_vs)/rho_vs in the runs", "units": "1", "value": 1e-5}
  ],
  "given": [
    "Derivative(phi(x, y, z, t), t) + M*(6*gamma/eps*phi(x, y, z, t)*(1 - phi(x, y, z, t))*(1 - 2*phi(x, y, z, t)) - 6*gamma*eps*(Derivative(phi(x, y, z, t), (x, 2)) + Derivative(phi(x, y, z, t), (y, 2)) + Derivative(phi(x, y, z, t), (z, 2)))) - alpha(T(x, y, z, t))*phi(x, y, z, t)**2*(1 - phi(x, y, z, t))**2*(rho_v(x, y, z, t) - rho_vs(T(x, y, z, t)))",
    "C(phi(x, y, z, t), T(x, y, z, t))*Derivative(T(x, y, z, t), t) - Derivative(k(phi(x, y, z, t), T(x, y, z, t))*Derivative(T(x, y, z, t), x), x) - Derivative(k(phi(x, y, z, t), T(x, y, z, t))*Derivative(T(x, y, z, t), y), y) - Derivative(k(phi(x, y, z, t), T(x, y, z, t))*Derivative(T(x, y, z, t), z), z) - rho_i*L(T(x, y, z, t))*Derivative(phi(x, y, z, t), t)",
    "(1 - phi(x, y, z, t))*Derivative(rho_v(x, y, z, t), t) - Derivative(D_v(T(x, y, z, t))*(1 - phi(x, y, z, t))*Derivative(rho_v(x, y, z, t), x), x) - Derivative(D_v(T(x, y, z, t))*(1 - phi(x, y, z, t))*Derivative(rho_v(x, y, z, t), y), y) - Derivative(D_v(T(x, y, z, t))*(1 - phi(x, y, z, t))*Derivative(rho_v(x, y, z, t), z), z) + (rho_i - rho_v(x, y, z, t))*Derivative(phi(x, y, z, t), t)"
  ],
  "assumptions": [
    "the three given residuals are exact: R_phi = 0, R_T = 0 and R_rho = 0 hold together",
    "M, gamma, eps and rho_i are constant; alpha, C, k, L, D_v and rho_vs are general positive functions of their arguments",
    "the values given are those of the current model at -20 C (the reference state); a symbol with no value is a function"
  ],
  "deliverable": {
    "kind": "sympy_system",
    "allowed_names": ["X", "Y", "Z", "s", "phi", "theta", "u", "Derivative"],
    "notes": "result_expr is the tuple (R_phi, R_T, R_u). It may use only X, Y, Z, s, the fields phi(X, Y, Z, s), theta(X, Y, Z, s) and u(X, Y, Z, s), Derivative, numbers, + - * / ** and parentheses, and the dimensionless groups and functions you declare; no dimensional symbol of the work order may appear in it. Declare in your symbols every scale (ell, t_c, T_r, DT, and those of your vapor map, with the map rho_v itself as the symbol rho_v_map), every group and every dimensionless function, each with a definition: a SymPy expression in the work order's symbols (for a dimensionless function, for example khat(phi, theta) with definition k(phi, T_r + DT*theta)/k_i). Give each group and each time scale its typical value at the reference state."
  },
  "hand_estimate": null
}
```

Left out of scope, so that T08 stays small: boundary conditions, and the reduced (leading-order) system.
Dropping terms is a "reduce" task. The natural place for it is T09, which needs the expansion anyway.

## 5. The grader: checks with no answer key

The given system is the reference. The agent sees it too, so the grader holds nothing the agent lacks.
What the grader checks is that the agent's transformation is exact and that its bookkeeping is honest.
All of G1–G6 must pass.

| Check | What it does | What it catches |
|---|---|---|
| **G1 names** | `result_expr` parses as a 3-tuple. It uses only $X, Y, Z, s$, $\phi, \theta, u$ and the agent's declared dimensionless names. No dimensional symbol of the work order appears. | A dimensional parameter left in an equation; a missing equation. |
| **G2 units** | Each scale has the units of what it scales ($\ell$: m, $t_c$: s, $T_r$ and $\Delta T$: K, the vapor map: kg m$^{-3}$). Each group and each dimensionless function's definition is dimensionless (pint). | A group that is not dimensionless, for example $\ell$ where $\ell^2$ belongs. |
| **G3 identity** (the main check) | For each equation $k$: put the agent's map and definitions into the given $R_k$, with random $O(1)$ parameter values (not physical ones), random smooth positive material functions and three random polynomial test fields. Require $\hat R_k = c_k R_k$ with one nonzero $c_k$ for every field and every point of a draw. 3 draws. | A wrong chain-rule factor or group; any dropped, doubled or sign-flipped term, of any physical size; an equation divided by something that depends on the fields; different time scales in different equations. |
| **G4 invertible map** | $\ell$, $t_c$, $\Delta T > 0$ at the reference values, and $\partial\rho_v/\partial u \neq 0$. | A map that loses information. |
| **G5 independent groups** | The Jacobian of ($\ln$ of each group) with respect to ($\ln$ of each parameter), at the reference values, has full row rank. The number of groups is at most the number of parameters minus the rank of their dimension matrix (Buckingham). | A redundant group, for example $\Pi_3 = \Pi_1\Pi_2$ declared separately. |
| **G6 sizes** | Each group's and each time scale's stated typical value equals its definition, evaluated at the work order's values, within 5 %. | Arithmetic errors in the size estimates that the user and T09 will rely on. |

**[checked]** I prototyped G3 on the vapor equation (1D, in the session scratchpad, not committed).
A hand-written correct nondimensionalization gives a constant ratio over three test fields on every draw.
The same answer with the exchange counted twice ($(1/\beta - 2u)\,\partial_s\phi$) gives a ratio that
changes from field to field, so it fails.

**Reported for you, not pass/fail:**

- the groups and their sizes;
- the agent's time scales, against the code's printout (`tau_sub`, $\rho_i/\rho_{vs}$);
- the terms the agent calls small, and the group that makes each one small;
- the coefficient of each equation's own rate term;
- the vapor-to-interface time ratio next to the code's $\xi_v = 10^{-3}$.

**Grader self-test** (`key_check` and `tests/test_grade_T08.py`, the same pattern as T07). I write one
reference nondimensionalization by hand. It must pass G1–G6. Each of these planted errors must fail at
least one check:

1. exchange counted twice, $(1/\beta - 2u)$ (G3);
2. $u\,\partial_s\phi$ dropped as small (G3; the system must stay exact);
3. diffusion group with $\ell$ for $\ell^2$ (G2);
4. vapor equation divided by $(1-\phi)$ (G3);
5. $\rho_i$ left in $\hat R_u$ (G1);
6. a redundant group (G5);
7. one typical value off by a factor of 10 (G6);
8. latent-heat sign flipped (G3);
9. two equations with different time scales (G3);
10. curvature coefficient wrong by a factor of 2 (G3).

## 6. How T08 settles the held T05 test

T05 is dimensional, and the grader evaluates it at physical sizes. There, an error of size
$\rho_v/\rho_i$ is indistinguishable from the approximation an agent may legitimately make. In T08, the
same ratio is the group $\beta$. G3 samples $\beta$ at random $O(1)$ values and requires an exact
identity, so the error is caught outright (prototype above). Dropping the term is not allowed in the
graded system; it may only be reported as small.

This supports **option (b)** of doc 30 §1. The T05 test stays, inverted: it asserts that the T05 grader
passes this error. Its docstring says why, and names T08's G3 and planted error 1 as the check that
catches it. Your decision; see §7.

## 7. Decisions for you

1. **Schema changes** (§4b): the optional `definition` on answer symbols, the `sympy_system` kind, the
   optional `value` on work-order symbols, and one line in the definition. These are edits to the agent's
   folder. *Feeds:* the grader cannot check a system without them. My recommendation: approve.
2. **Renderer** (`scripts/agentkit/render_derivation.py`): show the definition column and the three
   equations. The file is in the barrier fingerprint, so the launcher would require a red-team rerun
   (one launcher call, about a minute) before the first T08 call. *Feeds:* readability of
   `derivation.md` only; G1–G6 do not need it. Optional.
3. **The given system** (§2): the code's equations with general $C$, $L(T)$ and $\alpha(T)$.
   *Feeds:* doc 30 §5 and T09 stay open. My recommendation: approve as written.
4. **Length scale:** the agent's choice (recommended; it is part of the task), or prescribe
   $\ell = R_g$.
5. **The held T05 test:** option (b) as in §6, or another option from doc 30 §1.
6. **The reference values** (§3): check the table, in particular $\sigma_0 = 10^{-5}$. I inferred it from
   the experiment names; it is not a model input.

Still open from doc 30, not needed for T08: the mixture heat capacity (§5) and removing the clamps
(smooth floors). Both come due at step 3, when the agent writes C code. The temperature-gradient
comparison stays on hold.

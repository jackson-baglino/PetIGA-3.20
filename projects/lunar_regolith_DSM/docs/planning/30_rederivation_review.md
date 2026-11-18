# 30. Review of the blind rederivation (T04–T07): answers, the code comparison, and the derivation format

Date: 2026-11-18. This is the long form of the answers to your review of T04–T07. The comparison of the
derived model against the code is in `studies/derived_model/README.md`.

Tags: **[checked]** means verified by a script or in a file, with the place given. **[inferred]** means my
reasoning, not verified. **[unverified]** means it needs a source I have not read.

---

## 1. The held T05 test, in plain terms

**What the grader does.** The T05 grader has no answer key. It runs physical checks on the agent's vapor
equation. One of them, the "exchange" check, asks: when $\phi$ changes, does the equation move water between
ice and vapor at one fixed ice density $\rho^*$? It computes the $\rho^*$ the equation implies at three
states and requires them to agree to a relative $10^{-6}$.

**What the tests of the grader do.** `tests/test_grade_T05.py` feeds the grader equations that are wrong on
purpose ("planted errors") and asserts that it fails each one. This is how we know the grader can catch
mistakes.

**The planted error that is not caught.** `test_exchange_counted_twice` writes the storage term as
$\partial_t[(1-\phi)\rho_v]$ and also keeps the exchange term $(\rho_i-\rho_v)\,\partial_t\phi$. Expanded,
that is

$$(1-\phi)\,\partial_t\rho_v + (\rho_i - 2\rho_v)\,\partial_t\phi,$$

so the $\rho_v\,\partial_t\phi$ term is counted twice. The error is $\rho_v\,\partial_t\phi$, and next to
$\rho_i\,\partial_t\phi$ it is $\rho_v/\rho_i \approx 10^{-6}/920 \approx 10^{-9}$. That is below the
$10^{-6}$ tolerance, so the grader passes the equation and the test fails.

**Why the tolerance is $10^{-6}$ and not 0.** An agent may legitimately drop $\rho_v\,\partial_t\phi$
because it is $10^{-9}$ of the term next to it: then it writes $\rho_i\,\partial_t\phi$. That approximation
is the same size as the error, so no tolerance can pass the first and fail the second. **The grader is
blind to anything of size $\rho_v/\rho_i$ in that term.** This blind spot is real.

**The three options:**

- **(a) Delete the test.** The blind spot is then recorded nowhere.
- **(b) Keep it, inverted (my recommendation).** It asserts that the grader *passes* this error, and its
  docstring says why. That records the blind spot where a future reader of the grader will see it. What
  catches this error is the side-by-side comparison with the code, plus your review. If someone later makes
  the grader stricter, this test fails and reminds them that the approximation must stay allowed.
- **(c) Make the tolerance exact.** The grader then fails an agent that correctly drops a $10^{-9}$ term.

## 2. Step S3 of T07: the driving force

**Yes, this is the phase-change term.** $\Delta g$ is the energy per unit volume of ice released when vapor
at density $\rho_v$ deposits as ice. The derivation, step by step:

1. **Coexistence.** At a flat interface at temperature $T$, ice and vapor coexist when
   $\rho_v = \rho_{vs}(T)$. There the chemical potentials (Gibbs energy per unit mass) are equal:
   $\mu_{\text{ice}}(T) = \mu_{\text{vap}}(T, p_{s})$.
2. **Ideal vapor.** At fixed $T$, $d\mu = dp/\rho = R_vT\,dp/p$, so
   $\mu_{\text{vap}}(T,p) = \mu_{\text{vap}}(T,p_s) + R_vT\ln(p/p_s)$. With $p = \rho_vR_vT$,
   $p/p_s = \rho_v/\rho_{vs}$.
3. **Incompressible ice.** $\mu_{\text{ice}}$ does not change with the small pressure difference, so it
   stays equal to $\mu_{\text{vap}}(T, p_s)$.
4. **Energy released per kilogram deposited:**
   $\mu_{\text{vap}} - \mu_{\text{ice}} = R_vT\ln(\rho_v/\rho_{vs})$. **Per unit volume of ice**, multiply by
   $\rho_i$:

$$\Delta g = \rho_iR_vT\ln\frac{\rho_v}{\rho_{vs}(T)} \approx \rho_iR_vT\,\frac{\rho_v-\rho_{vs}(T)}{\rho_{vs}(T)}.$$

5. **Linearize.** The last step uses $\ln(1+\sigma)\approx\sigma$, where
   $\sigma = (\rho_v-\rho_{vs})/\rho_{vs}$ is the supersaturation. The relative error is about $\sigma/2$.
   Our runs have $\sigma$ of order $10^{-5}$ (for example the `sigp1e-5` experiments), so the error is
   negligible. [inferred from the experiment names; check against your runs]

**How it enters the equation.** $G = F - \int\Delta g\,h(\phi)\,dV$. With $\Delta g > 0$ (supersaturated),
raising $\phi$ where $h' > 0$ lowers $G$, so the interface advances. With $\Delta g < 0$, it sublimates.

**A check you can do by hand [checked].** Set $v_n = 0$ in S8 for a curved interface, $\Delta g = \gamma K$.
That gives

$$\sigma = \frac{\gamma}{\rho_iR_vT}\,K = d_0K,\qquad d_0 = \frac{\gamma}{\rho_iR_vT} = \frac{\gamma\,(m_w/\rho_i)}{k_BT}.$$

This is the linearized Kelvin equation, and $d_0$ is exactly the code's capillary length (`D0_CAPILLARY`,
`include/NASA_types.h:15`, which cites Kaempfer & Plapp Eq. 13). So the agent's driving force reproduces
K&P's capillary length.

**How K&P get their term [inferred from the code, unverified against the paper].** The code builds the
coupling as $\lambda = a_1\varepsilon/d_{0,\text{sub}}$ and $\alpha = \lambda/\tau$
(`src/lunar_main.c`, kinetic setup, near `lambda_sub = a1 * user.eps / d0_sub`). That is the Karma–Rappel
construction: the coupling is written as $\lambda\,g'(\phi)\,u$, with $u$ a dimensionless supersaturation,
and $\lambda$ is then fixed by a thin-interface (matched-asymptotic) analysis so that the sharp-interface
limit gives Gibbs–Thomson with the right $d_0$. It is fixed afterwards, not derived from a free energy. To
confirm, read K&P Sec. II: their $\phi$ equation, their Eq. 1 for $u$, and the $\lambda$–$d_0$ relation.

**The two routes give the same coupling [checked].** `studies/derived_model/verification/coupling_map.py`
proves with SymPy that, with $M = M_{\text{code}}/(2\gamma)$ (where $M_{\text{code}}$ is `mob_sub`), the
agent's coupling divided by the code's $\alpha/\rho_i$ is

$$\frac{5}{a_1}\,\frac{T}{T_0}\,\frac{\rho_{vs}(T_0)}{\rho_{vs}(T)}.$$

The code uses $a_1 = 5$, so this is exactly 1 at $T = T_0$. The only difference: the code freezes
$T/\rho_{vs}(T)$ at the initial temperature, while the agent's term uses the local $T$.

## 3. Your relation for S6

You proposed $\alpha_{\text{sub}} = 30M\rho_iR_vT$. The checked relation is

$$\alpha_{\text{sub}} = \frac{30\,M\,\rho_i^2R_vT_0}{\rho_{vs}(T_0)}.$$

The extra factor $\rho_i/\rho_{vs}$ appears because the code writes its term as
$(\alpha/\rho_i)(\rho_v-\rho_{vs})$, while the agent writes it as $(\rho_v-\rho_{vs})/\rho_{vs}$.

On $h(\phi)$ and scaling $M$: agreed. One $M$ sets both relaxation and growth, so scaling $M$ scales the
whole phase-change rate. The code does the same through $\tau$, which sets both `mob_sub` and `alph_sub`
(`lunar_main.c`, comment above the `-beta_sub0` override).

## 4. S8 values

Noted: $\alpha_c = 1$ is never used, 200 K is too cold for our runs, and $\varepsilon = 1\,\mu$m is at the
upper end. The estimate only illustrates the formula. The working value of $M$ comes from the code's
calibration ($M = M_{\text{code}}/(2\gamma) = 2.24\times10^{-10}$ m$^3$ J$^{-1}$ s$^{-1}$ at
$-20\,^\circ$C, $\alpha_c = 10^{-3}$, printed by the D1 runs).

## 5. The heat capacity in T04: why the mixture rule

**It was not a choice; it follows from the definition of energy [checked in the T04 derivation, S1 and
S7].** The enthalpy per unit volume of a mixture is the sum over phases of volume fraction × density ×
specific enthalpy, $H = \phi\rho_ih_i + (1-\phi)\rho_ah_a + \dots$. Its derivative in $T$ gives

$$C = \frac{\partial H}{\partial T} = \phi\rho_ic_i + (1-\phi)\rho_ac_a.$$

**What the code does instead.** It computes $\rho(\phi)\,c_p(\phi)$, with both $\rho$ and $c_p$
volume-averaged (`material_properties.c`, `HeatCap` and `Density`). A specific heat (per kilogram) must be
averaged by **mass**, not by volume: $c_{\text{mix}} = \sum_k\phi_k\rho_kc_k/\sum_k\phi_k\rho_k$. Then
$\rho\,c_{\text{mix}}$ is the mixture rule above. The two differ by

$$\rho(\phi)c_p(\phi) - C = -\phi(1-\phi)(\rho_i-\rho_a)(c_i-c_a) \;<\; 0.$$

At $\phi = 1/2$, with the code's values, A1's value is 0.767 of the mixture value (691 176 against
901 320 J m$^{-3}$ K$^{-1}$; `coupling_map.py`).

**What the mixture rule buys:**

1. **Energy is booked correctly.** The heat stored in the interface band is the heat that is there. The
   product rule stores up to 23 % too little at $\phi = 1/2$.
2. **Consistency.** $C$ and the latent term $\rho_iL\,\partial_t\phi$ come from the same energy balance, so
   the continuous equations conserve energy. The product rule has no derivation to stand on.
3. **It costs nothing:** the same number of operations.

**What it does not change.** Both rules agree in pure ice and pure gas. They differ only inside the diffuse
band, a few $\varepsilon$ wide, so the effect shrinks as $\varepsilon \to 0$ and is small when $T$ changes
slowly. The run comparison in `studies/derived_model/README.md` measures it on one case. My recommendation:
use the mixture rule in the new model.

## 6. A fixed format for `derivation.md`

**The problem.** The three files you read have three styles: LaTeX in T07, ASCII in T04, and symbols that
are not defined. The agent writes the file freehand, and the definition only says "one displayed equation
per step".

**Proposal (recommended): the agent fills structured steps, and a script writes `derivation.md`.**
Consistency then comes by construction, not by asking.

1. **Schema** (`result.schema.json`). Each entry of `steps` gets fields:
   - `id`;
   - `title`;
   - `why`: one or two sentences saying what this step does and why it is the next step;
   - `from`: the steps, assumptions or given items it uses;
   - `operation`: for example "integrate by parts" or "substitute S3 into S5";
   - `latex`: the equation, in LaTeX;
   - `check`: the check function in `scratch/check.py` and what it printed.

   Each entry of `symbols` gets a `latex` field (for example `\rho_{vs}(T)`), next to its meaning, units and
   typical value. `assumptions` become numbered A1…An.
2. **Renderer** (`scripts/agentkit/render_derivation.py`, run by the launcher after the call). It writes
   `derivation.md` in one fixed layout:
   - problem;
   - symbol table (LaTeX, meaning, units, typical value);
   - assumptions;
   - steps, each as: the words, then the displayed equation, then "uses: …", then "checked by: …";
   - result;
   - neglected effects, with reasons.
3. **Lint**, which the grader reports as a format check:
   - every step has words and an equation;
   - every LaTeX symbol in an equation is in the symbol table;
   - no Unicode Greek letters and no ASCII names in prose.
4. **Definition.** Step 4 changes from "write `derivation.md`" to "fill each step's fields; the
   derivation file is generated from them". This needs your OK, because it edits the definition.

**Alternative:** keep freehand writing, add a template and the lint. It is cheaper, but style still drifts.

Once you approve, T04–T07 can be rerun to regenerate the four files in the new format. That is about
5 minutes of agent time.

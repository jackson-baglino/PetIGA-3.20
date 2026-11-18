# The derived model (D1) against the current model (A1)

**Question.** If we write the residual and Jacobian from the equations `beth-mathematician` derived blind
(tasks T04 to T07, 2026-11-18), do we get the current model back? Where the two differ, how much does it
change a run?

**What D1 is.** It is `src/assembly_derived.c`, selected with `-residual_model derived` (default `A1`). It was written
from the agent's `answer.json` results, with the code's weak-form conventions:

$$\partial_t\phi + M\Big[\tfrac{6\gamma}{\varepsilon}\phi(1-\phi)(1-2\phi) - 6\gamma\varepsilon\nabla^2\phi\Big]
 - 30M\rho_i R_v T\,\phi^2(1-\phi)^2\,\frac{\rho_v-\rho_{vs}(T)}{\rho_{vs}(T)} = 0 \qquad \text{(T07)}$$

$$C(\phi)\,\partial_t T - \nabla\cdot\big(k(\phi)\nabla T\big) - \rho_i L\,\partial_t\phi = 0,\qquad
C(\phi) = \phi\rho_i c_i + (1-\phi)\rho_a c_a \qquad \text{(T04)}$$

$$(1-\phi)\,\partial_t\rho_v - \nabla\cdot\big(D(T)(1-\phi)\nabla\rho_v\big) + (\rho_i-\rho_v)\,\partial_t\phi = 0
\qquad \text{(T05)}$$

The three runs are `20261118-124517` (T07), `20261118-120651` (T04) and `20261118-123403` (T05), under
`runs/model-development/beth-mathematician/`.

**Not from the derivation.** These are kept identical to A1 so the comparison isolates the physics:

- the time-scaling factors $\xi_T$ and $\xi_v$;
- the floor $\max(1-\phi, \phi_{a,\min})$ in the vapor row and the clamp of $\phi$ to $[0,1]$ in the
  coefficients;
- the material laws $k(\phi)$, $D(T)$, $\rho_{vs}(T)$ and the constant $L$;
- $\rho_v c_v$, left out of $C$ (it is about $10^{-6}$ of $\rho_a c_a$); T04's $L(T)$ is taken as constant.

D1 has no wall, open-face or axisymmetric terms, and the solver refuses those settings with it.

## The parameter map (checked)

`verification/coupling_map.py` (SymPy) proves the map. With $M = M_{\text{code}}/(2\gamma)$, where
$M_{\text{code}}$ is `mob_sub`, D1's gradient and well coefficients equal A1's $3M_{\text{code}}\varepsilon$
and $3M_{\text{code}}/\varepsilon$. The coupling coefficients are related by

$$\frac{30M\rho_iR_vT/\rho_{vs}(T)}{\alpha/\rho_i} = \frac{5}{a_1}\,\frac{T}{T_0}\,\frac{\rho_{vs}(T_0)}{\rho_{vs}(T)},$$

where $\alpha$ is `alph_sub`. The code's $a_1 = 5$ makes this exactly 1 at $T = T_0$. So the coupling that
the code got from the thin-interface analysis is the one the agent got from the free energy. The only
difference is that the code freezes $T/\rho_{vs}(T)$ at the initial temperature $T_0$. In terms of $M$,
$\alpha = 30M\rho_i^2R_vT_0/\rho_{vs}(T_0)$.

Two options let D1 reproduce A1:

- `-derived_frozen 1` evaluates $T/\rho_{vs}(T)$ at $T_0$;
- `-derived_cp_mix 0` uses A1's $\rho(\phi)\,c_p(\phi)$.

At $\phi = 1/2$, A1's heat capacity is $0.767$ of the mixture value: the difference is
$-\phi(1-\phi)(\rho_i-\rho_a)(c_i-c_a)$.

## Checks

`-test_derived 1` evaluates both residuals on one state, prints the numbers and exits. The state is the
initial condition with $\phi$ mapped into $[0.02, 0.98]$, $T \pm 1$ K and $\rho_v \pm 10\,\%$ at random,
and a random $\partial_t U$. It then checks each Jacobian against central finite differences in both $U$
and $\partial_t U$. The case is the sintering pair (`sinterpair_2D_L120um_eps0.86um`), on 1 rank.

| Setting | $\phi$ residual, relative difference | $T$ | $\rho_v$ | Run folder (`~/SimulationResults/lunar_regolith_DSM/scratch/sinterpair_2D_L120um_eps0.86um/`) |
|---|---|---|---|---|
| matched (`-derived_cp_mix 0 -derived_frozen 1`) | $1.8\times10^{-16}$ | 0 | 0 | `2026-11-18__15.17.51_..._testD1_matched` |
| derived defaults | $4.0\times10^{-2}$ | $5.4\times10^{-12}$ | 0 | `2026-11-18__15.18.44_..._testD1_default` |

The matched row says D1 equals A1 to round-off. In the defaults row, the $\phi$ difference comes from the local
$T$ in the coupling, which is large here because of the random $\pm1$ K; the $T$ difference comes from the heat
capacity, which multiplies only $\partial_t T$.

Jacobian, worst $\|Jv - \text{FD}\|/\|Jv\|$ over 3 random directions (both settings give the same values to
three digits):

| | $\phi$ | $T$ | $\rho_v$ |
|---|---|---|---|
| A1 | $9.8\times10^{-12}$ | $1.7\times10^{-8}$ | $6.9\times10^{-10}$ |
| D1 | $1.0\times10^{-11}$ | $1.7\times10^{-8}$ | $6.9\times10^{-10}$ |

**Limits of this check.**

- With $\phi$ at exactly 0 or 1 (the raw initial condition), both Jacobians miss the finite difference by
  up to 11 % in the $T$ row. That comes from the clamp kinks, not from either Jacobian. The test therefore
  maps $\phi$ inside.
- D1's Jacobian has the term $\partial_\phi C\,\partial_tT$, which A1's leaves out. At these magnitudes it
  is about $10^{-10}$ of the row, so this test cannot see it either way.

## The time-dependent comparison

The case is the sintering pair at $-20\,^\circ$C, `dtlad60_T-20_sealed_ac1e-3_dt2tau` with
`-t_final 1.0e6`, on 2 ranks. It was run once with A1 and once with `-residual_model derived` (defaults). It
is compared with `verification/compare_runs.py`.

The run folders are under `~/SimulationResults/lunar_regolith_DSM/scratch/sinterpair_2D_L120um_eps0.86um/`:
`2026-11-18__15.19.22_..._cmpA1` (A1) and `2026-11-18__15.21.44_..._cmpD1` (derived). Each took about
2.5 minutes. Both took the same 147 steps to $t = 1.001\times10^6$ s.

| Quantity | max $\lvert\text{derived} - \text{A1}\rvert$ | Scale |
|---|---|---|
| interface measure (SSA_evo col. 0) over time | $2.2\times10^{-7}$ relative | falls 18 % (4.505e-5 to 3.703e-5): the neck grows |
| ice volume over time | 0 | constant, 3.549523e-9 |
| vapor content over time | $1.6\times10^{-7}$ relative | |
| total water mass over time | $3.1\times10^{-7}$ relative | |
| $\phi$ at the end, pointwise | $1.1\times10^{-6}$ | 0 to 1 |
| $T$ at the end, pointwise | $1.2\times10^{-6}$ K | uniform, $-20.00005\,^\circ$C |
| $\rho_v$ at the end, pointwise | $1.0\times10^{-10}$ kg m$^{-3}$ | $8.487\times10^{-4}$ |

**Reading.** On this case the two models give the same result to about $10^{-6}$. That is the expected
answer: $T$ stays within $5\times10^{-5}$ K of $T_0$, so the local-$T$ coupling differs from the frozen
one by about $\partial_T\ln\rho_{vs}\cdot\Delta T \approx 0.1\ \text{K}^{-1}\times 5\times10^{-5}\ \text{K}
= 5\times10^{-6}$ [inferred]. The heat capacity multiplies $\partial_tT \approx 0$.

**This case cannot separate the models.** Telling them apart needs a temperature gradient. A run with the
`tgrad` experiments, which use a wedge geometry, is the next comparison. Two cautions:

- I have not checked those geometries for wall terms, which D1 refuses.
- A $\pm 5$ K span changes $\rho_{vs}$ by a factor of about 2.6, so the coupling would differ by tens of
  percent from one end to the other [inferred].

# Sizing the model from a constant $\alpha_c$ inside the literature band

The companion study, `studies/libbrecht_kinetics/`, argues that Libbrecht's
$\sigma_0(T)$ cannot set $\alpha_c$: no single reference $\sigma$ keeps
$\alpha_c = \exp(-\sigma_0/\sigma)$ inside the literature band across
−40…−1 °C, and following it down in temperature demands an interface width
below the size of a water molecule.

This study assumes that conclusion and asks the practical question instead. If
$\alpha_c$ is simply a **constant chosen inside the band the literature
supports** — $10^{-3} < \alpha_c < 10^{-1}$ (Libbrecht 2017; Braun, Fourteau &
Löwe, *The Cryosphere* **18**, 1653, 2024) — what does each choice cost?

$$\alpha_c \;\to\; \beta_{sub} \propto 1/\alpha_c \;\to\; \text{the K\&P bounds}
\;\to\; \epsilon \;\to\; N_x$$

No $\sigma_0$ table, no nucleation law, no extrapolation — three numbers the
literature already endorses, pushed through `comp_eps.compute_eps()`.

## Regenerate

```bash
python preprocess/plot_alpha_constant_constraints.py --out studies/alpha_c_sizing
```

Reference case (needed because Eq. 45 wants a front velocity and $N_x$ wants a
domain): the Molaro dom2 production domain, $L_x \times L_y = 450 \times
225\ \mu$m, $R_{ave} = 86.75\ \mu$m, at the measured neck rate
$v_n = 3.416\times10^{-9}$ m/s. $L_x$ enters only as
$N_x = \lceil\sqrt{2}L_x/\epsilon\rceil$, i.e. *linearly*.

## Figures

| file | what it shows |
|---|---|
| `fig1_beta_sub_vs_T.png` | $\beta_{sub}(T)$ at each $\alpha_c$, against M&F's Table S1 range |
| `fig2_eps_vs_T.png` | $\epsilon(T)$, shaded by which bound binds |
| `fig3_Nx_vs_T.png` | $N_x(T)$, against the practical ceiling and the run we do |
| `fig4_bounds_vs_alpha.png` | the four bounds against $\alpha_c$ — **why** there is an optimum |
| `fig5_Nx_vs_alpha.png` | $N_x$ against $\alpha_c$ at three temperatures — **where** it is |
| `fig6_fkin_vs_alpha.png` | is the sharp-interface mapping still meaningful? |
| `fig0_master.png` | all six, 3 × 2, one shared legend — 10 × 6.6 in |
| `<panel>_legend.png` | each panel's own key, 3 in wide |
| `fig7_legend.png` | one key for the whole set, 10 in wide |
| `alpha_c_sizing.csv` | the numbers |

Format and policy are identical to the Libbrecht set — panels 3 × 2.4 in,
master 10 × 6.6 in, nothing below 10 pt, 10 in a hard width ceiling — because
both import `preprocess/figstyle.py`. Neither set can drift from the other.

## The result: $\alpha_c$ is not monotone in cost

It is tempting to read "$\alpha_c$ is a free parameter, so pick a small one —
the kinetics are slow but the mesh is cheap". The bounds say otherwise, because
two of them move in **opposite** directions:

| bound | scales as | as $\alpha_c$ falls |
|---|---|---|
| B-HEAT, B-VAPOR | $\epsilon \le s\,D\,\beta_{HK} \propto 1/\alpha_c$ | **loosens** |
| B-KINETIC | $\epsilon \le s\,d_0/(\beta_{sub}v_n) \propto \alpha_c$ | **tightens** |

So $N_x(\alpha_c)$ is V-shaped and there is a genuine optimum. Panel 4 is the
mechanism — the two lines crossing — and panel 5 is where the window lands:

| $T$ | cheapest $\alpha_c$ | $N_x$ there | $N_x$ at $\alpha_c = 10^{-4}$ | at $\alpha_c = 1$ |
|---|---|---|---|---|
| −40 °C | 0.11 | 2 339 | 2.6e6 | 2.1e4 |
| −20 °C | 0.039 | 867 | 3.4e5 | 2.2e4 |
| −5 °C | 0.020 | 449 | 9.0e4 | 2.2e4 |

Both ends are expensive, and both exceed the $N_x = 10^4$ practical ceiling.
The optimum sits at the **top** of the literature band, not the bottom.

## What each of the three costs

| $T$ | $\alpha_c$ | $\beta_{sub}$ [s/m] | $\epsilon$ | binds | $N_x$ | $f_{kin}$ |
|---|---|---|---|---|---|---|
| −40 °C | $10^{-1}$ | 6.7e5 | 242 nm | B-KINETIC | 2 626 | 0.68 |
| | $10^{-2}$ | 6.7e6 | 24 nm | B-KINETIC | 26 254 | 1.00 |
| | $10^{-3}$ | 6.7e7 | 2.4 nm | B-KINETIC | 262 531 | 1.00 |
| −20 °C | $10^{-1}$ | 7.9e4 | 296 nm | B-HEAT | 2 153 | 0.64 |
| | $10^{-2}$ | 7.9e5 | 188 nm | B-KINETIC | 3 389 | 0.97 |
| | $10^{-3}$ | 7.9e6 | 19 nm | B-KINETIC | 33 887 | 1.00 |
| −5 °C | $10^{-1}$ | 2.0e4 | 287 nm | B-HEAT | 2 216 | 0.65 |
| | $10^{-2}$ | 2.0e5 | 709 nm | B-KINETIC | 898 | 0.88 |
| | $10^{-3}$ | 2.0e6 | 71 nm | B-KINETIC | 8 977 | 1.00 |

Read across: **$\alpha_c = 10^{-3}$ is the expensive choice**, costing
$3.4\times10^4$ nodes across $x$ at −20 °C and $2.6\times10^5$ at −40 °C — two
to twenty-six times the practical ceiling. $\alpha_c = 10^{-1}$ and $10^{-2}$
both sit comfortably under it everywhere in −40…−1 °C.

## The catch at the top of the band

Raising $\alpha_c$ shrinks $\beta_{HK}$, which shrinks the binding $\epsilon$,
which makes the thin-interface corrections a *larger* share of $\tau_{sub}$.
Panel 6 plots that share, $f_{kin}$: it falls from 1.00 at $\alpha_c \le
10^{-2}$ to about 0.64 at $\alpha_c = 10^{-1}$ and stays there. `comp_eps`
flags below 0.5 and warns below 0.1, so 0.64 is *acceptable but no longer
generous* — the O($\delta^2$) residual in the Eq. (9) mapping is real at that
end even though the mesh is cheapest there.

That is the whole trade, in one line: the mesh wants $\alpha_c$ high, the
asymptotics want it low, and $10^{-2}$ is the value that is comfortable on both
counts across the full temperature range.

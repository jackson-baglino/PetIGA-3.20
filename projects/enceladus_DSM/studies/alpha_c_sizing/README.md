# $\alpha_c$ as a knob: sweeping it over $[10^{-4}, 1]$

The companion study, `studies/libbrecht_kinetics/`, argues that Libbrecht's
$\sigma_0(T)$ cannot set $\alpha_c$. This one takes $\alpha_c$ as what it
actually is in the solver — a constant we choose — and puts it on the x-axis,
so every panel answers *"if I pick this $\alpha_c$, what do I get?"*

$$\alpha_c \;\to\; \beta_{sub} \propto 1/\alpha_c \;\to\;
\text{the K\&P bounds} \;\to\; \epsilon \;\to\; N_x$$

The band the literature supports, $10^{-3} < \alpha_c < 10^{-1}$ (Libbrecht
2017; Braun, Fourteau & Löwe, *The Cryosphere* **18**, 1653, 2024), is shaded
on **every** panel: it is the only region a choice is allowed to come from. The
sweep runs a decade past it on each side so the shape of the trade is visible.

## Regenerate

```bash
python preprocess/plot_alpha_c_sweep.py --out studies/alpha_c_sizing
```

Reference case (needed because Eq. 45 wants a front velocity and $N_x$ wants a
domain): the Molaro dom2 production domain, $L_x \times L_y = 450 \times
225\ \mu$m, $R_{ave} = 86.75\ \mu$m, at the measured neck rate
$v_n = 3.416\times10^{-9}$ m/s. $L_x$ enters only as
$N_x = \lceil\sqrt{2}L_x/\epsilon\rceil$, i.e. *linearly*.

## Figures

| file | what it shows |
|---|---|
| `fig1_beta_vs_alpha.png` | $\beta_{sub}$ and $\beta_{HK}$ against $\alpha_c$ |
| `fig2_bounds_vs_alpha.png` | the four K&P bounds — **why** there is a tent |
| `fig3_eps_vs_alpha.png` | the interface width that survives them |
| `fig4_Nx_vs_alpha.png` | the mesh it demands |
| `fig5_fkin_vs_alpha.png` | is the sharp-interface mapping still meaningful? |
| `fig6_Lstar_vs_alpha.png` | when does $\alpha_c$ stop mattering at all? |
| `fig0_master.png` | all six, 3 × 2, one shared legend — 10 × 6.6 in |
| `<panel>_legend.png` | each panel's own key, 3 in wide |
| `fig7_legend.png` | one key for the whole set, 10 in wide |
| `alpha_c_sweep.csv` | the numbers |

Format and policy are identical to the Libbrecht set — panels 3 × 2.4 in,
master 10 × 6.6 in, nothing below 10 pt, 10 in a hard width ceiling — because
both import `preprocess/figstyle.py`. Neither set can drift from the other.

## Why $\beta_{sub}$ carries a temperature dependence at all

There are two coefficients, and only one of them really moves with $T$:

$$\beta_{HK}(T) = \frac{1}{\alpha_c}\sqrt{\frac{2\pi m}{k_B T}}
\quad\text{(scaled, K\&P }\beta'\text{)}, \qquad
\beta_{sub}(T) = \frac{\beta_{HK}}{\rho_{vs}(T)/\rho_i}
\quad\text{(unscaled, K\&P }\beta_0\text{)}$$

$\beta_{HK}$ is the physical Hertz–Knudsen coefficient, and its only
temperature dependence is the $\sqrt{T}$ in the mean thermal speed — worth
**8%** over −40…−1 °C, which is nothing. All of the spread in $\beta_{sub}$ is
the $\rho_{vs}(T)/\rho_i$ factor that converts between the two, and $\rho_{vs}$
falls **44×** over that same range.

So panel 1's three curves are not three different kinetics. They are *one*
kinetics seen through a saturation vapour density that collapses as it gets
cold — which is why the panel draws $\beta_{HK}$ underneath them as a single
line: the ~$10^6$ gap between the families is exactly $\rho_{vs}/\rho_i$. The
solver is passed `-beta_sub0` (unscaled) and rescales internally, which is why
$\beta_{sub}$ is the one plotted.

## The result: $\alpha_c$ is not monotone in cost

It is tempting to read "$\alpha_c$ is free, so pick a small one — the kinetics
are slow but the mesh is cheap". The bounds say otherwise, because two of them
move in **opposite** directions:

| bound | scales as | as $\alpha_c$ falls |
|---|---|---|
| B-HEAT, B-VAPOR | $\epsilon \le s\,D\,\beta_{HK} \propto 1/\alpha_c$ | **loosens** |
| B-KINETIC | $\epsilon \le s\,d_0/(\beta_{sub}v_n) \propto \alpha_c$ | **tightens** |

Panel 2 is that crossing. $\epsilon(\alpha_c)$ is therefore a tent (panel 3)
and $N_x(\alpha_c)$ a V (panel 4), with a genuine optimum:

| $T$ | cheapest $\alpha_c$ | $\epsilon$, $N_x$ there | $N_x$ at $10^{-4}$ | at $1$ | $N_x \le 10^4$ for |
|---|---|---|---|---|---|
| −40 °C | 0.112 | 272 nm, 2 339 | 2.6e6 | 2.1e4 | 0.026 – 0.48 |
| −20 °C | 0.039 | 737 nm, 864 | 3.4e5 | 2.2e4 | 0.0035 – 0.46 |
| −5 °C | 0.020 | 1 419 nm, 449 | 9.0e4 | 2.2e4 | 0.0009 – 0.45 |

Both ends are expensive, and both exceed the $N_x = 10^4$ practical ceiling.
**The optimum sits at the top of the literature band, not the bottom** — and
the affordable window narrows as it gets colder: at −5 °C almost the whole band
is affordable, at −40 °C only $\alpha_c \gtrsim 0.026$ is.

## Two things that bound the top end

**Asymptotic validity (panel 5).** Raising $\alpha_c$ shrinks $\beta_{HK}$,
which shrinks the binding $\epsilon$, which makes the thin-interface
corrections a *larger* share of $\tau_{sub}$. That share, $f_{kin}$, falls from
1.00 at $\alpha_c \le 10^{-2}$ to about 0.64 at $\alpha_c \ge 10^{-1}$.
`comp_eps` flags below 0.5 and warns below 0.1, so 0.64 is acceptable but no
longer generous: the O($\delta^2$) residual in the M&F Eq. (9) mapping is real
at that end.

**Diminishing returns (panel 6).** $L^* = \beta_{HK}D_v$ is the crossover
length: features smaller than $L^*$ are attachment-limited, so $\alpha_c$
controls them; features larger are vapour-diffusion-limited, and $\alpha_c$
stops mattering. $L^*$ drops below the grain radius at $\alpha_c \approx
1.5\times10^{-3}$ and below the neck radius at $\approx 10^{-2}$. Above that,
raising $\alpha_c$ buys progressively less real speed-up — which is
independently what `studies/molaro_2019/alpha_c_estimate.csv` found, that
$\alpha_c$ is exhausted as a knob near 0.1 because the process is already
transport-limited there.

$L^*$ is also nearly temperature-independent — the three curves lie on top of
one another — for the same reason $\beta_{HK}$ is.

## The trade, in one line

The mesh wants $\alpha_c$ high, the asymptotics want it low, and the transport
regime says the high end buys less than it looks like it should.
$\alpha_c \sim 10^{-2}$ is the value that is comfortable on all three counts
across the full temperature range; $\alpha_c = 0.1$, which the Molaro campaign
runs, is at the cheap end for mesh and the marginal end for $f_{kin}$.

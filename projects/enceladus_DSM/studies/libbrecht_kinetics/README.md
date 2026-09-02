# Libbrecht (2017) $\sigma_0(T)$ as a source for $\alpha_c$ — why we do not use it

Libbrecht, *Annu. Rev. Mater. Res.* **47**, 271 (2017) gives a critical
supersaturation $\sigma_0(T)$ for ice crystal growth together with the
nucleation form for the attachment coefficient

$$\alpha_c(T,\sigma) = \exp\!\left[-\sigma_0(T)/\sigma\right].$$

Using it would turn $\alpha_c$ — the model's only genuinely free parameter —
into a function of the local state. This directory is the record of why that
does not work here, regenerated against the *current* `preprocess/comp_eps.py`
(the earlier version of these figures predates several fixes to the sizer).

## Regenerate

```bash
python preprocess/plot_libbrecht_constraints.py --out studies/libbrecht_kinetics
```

The **argument** uses only Libbrecht's own chamber supersaturations. The
**reference case** — needed because Eq. (45) wants a front velocity and $N_x$
wants a domain — defaults to a $450 \times 225\ \mu$m grain pair
($R_{ave} = 86.75\ \mu$m) at the measured front velocity
$v_n = 3.416\times10^{-9}$ m/s, the integrated Fig. 11 neck-growth rate from
`studies/molaro_2019/alpha_c_estimate.csv`. It sets the scale of the vertical
axes, not the conclusion. Every $\epsilon$ and every mesh count comes from
`comp_eps.compute_eps()` itself, so these plots cannot drift from what the
sizer does.

## Figures

| file | what it shows |
|---|---|
| `fig1_sigma0_vs_T.png` | $\sigma_0(T)$: Libbrecht's table and the log–log interpolant `comp_eps.sigma0()` builds from it |
| `fig2_alpha_vs_sigma0.png` | $\alpha_c$ against $\sigma_0$, at each chamber $\sigma$ |
| `fig3_alpha_vs_T.png` | the same, in temperature |
| `fig4_beta_sub_vs_T.png` | $\beta_{sub}(T) \propto 1/\alpha_c$, with the defining equation on the axes |
| `fig5_eps_vs_T.png` | $\epsilon(T)$ from the four K&P bounds, shaded by which bound binds |
| `fig6_Nx_vs_T.png` | the $N_x(T)$ that follows, against the $N_x = 10^4$ practical ceiling |
| `fig0_master.png` | all six at once, 3 × 2, one shared legend — 10 × 6.6 in |
| `fig7_legend.png` | the shared legend alone, for laying the small panels out by hand |
| `libbrecht_constraints.csv` | the numbers |

**Axis windows.** Plotted full-range these quantities span sixty decades and
the usable band is a hairline, so every axis is clipped to the *literature*
band opened up by two decades on each side — α_c to [1e-5, 2.5], β_sub to
[1e3, 2e8] around M&F's [2e4, 2e6], and N_x to [4e2, 1e6] around the
N_x = 10⁴ practical ceiling. Panel 5's ε window is the exact reciprocal of
panel 6's N_x window, so the two read against each other.

**Format.** The six panels are **3.00 × 2.40 in** at 200 dpi, sized to sit
several-to-a-slide, and `fig0_master.png` is all six inside **10.00 × 6.60 in**.
Nothing is set below 10 pt anywhere, and `_save()` asserts the 10 in width
ceiling on every figure rather than trusting a figsize.

At 3 in a panel cannot hold a title, prose annotations *and* a legend, so below
`--panel_width 6` they are drawn in **compact mode**: short title, short y-label
where the full one steals width, no prose annotations, no per-panel legend.
`fig7_legend.png` (10 × 1.15 in) is that legend on its own — place it once under
a row of panels. Compact mode is a flag the same panel functions consult, not a
second set of drawing code, so the small panels and the master cannot disagree.
`--panel_width 10` gives the fully annotated panels back.

> If PowerPoint inserts a figure larger than its stated size, it ignored the
> PNG's DPI tag and assumed 96 dpi. Set the width box back — nothing is lost,
> the image is simply 200 dpi at that size. `--dpi 96` produces drop-in sizing
> instead, at some cost in sharpness on a projector.

**Only Libbrecht's own chamber conditions are plotted** — $\sigma = 10^{-1}$
(solid) and $\sigma = 10^{-2}$ (dashed). Nothing here depends on the conditions
of one of our simulations, so nothing here can be answered with "your boundary
condition is wrong". See *What is deliberately not plotted*, below.

## The argument, in three steps

**1. The law cannot sit inside the literature band at any one σ.** The band the
literature supports is $10^{-3} < \alpha_c < 10^{-1}$ (Libbrecht 2017; Braun,
Fourteau & Löwe 2024). Evaluate $\alpha_c = \exp(-\sigma_0/\sigma)$ at
Libbrecht's *own upper* chamber condition and it overshoots; at his *own lower*
one it undershoots, once it is cold enough:

| $\sigma$ | $\alpha_c$ at −1 °C | at −20 °C | at −40 °C | vs. the band |
|---|---|---|---|---|
| $10^{-1}$ | 0.96 | 0.71 | 0.33 | 3–10× **above** the ceiling, everywhere |
| $10^{-2}$ | 0.67 | 3.1e-2 | 1.7e-5 | **below** the floor colder than −29.7 °C |

That is structural, not a tuning problem. $\sigma_0$ spans a factor 27 over
−2 to −40 °C and it sits in the *exponent*, so no single reference $\sigma$
holds $\alpha_c$ inside one decade across the range.

**2. K&P Eq. (45) turns that into mesh.** $\beta_{sub} \propto 1/\alpha_c$, and
the kinetic bound is $\epsilon \le s\,d_0/(\beta_{sub}v_n)$, so a falling
$\alpha_c$ is not just slow physics we can wait out. At $\sigma = 10^{-2}$ —
Libbrecht's own condition:

| $T$ | $\alpha_c$ | $\beta_{sub}$ [s/m] | $\epsilon$ | binds | $N_x$ |
|---|---|---|---|---|---|
| −20 °C | 3.1e-2 | 2.6e5 | 0.58 µm | B-KINETIC | 1.1e3 |
| −30 °C | 9.1e-4 | 2.4e7 | 6.4 nm | B-KINETIC | 9.9e4 |
| −40 °C | 1.7e-5 | 4.0e9 | **0.40 Å** | B-KINETIC | **1.6e7** |

For scale: the reference run is $\epsilon = 0.118\ \mu$m, $N_x = 5394$, 14.5 M
nodes; one water molecule is $a = (m/\rho_i)^{1/3} = 3.19$ Å, so $L_x/a =
1.4\times10^6$ is *one node per molecule across the domain*. At −40 °C the law
asks for an $\epsilon$ an **eighth of a water molecule** and eleven times more
nodes than the domain has molecules. It has passed $N_x = 10^4$ by −25 °C.
That is not an expensive mesh; it is a category error, because the continuum
phase field has no meaning below $a$.

**3. It is the wrong experiment for our problem.** Libbrecht grew single ice
crystals from vapour in a near-vacuum chamber. That isolates any dependence
$\alpha_c$ has on supersaturation — which is what makes the data valuable — but
it also means $\sigma$ was *imposed externally*. Our simulations are not in
vacuum: there is a continuous vapour field between the grains, so $\sigma$ is a
**solution variable**, set by the local curvature and the wall humidity. And it
is a nine-point table, non-monotonic (the kink at −6/−7 °C is a digitisation
artifact), which `comp_eps.libbrecht2_params()` can only fit by freeing the
prefactor and rescaling — after which $A$ pins to the clamp ceiling.

## What is deliberately not plotted

Sintering is capillarity-driven, so its supersaturation is $d_0/\rho_f$ — lower
again than either chamber value, for any micron-scale fillet. That only makes
the result worse, and earlier versions of this figure drew it for a particular
geometry (an imposed wall undersaturation, and a neck-fillet $d_0/\rho$). Those
curves are gone, for two reasons:

- they are specific to a **single digital experiment**, and the case being made
  is general;
- at those $\sigma$ the law returns an $\alpha_c$ decades below the literature
  band, so the $\beta_{sub}$ it implies is one **we would never run** — drawing
  it invites the reading that we do.

A reference case is still unavoidable: Eq. (45) needs a front velocity and
$N_x$ needs a domain. The defaults are a 450 × 225 µm grain pair and the
measured integrated neck rate 3.416e-9 m/s from `studies/molaro_2019/`. They
set the **scale** of the vertical axes; they do not carry the argument, which
is about where $\alpha_c(T)$ lands relative to a band that has nothing to do
with them.

## The escape route, and why it is circular

Passing `--vn_feature R_feat` to `comp_eps.py` feeds $v_n$ back
self-consistently from the kinetics, $v_n = d_0/(\beta_{sub} R_{feat})$. Then
Eq. (45) collapses to $\epsilon \le s\,R_{feat}$ — independent of
$\beta_{sub}$, hence of $\alpha_c$ — and the mesh becomes affordable. But it is
affordable only because the constraint has been made insensitive to the
quantity under test. The cost moves from resolution to runtime: physical
$t_{final} \propto 1/\alpha_c$, so a 2.5 h experiment needs ~65 h simulated at
$\alpha_c = 10^{-4}$ and ~26 days at $10^{-5}$.

Fixing $v_n$ to a *measured* value is what makes Eq. (45) a real constraint,
and it is why the default here is Molaro's own neck rate rather than
`comp_eps.py`'s generic 1e-9 m/s.

## What we use instead

A constant $\alpha_c$ inside the literature band
$10^{-3} \le \alpha_c \le 10^{-1}$ (Libbrecht 2017; Braun, Fourteau & Löwe,
*The Cryosphere* **18**, 1653, 2024), currently $\alpha_c = 0.1$ for the Molaro
campaign — a value `studies/molaro_2019/alpha_c_estimate.csv` shows is already
exhausted as a knob, because the process is transport-limited there
($\Lambda = 0.29$). Libbrecht's data is kept as what it can honestly support:
the **sign** of the temperature dependence, and a plausible band.

See `preprocess/plot_alpha_models.py` for the alternative laws that were
considered (Arrhenius in $T$; the two-parameter Libbrecht fit). The constraints
themselves are written out as LaTeX in `docs/tex/` —
`constraints_iguanatex.txt` for copy-pasting single equations into PowerPoint
via IguanaTeX, the `.tex` fragments for a written document.

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

Defaults are the Molaro dom2 production geometry ($450 \times 225\ \mu$m,
$R_{ave} = 86.75\ \mu$m) and the **measured** front velocity
$v_n = 3.416\times10^{-9}$ m/s — the integrated Fig. 11 neck-growth rate from
`studies/molaro_2019/alpha_c_estimate.csv`. Every $\epsilon$ and every mesh
count comes from `comp_eps.compute_eps()` itself, so these plots cannot drift
from what the sizer does.

## Figures

| file | what it shows |
|---|---|
| `fig1_sigma0_vs_T.png` | $\sigma_0(T)$: Libbrecht's table and the log–log interpolant `comp_eps.sigma0()` builds from it |
| `fig2_alpha_vs_sigma0.png` | $\alpha_c$ against $\sigma_0$, at each $\sigma$ in play |
| `fig3_alpha_vs_T.png` | the same, in temperature |
| `fig4_beta_sub_vs_T.png` | $\beta_{sub}(T) \propto 1/\alpha_c$, with the defining equation on the axes |
| `fig5_eps_vs_T.png` | $\epsilon(T)$ from the four K&P bounds, shaded by which bound binds |
| `fig6_Nx_vs_T.png` | the $N_x(T)$ that follows, against the production mesh and the molecular limit |
| `fig0_overview.png` | all six on one slide |
| `libbrecht_constraints.csv` | the numbers |

**Axis windows.** Plotted full-range these quantities span sixty decades and
the usable band is a hairline, so every axis is clipped to the *literature*
band opened up by two decades on each side — α_c to [1e-5, 30], β_sub to
[2e2, 2e8] around M&F's [2e4, 2e6], and N_x to [3e1, 1e6] around the
N_x = 10⁴ practical ceiling. Panel 5's ε window is the exact reciprocal of
panel 6's N_x window, so the two read against each other.

Each of panels 2, 3, 5 and 6 reserves a band that *no data can reach* and puts
its legend there, rather than laying it over the curves: above α_c = 1 in
panels 2–3 (every impinging molecule sticks — nothing is above it), above the
largest resolvable ε in panel 5, and below the smallest N_x in panel 6.

Where a curve leaves the frame it is marked on the edge with the extreme it
reaches, so nothing is hidden; the full range is in
`libbrecht_constraints.csv` and in the table below.

The four supersaturations plotted, and where each comes from:

- $\sigma = 10^{-1},\ 10^{-2}$ — **Libbrecht's own chamber conditions** (dashed)
- $\sigma = 2.7\times10^{-3}$ — **our** imposed Molaro wall undersaturation
  $1-h$, from `studies/molaro_2019/vapor_bc_estimate.csv` (solid)
- $\sigma = 4.5\times10^{-4}$ — **our** Gibbs–Thomson supersaturation at a neck
  fillet, $d_0/\rho_f$ (solid)

## The argument, in three steps

**1. Tractability.** $\beta_{sub} \propto 1/\alpha_c$, and K&P Eq. (45) is
$\epsilon \le s\,d_0/(\beta_{sub}v_n)$, so the mesh is *exponential* in
$\sigma_0/\sigma$. At $T = -20\ ^\circ$C:

| $\sigma$ | $\alpha_c$ | $\beta_{sub}$ [s/m] | $\epsilon$ [m] | binds | $N_x$ | 2D nodes |
|---|---|---|---|---|---|---|
| $10^{-1}$ (Libbrecht) | 7.1e-1 | 1.1e4 | 4.2e-8 | B-HEAT | 1.5e4 | 1.2e8 |
| $10^{-2}$ (Libbrecht) | 3.1e-2 | 2.6e5 | 5.8e-7 | B-KINETIC | 1.1e3 | 6.1e5 |
| $2.7\times10^{-3}$ (ours) | 2.5e-6 | 3.2e9 | 4.7e-11 | B-KINETIC | 1.4e7 | 9.2e13 |
| $4.5\times10^{-4}$ (ours) | 1e-30 (floor) | 7.9e33 | 1.9e-35 | B-KINETIC | 3.4e31 | 5.7e62 |

For scale: the production run is $\epsilon = 0.118\ \mu$m, $N_x = 5394$,
14.5 M nodes; one water molecule is $a = (m/\rho_i)^{1/3} = 3.19$ Å, so
$L_x/a = 1.4\times10^6$ is *one node per molecule across the domain*. At our
own wall BC the law already demands ~10× more nodes across $x$ than that, with
$\epsilon$ half an Ångström. At the fillet supersaturation it demands
$\epsilon$ twenty-five orders of magnitude below a water molecule. That is not
an expensive mesh; it is a category error, because the continuum phase field
has no meaning below $a$.

**2. It is the wrong experiment.** Libbrecht grew single ice crystals from
vapour in a near-vacuum chamber. That isolates any dependence $\alpha_c$ has on
supersaturation — which is exactly what makes the data valuable — but it also
means $\sigma$ was *imposed externally* and held one to two decades above ours.
Our simulations are not in vacuum: there is a continuous vapour field between
the grains, so $\sigma$ is a **solution variable**, set by the local curvature
and the wall humidity. Applying a nucleation law two decades below where it was
calibrated puts the entire extrapolation inside an exponential.

**3. The table cannot be fitted anyway.** Nine points; non-monotonic (the kink
at $-6/-7\ ^\circ$C is a digitisation artifact); and $\sigma_0$ spans a factor
27 over $-2$ to $-40\ ^\circ$C. Since $\ln\alpha_c = -\sigma_0/\sigma$, no
single reference $\sigma$ holds $\alpha_c$ inside one decade across that range.
`comp_eps.libbrecht2_params()` frees the prefactor and rescales the table by
one factor $f$ to make a fit possible at all; the result is that $A$ pins to
the clamp ceiling and the $\sigma$-response survives only below
$\sigma \sim 10^{-3}$.

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

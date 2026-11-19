# T07_phase_field

Run `20261118-173013_T07_phase_field` · status **done** · written by `scripts/agentkit/render_derivation.py` from `answer.json`; do not edit by hand.

## Work order (verbatim)

> We are modeling the pore space of a granular material (lunar regolith) in which water ice sublimates and deposits. Treat the pore space as a three-dimensional region filled with solid ice and gas (air that contains water vapor). The local volume fraction of ice is the phase field phi(x, y, z, t), between 0 and 1: phi = 1 in ice and phi = 0 in gas. At phase coexistence the ice-gas interface is described by the free energy F[phi] given below, which you derived in an earlier task. Away from coexistence, ice deposits where the vapor density rho_v exceeds the saturation vapor density rho_vs(T) at the local temperature, and sublimates where rho_v is below it. There is no nucleation: ice grows or shrinks only where ice already meets gas. In this task rho_v and T are known fields; their equations are separate tasks. Derive the equation for how phi evolves in time, combining the relaxation of the interface toward the equilibrium described by F with the growth or sublimation driven by the departure of rho_v from rho_vs(T). State every assumption, and list every effect you leave out with the reason. Keep the coefficients general where you can, and say what each one means physically. Write the result as R = 0 in 3D Cartesian coordinates and return R.

## Symbols

| Symbol | Meaning | Units | Typical value (SI) | Source |
|---|---|---|---|---|
| $x$ | first position coordinate | $\mathrm{m}$ |  | work order |
| $y$ | second position coordinate | $\mathrm{m}$ |  | work order |
| $z$ | third position coordinate | $\mathrm{m}$ |  | work order |
| $t$ | time | $\mathrm{s}$ |  | work order |
| $\phi$ | local volume fraction of ice (phase field), between 0 and 1 | 1 |  | work order |
| $T$ | temperature | $\mathrm{K}$ |  | work order |
| $\rho_v$ | vapor density: mass of water vapor per unit volume of gas | $\mathrm{kg\,m^{-3}}$ |  | work order |
| $\rho_{vs}(T)$ | saturation vapor density over a flat ice surface at temperature T | $\mathrm{kg\,m^{-3}}$ |  | work order |
| $\gamma$ | energy per unit area of the ice-gas interface, constant | $\mathrm{J\,m^{-2}}$ |  | work order |
| $\varepsilon$ | the length that sets the width of the diffuse interface, constant | $\mathrm{m}$ |  | work order |
| $M$ | interface mobility (Allen-Cahn kinetic coefficient): how fast $\phi$ relaxes down the gradient of $F$; with it a curved interface moves at $v_n=-6\gamma\varepsilon M\kappa$ | $\mathrm{m^{3}\,J^{-1}\,s^{-1}}$ | 1e-06 | introduced |
| $\alpha$ | deposition/sublimation rate coefficient; a flat interface advances into the gas at $v_n = 6\varepsilon\alpha(\rho_v-\rho_{vs})$, so $6\varepsilon\alpha = k/\rho_i$ with $k$ a kinetic mass-transfer velocity (Hertz-Knudsen: $k=a_s\sqrt{k_BT/2\pi m}$) and $\rho_i$ the ice density | $\mathrm{m^{3}\,kg^{-1}\,s^{-1}}$ | 1e+03 | introduced |
| $p(\phi)$ | interface-localizing function $30\phi^2(1-\phi)^2$, zero in bulk ice and gas, $\int_0^1 p\,d\phi = 1$ | 1 | 1 | introduced |
| $\mu$ | variational derivative $\delta F/\delta\phi$ | $\mathrm{J\,m^{-3}}$ | 1e+03 | introduced |
| $\phi_0(n)$ | equilibrium flat-interface profile along the normal coordinate $n$ | 1 | 0.5 | introduced |
| $n$ | coordinate normal to a flat interface, pointing into the ice | $\mathrm{m}$ | 1e-06 | introduced |
| $v_n$ | normal velocity of the interface, positive when ice grows | $\mathrm{m\,s^{-1}}$ | 1e-09 | introduced |
| $\kappa$ | sum of principal curvatures of the interface (positive for a convex ice grain) | $\mathrm{m^{-1}}$ | 1e+04 | introduced |
| $\Delta$ | supersaturation $\rho_v-\rho_{vs}(T)$ | $\mathrm{kg\,m^{-3}}$ | 1e-08 | introduced |

## Given

- `F = Integral(3*gamma/eps*phi**2*(1 - phi)**2 + 3*gamma*eps*(Derivative(phi, x)**2 + Derivative(phi, y)**2 + Derivative(phi, z)**2), over the volume)`

## Assumptions

- **A1.** A1 (given): $\rho_v$ and $T$ are known, smooth fields.
- **A2.** A2 (given): the flat interface described by $F$ has energy $\gamma$ per unit area and a 10-90 % width of $2\ln 9\,\varepsilon$.
- **A3.** A3: Near equilibrium, the interface relaxes by linear, non-conserved (Allen-Cahn) gradient flow $\partial_t\phi=-M\,\delta F/\delta\phi$ with a positive scalar mobility $M$ (isotropic, taken constant).
- **A4.** A4: The growth rate depends linearly on the supersaturation $\rho_v-\rho_{vs}(T)$ (linear kinetics, as in Hertz-Knudsen), with a constant positive coefficient $\alpha$.
- **A5.** A5: No nucleation: the source is multiplied by $p(\phi)=30\phi^2(1-\phi)^2$, which is zero in bulk ice and gas, so ice changes only at the existing interface. Any other localized $p$ with $\int_0^1 p\,d\phi=1$ would give the same sharp-interface limit.
- **A6.** A6: Relaxation and phase change act in parallel, so their rates add.

## Derivation

### S1. Variational derivative of F

The interface relaxes down the gradient of $F$, so we need $\delta F/\delta\phi$. The Euler-Lagrange operator applied to the integrand gives it; the boundary term is dropped with no-flux (natural) boundary conditions.

$$
\mu=\frac{\delta F}{\delta\phi}=\frac{6\gamma}{\varepsilon}\phi(1-\phi)(1-2\phi)-6\gamma\varepsilon\nabla^2\phi
$$

*Uses:* given. *Operation:* Euler-Lagrange derivative of the integrand of F. *Checked by:* `check_s1` in `scratch/check.py`: S1 dF/dphi matches: True

### S2. Relaxation part (Allen-Cahn)

$\phi$ is not conserved (ice turns into vapor), so the simplest relaxation that always lowers $F$ is the non-conserved gradient flow, with a positive mobility $M$. Then $dF/dt=-\int M\mu^2\,dV\le 0$.

$$
\partial_t\phi\big|_{\rm relax}=-M\mu
$$

*Uses:* S1, A3. *Operation:* postulate linear non-conserved gradient flow. *Checked by:* none (Constitutive choice; units checked in check_units: M*gam/eps -> 1/s.)

### S3. Phase-change source term

Deposition must raise $\phi$ where $\rho_v>\rho_{vs}(T)$ and lower it where $\rho_v<\rho_{vs}$, to linear order in the departure. With no nucleation the source must vanish in pure gas and pure ice, so we multiply by $p(\phi)=30\phi^2(1-\phi)^2$, which is zero at $\phi=0,1$ with zero slope there (bulk phases stay fixed points) and is normalized to $\int_0^1p\,d\phi=1$.

$$
\partial_t\phi\big|_{\rm phase}=\alpha\,p(\phi)\left(\rho_v-\rho_{vs}(T)\right),\quad p(\phi)=30\phi^2(1-\phi)^2
$$

*Uses:* A4, A5, A6. *Operation:* postulate linear driving term localized to the interface. *Checked by:* `check_s3` in `scratch/check.py`: S3 p(0), p(1): 0 0  integral 0..1: 1

### S4. Evolution equation

The two mechanisms act in parallel, so their rates add. Moving everything to one side gives the residual $R$.

$$
R=\partial_t\phi+M\left[\frac{6\gamma}{\varepsilon}\phi(1-\phi)(1-2\phi)-6\gamma\varepsilon\nabla^2\phi\right]-30\alpha\phi^2(1-\phi)^2(\rho_v-\rho_{vs}(T))=0
$$

*Uses:* S2, S3. *Operation:* add S2 and S3, substitute S1, move to left side. *Checked by:* `check_units` in `scratch/check.py`: every term has units 1/second

### S5. Flat equilibrium profile is consistent with A2

To read off what $M$ and $\alpha$ mean, we need the flat profile with $\mu=0$. The logistic profile solves it, has 10-90% width $2\ln 9\,\varepsilon$ and energy $\gamma$ per area, as A2 states.

$$
\phi_0(n)=\frac{1}{1+e^{-n/\varepsilon}},\quad \int_{-\infty}^{\infty}(\phi_0')^2dn=\frac{1}{6\varepsilon}
$$

*Uses:* S1, A2. *Operation:* solve mu=0 in 1D, integrate. *Checked by:* `check_s5` in `scratch/check.py`: profile solves dF/dphi=0: True; 10-90 width: 4*epsilon*log(3); energy per area: gamma; int phi0_n^2 dn: 1/(6*epsilon)

### S6. Meaning of the coefficients

Put a front $\phi=\phi_0(n+v_nt)$ into S4, multiply by $\phi_0'$ and integrate across the interface; the $\mu$ term vanishes for the flat profile, and the Laplacian term gives the usual Allen-Cahn curvature law. This shows that $\alpha$ sets the deposition speed and $M$ sets the curvature-driven (Gibbs-Thomson-like) motion.

$$
v_n=6\varepsilon\alpha\,(\rho_v-\rho_{vs}(T))-6\gamma\varepsilon M\kappa
$$

*Uses:* S4, S5. *Operation:* solvability: multiply by phi_0' and integrate over n. *Checked by:* `check_s6` in `scratch/check.py`: S6 V = [6*Delta*alpha*epsilon]; diffusion coefficient on laplacian: 6*M*epsilon*gamma

## Result

$$
R = \partial_t\phi + M\left[\frac{6\gamma}{\varepsilon}\phi(1-\phi)(1-2\phi) - 6\gamma\varepsilon\left(\partial_{xx}\phi+\partial_{yy}\phi+\partial_{zz}\phi\right)\right] - 30\,\alpha\,\phi^2(1-\phi)^2\left(\rho_v-\rho_{vs}(T)\right) = 0
$$

Graded SymPy form:

```
Derivative(phi(x, y, z, t), t) + M*(6*gamma/eps*phi(x, y, z, t)*(1 - phi(x, y, z, t))*(1 - 2*phi(x, y, z, t)) - 6*gamma*eps*(Derivative(phi(x, y, z, t), (x, 2)) + Derivative(phi(x, y, z, t), (y, 2)) + Derivative(phi(x, y, z, t), (z, 2)))) - alpha*30*phi(x, y, z, t)**2*(1 - phi(x, y, z, t))**2*(rho_v(x, y, z, t) - rho_vs(T(x, y, z, t)))
```

## Neglected effects

| Effect | Why it is left out |
|---|---|
| Anisotropy of ice surface energy and kinetics (faceting) | Taken as isotropic: $\gamma$, $M$ and $\alpha$ are scalars. This is acceptable for rounded pore-scale grains; anisotropy would add orientation-dependent factors. |
| Temperature dependence of $M$ and $\alpha$ | Kept constant for simplicity. Their Arrhenius or $\sqrt{T}$ variation can be restored by writing $M(T)$ and $\alpha(T)$ without changing the form of the equation. |
| Nonlinear kinetics in large supersaturation | The pore vapor is close to saturation ($\|\rho_v-\rho_{vs}\|/\rho_{vs}\ll1$), so the linear term is the leading order and the dropped terms are $O((\Delta/\rho_{vs})^2)$. |
| Curvature dependence of the saturation density (explicit Gibbs-Thomson term in $\rho_{vs}$) | It is represented by the $M$-term ($v_n\propto-\kappa$). Matching $6\gamma\varepsilon M=6\varepsilon\alpha\rho_{vs}\gamma/(\rho_i R_vT)\cdot$(molar factors) would make it quantitative; that choice is left to the calibration of $M$. The size is $\rho_{vs}\gamma/(\rho_iR_vT\,r)\sim10^{-4}\rho_{vs}$ for $r\sim1\,\mu$m. |
| Mechanical stress, elasticity and grain motion of ice and regolith | The regolith grains and ice are treated as rigid; stresses are outside this model. |
| Thermal noise / fluctuations | Deterministic model; noise is needed only for nucleation, which is excluded. |
| Wetting of regolith grains by ice (contact angle at grain walls) | Boundary conditions at the grain surfaces are taken as natural (no-flux), as needed to drop the boundary term in S1; a wetting energy would add a surface term. |

## Questions

- The interface function $p(\phi)=30\phi^2(1-\phi)^2$ is my choice. Any localized $p$ with $\int_0^1p\,d\phi=1$ gives the same sharp-interface speed. If your model uses a different one, such as $\phi^2(1-\phi)^2$ without the 30, $\alpha$ has to be rescaled to match.
- Should $M$ be tied to $\alpha$ through the Gibbs-Thomson relation, so that curvature-driven sublimation agrees with the Kelvin equation? Or should it stay a free parameter?


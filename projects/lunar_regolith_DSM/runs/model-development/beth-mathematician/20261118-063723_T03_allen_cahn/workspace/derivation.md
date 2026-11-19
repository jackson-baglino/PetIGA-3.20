# T03: Allen-Cahn equation from the Ginzburg-Landau free energy

Notation: $f(\phi) = W\phi^2(1-\phi)^2$ (J m^-3), $\nabla = (\partial_x, \partial_y)$, $\Omega\subset\mathbb{R}^2$.

**S1 (given).** Free energy (units J per meter of depth in 2D, i.e., J m^-3 times m^2):
$$F[\phi] = \int_\Omega \Big[\tfrac{\kappa}{2}\big(\phi_x^2 + \phi_y^2\big) + W\phi^2(1-\phi)^2\Big]\,dA .$$
Units check: $\kappa\,\phi_x^2$ = J m^-1 · m^-2 = J m^-3; $W\phi^2(1-\phi)^2$ = J m^-3.

**S2 (first variation: perturb $\phi\to\phi+\varepsilon\eta$, differentiate in $\varepsilon$ at $\varepsilon=0$).**
$$\delta F = \frac{d}{d\varepsilon}F[\phi+\varepsilon\eta]\Big|_{\varepsilon=0} = \int_\Omega \Big[\kappa\,(\phi_x\eta_x + \phi_y\eta_y) + f'(\phi)\,\eta\Big]\,dA .$$

**S3 (derivative of the double well, chain rule).**
$$f'(\phi) = 2W\phi(1-\phi)^2 - 2W\phi^2(1-\phi) = 2W\phi(1-\phi)(1-2\phi)\quad[\text{J m}^{-3}].$$

**S4 (integration by parts / divergence theorem on the gradient term).**
$$\int_\Omega \kappa\,\nabla\phi\cdot\nabla\eta\,dA = \oint_{\partial\Omega}\kappa\,\eta\,\nabla\phi\cdot\mathbf n\,ds - \int_\Omega \kappa\,\eta\,(\phi_{xx}+\phi_{yy})\,dA .$$

**S5 (boundary term vanishes since $\eta=0$ on $\partial\Omega$; substitute S3, S4 into S2).**
$$\delta F = \int_\Omega \Big[2W\phi(1-\phi)(1-2\phi) - \kappa(\phi_{xx}+\phi_{yy})\Big]\eta\,dA .$$

**S6 (read off the variational derivative, $\delta F = \int_\Omega \mu\,\eta\,dA$ for all admissible $\eta$).**
$$\mu = \frac{\delta F}{\delta\phi} = 2W\phi(1-\phi)(1-2\phi) - \kappa(\phi_{xx}+\phi_{yy})\quad[\text{J m}^{-3}].$$
Units: $\kappa\phi_{xx}$ = J m^-1 · m^-2 = J m^-3. Consistent.

**S7 (Allen-Cahn, given form).**
$$\phi_t = -M\mu .$$
Units: $M\mu$ = m^3 J^-1 s^-1 · J m^-3 = s^-1 = units of $\phi_t$.

**S8 (add $M\mu$ to both sides, substitute S6).**
$$R \equiv \phi_t + M\Big[2W\phi(1-\phi)(1-2\phi) - \kappa(\phi_{xx}+\phi_{yy})\Big] = 0 .$$
The coefficient of $\phi_t$ is 1. No approximation was made.

SymPy check (`scratch/check.py`): Euler-Lagrange expression of the integrand minus S6 simplifies to 0; derivative of $f$ minus S3 simplifies to 0.

# T06: Free energy of the ice–gas phase field at coexistence

Checks: `scratch/check.py` (SymPy for the algebra and calculus, pint for the units).

**S1 (form chosen).** A Ginzburg–Landau functional with a bulk part and a gradient part:
$$F[\phi]=\int_V f\,dV,\qquad f = W(\phi)+\tfrac{\kappa}{2}|\nabla\phi|^2,\qquad W(\phi)=H\,\phi^2(1-\phi)^2 .$$
H [J m^-3] is the barrier height and κ [J m^-1] is the gradient coefficient.
- W is a double well. Its minima make bulk ice (φ = 1) and bulk gas (φ = 0) the stable states. W(0) = W(1) = 0, so the two wells are equally deep. That is the coexistence assumption: there is no driving force, so there is no term tilting one well against the other.
- The gradient term penalizes sharp changes in φ. Without it the interface would collapse to zero width. With it, the balance between the gradient term and W fixes a finite width and a finite energy per unit area.

**S2 (wells; checked).** W(0) = W(1) = 0, W′(0) = W′(1) = 0, W″(0) = W″(1) = 2H > 0, and the barrier is W(1/2) = H/16.

**S3 (Euler–Lagrange equation for a flat interface along x; δF/δφ = 0).**
$$\kappa\,\phi'' = W'(\phi).$$
Multiply by φ′ and integrate, using φ → 0 or 1 and φ′ → 0 far from the interface. This gives the first integral
$$\tfrac{\kappa}{2}\phi'^2 = W(\phi).$$

**S4 (profile; checked: both residuals are 0).**
$$\phi(x)=\tfrac12\left[1+\tanh\frac{x}{2\delta}\right],\qquad \delta=\sqrt{\kappa/(2H)} .$$

**S5 (energy per unit area).** By S3 the two terms of f are equal everywhere on the profile, so the integrand is 2W. Change the variable from x to φ, using dx = dφ/φ′ = dφ·sqrt(κ/(2W)):
$$\gamma=\int_{-\infty}^{\infty}2W\,dx=\int_0^1\sqrt{2\kappa W}\,d\phi=\frac{\sqrt{2\kappa H}}{6}.$$

**S6 (choose δ = eps and impose the energy γ; solved with SymPy).** Set sqrt(κ/(2H)) = eps and sqrt(2κH)/6 = γ. Then
$$H=\frac{3\gamma}{\varepsilon},\qquad \kappa=6\gamma\varepsilon .$$

**S7 (result; checked: integrating f over the tanh profile in x gives exactly γ; pint: every term is in J m^-3).**
$$f=\frac{3\gamma}{\varepsilon}\,\phi^2(1-\phi)^2+3\gamma\varepsilon\,|\nabla\phi|^2 .$$
The equilibrium profile is φ = ½[1 + tanh(x/(2ε))].

**Width.** I define the width as w = 1/max|φ′| = (Δφ = 1)/(slope at the center). The slope at the center is 1/(4ε), so w = **4ε**. Over this width φ goes from 0.12 to 0.88. For comparison, the 10–90 % width is 4 ln 3 · ε ≈ 4.39 ε.

## Assumptions
- Temperature is uniform, constant and at coexistence (given). So the wells are equally deep and there is no tilt.
- The interface energy is isotropic (scalar γ and scalar κ).
- γ and eps are constants.
- Only the free energy of the ice–gas pair is modeled here. The grains are not included.

## Effects left out
- **Driving force / tilt term** (for example, a term in (T − T_m) or in the vapor supersaturation times p(φ)): it is zero at coexistence by assumption.
- **Anisotropy of the ice surface energy:** the hexagonal anisotropy of ice γ is a few percent. It is not needed for the coexistence functional.
- **Contact with the regolith grains (wetting, contact angle):** these are boundary terms. They are outside this bulk density.
- **Elastic strain energy and curvature-dependent corrections:** they are small next to γ/ε for pore-scale interfaces.
- **Thermal fluctuations:** kT/(γ ε²) is about 1e-3 or less for ε of 10 nm or more, with γ ≈ 0.1 J m^-2 at 200 K.
- **Vapor density and temperature fields:** they belong to the coupled model, not to this isothermal coexistence energy.

# T06: free energy density for the ice–gas phase field

Checks: `scratch/check.py` (SymPy + pint).

**S1 (ansatz).** Ginzburg–Landau form, at coexistence (no tilt term):
$$F[\phi]=\int_V f\,dV,\qquad f = A\,W(\phi)+\frac{\kappa}{2}|\nabla\phi|^2,\qquad W=\phi^2(1-\phi)^2 .$$
A [J m^-3] is the well height (cost of mixed states); kappa [J m^-1] penalizes gradients. Units: both terms J m^-3.
- The double well makes phi = 0 (gas) and phi = 1 (ice) the stable bulk states.
- Equal minima W(0) = W(1) = 0 encode coexistence (no driving force; a phi^3-type tilt term would represent supersaturation and is left out).
- The gradient term makes sharp jumps costly, giving a finite width and an energy per area. Without it the interface would be sharp with zero energy; without W it would spread without bound.

**S2 (stationary points of W).** W'(phi) = 2phi(1-phi)(1-2phi): zeros at 0, 1/2, 1; W''(0) = W''(1) = 2 > 0 (minima), 1/2 is a maximum. (check.py; the solver was told phi > 0 so it lists 1/2, 1; phi = 0 is a zero by inspection.)

**S3 (Euler–Lagrange, flat interface along x).** delta F/delta phi = 0:
$$\kappa\,\phi'' = A\,W'(\phi).$$
Multiply by phi' and integrate from -infinity with phi -> 0, phi' -> 0:
$$\tfrac{\kappa}{2}\phi'^2 = A\,W(\phi)\;\Rightarrow\;\phi' = \sqrt{2A/\kappa}\;\phi(1-\phi).$$

**S4 (profile).** Integrate S3 with phi(0) = 1/2, delta = sqrt(kappa/(2A)) [m]:
$$\phi_0(x)=\tfrac12\left[1+\tanh\!\frac{x}{2\delta}\right].$$
check.py: Euler–Lagrange residual 0, first-integral residual 0.

**S5 (energy per area).** Using equipartition from S3, then change of variables x -> phi:
$$\gamma=\int_{-\infty}^{\infty}\Big[AW+\tfrac{\kappa}{2}\phi'^2\Big]dx=\int\kappa\phi'^2dx=\int_0^1\sqrt{2A\kappa}\,\phi(1-\phi)\,d\phi=\frac{\sqrt{2A\kappa}}{6}.$$
check.py: integral = sqrt(2) sqrt(A) sqrt(kappa)/6.

**S6 (solve for coefficients).** Impose delta = eps and S5 = gamma:
$$A=\frac{3\gamma}{\varepsilon},\qquad \kappa=6\gamma\varepsilon .$$
check.py: SymPy solve; direct integration of f over phi_0 with gamma = 1, eps = 1/3 gives 1. Units: gamma/eps = J m^-3, gamma*eps*|grad phi|^2 = J m^-3 (pint).

**S7 (result).**
$$f=\frac{3\gamma}{\varepsilon}\,\phi^2(1-\phi)^2+3\gamma\varepsilon\,|\nabla\phi|^2 .$$

**S8 (width).** Define the width as the distance over which the tangent at the midpoint goes from 0 to 1, w = 1/phi_0'(0) = **4 eps** (check.py). With the alternative 10%–90% definition, w = 4 ln 3 eps ≈ 4.39 eps.

## Assumptions
- Uniform, constant temperature at coexistence (given): equal well depths, no bulk driving force.
- Isotropic interface energy (gamma independent of interface orientation).
- Constant gamma and eps; scalar non-conserved order parameter; the 1D flat-interface solution sets the coefficients.
- Polynomial double well phi^2(1-phi)^2 (a modeling choice; other wells change the numerical factors).

## Left out
- Bulk driving force (supersaturation / latent heat tilt): zero at coexistence by assumption.
- Vapor density / mass and temperature fields in f: coupling handled elsewhere; not needed at uniform T and coexistence.
- Grain (regolith) phase and ice–grain, gas–grain wetting energies: the work order treats only the ice–gas pore space.
- Anisotropy of ice surface energy: assumed isotropic.
- Elastic strain energy in ice, curvature/Gibbs–Thomson corrections: no stress model; curvature effect emerges dynamically, not in f.
- Higher-gradient terms (second derivatives): excluded by the deliverable; they correct the energy at order (eps/R)^2 for interfaces of curvature radius R.

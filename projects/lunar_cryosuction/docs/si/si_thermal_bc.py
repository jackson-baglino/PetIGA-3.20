"""Supplementary methods: thermal boundary conditions for heating of the regolith walls.

    ../lunar_regolith_DSM/venv_docparse/bin/python docs/si/si_thermal_bc.py

Formulation: src/assembly.c (commit 4d7e0a61). Verification numbers:
studies/heatflux_bc/verification/README.md and verify_heatflux.csv.
"""
from si_common import EQ, TBL, H, P, save, title

title("Supplementary Methods: thermal boundary conditions for heating of the regolith",
      "Draft, 13 November 2026. Implementation: src/assembly.c (commit 4d7e0a61); checks: "
      "studies/heatflux_bc/verification/.")

H("Purpose", 1)
P("Heating of the regolith is represented by conditions on the pore walls: either a prescribed heat flux or "
  "a wall temperature that follows one-dimensional heat diffusion through the regolith. Both are optional and "
  "apply only to the faces named; all other faces keep their previous conditions.")

H("Energy equation and its boundary term", 1)
P("The temperature equation is $\\rho c\\,\\partial_t T = \\xi_T\\,\\nabla\\cdot(k\\nabla T) + "
  "\\xi_T\\,\\rho_{\\mathrm{ice}} L\\,\\partial_t\\phi$, with the conductivity $k(\\phi)$ of the local phase "
  "mixture and $\\xi_T = 1$. Its weak form with basis function $N_a$ contains the boundary integral")
EQ(r"-\,\xi_T\int_{\Gamma} N_a\,k\,\frac{\partial T}{\partial n}\,d\Gamma")
P("which is zero on insulated faces. The two conditions below set this integral.")

H("Prescribed heat flux", 1)
P("On the flux faces $\\Gamma_q$ the condition $k\\,\\partial T/\\partial n = q$ (with $q > 0$ heating the "
  "domain, switched on at $t \\ge t_0$) replaces the integral by a known term:")
EQ(r"R_a \leftarrow R_a - \xi_T\int_{\Gamma_q} q\,N_a\,d\Gamma")
P("The term does not depend on the solution, so the Jacobian is unchanged.")

H("Wall temperature from one-dimensional diffusion", 1)
P("For a regolith wall heated from one end, with $s$ the distance along the wall from the heated end, "
  "$\\kappa_r$ the thermal diffusivity of the regolith and a temperature step $\\Delta T$ at $s = 0$ from "
  "time $t_0$, the wall temperature is")
EQ(r"T_w(s,t) = T_{\mathrm{init}} + \Delta T\,\operatorname{erfc}\!\left(\frac{s}{2\sqrt{\kappa_r\,(t - t_0)}}\right)")
P("and $T_w = T_{\\mathrm{init}}$ before $t_0$. Two simpler forms are available: the same profile held at a "
  "fixed time $t_{\\mathrm{ref}}$, and a linear profile $T_w = T_{\\mathrm{init}} + \\Delta T\\,(1 - s/L)$ "
  "for $s < L$. Heat relaxes across a pore of width $\\ell$ in a time of order $\\ell^2/\\kappa$, much shorter "
  "than the time over which the ice changes, while the heated zone along the wall advances over times of "
  "order $s^2/\\kappa_r$; the fixed-time profile is therefore a reasonable simplification when the wall "
  "temperature changes little during a run.")

H("Weak (Nitsche) imposition of the wall temperature", 1)
P("The wall temperature is imposed weakly by the symmetric Nitsche method. On the wall faces $\\Gamma_w$ "
  "the residual of $N_a$ receives")
EQ(r"R_a \leftarrow R_a + \xi_T\int_{\Gamma_w}\Bigl[-N_a\,k\,\frac{\partial T}{\partial n}"
   r" - k\,\frac{\partial N_a}{\partial n}\,(T - T_w) + \beta\,N_a\,(T - T_w)\Bigr]d\Gamma,"
   r"\qquad \beta = \frac{C\,k(\phi)}{h_n}")
P("with $h_n$ the element size normal to the face and $C = 10$. The first term is the boundary integral of "
  "the weak form; the second makes the formulation symmetric; the third enforces $T = T_w$. The method is "
  "consistent, so a temperature field that lies in the discrete space is reproduced exactly (check c2). "
  "We use it rather than fixed nodal values because the solver's fixed boundary values are constant per "
  "face, quadratic spline control values do not interpolate $T_w$ along a face, and the time integrator "
  "evaluates the residual at intermediate times, where a weak term is evaluated consistently with the "
  "interior terms. The Jacobian includes the derivatives of these terms with respect to $T$ and $\\phi$.")

H("Verification", 1)
P("Each check compares a simulation with an exact solution or a conservation law (deterministic, scripted: "
  "verify_heatflux.sh). Every field is a single phase (all air or all ice), so phase change does not enter.",
  after=3)
TBL(["Check", "Quantity", "Measured", "Tolerance"], [
    ["Steady flux, $T$ fixed on the far face: linear $T$, slope $q/k$ (two meshes)",
     "max error / $(qL/k)$", "1.1e-14", "1e-8"],
    ["Flux into a deep slab: $T(0,t) - T_0 = (2q/k)\\sqrt{\\alpha t/\\pi}$ (air and ice, four times)",
     "relative error", "≤ 3.4e-5", "1e-3"],
    ["erfc wall temperature reproduced in the interior (two meshes)", "max error / $\\Delta T$", "≤ 4.0e-4",
     "1e-3"],
    ["Same, order of convergence in space", "order", "1.98", "> 1.5"],
    ["Linear wall temperature, steady (consistency)", "max error / $\\Delta T$", "7.1e-15", "1e-8"],
    ["Energy budget, flux runs", "stored minus input / input", "1.0e-5", "1e-4"],
    ["Boundary Jacobian against finite differences", "relative error", "2.8e-8", "1e-6"],
    ["Runs without thermal faces unchanged (with and without contact-angle walls)", "SSA_evo.dat",
     "identical", "exact"],
], [3.3, 1.4, 0.85, 0.75], size=8.5)

H("Not covered", 1)
P("Latent heat in the energy budget (no interface in the checks), mapped (-geom_file) meshes, axisymmetric "
  "runs and runs on several processes have not been checked for these conditions.")

save("SI_thermal_BC.docx")

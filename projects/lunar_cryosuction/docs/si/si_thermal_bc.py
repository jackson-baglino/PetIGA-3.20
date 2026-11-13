"""Supplementary methods: thermal boundary conditions for heating of the regolith walls.

    venv_cryosuction/bin/python docs/si/si_figures.py docs/si/fig
    ../lunar_regolith_DSM/venv_docparse/bin/python docs/si/si_thermal_bc.py

Formulation: src/assembly.c (commit 4d7e0a61); material constants: src/cryosuction_main.c lines
171-178; mixtures: src/material_properties.c. Verification: studies/heatflux_bc/verification/
(verify_heatflux.csv, rerun on main after the merge, 42 PASS). Nitsche (1971) DOI checked on Crossref.
"""
from pathlib import Path

from si_common import EQ, TBL, H, P, save, title  # noqa: I001  (sets the path for docx_helpers)
from docx_helpers import figure

title("Supplementary Methods: thermal boundary conditions for heating of the regolith",
      "Draft 3.1, 13 November 2026.")

H("Purpose and setting", 1)
P("The computational domain is a single pore (ice and gas, 2D); the regolith around it is not meshed "
  "(Fig. S2). The model solves ice sublimation and deposition with a phase field, in air at 1 atm and −20 °C "
  "(a laboratory analog). The pore is oriented along the direction of heating: its axis $x$ points away from "
  "the heated surface (to the right in Fig. S2), its end nearest the surface lies at depth $z_0$, and the walls "
  "are the pore's sides. To represent warming of the regolith, two optional conditions "
  "can be set on the walls: a prescribed heat flux, or a wall temperature equal to that of the regolith "
  "computed by one-dimensional heat diffusion from the surface. The two pore ends keep the condition chosen "
  "for the run (a fixed temperature or insulation).")
figure(str(Path(__file__).parent / "fig/fig_thermal.png"),
       "Figure S2. Left: the pore is the mesh; its walls (red) carry the thermal condition, $n$ is the outward "
       "normal and $s$ the depth below the heated surface measured along the wall; the surface lies a distance "
       "$z_0$ beyond the left end. Right: the wall temperature of Eq. (S4) against depth at three times; "
       "$\\ell_s$ is any reference length.", width=6.3)

H("Symbols", 1)
TBL(["Symbol", "Meaning"], [
    ["$x$, $L_x$, $\\ell$", "position along the pore axis; pore length; pore width (gap)"],
    ["$\\Gamma$, $\\Gamma_q$, $\\Gamma_w$", "boundary of the pore; faces with a prescribed flux; faces with a "
     "prescribed wall temperature"],
    ["$n$", "outward unit normal; $\\partial T/\\partial n = \\nabla T\\cdot n$"],
    ["$s$", "depth below the heated surface along the wall, $s = x + z_0$, with $z_0$ the depth of the pore's "
     "upper end"],
    ["$\\phi$", "ice phase field: 1 in ice, 0 in gas"],
    ["$k(\\phi)$, $\\rho c(\\phi)$", "conductivity and volumetric heat capacity, linear in $\\phi$: "
     "$k$ = 2.29 (ice), 0.02 (air) $\\mathrm{W/(m\\,K)}$; $\\rho c = 919\\times1960$ (ice), $1.341\\times1044$ (air) $\\mathrm{J/(m^3\\,K)}$"],
    ["$\\kappa = k/\\rho c$", "diffusivity in the pore: $1.3\\times10^{-6}$ (ice), $1.4\\times10^{-5}$ (air) "
     "$\\mathrm{m^2/s}$"],
    ["$\\rho_{\\mathrm{ice}}$, $L_s$", "ice density, 919 $\\mathrm{kg/m^3}$; latent heat of sublimation, "
     "$2.83\\times10^{6}$ J/kg"],
    ["$t_{\\mathrm{run}}$, $\\ell_s$", "length of a run; a reference length used only in Fig. S2"],
    ["$q$", "prescribed heat flux into the pore ($\\mathrm{W/m^2}$)"],
    ["$T_{\\mathrm{init}}$, $\\Delta T$, $t_0$", "regolith temperature before heating; temperature step at the "
     "surface; time heating starts"],
    ["$\\kappa_r$", "thermal diffusivity of the regolith (an input of the run)"],
    ["$t_{\\mathrm{ref}}$, $L_h$", "fixed time of the frozen profile; length of the linear profile"],
    ["$N_a$, $R_a$", "quadratic B-spline basis function $a$; the Galerkin residual for $N_a$ (below)"],
    ["$h_n$, $C$", "element size normal to a face; Nitsche constant"],
], [1.5, 5.0], size=8.5)

H("Discretization", 1)
P("The heat equation is")
EQ(r"\rho c\,\frac{\partial T}{\partial t} = \nabla\cdot(k\nabla T) + "
   r"\rho_{\mathrm{ice}}L_s\,\frac{\partial\phi}{\partial t}")
P("(The model carries a time-scaling factor on the conduction and latent terms; it equals 1 in every run and "
  "is omitted here.) It is solved with quadratic B-splines (isogeometric Galerkin method) and the generalized-$\\alpha$ time integrator. "
  "Multiplying by $N_a$, integrating over the pore and integrating the conduction term by parts gives the "
  "residual, which is set to zero at every step:")
EQ(r"R_a = \int\Bigl[\rho c\,\frac{\partial T}{\partial t}N_a + k\,\nabla T\cdot\nabla N_a - "
   r"\rho_{\mathrm{ice}}L_s\frac{\partial\phi}{\partial t}N_a\Bigr]dV - "
   r"\int_{\Gamma} N_a\,k\,\frac{\partial T}{\partial n}\,d\Gamma")
P("The last term is the heat entering through the boundary. It is zero on an insulated face; the two "
  "conditions below specify it.")

H("Prescribed heat flux", 1)
P("On $\\Gamma_q$ the heat entering the pore is $k\\,\\partial T/\\partial n = q$ from time $t_0$ on, so the "
  "boundary term becomes the known quantity")
EQ(r"-\int_{\Gamma_q} q\,N_a\,d\Gamma")
P("It does not depend on $T$ or $\\phi$, so the Jacobian of the residual is unchanged.")

H("Wall temperature from one-dimensional diffusion", 1)
P("The regolith along the wall is treated as a semi-infinite conductor whose surface ($s = 0$) is raised by "
  "$\\Delta T$ at $t_0$. Its temperature at depth $s$ is imposed on $\\Gamma_w$:")
EQ(r"T_w(s,t) = T_{\mathrm{init}} + \Delta T\,\mathrm{erfc}\!\left(\frac{s}{2\sqrt{\kappa_r\,(t - t_0)}}\right),"
   r"\qquad T_w = T_{\mathrm{init}}\ \ (t < t_0)")
P("Two simpler forms are available: the same profile held at $t = t_0 + t_{\\mathrm{ref}}$, and a linear "
  "profile $T_w = T_{\\mathrm{init}} + \\Delta T\\,(1 - s/L_h)$ for $s < L_h$, $T_w = T_{\\mathrm{init}}$ beyond. "
  "The coupling is one-way: heat drawn by the pore, including latent heat, does not change $T_w$; this "
  "assumes the pore draws little heat compared with what the regolith conducts.")
P("**Scales.** Heat crosses a pore $\\ell$ = 150 µm wide in $\\ell^2/\\kappa$ ≈ 2 ms (air) to 20 ms (ice), "
  "while the ice changes over days, so the pore temperature follows the wall without lag. The wall "
  "temperature at depth $z_0$ changes over about $z_0^2/\\kappa_r$ (for example 12 days for $z_0$ = 10 cm and "
  "$\\kappa_r = 10^{-8}$ $\\mathrm{m^2/s}$). Along the wall it changes at "
  "$\\partial T_w/\\partial s = -\\Delta T\\,e^{-s^2/(4\\kappa_r(t - t_0))}/\\sqrt{\\pi\\kappa_r(t - t_0)}$; this "
  "along-wall gradient is what the pore feels, and it varies little across a pore a few millimeters long. The frozen profile "
  "replaces the time-dependent one when the run is much shorter than $t_{\\mathrm{ref}}$: the profile changes "
  "with $\\sqrt{t - t_0}$, so over a run of length $t_{\\mathrm{run}}$ its argument changes by about "
  "$t_{\\mathrm{run}}/(2t_{\\mathrm{ref}})$.")

H("Imposing the wall temperature weakly (Nitsche)", 1)
P("The wall temperature is imposed with the symmetric Nitsche method (Nitsche, 1971). On $\\Gamma_w$ the "
  "boundary term of $R_a$ is replaced by")
EQ(r"\int_{\Gamma_w}\Bigl[-N_a\,k\,\frac{\partial T}{\partial n}"
   r" - k\,\frac{\partial N_a}{\partial n}\,(T - T_w) + \beta\,N_a\,(T - T_w)\Bigr]d\Gamma,"
   r"\qquad \beta = \frac{C\,k(\phi)}{h_n},\quad C = 10")
P("The first term is the heat entering through the wall; the second keeps the discrete problem symmetric; the "
  "third drives $T$ toward $T_w$, with a strength that scales with the conductivity and the inverse element "
  "size. The method is consistent: the exact solution of the boundary-value problem satisfies these equations, "
  "so, for example, a linear temperature field is reproduced to round-off (check H5). Stability requires $C$ "
  "above a constant that depends on the spline degree; $C = 10$ was chosen and passed all checks, but its "
  "sensitivity was not tested. The terms are evaluated at the same intermediate time of each generalized-$\\alpha$ "
  "step as the rest of the residual, and their derivatives with respect to $T$ and $\\phi$ are included in "
  "the Jacobian. Fixing nodal values instead would not work here: spline coefficients on a face are not "
  "temperatures on the face, the solver's fixed values are one constant per face, and fixed values have zero "
  "rate of change in the time integrator.")

H("Verification", 1)
P("Each check runs a simulation and compares it with an exact solution or a conservation law. Every field is "
  "a single phase (all air or all ice). Errors are maximum differences over the domain, divided by the scale "
  "given.", after=3)
TBL(["ID", "Setup", "Quantity", "Measured", "Tolerance"], [
    ["H1", "Air box $200\\times20$ µm; flux $q$ = 100 $\\mathrm{W/m^2}$ on one end, $T$ fixed on the other, walls insulated; "
     "steady; element sizes 5 and 2.5 µm", "error in $T$ / $(qL_x/k)$", "1.1e-14", "1e-8"],
    ["H2", "Box 600 µm long, flux on one end; air ($q$ = 100 $\\mathrm{W/m^2}$, $t$ = 0.25–2 ms) and ice ($q$ = 1000 $\\mathrm{W/m^2}$, "
     "$t$ = 2.5–20 ms): heated depth ≤ 170 µm", "error / surface rise $(2q/k)\\sqrt{\\kappa t/\\pi}$", "≤ 3.4e-5",
     "1e-3"],
    ["H3", "Air box; erfc wall temperature on all faces with $\\kappa_r = \\kappa_{\\mathrm{air}}$; element "
     "sizes 8 and 4 µm; four times", "error in $T$ / $\\Delta T$", "≤ 4.0e-4", "1e-3"],
    ["H4", "H3 at the time where a run with $\\Delta t/4$ changes the error by 0.2 % (spatial error dominant)",
     "order, element size halved", "1.98", "> 1.5"],
    ["H5", "Linear wall temperature on all faces; steady", "error in $T$ / $\\Delta T$", "7.1e-15", "1e-8"],
    ["H6", "H2: stored heat against heat supplied, all output times", "difference / heat supplied",
     "1.0e-5", "1e-4"],
    ["H7", "Nitsche terms of the Jacobian against finite differences", "relative error", "2.8e-8", "1e-6"],
    ["H8", "Runs without thermal faces, with and without contact-angle walls",
     "output against the previous code", "identical bit for bit", "exact"],
], [0.35, 2.85, 1.45, 0.85, 0.75], size=8.5)
P("H6: the difference arises in the first time step, where the integrator starts from a zero heating rate "
  "while the flux switches on, and is the same at every later output time. H4: quadratic splines can reach "
  "third order in this norm for smooth solutions; the measured second order has not been traced to a cause "
  "(the erfc profile has a step at $s = 0$, $t = t_0$).", after=4)

H("Not covered", 1)
P("Ice–gas interfaces at a heated wall (latent heat in the budget, the jump of $k$ inside $\\beta$ where an "
  "interface meets the wall, the $\\phi$ part of the Jacobian with $\\phi$ varying), curved-wall meshes, "
  "axisymmetric runs and runs on several processes.")

H("Reference", 2)
P("Nitsche, J. (1971). Über ein Variationsprinzip zur Lösung von Dirichlet-Problemen bei Verwendung von "
  "Teilräumen, die keinen Randbedingungen unterworfen sind. Abh. Math. Sem. Univ. Hamburg 36, 9–15. "
  "doi:10.1007/BF02995904", size=8.5)
P("Code: src/assembly.c and studies/heatflux_bc/verification/ (commit 4d7e0a61).", size=8.5)

save("SI_thermal_BC.docx")

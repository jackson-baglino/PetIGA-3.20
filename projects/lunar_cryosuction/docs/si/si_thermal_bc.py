"""Supplementary methods: thermal boundary conditions for heating of the regolith walls.

    venv_cryosuction/bin/python docs/si/si_figures.py docs/si/fig
    ../lunar_regolith_DSM/venv_docparse/bin/python docs/si/si_thermal_bc.py

Formulation: src/assembly.c (commit 4d7e0a61); property mixtures: src/material_properties.c.
Verification numbers: studies/heatflux_bc/verification/README.md and verify_heatflux.csv (rerun on
main after the merge, 42 PASS).
"""
from pathlib import Path

from si_common import EQ, TBL, H, P, save, title  # noqa: I001  (sets the path for docx_helpers)
from docx_helpers import figure

title("Supplementary Methods: thermal boundary conditions for heating of the regolith",
      "Draft 2, 13 November 2026. Implementation: src/assembly.c (commit 4d7e0a61); checks: "
      "studies/heatflux_bc/verification/.")

H("Purpose", 1)
P("The computational domain is the pore (ice and gas); the regolith around it is not meshed (Fig. S2). To "
  "represent heating of the regolith, for example warming that starts at the surface and diffuses downward, "
  "two optional conditions can be set on the pore walls: a prescribed heat flux, or a wall temperature that "
  "follows one-dimensional heat diffusion through the regolith. Faces not named keep the condition chosen for "
  "the run (a fixed temperature on the pore ends, or insulation).")
figure(str(Path(__file__).parent / "fig/fig_thermal.png"),
       "Figure S2. Left: the pore is the mesh; its walls (red) are faces of the mesh, $n$ is the outward normal "
       "and $s$ the distance along the wall from the heated end. Right: the wall temperature of Eq. (S4) at "
       "three times; $\\ell_s$ is any reference length along the wall.", width=6.3)

H("Symbols", 1)
TBL(["Symbol", "Meaning"], [
    ["$\\phi$", "ice phase field: 1 in ice, 0 in gas"],
    ["$k(\\phi)$, $\\rho c(\\phi)$", "thermal conductivity and volumetric heat capacity, linear in $\\phi$ between "
     "the values for ice and air ($k_{\\mathrm{ice}}/k_{\\mathrm{air}} = 115$)"],
    ["$\\rho_{\\mathrm{ice}}$, $L_s$", "ice density; latent heat of sublimation"],
    ["$\\xi_T$", "the model's scaling factor for the heat equation; $\\xi_T = 1$ in every run"],
    ["$n$", "outward unit normal of a face; $\\partial T/\\partial n = \\nabla T\\cdot n$"],
    ["$N_a$, $R_a$", "spline basis function $a$ (quadratic); the residual of the heat equation for $N_a$"],
    ["$q$", "prescribed heat flux into the pore [$\\mathrm{W/m^2}$]; positive heats the pore"],
    ["$T_{\\mathrm{init}}$, $\\Delta T$, $t_0$", "wall temperature before heating; temperature step at the "
     "heated end; time heating starts"],
    ["$\\kappa_r$", "thermal diffusivity of the regolith [$\\mathrm{m^2/s}$]"],
    ["$\\kappa$, $\\ell$", "thermal diffusivity in the pore; pore width"],
], [1.6, 4.9], size=8.5)

H("Heat equation and its boundary term", 1)
P("The temperature satisfies")
EQ(r"\rho c\,\frac{\partial T}{\partial t} = \xi_T\,\nabla\cdot(k\nabla T) + "
   r"\xi_T\,\rho_{\mathrm{ice}}L_s\,\frac{\partial\phi}{\partial t}")
P("Multiplying by $N_a$ and integrating by parts gives, besides the volume terms, the boundary term")
EQ(r"-\,\xi_T\int_{\Gamma} N_a\,k\,\frac{\partial T}{\partial n}\,d\Gamma")
P("over the boundary $\\Gamma$ of the pore. On an insulated face the term is zero. The two conditions below "
  "specify it on the faces where they are set.")

H("Prescribed heat flux", 1)
P("On the faces $\\Gamma_q$ the heat flowing into the pore is $k\\,\\partial T/\\partial n = q$ from time "
  "$t_0$ on. The boundary term becomes a known quantity and is added to the residual:")
EQ(r"R_a \leftarrow R_a - \xi_T\int_{\Gamma_q} q\,N_a\,d\Gamma")
P("It does not depend on $T$ or $\\phi$, so the Jacobian is unchanged.")

H("Wall temperature from one-dimensional diffusion", 1)
P("The regolith along a wall is treated as a semi-infinite one-dimensional conductor whose end ($s = 0$) is "
  "raised by $\\Delta T$ at time $t_0$. Its temperature, imposed on the wall faces $\\Gamma_w$, is")
EQ(r"T_w(s,t) = T_{\mathrm{init}} + \Delta T\,\mathrm{erfc}\!\left(\frac{s}{2\sqrt{\kappa_r\,(t - t_0)}}\right),"
   r"\qquad T_w = T_{\mathrm{init}}\ \ (t < t_0)")
P("Two simpler forms are available: the same profile held at one time $t_0 + t_{\\mathrm{ref}}$, and a linear "
  "profile $T_w = T_{\\mathrm{init}} + \\Delta T\\,(1 - s/L_h)$ for $s < L_h$ and $T_w = T_{\\mathrm{init}}$ "
  "beyond. The coupling is one-way: the pore does not change $T_w$, so heat drawn by the pore, including latent "
  "heat of deposition at the wall, is supplied by the wall without limit.")
P("Three times set which form applies. Heat relaxes across the pore in about $\\ell^2/\\kappa$; the ice "
  "changes over days or longer; the heated zone along the wall advances over about $s^2/\\kappa_r$. When "
  "$\\ell^2/\\kappa$ is much shorter than the ice time, the pore temperature follows the wall without lag. The "
  "profile held at $t_{\\mathrm{ref}}$ is a valid replacement for the time-dependent one when the run is much "
  "shorter than $t_{\\mathrm{ref}}$, because the erfc profile changes with $\\sqrt{t - t_0}$.")

H("Imposing the wall temperature weakly (Nitsche)", 1)
P("On $\\Gamma_w$ the residual of $N_a$ receives")
EQ(r"R_a \leftarrow R_a + \xi_T\int_{\Gamma_w}\Bigl[-N_a\,k\,\frac{\partial T}{\partial n}"
   r" - k\,\frac{\partial N_a}{\partial n}\,(T - T_w) + \beta\,N_a\,(T - T_w)\Bigr]d\Gamma,"
   r"\qquad \beta = \frac{C\,k(\phi)}{h_n}")
P("with $h_n$ the physical size of the boundary element normal to the face and $C = 10$. The first term is the "
  "boundary term of the weak form, the second makes the formulation symmetric, and the third drives $T$ to "
  "$T_w$. The method is consistent: an exact solution that lies in the spline space satisfies these equations "
  "exactly (check H5). The terms are evaluated at the same intermediate time as the rest of the residual in "
  "each step of the time integrator, and their Jacobian with respect to $T$ and $\\phi$ is included. We use "
  "this rather than fixing nodal values because the solver's fixed values are one constant per face, "
  "quadratic spline coefficients do not equal the temperature on the face, and the time integrator sets the "
  "rate of change of a fixed value to zero, which would drop the heat-capacity term of a time-varying wall "
  "temperature. $C = 10$ is a choice; its sensitivity was not tested.")

H("Verification", 1)
P("Each check compares a scripted simulation (verify_heatflux.sh) with an exact solution or a conservation law. "
  "Every field is a single phase (all air or all ice); errors are maximum differences over the domain.",
  after=3)
TBL(["ID", "Setup", "Quantity", "Measured", "Tolerance"], [
    ["H1", "Box; flux $q$ on one end, $T$ fixed on the other, insulated walls; steady; two meshes",
     "error in $T$ / $(q L_x/k)$", "1.1e-14", "1e-8"],
    ["H2", "Long box (600 µm), flux $q$ on one end; air and ice; four output times",
     "error vs $(2q/k)\\sqrt{\\kappa t/\\pi}$", "≤ 3.4e-5", "1e-3"],
    ["H3", "Box with the erfc wall temperature on all faces and $\\kappa = \\kappa_r$; two meshes",
     "error in $T$ / $\\Delta T$", "≤ 4.0e-4", "1e-3"],
    ["H4", "H3 at $t_f/4$, mesh halved", "order of convergence", "1.98", "> 1.5"],
    ["H5", "Linear wall temperature on all faces; steady", "error in $T$ / $\\Delta T$", "7.1e-15", "1e-8"],
    ["H6", "H2: stored heat against heat supplied", "difference / supplied", "1.0e-5", "1e-4"],
    ["H7", "Nitsche terms of the Jacobian against finite differences", "relative error", "2.8e-8", "1e-6"],
    ["H8", "Runs that name no thermal face, with and without contact-angle walls",
     "output against the previous code", "identical bit for bit", "exact"],
], [0.35, 2.75, 1.45, 0.95, 0.75], size=8.5)
P("The tolerances lie well above the round-off or discretization error expected for each check and well below "
  "an error a wrong sign or factor would give. The difference in H6 arises in the first step only: the time "
  "integrator starts from a zero heating rate while the flux switches on at $t = 0$, and the offset is constant "
  "afterwards. The order in H4 is 2 although quadratic splines can give 3 for smooth solutions; the erfc "
  "solution has a step at $s = 0$ and $t = t_0$, and the cause of the lower order has not been isolated.",
  after=4)

H("Not covered", 1)
P("Ice–gas interfaces at a heated wall (latent heat in the budget, the jump of $k$ inside $\\beta$ at a contact "
  "line, the $\\phi$ part of the Jacobian with $\\phi$ varying), mapped (curved-wall) meshes, axisymmetric runs, "
  "and runs on several processes have not been checked for these conditions.")

save("SI_thermal_BC.docx")

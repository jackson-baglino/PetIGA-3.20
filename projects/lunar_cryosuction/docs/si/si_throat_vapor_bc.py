"""Supplementary methods: vapor boundary conditions for the throat (type 5) runs.

    ../lunar_regolith_DSM/venv_docparse/bin/python docs/si/si_throat_vapor_bc.py

Numbers: studies/cryosuction/scripts/plan_estimates.py --wt 4e-5 8e-5 (section 8c,
commit 96c6a2ab), copied into the table below with the command that produced them.
"""
from si_common import EQ, TBL, H, P, save, title

title("Supplementary Methods: vapor boundary conditions for the throat simulations",
      "Draft, 13 November 2026. Model estimates used to choose boundary conditions; "
      "source: studies/cryosuction/scripts/plan_estimates.py, sections 8 to 8c.")

H("Problem", 1)
P("A lens that spans the pore changes only through its two faces. With no other ice in the pore, each face "
  "exchanges vapor only through the gas between it and the end of the pore. The vapor condition at the ends "
  "therefore sets whether the lens grows or sublimates, and it can bias the motion of the lens. We use a "
  "periodic cell for the runs without a temperature gradient and fixed end humidities, chosen below, for the "
  "runs with a gradient.")

H("Face equilibrium", 1)
P("Let $g(x)$ be the local gap, $\\theta$ the contact angle (60°), and $\\omega_i$ the angle by which the "
  "walls open toward the gas at face $i$ ($i = w$ for the warm face, $c$ for the cold face). The face is an "
  "arc of curvature")
EQ(r"\chi_i = \frac{2\cos(\theta + \omega_i)}{g_i}")
P("and, with the Gibbs–Thomson relation linearized in the capillary length $d_0$, its equilibrium vapor "
  "density relative to $\\rho_{vs}(T)$ is")
EQ(r"h_i = 1 - d_0\,\chi_i")

H("Flux at each face", 1)
P("Per unit depth, face $i$ gains ice area at")
EQ(r"Q_i = \frac{\Delta_i}{R_i},\qquad R_i = \frac{\beta_0}{g_i} + K\int_{\mathrm{column}\ i}\frac{dx}{g(x)},"
   r"\qquad K = \frac{\rho_{\mathrm{ice}}}{D_v\,\rho_{vs}}")
P("$R_i$ adds the attachment resistance of the face ($\\beta_0$) and the diffusion resistance of the gas "
  "column between the face and its end. For a uniform gap $g$ and column length $l$, $Q_i/g$ reduces to "
  "$\\Delta/(\\beta_0 + K l)$, the face speed that the straight-pore runs followed within 1 %. With the ends "
  "held at $\\rho_v = h\\,\\rho_{vs}(T_{\\mathrm{end}})$, an isothermal lens, and "
  "$s = d\\ln\\rho_{vs}/dT$, the supersaturations seen by the faces are")
EQ(r"\Delta_w = h - h_w + s\,\Delta T_w,\qquad \Delta_c = h - h_c - s\,\Delta T_c")
P("where $\\Delta T_w$ and $\\Delta T_c$ are the temperature drops across the warm and cold gas columns "
  "(both zero without a gradient).")

H("Growth rate of the lens", 1)
P("Adding the two faces gives")
EQ(r"\frac{dA}{dt} = Q_w + Q_c = a\,\bigl(h - h_0\bigr),\qquad a = \frac{1}{R_w} + \frac{1}{R_c}")
EQ(r"h_0 = \frac{h_w/R_w + h_c/R_c + s\,(\Delta T_c/R_c - \Delta T_w/R_w)}{a}")
P("$h_0$ is the end humidity at which the lens neither grows nor sublimates at that instant. It is a mean of "
  "the two face equilibria weighted by $1/R_i$, so the face with the shorter gas column counts more. "
  "$h_0$ changes as the lens moves, so a constant $h$ can cancel growth only on average.")

H("Choice of a constant end humidity", 1)
P("We choose $h$ to minimize the squared growth rate over the run:")
EQ(r"J(h) = \int_0^{t_f}\Bigl(\frac{dA}{dt}\Bigr)^2 dt = \int_0^{t_f} a^2\,(h - h_0)^2\,dt")
EQ(r"\frac{dJ}{dh} = 2\int_0^{t_f} a^2\,(h - h_0)\,dt = 0\;\;\Rightarrow\;\;"
   r"h_{LS} = \frac{\int_0^{t_f} a^2 h_0\,dt}{\int_0^{t_f} a^2\,dt},\qquad "
   r"\frac{d^2J}{dh^2} = 2\int_0^{t_f} a^2\,dt > 0")
P("so $h_{LS}$ is the unique minimum. Requiring instead zero net change of area over the run gives")
EQ(r"h_{\mathrm{net}} = \frac{\int_0^{t_f} a\,h_0\,dt}{\int_0^{t_f} a\,dt}")
P("Both need $a(t)$ and $h_0(t)$ along the path of the lens. We take the path from the same model, with "
  "the lens center moving toward $x = 0$ at")
EQ(r"v = \tfrac12\Bigl(\frac{Q_w}{g_w} - \frac{Q_c}{g_c}\Bigr)")
P("and the lens length set by its area; the path depends weakly on $h$, so path and $h$ are iterated three "
  "times. Table S1 compares these choices with the unweighted mean curvature, "
  "$h = 1 - d_0\\langle\\chi\\rangle$, and with saturation, $h = 1$.")

P("**Table S1.** Change of lens area after 300 days for each end humidity. Pore 400 µm long with walls "
  "$g = \\sqrt{w_t^2 + c^2\\xi^2}$ (150 µm at the ends, no flat sections), lens 100 µm long started with its "
  "center 70 µm on the cold side of the throat, $\\theta$ = 60°. Command: plan_estimates.py --wt 4e-5 8e-5.",
  size=8.5, after=2)
TBL(["Throat, G", "$h_{LS}$", "area at $h_{LS}$", "area at $h_{\\mathrm{net}}$",
   "area at $1 - d_0\\langle\\chi\\rangle$", "area at $h = 1$"], [
    ["40 µm, 0 K/m", "0.999 991 34", "+0.35 %", "+0.08 %", "−30.2 %", "+86.5 %"],
    ["80 µm, 0 K/m", "0.999 992 72", "+0.03 %", "0.00 %", "−11.1 %", "+68.5 %"],
    ["40 µm, 1.56 K/m", "0.999 990 52", "+0.39 %", "+0.09 %", "−14.3 %", "+83.3 %"],
    ["80 µm, 1.47 K/m", "0.999 992 06", "+0.58 %", "+0.06 %", "−7.0 %", "+59.0 %"],
], [1.15, 1.05, 0.95, 1.05, 1.15, 0.95], size=8.5)

H("Periodic cell without a gradient", 1)
P("Even at $h_{LS}$ the ends bias the motion: in the pore of Table S1 with a 40 µm throat and no gradient, "
  "the lens passes the throat by 37 µm within 300 days, because the face nearer an end has the smaller "
  "$R_i$ and its mismatch $h - h_i$ drives a flux comparable to the curvature difference between the faces. "
  "For the runs without a gradient we therefore use a cell that is periodic in $x$, with the throat at its "
  "center. The two faces are then joined by one gas column, and")
EQ(r"Q_w = -Q_c = \frac{h_c - h_w}{\beta_0/g_w + \beta_0/g_c + K\int_{\mathrm{gas}} dx/g(x)}")
P("The ice area is conserved exactly, the lens moves toward the face of larger curvature, and it comes to "
  "rest where $\\chi_w = \\chi_c$, centered on the throat of a symmetric cell. No vapor parameter is set.")

H("Runs with a gradient", 1)
P("A linear temperature profile cannot be periodic, so the runs with a gradient keep fixed ends at "
  "$\\rho_v = h_{\\mathrm{net}}\\,\\rho_{vs}(T_{\\mathrm{end}})$, with $h_{\\mathrm{net}}$ computed for each "
  "run. A straight pore with the same ends, lens and gradient serves as a control for the drift caused by "
  "the ends.")

H("Assumptions", 1)
P("Isothermal ice; the static contact angle holds at moving contact lines; quasi-steady, one-dimensional "
  "gas columns; lens length from its area. The estimates choose boundary conditions only; the simulations "
  "test them.")

save("SI_throat_vapor_BC.docx")

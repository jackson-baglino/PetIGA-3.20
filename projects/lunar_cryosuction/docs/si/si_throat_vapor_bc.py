"""Supplementary methods: vapor boundary conditions for the throat (type 5) runs.

    venv_cryosuction/bin/python docs/si/si_figures.py docs/si/fig
    ../lunar_regolith_DSM/venv_docparse/bin/python docs/si/si_throat_vapor_bc.py

Numbers: studies/cryosuction/scripts/plan_estimates.py --wt 4e-5 8e-5 (section 8c, commit
3b544302) and its printed constants at -20 C.
"""
from si_common import EQ, TBL, H, P, save, title
from pathlib import Path

from docx_helpers import figure

title("Supplementary Methods: vapor boundary conditions for the throat simulations",
      "Draft 2, 13 November 2026. Estimates used to choose the boundary conditions; source: "
      "studies/cryosuction/scripts/plan_estimates.py, section 8c.")

H("Purpose", 1)
P("In the throat simulations an ice lens spans a pore that narrows to a throat (Fig. S1). The lens can move "
  "only by sublimation at one face and deposition at the other; the vapor travels through the gas between each "
  "face and the end of the pore. The vapor density held at the two ends therefore decides whether the lens "
  "also grows or shrinks while it moves. This section derives the end condition that keeps the amount of ice "
  "constant, and explains why the runs without a temperature gradient use a periodic pore instead.")
figure(str(Path(__file__).parent / "fig/fig_throat.png"), "Figure S1. Pore with a throat and an ice lens (2D, per unit depth). Each face "
       "exchanges vapor with the end of the pore on its side through the gas column between them. The ends "
       "hold the vapor density at a fraction h of saturation at the local temperature.", width=6.3)

H("Symbols", 1)
TBL(["Symbol", "Meaning", "Value at −20 °C"], [
    ["$g(x)$, $w_t$", "gap between the walls; its smallest value (the throat)", ""],
    ["$\\theta$", "contact angle, measured through the ice", "60°"],
    ["$\\omega_i$", "angle between the walls at face $i$ and the pore axis, positive when the gap widens "
     "into the face's gas, negative when it narrows; $\\tan\\omega_i = \\pm\\tfrac12\\,dg/dx$", ""],
    ["$\\chi_i$", "curvature of face $i$, positive when the face is concave seen from the gas", ""],
    ["$\\rho_{vs}(T)$", "saturation vapor density over flat ice", "$8.49\\times10^{-4}\\ \\mathrm{kg/m^3}$"],
    ["$h$", "vapor density held at the ends, as a fraction of $\\rho_{vs}$ at the end's temperature", ""],
    ["$d_0$", "capillary length, $\\gamma a^3/(k_B T)$", "$1.015\\times10^{-9}$ m"],
    ["$\\beta_0$", "kinetic coefficient: a face moves at $v = \\Delta/\\beta_0$ for supersaturation "
     "$\\Delta$", "$7.94\\times10^{5}$ s/m"],
    ["$D_v$, $\\rho_{\\mathrm{ice}}$", "vapor diffusivity in the gas; ice density",
     "$1.90\\times10^{-5}\\ \\mathrm{m^2/s}$; $919\\ \\mathrm{kg/m^3}$"],
    ["$K$", "$\\rho_{\\mathrm{ice}}/(D_v\\,\\rho_{vs})$, diffusion resistance per length", "$5.71\\times10^{10}\\ \\mathrm{s/m^2}$"],
    ["$s$", "$d\\ln\\rho_{vs}/dT$", "0.0961 1/K"],
    ["$A$", "lens area (per unit depth)", ""],
    ["$G$", "imposed temperature gradient, warm end at $x = 0$", ""],
], [1.05, 3.85, 1.6], size=8.5)

H("Step 1: equilibrium of each face", 1)
P("A face spans the local gap $g_i$ and meets both walls at the contact angle; with the walls inclined by "
  "$\\omega_i$ it is a circular arc of curvature")
EQ(r"\chi_i = \frac{2\cos(\theta + \omega_i)}{g_i}")
P("(the arc radius is $g_i/[2\\cos(\\theta+\\omega_i)]$). By the Gibbs–Thomson relation, linearized in "
  "$d_0\\chi_i \\ll 1$, the vapor density in equilibrium with the face, as a fraction of $\\rho_{vs}$, is")
EQ(r"h_i = 1 - d_0\,\chi_i")
P("A concave face ($\\chi_i > 0$) is in equilibrium below saturation; a narrower gap gives a lower $h_i$.")

H("Step 2: rate at which each face gains ice", 1)
P("Vapor reaches face $i$ through two resistances in series: attachment at the face and diffusion along the "
  "gas column to the end (Fig. S1). Per unit depth, the face gains ice area at")
EQ(r"Q_i = \frac{\Delta_i}{R_i},\qquad R_i = \frac{\beta_0}{g_i} + K\int_{\mathrm{column}\ i}\frac{dx}{g(x)}")
P("where $\\Delta_i$ is the supersaturation between the end and the face, and the integral runs from the face "
  "to its end. For a uniform gap $g$ and a column of length $l$ this gives the face speed "
  "$Q_i/g = \\Delta_i/(\\beta_0 + K l)$, which the straight-pore (type 3) simulations of this study followed "
  "within 1 %.")
P("The ends hold $\\rho_v = h\\,\\rho_{vs}(T_{\\mathrm{end}})$. Ice conducts heat 115 times better than air, so "
  "the lens is nearly isothermal and the imposed temperature difference $G L$ falls across the two gas "
  "columns, split in proportion to their thermal resistances, $\\Delta T_i = G L\\,I_i/(I_w + I_c)$ with "
  "$I_i = \\int_{\\mathrm{column}\\ i} dx/g$. Linearizing $\\rho_{vs}(T_{\\mathrm{end}})/\\rho_{vs}"
  "(T_{\\mathrm{lens}}) = 1 \\pm s\\,\\Delta T_i$, the faces see")
EQ(r"\Delta_w = h - h_w + s\,\Delta T_w,\qquad \Delta_c = h - h_c - s\,\Delta T_c")
P("($w$: warm face, toward $x = 0$; $c$: cold face; both $\\Delta T_i \\ge 0$; both are zero without a gradient).")

H("Step 3: growth of the lens and the growth-free end humidity", 1)
P("The lens area changes at the sum of the two face rates:")
EQ(r"\frac{dA}{dt} = Q_w + Q_c = a\,\bigl(h - h_0\bigr),\qquad a = \frac{1}{R_w} + \frac{1}{R_c},\qquad "
   r"h_0 = \frac{h_w/R_w + h_c/R_c + s\,(\Delta T_c/R_c - \Delta T_w/R_w)}{a}")
P("$h_0$ is the end humidity at which the lens neither grows nor shrinks at a given moment: the two face "
  "equilibria averaged with weights $1/R_i$, plus the effect of the gradient. As the lens moves, $\\chi_i$ and "
  "$R_i$ change, so $h_0$ changes, and a constant $h$ can keep the amount of ice constant only over the run as "
  "a whole.")

H("Step 4: choice of a constant end humidity", 1)
P("We require that the lens area be the same at the end of the run (time $t_f$) as at the start. Integrating "
  "Step 3 and setting the result to zero gives")
EQ(r"\int_0^{t_f}\frac{dA}{dt}\,dt = 0\;\;\Rightarrow\;\; h_{\mathrm{net}} = "
   r"\frac{\int_0^{t_f} a\,h_0\,dt}{\int_0^{t_f} a\,dt}")
P("As a check that this also keeps the lens close to constant during the run, the $h$ that minimizes the "
  "squared growth rate, $J(h) = \\int_0^{t_f}(dA/dt)^2\\,dt$, is")
EQ(r"\frac{dJ}{dh} = 2\int_0^{t_f} a^2\,(h - h_0)\,dt = 0\;\;\Rightarrow\;\; h_{LS} = "
   r"\frac{\int_0^{t_f} a^2 h_0\,dt}{\int_0^{t_f} a^2\,dt}\qquad\Bigl(\frac{d^2J}{dh^2} = "
   r"2\int_0^{t_f} a^2\,dt > 0\Bigr)")
P("Both integrals need $a(t)$ and $h_0(t)$ along the path of the lens. We compute the path with the same "
  "model: each face moves at $Q_i/g_i$, the lens center at their mean, and the lens length follows from its "
  "area ($A = \\int g\\,dx$ over the lens). The path depends on $h$, so $h$ is updated from the path until it "
  "changes by less than $10^{-12}$ (7 to 11 passes).")
P("**Table S1.** Lens area after 300 days for four end humidities. Pore 400 µm long with gap "
  "$g = \\sqrt{w_t^2 + c^2 (x - 200\\ \\mu\\mathrm{m})^2}$ and $g$ = 150 µm at both ends; lens 100 µm long, "
  "started with its center 70 µm on the cold side of the throat; $\\theta$ = 60°. "
  "$1 - d_0\\langle\\chi\\rangle$ uses the mean curvature of the two faces over the run, unweighted.",
  size=8.5, after=2)
TBL(["$w_t$, $G$", "$h_{\\mathrm{net}}$", "$h_{LS}$", "area at $h_{\\mathrm{net}}$", "at $h_{LS}$",
     "at $1 - d_0\\langle\\chi\\rangle$", "at $h = 1$"], [
    ["40 µm, 0 K/m", "0.999 990 87", "0.999 990 87", "0.00 %", "−0.03 %", "−28.8 %", "+88.7 %"],
    ["40 µm, 1 K/m", "0.999 991 47", "0.999 991 44", "0.00 %", "−0.27 %", "−11.1 %", "+77.4 %"],
    ["40 µm, 2 K/m", "0.999 991 69", "0.999 991 71", "0.00 %", "+0.10 %", "−26.8 %", "+78.3 %"],
    ["80 µm, 0 K/m", "0.999 992 27", "0.999 992 27", "0.00 %", "+0.02 %", "−12.8 %", "+70.2 %"],
    ["80 µm, 1 K/m", "0.999 991 15", "0.999 991 16", "0.00 %", "+0.12 %", "−2.4 %", "+68.3 %"],
    ["80 µm, 2 K/m", "0.999 991 57", "0.999 991 59", "0.00 %", "+0.10 %", "−4.6 %", "+63.0 %"],
], [1.05, 1.05, 1.05, 0.95, 0.75, 0.95, 0.75], size=8.5)
P("$h_{\\mathrm{net}}$ and $h_{LS}$ agree to $3\\times10^{-8}$, so the zero-net choice also keeps the growth "
  "rate near its smallest. The unweighted mean curvature leaves the lens 2–29 % smaller, because the face with "
  "the shorter gas column (smaller $R_i$) dominates $h_0$. These values are model estimates: a change of "
  "$10^{-7}$ in $h$ changes the area by about 1 % over 300 days, which is below the error of the "
  "approximations listed at the end; the area change in each simulation is therefore checked directly.",
  after=4)

H("Runs without a temperature gradient: periodic pore", 1)
P("Even at $h_{\\mathrm{net}}$ the ends bias the motion. In Table S1 ($w_t$ = 40 µm, no gradient) the lens does "
  "not stop at the throat: its center ends 34 µm past it. The face nearer an end has the smaller $R_i$, so the "
  "mismatch $h - h_i$ at that face drives a flux as large as the curvature difference between the two faces. "
  "For the runs without a gradient we therefore use a pore that is periodic in $x$: one cell of an infinite "
  "chain of identical throats, with the throat at the cell center and gap "
  "$g = \\sqrt{w_t^2 + (cL/\\pi)^2\\sin^2(\\pi\\xi/L)}$, where $\\xi$ is the distance from the throat; this "
  "equals the throat profile of Table S1 near the throat and joins smoothly at the cell edges. The two faces "
  "are then connected by one gas column through the neighboring cell, and")
EQ(r"Q_w = -Q_c = \frac{h_c - h_w}{\beta_0/g_w + \beta_0/g_c + K\int_{\mathrm{gas}} dx/g(x)}")
P("The ice area is conserved exactly. Vapor moves from the face with the higher equilibrium humidity (smaller "
  "curvature) to the other, so the lens moves toward the narrower gap and stops where the two faces have equal "
  "curvature: centered on the throat (stable) or centered on the widest point (unstable). No end humidity is "
  "needed.")

H("Runs with a temperature gradient", 1)
P("A linear temperature profile cannot be periodic, so these runs keep the ends at "
  "$\\rho_v = h_{\\mathrm{net}}\\,\\rho_{vs}(T_{\\mathrm{end}})$, with $h_{\\mathrm{net}}$ computed as above for "
  "each run. A straight pore with the same lens, gradient and its own $h_{\\mathrm{net}}$ gives the drift caused "
  "by the ends alone, which is subtracted from the drift in the throat pore.")

H("Approximations", 1)
P("Isothermal ice; one-dimensional, quasi-steady gas columns; the face length taken as the gap $g_i$ (the arc "
  "is about 5 % longer at 60°); the static contact angle at moving contact lines; Gibbs–Thomson linearized in "
  "$d_0\\chi$. The estimates set the boundary conditions; the simulations do not use them otherwise.")

save("SI_throat_vapor_BC.docx")

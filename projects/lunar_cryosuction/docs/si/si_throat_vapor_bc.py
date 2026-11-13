"""Supplementary methods: vapor boundary conditions for the throat (type 5) runs.

    venv_cryosuction/bin/python docs/si/si_figures.py docs/si/fig
    ../lunar_regolith_DSM/venv_docparse/bin/python docs/si/si_throat_vapor_bc.py

Numbers: studies/cryosuction/scripts/plan_estimates.py --wt 4e-5 8e-5 (section 8c) and its printed
constants at -20 C; rho_vs: src/material_properties.c RhoVS_I (Wexler 1977 coefficients, DOI checked
on Crossref); D_v, beta0, d0: preprocess/mesh_params/physics.py.
"""
from pathlib import Path

from si_common import EQ, TBL, H, P, save, title  # noqa: I001  (sets the path for docx_helpers)
from docx_helpers import figure

title("Supplementary Methods: vapor boundary conditions for the throat simulations",
      "Draft 3.1, 13 November 2026.")

H("Purpose", 1)
P("The simulations solve a phase-field model of ice sublimation and deposition in 2D pores at −20 °C "
  "in air at 1 atm, a laboratory analog of ice in regolith pores. In the throat simulations an ice lens spans "
  "a pore that narrows to a throat (Fig. S1). The lens moves only "
  "by sublimation at one face and deposition at the other; the vapor travels through the gas between each face "
  "and the end of the pore on that side. The vapor density held at the two ends therefore decides whether the "
  "lens also grows or shrinks while it moves. This section derives the end condition that returns the lens to "
  "its initial area at the end of a run, and explains why the runs without a temperature gradient use a "
  "periodic pore instead. The effects are small and of similar size: curvature shifts a face's equilibrium "
  "humidity by $d_0\\chi \\approx 2.5\\times10^{-5}$ (40 µm gap), a gradient of 1 K/m shifts the ends by "
  "$sGL \\approx 4\\times10^{-5}$, and an error of $10^{-7}$ in the end humidity changes the lens area by "
  "about 1 % over 300 days. This is why the end humidity is given to seven digits.")
figure(str(Path(__file__).parent / "fig/fig_throat.png"),
       "Figure S1. Pore with a throat and an ice lens (2D, per unit depth). Each face exchanges vapor with the "
       "end of the pore on its side through the gas column between them. The ends hold the vapor density at a "
       "fraction $h$ of saturation at the local temperature. The walls are impermeable and insulating. Schematic: "
       "the lens of Table S1 starts 70 µm to the right of the throat.",
       width=6.3)

H("Setting and symbols", 1)
P("The pore runs from $x = 0$ to $x = L$ = 400 µm. Its gap is $g(x) = \\sqrt{w_t^2 + c^2(x - L/2)^2}$, with the "
  "throat width $w_t$ at $x = L/2$ and $c = \\sqrt{W^2 - w_t^2}/(L/2)$ so that $g = W$ = 150 µm at both ends. "
  "The walls carry no flux of heat or vapor. A temperature gradient $G$ is imposed by fixing the temperature "
  "at the two ends, warm at $x = 0$. The gas is air at 1 atm (a laboratory analog). The lens has a warm face "
  "$w$ (toward $x = 0$) and a cold face $c$; the labels only name the two faces when $G = 0$.", after=3)
TBL(["Symbol", "Meaning", "Value at −20 °C"], [
    ["$\\theta$", "contact angle, measured through the ice", "60°"],
    ["$\\omega_i$", "angle between the walls at face $i$ and the pore axis, positive when the gap widens into "
     "the face's gas: $\\tan\\omega_w = -\\tfrac12 g'(x_w)$, $\\tan\\omega_c = +\\tfrac12 g'(x_c)$", ""],
    ["$\\chi_i$", "curvature of face $i$; positive when the face is concave seen from the gas", ""],
    ["$\\rho_{vs}(T)$", "saturation vapor density over ice in the model: $0.62\\,p_{vs}/(p_{\\mathrm{atm}} - "
     "p_{vs})$ times a fixed air density 1.341 $\\mathrm{kg/m^3}$, with $p_{vs}$ from Wexler (1977)",
     "$8.49\\times10^{-4}\\ \\mathrm{kg/m^3}$"],
    ["$s$", "$d\\ln\\rho_{vs}/dT$ of that expression", "0.0961 1/K"],
    ["$h$", "vapor density held at the ends, as a fraction of $\\rho_{vs}$ at the end's temperature", ""],
    ["$h_i^{\\mathrm{eq}}$", "vapor density in equilibrium with face $i$, as a fraction of $\\rho_{vs}$", ""],
    ["$\\Delta_i$", "supersaturation that drives face $i$: (vapor density at the end − equilibrium density of the "
     "face) / $\\rho_{vs}(T_{\\mathrm{lens}})$", ""],
    ["$d_0$", "capillary length $\\gamma v_m/(k_B T)$; $\\gamma$ = 0.109 $\\mathrm{J/m^2}$, $v_m$ = volume of one water "
     "molecule in ice, $k_B$ Boltzmann's constant", "$1.015\\times10^{-9}$ m"],
    ["$\\beta_0$", "attachment coefficient, $\\beta_0 = \\alpha^{-1}\\sqrt{2\\pi m/(k_B T)}\\,"
     "\\rho_{\\mathrm{ice}}/\\rho_{vs}$ (Hertz–Knudsen) with deposition coefficient $\\alpha$ = 0.01 and "
     "$m$ the mass of a water molecule; a face moves at $\\Delta/\\beta_0$ when attachment alone limits it",
     "$7.94\\times10^{5}$ s/m"],
    ["$D_v$", "vapor diffusivity in air at 1 atm, $2.178\\times10^{-5}(T/273.15)^{1.81}$ $\\mathrm{m^2/s}$",
     "$1.90\\times10^{-5}\\ \\mathrm{m^2/s}$"],
    ["$K$", "$\\rho_{\\mathrm{ice}}/(D_v\\,\\rho_{vs})$, with $\\rho_{\\mathrm{ice}}$ = 919 $\\mathrm{kg/m^3}$",
     "$5.70\\times10^{10}\\ \\mathrm{s/m^2}$"],
    ["$A$", "lens area per unit depth", ""],
], [0.85, 4.15, 1.5], size=8.5)

H("Step 1: equilibrium of each face", 1)
P("A face spans the local gap $g_i$ and meets each wall at the contact angle; with the walls inclined by "
  "$\\omega_i$ it is a circular arc of radius $g_i/[2\\cos(\\theta + \\omega_i)]$, so its curvature is")
EQ(r"\chi_i = \frac{2\cos(\theta + \omega_i)}{g_i}")
P("By the Gibbs–Thomson relation, linearized because $d_0\\chi_i \\sim 10^{-5}$,")
EQ(r"h_i^{\mathrm{eq}} = 1 - d_0\,\chi_i")
P("A concave face is in equilibrium below saturation, and a narrower gap gives a lower $h_i^{\\mathrm{eq}}$.")

H("Step 2: rate at which each face gains ice", 1)
P("Vapor reaches face $i$ in two steps in series. Attachment: the face advances at "
  "$v_i = \\Delta_i^{\\mathrm{face}}/\\beta_0$, so the area gained per time is $Q_i = g_i v_i$ and the "
  "supersaturation used up at the face is $\\beta_0 Q_i/g_i$. Diffusion: in the quasi-steady gas column the "
  "same mass flux $\\rho_{\\mathrm{ice}}Q_i$ crosses every section of width $g(x)$, so Fick's law gives a drop "
  "of $K Q_i\\,dx/g(x)$ in the saturation fraction over a length $dx$. Adding the two drops between the end "
  "and the face:")
EQ(r"Q_i = \frac{\Delta_i}{R_i},\qquad R_i = \frac{\beta_0}{g_i} + K\int_{\mathrm{column}\ i}\frac{dx}{g(x)}")
P("where the integral runs from the face (its apex on the axis) to the end. For a uniform gap $g$ and a column "
  "of length $l$ the face speed is $\\Delta_i/(\\beta_0 + K l)$; simulations of a lens in a straight pore "
  "(this study) followed that law within 1 %.")

H("Step 3: what each face sees", 1)
P("Because the walls are insulating and ice conducts heat 115 times better than air (2.29 against "
  "0.02 $\\mathrm{W\\,m^{-1}\\,K^{-1}}$), the lens is nearly isothermal and the imposed temperature difference $GL$ falls across the "
  "two gas columns. Heat crosses the same varying gap as vapor, so it splits as "
  "$\\Delta T_i = G L\\,I_i/(I_w + I_c)$ with $I_i = \\int_{\\mathrm{column}\\ i}dx/g$; both "
  "$\\Delta T_i \\ge 0$ (warm end above the lens, lens above the cold end). The ends hold "
  "$\\rho_v = h\\,\\rho_{vs}(T_{\\mathrm{end}}) \\approx h\\,\\rho_{vs}(T_{\\mathrm{lens}})(1 \\pm s\\,\\Delta T_i)$; "
  "dropping the product $(h - 1)s\\Delta T_i$, which is smaller by $10^{-5}$,")
EQ(r"\Delta_w = h - h_w^{\mathrm{eq}} + s\,\Delta T_w,\qquad \Delta_c = h - h_c^{\mathrm{eq}} - s\,\Delta T_c")

H("Step 4: growth of the lens", 1)
P("The lens area changes at the sum of the two face rates:")
EQ(r"\frac{dA}{dt} = Q_w + Q_c = a\,\bigl(h - h_0\bigr),\qquad a = \frac{1}{R_w} + \frac{1}{R_c},\qquad "
   r"h_0 = \frac{h_w^{\mathrm{eq}}/R_w + h_c^{\mathrm{eq}}/R_c + s\,(\Delta T_c/R_c - \Delta T_w/R_w)}{a}")
P("$h_0$ is the end humidity at which the lens neither grows nor shrinks at that moment: the two face equilibria "
  "averaged with weights $1/R_i$, plus the gradient's effect. It changes as the lens moves.")

H("Step 5: choice of a constant end humidity", 1)
P("We take the constant $h$ for which the lens ends the run (time $t_f$) with its initial area. Integrating "
  "Step 4 and setting the result to zero:")
EQ(r"h_{\mathrm{net}} = \frac{\int_0^{t_f} a\,h_0\,dt}{\int_0^{t_f} a\,dt}")
P("The integrals are taken along the path of the lens, computed with the same model: the lens center moves at "
  "the mean of the two face speeds $Q_w/g_w$ (toward $x = 0$) and $-Q_c/g_c$, the area changes at $dA/dt$, and "
  "the lens length is the one that holds that area between the walls. The path depends on $h$, so $h$ is "
  "recomputed from the path until it changes by less than $10^{-12}$ (7 to 11 passes). On the converged path, "
  "the $h$ that minimizes $\\int (dA/dt)^2dt$ is $h_{LS} = \\int a^2 h_0\\,dt/\\int a^2\\,dt$; it differs from "
  "$h_{\\mathrm{net}}$ by less than $3\\times10^{-8}$.")
P("**Table S1.** Pore and lens as above; lens 100 µm long, started with its center 70 µm on the cold side of the "
  "throat; 300 days. Area changes at the end of the run (largest change during the run in parentheses). The "
  "end change at $h_{\\mathrm{net}}$ is zero by construction. $1 - d_0\\langle\\chi\\rangle$ uses the curvature "
  "of the two faces averaged over faces and time without weights.", size=8.5, after=2)
TBL(["$w_t$, $G$", "$h_{\\mathrm{net}}$", "area at $h_{\\mathrm{net}}$", "at $1 - d_0\\langle\\chi\\rangle$",
     "at $h = 1$", "lens center at the end"], [
    ["40 µm, 0 K/m", "0.999 990 9", "0 (0.8 %)", "−28.8 %", "+88.7 %", "34 µm past the throat, warm side"],
    ["40 µm, 1 K/m", "0.999 991 5", "0 (3.0 %)", "−11.1 %", "+77.4 %", "28 µm past, warm side"],
    ["40 µm, 2 K/m", "0.999 991 7", "0 (3.2 %)", "−26.8 %", "+78.3 %", "81 µm past, warm side"],
    ["80 µm, 0 K/m", "0.999 992 3", "0 (0.5 %)", "−12.8 %", "+70.2 %", "58 µm past, warm side"],
    ["80 µm, 1 K/m", "0.999 991 2", "0 (0.7 %)", "−2.4 %", "+68.3 %", "3 µm past, warm side"],
    ["80 µm, 2 K/m", "0.999 991 6", "0 (3.1 %)", "−4.6 %", "+63.0 %", "63 µm past, warm side"],
], [1.05, 1.0, 1.05, 1.05, 0.8, 1.55], size=8.5)
P("Weighting matters: the unweighted mean curvature leaves the lens 2–29 % smaller. A change of $10^{-7}$ in "
  "$h$ changes the area by about 1 % over 300 days, and the approximations listed below can shift $h_0$ by "
  "a few times $10^{-7}$, so $h_{\\mathrm{net}}$ is a starting value and the lens area is monitored in each run. "
  "The final lens positions are not monotonic in $G$ (80 µm throat: 58, 3 and 63 µm for $G$ = 0, 1, 2 K/m); we "
  "have not traced why, and only the $G = 0$ rows are used below.", after=4)

H("Runs without a temperature gradient: periodic pore", 1)
P("Table S1 shows that even at $h_{\\mathrm{net}}$ the ends move the lens: with no gradient it does not stop at "
  "the throat but ends 34 µm (40 µm throat) to 58 µm (80 µm throat) past it on the warm side, because the face "
  "nearer an end has the smaller $R_i$ and its offset $h - h_i^{\\mathrm{eq}}$ drives a flux comparable to the "
  "difference between the two face equilibria. The runs without a gradient therefore use a pore that is "
  "periodic in $x$: one cell of an infinite chain of identical throats, 400 µm long, throat at the center, "
  "$g = \\sqrt{w_t^2 + (c_pL/\\pi)^2\\sin^2(\\pi\\xi/L)}$ with $\\xi = x - L/2$ and "
  "$c_p = \\pi\\sqrt{W^2 - w_t^2}/L$, so that $g = W$ at the cell edges. It has the same throat and widest "
  "gaps as the pore above, but near the throat its walls are $\\pi/2$ times steeper, and the wall slope is zero "
  "at the cell edges, where the cells join. The two faces are "
  "connected by one gas column through the neighboring cell, so")
EQ(r"Q_1 = -Q_2 = \frac{h_2^{\mathrm{eq}} - h_1^{\mathrm{eq}}}{\beta_0/g_1 + \beta_0/g_2 + "
   r"K\int_{\mathrm{gas}} dx/g(x)}")
P("for faces 1 and 2. The ice area is conserved (apart from the small vapor content of the gas). Vapor moves "
  "from the face with the higher equilibrium humidity (smaller curvature) to the other, so the lens moves toward "
  "the narrower gap and stops where the faces have equal curvature: centered on the throat (stable) or centered "
  "on the widest point (unstable). No end humidity is needed.")

H("Runs with a temperature gradient", 1)
P("A linear temperature profile cannot be periodic, so these runs use the pore above with the ends at "
  "$\\rho_v = h_{\\mathrm{net}}\\,\\rho_{vs}(T_{\\mathrm{end}})$, $h_{\\mathrm{net}}$ computed for each run. A "
  "straight pore of width 150 µm with the same lens, gradient and its own $h_{\\mathrm{net}}$ shows how far "
  "the ends and the gradient alone move the lens; the difference between the two pores is attributed to the "
  "throat, to first order.")

H("Approximations", 1)
P("Isothermal ice; quasi-steady, one-dimensional gas columns; the face length taken as the gap $g_i$ (the arc is "
  "about 5 % longer at 60°); the static contact angle at moving contact lines (the moving angle overshoots it); "
  "Gibbs–Thomson linearized in $d_0\\chi$. These estimates set the boundary conditions only; the simulations do "
  "not use them otherwise.")

H("Reference", 2)
P("Wexler, A. (1977). Vapor pressure formulation for ice. J. Res. Natl. Bur. Stand. 81A, 5. "
  "doi:10.6028/jres.081A.003", size=8.5)

save("SI_throat_vapor_BC.docx")

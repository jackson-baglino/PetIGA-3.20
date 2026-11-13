"""Build the concise run plan (2026-11-12): physics, then one section per
simulation type with its schematic, setup, varied parameters, run table and
analysis.

    ../lunar_regolith_DSM/venv_docparse/bin/python docs/run_plan/build_runplan2.py <fig_dir> <out.docx>

Figures: docs/run_plan/sketches.py <fig_dir> (no scratch geometries).
Numbers: studies/cryosuction/scripts/plan_estimates.py (sections 2, 3, 5, 6b,
7), the tests files in studies/cryosuction/ and RUNS.md. The gradients are
set once below (G_PROD, G_HIGH); every table and sentence uses them.
Check the result with audit.py.
"""
import math
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
FIGDIR, OUTPATH = sys.argv[1], sys.argv[2]
sys.argv = [sys.argv[0], FIGDIR, "", OUTPATH]          # docx_helpers reads nothing from argv
os.environ["RUNPLAN_DATE"] = "2026-11-12"              # footer date (docx_helpers)
from docx_helpers import *  # noqa: E402,F401,F403
from docx_helpers import GRAYT, H, Inches, P, Pt, T, WD_ALIGN_PARAGRAPH, doc, rich  # noqa: E402

# ------------------------------------------------------------------ gradients (one place)
G_PROD = (10, 20)          # production gradients [K/m]; chosen by the user 2026-11-12
G_HIGH = 50                # a few runs above the face-stability estimate of a 150 um pore (user)
G_ONE = G_PROD[1]          # gradient of the single-G series (position, size, vapor set)
S_LN = 9.60556e-2          # d ln(rho_vs)/dT at -20 C [1/K] (plan_estimates.py, inputs)

GP = ", ".join(str(g) for g in G_PROD)
GA = f"{GP}, {G_HIGH}"


def dz_mm(h, G):
    """Harvest cutoff distance (plan_estimates.py section 5)."""
    return math.log(1 / h) / (S_LN * G) * 1e3


for sec in doc.sections:
    sec.top_margin = sec.bottom_margin = Inches(0.7)


def fig(name, caption, width=6.5):
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.keep_with_next = True
    p.paragraph_format.space_after = Pt(0)
    p.add_run().add_picture(f"{FIGDIR}/{name}.png", width=Inches(width))
    c = doc.add_paragraph()
    c.paragraph_format.space_after = Pt(3)
    rich(c, caption, size=8, italic=True, color=GRAYT)


def para(label, text):
    p = P(f"**{label}** {text}", size=9.5)
    p.paragraph_format.space_after = Pt(3)
    return p


def runs(rows, var="Varied"):
    T(["ID", var, "G (K/m)", "Vapor", "Days", "State"], rows, [0.7, 2.45, 0.85, 0.9, 0.55, 1.15], size=8)


# ================================================================== title
doc.add_paragraph("Run Plan: Ice in Regolith Pores", style="Title")
P("Paper 2 and AGU 2026. 12 November 2026; replaces the run plan of 11 November. Run records: "
  "studies/cryosuction/RUNS.md.", size=9, color=GRAYT, after=6)

# ================================================================== 1
H("1. Physics and problems", 1)
P("All runs solve the two-phase phase-field model of ice sublimation and deposition (Kaempfer and Plapp "
  "2009) in 2D for the ice phase φ, the temperature T and the vapor density ρv (air is 1 − φ), at "
  "T0 = −20 °C. The pore walls are regolith: they carry a contact-angle condition set by the wall "
  "energies (θ = 60° unless stated) and no flux of heat or vapor. The pore ends hold fixed T and either "
  "fixed ρv (vapor reservoirs) or zero vapor flux; a temperature gradient G is imposed between the ends "
  "(warm end at x = 0) or the pore is isothermal. The runs measure how ice redistributes in a pore: the "
  "drift and growth of single bodies, the exchange of mass between bodies, the effect of pore shape "
  "(taper, throat, neck) on the motion of a lens, and the rate at which ice leaves a pore through a "
  "collector.")
P(f"**Gradients.** Production gradients are G = {GP} K/m (chosen 12 November), with a few runs at "
  f"G = {G_HIGH} K/m, above the face-stability estimate for a 150 µm pore (18.5 K/m in the gas, "
  "plan_estimates.py section 6b), for selected geometries and contact angles. The estimated migration speed is 1.5, 2.9 and 7.3 µm/d at G = 10, "
  "20 and 50 K/m (plan_estimates.py section 3).", size=9.5)
T(["Term", "Definition"], [
    ["closed", "zero vapor flux on both ends"],
    ["res", "ρv fixed on both ends at $\\rho_{vs}(T_0)(1 + d_0\\chi)$ for a reference curvature χ "
     "(h = 1: χ = 0)"],
    ["ls (local saturation)", "ρv fixed on each end at ρvs of that end's temperature"],
    ["dh", "isothermal; ρv fixed at the ls values of a gradient run (separates vapor supply from the T field)"],
    ["collector", "$x = 0$ at $h_{\\mathrm{open}}\\,\\rho_{vs}(T_0)$, $h_{\\mathrm{open}}$ = 0.99 unless "
     "stated; $x = L_x$ a sealed dead end"],
    ["numerics", "ε = 0.858 µm (0.800 µm for types 2 and 5, whose smallest feature is 20 µm); elements "
     "0.6 µm; $\\mathrm{dt}_{\\max} = 2\\tau_{\\mathrm{sub}}$; 50 000 DOF per core"],
    ["snapshots", "401 per run (-outp 0 -n_out 401) for movies; about 11 MB each on a 1.36 M DOF mesh "
     "(4.4 GB per type 6 run)"],
    ["states", "done: ran and analyzed; ran: analysis open; queued: inputs pushed, submission pending; "
     "ready: inputs built and preflight passed; to build: inputs not written"],
], [1.35, 5.25], size=8)

# ================================================================== type 1
H("2. Type 1: one wall droplet", 1)
fig("type1", "Figure 1. Type 1: 600 by 300 µm box, regolith wall at y = 0, one ice cap.")
para("Simulations.", "One ice cap (circular arc, 40° start) on a regolith wall in an imposed gradient. With "
     "closed vapor ends the ice mass is fixed and the cap can only drift; with fixed vapor at the ends it "
     "can also grow or shrink.")
para("Varied.", "G; vapor set (closed, res, ls, dh); cap position x; cap radius R.")
runs([
    ["1A-0", "control", "0", "closed", "10", "ready"],
    ["1A-1…3", "G", GA, "closed", "60, 60, 30", "to build"],
    ["1B-1…3", "x = 500, 300, 100 µm", str(G_ONE), "res", "60", "to build"],
    ["1B-4…6", "x = 500, 300, 100 µm", str(G_ONE), "ls", "60", "to build"],
    ["1B-7", "isothermal, open ends", "0", "res", "10", "ready"],
    ["1B-8", "R 50 µm", str(G_HIGH), "ls", "30", "to build"],
    ["1C-1, 2", "R = 25, 100 µm", str(G_ONE), "closed", "60", "to build"],
    ["1C-3, 4", "R = 25, 100 µm", str(G_ONE), "ls", "60", "to build"],
    ["1D-1, 2", "vapor twins of 1B-5 and 1B-8", "0", "dh", "60, 30", "to build"],
])
para("Analysis.", "ice_bodies.py: cap area A(t), center of mass $x_c(t)$ and drift speed, position and "
     "speed of each contact line. water_budget.py: ice, stored vapor and outflow at the ends (res, ls, dh). "
     "contact_angle.py: angle of each face against time. Reported against G, R, x and vapor set.")

# ================================================================== type 2
H("3. Type 2: several droplets", 1)
fig("type2", "Figure 2. Type 2: (a) the pair, 600 by 300 µm; (b) the five-cap chain, 900 by 300 µm. "
    "Other layouts: pair with the large cap on the cold side; three equal caps; three caps with a small "
    "center cap; two 40 µm caps 20, 80 and 200 µm apart.")
para("Simulations.", "Two or more caps on one wall exchange vapor. In a closed isothermal box the small "
     "cap of a pair loses mass to the large one; a gradient adds a transfer from the warm cap to the cold "
     "cap, which cancels the ripening drive of the pair at G* = 1.6 K/m (plan_estimates.py section 2).")
para("Varied.", "G (ladder around G*); which cap is warm; number and sizes of caps; cap spacing.")
runs([
    ["2A-1", "pair", "0", "closed", "60", "ready"],
    ["2A-2", "pair, open ends", "0", "res", "60", "ready"],
    ["2B-1…5", "pair, large cap warm; G", "0.8, 1.6, 3.2, 8, 32", "closed", "60", "ready"],
    ["2B-6", "pair, large cap cold", "3.2", "closed", "60", "ready"],
    ["2C-1…3", "chain: 3 equal, 3 small center, 5 random", str(G_ONE), "closed", "60", "to build"],
    ["2D-1…3", "gap 20, 80, 200 µm", str(G_ONE), "closed", "60", "to build"],
])
para("Analysis.", "ice_bodies.py: area of each cap $A_i(t)$; vanishing time of the smaller cap; the G at "
     "which $dA_{\\mathrm{small}}/dt$ changes sign (2B). section_flux.py: vapor flux through vertical lines "
     "between caps (2C, 2D). water_budget.py: mass conservation in the closed box.")

# ================================================================== type 3
H("4. Type 3: lens in a flat pore (done)", 1)
fig("type3", "Figure 3. Type 3: 500 µm channel, H = 125 µm (also 60 µm), 150 µm lens with flat faces at "
    "the start.")
para("Simulations.", "A lens spanning a straight channel, isothermal, with vapor reservoirs at both ends. "
     "The meniscus curvature $2\\cos\\theta/H$ sets the vapor density at each face.")
para("Varied.", "θ; H; vapor set.")
runs([
    ["3A-1…5", "θ = 30, 60, 90, 120, 150°; H 125 µm", "0", "res", "120", "done"],
    ["3B-1, 2", "θ = 60, 120°; H 60 µm", "0", "res", "120", "done"],
    ["3C-1", "closed control", "0", "closed", "60", "done"],
    ["3D-1, 2", "reservoir at x = 0 only; θ = 60, 120°", "0", "res / closed", "120", "done"],
])
para("Analysis.", "ice_bodies.py: face positions and speeds; analyze_type3.py and compare_velocity.py: "
     "face speed against $v = 2 d_0\\cos\\theta/H$ in series with diffusion to the reservoirs, with no fitted "
     "parameter; contact_angle.py: face angle. Measured: speed within 1 % of the series law "
     "(REPORT_type3).")

# ================================================================== type 4
H("5. Type 4: wedge (done)", 1)
fig("type4", "Figure 4. Type 4: wedge with wall tilt ω = 10° (also 2.5, 5, 20°), lens where the gap is "
    "100 µm.")
para("Simulations.", "A lens in a tapered channel, isothermal, with reservoirs at both ends. The two faces "
     "meet the walls at θ − ω and θ + ω, so their curvatures differ.")
para("Varied.", "ω; θ.")
runs([
    ["4A-0", "flat channel", "0", "res", "180", "done"],
    ["4A-1…4", "ω = 2.5, 5, 10, 20°", "0", "res", "180", "done"],
    ["4B-1, 2", "ω = 10°; θ = 90, 85°", "0", "res", "180", "done"],
])
para("Analysis.", "ice_bodies.py: lens center position and face speeds; analyze_type4.py: speeds against "
     "the series law with $\\theta \\mp \\omega$. Measured: agreement within 3–5 % (REPORT_type4).")

# ================================================================== type 5
H("6. Type 5: lens at a throat", 1)
fig("type5", "Figure 5. Type 5: 800 µm pore, 150 µm wide, throat at x = 400 µm, lens of 120 µm.")
para("Simulations.", "A lens spans a throat in a pore with local-saturation ends. With no gradient the "
     "lens rests where its two faces have equal curvature; a gradient displaces it toward the warm end. "
     "5B adds two caps and a grain near the throat in a closed pore; 5C starts the lens 100 µm from the "
     "throat.")
para("Varied.", "G; throat width $w_t$; start position; other ice bodies (5B).")
runs([
    ["5A-1…4", "$w_t$ = 40 µm; G", f"0, {GA}", "res (G 0), ls", "120", "ready (G 0); to build"],
    ["5A-5…7", "$w_t$ = 80 µm; G", f"0, {GP}", "res (G 0), ls", "120", "ready (G 0); to build"],
    ["5B-1, 2", "lens, 2 caps and a grain", f"0, {G_PROD[0]}", "closed", "120", "ready (G 0); to build"],
    ["5C-1, 2", "lens 100 µm toward $x = 0$, toward $x = L_x$", str(G_PROD[0]), "ls", "120", "to build"],
])
para("Analysis.", "ice_bodies.py: lens center $x_c(t)$, its rest position against G and $w_t$, face "
     "positions; contact_angle.py: face angles at rest; water_budget.py and ice_bodies.py: mass moved from "
     "the caps and grain to the lens (5B).")

# ================================================================== type 6
H("7. Type 6: lens through a constriction", 1)
fig("type6", "Figure 6. Type 6: 1100 µm pore, 150 µm wide, neck of width $w_c$ on x = 400–700 µm; the "
    "lens starts in the cold section (θ = 60°: concave faces; 6B: faces drawn as arcs of radius 7.5 mm, "
    "0.38 µm from flat; 6D: flat faces from a φ file).")
para("Simulations.", "A lens in the cold section moves toward the warm end under G and passes through a "
     "neck. $w_c$ = 150 µm is the straight-pore reference.")
para("Varied.", "$w_c$; G; θ.")
runs([
    ["6A-1…3", "$w_c$ = 40, 80, 150 µm; θ 60°", "100", "ls", "60", "ran"],
    ["6A-4…6", "$w_c$ = 40, 80, 150 µm; θ 60°", "300", "ls", "30", "ran; 6A-4 stopped at 9.6 d"],
    ["6B-1", "$w_c$ = 40 µm; θ 90°", "300", "ls", "30", "queued (resume from 11.5 d)"],
    ["6B-2", "$w_c$ = 80 µm; θ 90°", "300", "ls", "30", "queued"],
    ["6C-1…3", "$w_c$ = 40, 80, 150 µm; θ 60°", str(G_PROD[0]), "ls", "360", "to build"],
    ["6C-4…6", "$w_c$ = 40, 80, 150 µm; θ 60°", str(G_PROD[1]), "ls", "180", "to build"],
    ["6C-7…9", "$w_c$ = 40, 80, 150 µm; θ 60°", str(G_HIGH), "ls", "75", "to build"],
    ["6D-1…3", "straight, W 150 µm; θ 90°; flat + 0.5 µm seed", GA, "ls", "120, 120, 60", "to build"],
    ["6D-4", "straight, W 75 µm; θ 90°; flat + 0.5 µm seed", str(G_HIGH), "ls", "60", "to build"],
    ["6D-5", "straight, W 150 µm; θ 90°; flat, no seed", str(G_HIGH), "ls", "60", "to build"],
    ["6D-6", "straight, W 150 µm; θ 75°", str(G_HIGH), "ls", "60", "to build"],
], var="Varied")
para("6D, face stability.", "A depositing face in a pore of width W with no-flux walls is predicted "
     "unstable to the symmetric mode (axis and walls moving apart, wavenumber $k = 2\\pi/W$) when the "
     "gas gradient exceeds $G_s = (2\\pi/W)^2 d_0/s$, with growth rate "
     "$\\sigma = (D_v \\rho_{vs}/\\rho_i)\\,k\\,(s G_{\\mathrm{gas}} - d_0 k^2)$ (linear, flat base "
     "face; plan_estimates.py section 6b). $G_s$ is 18.5 K/m in the gas for W = 150 µm (16 K/m nominal in "
     "the straight pore) and 74 K/m for W = 75 µm, so 6D-1 is below, 6D-2 and 6D-3 above, and 6D-4 below "
     "it. The flat face with a cosine seed of 0.5 µm needs a flat-face feature in build_ice_layout.py "
     "(-ic_type phi_file; no solver change). 6D-5 shows whether the face deforms without a seed; 6D-6 "
     "and 6C-9 (θ 60°) give the effect of a curved starting face at the same G.")
para("Analysis.", "6D: amplitude A(t) = x(axis) − x(wall) of each face; σ from the slope of ln|A| "
     "against t, compared with the formula above; the G and W at which σ changes sign. 6A–6C: "
     "ice_bodies.py: x of the warm and cold faces on the axis and at each wall against t; "
     "times at which the lead face enters and leaves the neck; transit time against $w_c$. section_flux.py: "
     "vapor flux through a line across the neck against the face speed. water_budget.py: ice gained from "
     "the reservoirs. resolution.dat: Π. Estimated cost (plan_estimates.py section 7d): 22–108 h per 6C "
     "run; runs longer than 48 h continue with scripts/HPC/resume_run.sh.")

# ================================================================== type 7
H("8. Type 7: rough regolith pore with mixed ice", 1)
fig("type7", "Figure 7. Type 7, layout A: 2 mm pore with rough walls (gap 50–147 µm, 50 µm throat at "
    "x = 1 mm), three lenses at 60°, four caps and three grains at 40°; 24.6 % ice; collector at x = 0.")
para("Simulations.", "Ice leaves the pore through the collector at x = 0. 7A is isothermal; 7B insulates the "
     "dead end; 7C puts the warm end at the collector, 7D at the dead end.")
para("Varied.", "Thermal condition; G.")
runs([
    ["7A", "isothermal", "0", "collector", "60", "ready (pre-check)"],
    ["7B", "dead end insulated", "0", "collector", "60", "to build"],
    ["7C-1…3", "warm end at the collector", GA, "collector", "60", "to build"],
    ["7D-1…3", "warm end at the dead end", GA, "collector", "60", "to build"],
], var="Varied")
para("Analysis.", "ice_bodies.py: area of each body against t and the order in which bodies empty. "
     "water_budget.py with bflux.dat: harvest rate (vapor outflow at x = 0) and cumulative harvest. "
     "Recapture fraction: ice deposited on bodies divided by ice sublimated (ice_bodies.py areas). "
     "section_flux.py: flux through the 50 µm throat. Each run is 2.41 M DOF (48 cores).")

# ================================================================== S
H("9. Additional runs: S5 and S2", 1)
fig("S", "Figure 8. S5: 850 µm pore with throats of 30 and 90 µm, both plugged with ice; collector at x = 0, "
    "dead end at $x = L_x$.")
cut = [(G, h, dz_mm(h, G)) for G in (G_PROD[1], G_HIGH) for h in (0.999, 0.9995)]
para("S5, harvest cutoff.", "Collector at x = 0 on the warm end. Ice colder than the temperature at which "
     "ρvs equals the collector density does not lose mass to the collector; the cutoff lies "
     "$\\Delta z = \\ln(1/h_{\\mathrm{open}})/(s G)$ beyond the T0 point (s = dln ρvs/dT). Values: "
     + "; ".join(f"{d:.2f} mm (G = {G} K/m, $h_{{\\mathrm{{open}}}}$ = {h:g})" for G, h, d in cut) + ".")
para("S2, throat as an equivalent length.", "The S5 pore with one throat and the same pore without it, "
     "isothermal, collector at $h_{\\mathrm{open}}$ = 0.99. Measured: harvest time; the throat is expressed "
     "as the straight length that gives the same harvest time.")
T(["ID", "Varied", "G (K/m)", "Vapor", "Days", "State"], [
    ["S5-1…4", "$h_{\\mathrm{open}}$ = 0.999, 0.9995", f"{G_PROD[1]}, {G_HIGH}", "collector", "60", "to build"],
    ["S2-1, 2", "with and without the throat", "0", "collector", "60", "to build"],
], [0.7, 2.45, 0.85, 0.9, 0.55, 1.15], size=8)
para("Analysis.", "ice_bodies.py: area of each plug against t and the position beyond which ice stops "
     "losing mass (S5); water_budget.py: harvest against t (S2).")

H("References", 2)
P("Kaempfer, T. U. and Plapp, M. (2009). Phase-field modeling of dry snow metamorphism. Phys. Rev. E 79, "
  "031502. doi:10.1103/PhysRevE.79.031502", size=8.5)

doc.save(OUTPATH)
print("wrote", OUTPATH)

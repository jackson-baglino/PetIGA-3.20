"""Build the run plan (every simulation for paper 2 / AGU 2026, its setup and
its state) as a Word document.

    ../lunar_regolith_DSM/venv_docparse/bin/python docs/run_plan/build_runplan.py <fig_dir> <out.docx>

Figures come from docs/run_plan/sketches.py and t6_timelapse.py. Numbers in
the text come from studies/cryosuction/scripts/plan_estimates.py (sections 3,
5, 6), postprocess/run_health.py and postprocess/ice_bodies.py on the T6
batch, and the input files named in each table. Check the result with
audit.py (text outside equations).
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
FIGDIR, OUTPATH = sys.argv[1], sys.argv[2]
sys.argv = [sys.argv[0], FIGDIR, "", OUTPATH]          # docx_helpers reads nothing from argv
from docx_helpers import *  # noqa: E402,F401,F403
from docx_helpers import B, GRAYT, H, Inches, NAVY, P, T, WD_ALIGN_PARAGRAPH, box, doc, rich  # noqa: E402

GREEN, AMBER, BLUE = "E8F4EA", "FFF8E1", "EAF2FB"


def fig(name, caption, width=6.5):
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.keep_with_next = True
    p.add_run().add_picture(f"{FIGDIR}/{name}.png", width=Inches(width))
    c = doc.add_paragraph()
    rich(c, caption, size=8.5, italic=True, color=GRAYT)


def setup_table(rows):
    T(["Item", "Value"], rows, [1.6, 5.0], size=8.5, bold_first=True)


def runs_table(rows, extra="Variable"):
    T(["ID", extra, "G (K/m)", "Vapor BC", "Days", "State"], rows, [0.75, 2.35, 0.75, 0.95, 0.5, 1.3],
      size=8)


# ================================================================== TITLE
doc.add_paragraph("Run Plan: Ice in Regolith Pores", style="Title")
P("Every simulation for paper 2 and the AGU 2026 presentation: setup, purpose and state", size=12,
  color=NAVY, after=2)
P("11 November 2026; supersedes the run lists of SimulationPlan_v2 (6 November 2026); "
  "states from studies/cryosuction/RUNS.md", size=9, color=GRAYT, after=10)

box("How to read this document", [
    "Section 1 lists every simulation set and its state. Section 2 gives the model, the parameters "
    "shared by all runs and the boundary-condition sets. Section 3 reports what the type 6 runs showed "
    "about the temperature gradient and proposes the gradients for the remaining runs. Sections 4–10 "
    "describe each type: a schematic drawn from the actual mesh and initial ice, the setup, what the set "
    "tests, and the run matrix. Section 11 covers the additional runs, section 12 the items to "
    "investigate, and section 13 the order of work.",
    "**States.** Done: ran and analyzed. Ran: ran, analysis open. Ready: inputs built, preflight passed, "
    "not submitted. Revise: inputs built, but the gradient or duration changes under section 3. "
    "To build: inputs not written. Proposed: new in this document, not yet approved.",
], BLUE, "2E6DB4")

# ================================================================== 1
H("1. Status of every simulation set", 1)
T(["Set", "Question", "Runs", "State", "Next action"], [
    ["0, 0b, 0c", "Does the solver carry a harvesting run (dead-end channel, throat, rough pore)?",
     "about 20", "Done (5 Nov)", "none"],
    ["I1", "Is the 0.86 µm mesh accurate for the fastest fronts?", "1", "Done (7 Nov): within 3 %", "none"],
    ["1", "Drift and growth of one wall droplet vs gradient, vapor path, size", "19",
     "Ready; 11 Revise", "rebuild gradient runs after S7"],
    ["2", "Ripening of droplets vs gradient; shielding vs spacing", "14", "Ready; 6 Revise",
     "rebuild 2C, 2D after S7"],
    ["3", "Lens speed in a flat pore vs the series law", "10", "Done (7 Nov): law to 1 %", "none"],
    ["4", "Does a taper move a lens to the narrow end?", "7", "Done (7 Nov): yes", "optional: 4A-4 margin"],
    ["5", "Resting position of a lens at a throat vs gradient", "12", "Ready", "submit"],
    ["6", "Delay of a lens by a constriction", "7 ran; 7 proposed", "Ran (9–11 Nov) at G = 100, 300 K/m",
     "redesign at G ≤ 30 K/m after S7"],
    ["7", "Harvest order and rate in a rough pore with mixed ice", "6", "7A Ready; 7B–7D To build",
     "pre-check 7A; write 7B–7D"],
    ["S7", "Gradient at which a depositing lens face stays a meniscus", "4", "Proposed", "build and run first"],
    ["S2, S5", "Throat as an equivalent length; harvest cutoff vs gradient", "2 + 4", "To build",
     "after types 5 and 7"],
    ["8", "Lunar physics", "–", "Paper 3", "not in this model"],
], [0.55, 2.55, 0.65, 1.45, 1.4], size=8)

# ================================================================== 2
H("2. Model and shared setup", 1)
P("All runs use the two-phase phase-field model of sublimation and deposition (ice φ, temperature T, vapor "
  "density ρv; air is 1 − φ) with the regolith walls as boundaries that carry a contact-angle condition. "
  "The solver is unfrozen until production (FREEZE.md); every solver change is logged there, and the "
  "binary of each run is tied to its commit by the staged source copy in the run folder.")
setup_table([
    ["Temperature", "T0 = −20 °C for every run (eps and mesh are sized for it; the solver aborts at another "
     "temperature)"],
    ["Kinetics", "α_c = 10⁻², β0 = 7.94 × 10⁵ s/m; capillary length d0 = 1.015 × 10⁻⁹ m"],
    ["Vapor", "$\\rho_{vs}(T_0) = 8.487 \\times 10^{-4}\\ \\mathrm{kg/m^3}$; $D_v(T_0) = 1.898 \\times 10^{-5}\\ \\mathrm{m^2/s}$; "
     "$s = d\\ln\\rho_{vs}/dT = 0.0961\\ \\mathrm{K^{-1}}$"],
    ["Interface", "ε = 0.858 µm (types 1, 3, 4, 6, 7, S), 0.800 µm (types 2, 5); element size "
     "about 0.6 µm along the pore; the visible 5–95 % band is about 5.9ε ≈ 5 µm"],
    ["Time step", "adaptive, interface-CFL capped; dtmax = 2τ_sub = 1341 s (ε 0.858 µm) or 1154 s "
     "(ε 0.800 µm)"],
    ["Wall angle", "θ = 60° on every regolith wall (decision D1, to confirm), set by the wall energies "
     "(γ_is from Young)"],
    ["Initial ice", "lenses start at the wall angle (60°); droplets (caps) and grains start at 40°, a "
     "placeholder for liquid water on regolith (user, 9 Nov)"],
    ["Open faces", "φ: 90° face (zero normal gradient) on every non-wall face (phi_open_bc 0)"],
    ["Resources", "50 000 DOF per core; about 11 s per step at 28 cores and 1.4 M DOF (T6)"],
])
P("**Thermal conditions.** Isothermal: T = T0 on every face. Gradient: T fixed on both pore ends at "
  "T0 ∓ G·L/2, warm end at x = 0, walls insulating. A warm dead end (7D) puts the warm end at x = L_x.")
P("**Vapor sets** (fixed values are written as h = ρv/ρvs(T0)):")
T(["Set", "x = 0", "x = L_x", "Used by"], [
    ["closed", "zero flux", "zero flux", "1A, 1C, 2A-1, 2B–2D, 3C, 5B"],
    ["res", "h = 1 + d0χ (χ of a reference feature)", "same", "1B-1…3, 1B-7, 2A-2, 3A, 3B, 4, 5A (G = 0)"],
    ["ls (local saturation)", "ρvs(T(0))", "ρvs(T(L_x))", "1B-4…6, 1B-8, 1C-3/4, 5A, 5C, 6, S7"],
    ["dh", "h = 1 + Δh/2 (isothermal)", "h = 1 − Δh/2", "1D (vapor-only twin of an ls run)"],
    ["top", "h = 1", "zero flux", "3D"],
    ["collector", "$h = h_{\\mathrm{open}}$ (0.99 or as stated)", "zero flux (dead end)", "7, S2, S5"],
], [1.2, 2.0, 1.4, 2.0], size=8)

# ================================================================== 3
H("3. The temperature gradient", 1)
H("3.1 What the type 6 runs showed", 2)
P("The seven type 6 runs (batch 2026-11-09__15.09.21_type6) put a 150 µm long lens with 60° faces in a "
  "1100 µm dog-bone pore, local-saturation vapor at both ends and G = 100 or 300 K/m, warm end at x = 0. "
  "Five runs finished (27–35 h on 28 cores); the two G = 300 K/m runs with the 40 µm neck were still "
  "running at 46 h of a 48 h limit, at day 9.6 and 11.7 of 30.")
T(["Run", "Done", "Lead edge leaves neck (d)", "Ice at x = 0 (d)", "Peak ice", "Steps with Π > 1"], [
    ["W = 150 (no neck), G 100", "60 d", "14", "26", "+12.5 %", "0.5 %"],
    ["neck 80 µm, G 100", "60 d", "14", "28", "+11.2 %", "1.2 %"],
    ["neck 40 µm, G 100", "60 d", "18", "40", "+10.7 %", "4.1 %"],
    ["W = 150, G 300", "30 d", "5", "9.5", "+15.3 %", "1.6 %"],
    ["neck 80 µm, G 300", "30 d", "5", "10.5", "+13.7 %", "17 %"],
    ["neck 40 µm, G 300", "no (9.6 d)", "6.5", "–", "+12.3 %", "65 %"],
    ["neck 40 µm, θ 90°, G 300", "no (11.7 d)", "7.5", "–", "+13.7 %", "57 %"],
], [1.75, 0.75, 1.15, 0.95, 0.75, 1.15], size=8)
P("Times are read from snapshots 0.5–1 d apart (ice_bodies.py) and carry $\\pm 1$ d. Π = β0·v·ε/d0 is the "
  "resolution number of Kaempfer and Plapp (2009, Eq. 45); above 1 the front moves faster than the mesh "
  "resolves.")
fig("t6_wc80_G100", "Figure 1. Run with the 80 µm neck at G = 100 K/m, warm end at x = 0 (brown: walls from "
    "the mesh; blue: φ > 1/2). The warm face's contact lines run ahead along both walls while its center "
    "lags; the cold face flattens and retreats; at 15 d the lens opens along the axis and leaves one ice "
    "film on each wall.")
P("Three observations bear on the choice of gradient.")
B("**The lens does not move as a plug.** In every θ = 60° run, the face that gains ice (the warm face) "
  "deforms: the distance by which its wall contact lines lead its center on the axis grows from the initial "
  "20 µm (the 60° meniscus) to 75 µm at 3 d and 345 µm at 11 d (W = 150, G = 100 K/m), three times faster "
  "at G = 300 K/m. The face that loses ice (the cold face) flattens: its lead falls from 20 µm to 2–4 µm. "
  "The lens then splits into two wall films (day 4–15; day 37 for the 40 µm neck at G = 100 K/m), in the "
  "straight pore as well as with a neck, so the split is not caused by the neck.")
B("**A flat start stays flat.** The θ = 90° lens (flat faces, 90° wall angle) kept a flat warm face "
  "(lead below 0.4 µm) for 11.7 d at G = 300 K/m and did not split. The deformation grows from the "
  "curvature the 60° meniscus already has; it is amplified, not created, by the gradient.")
B("**Resolution.** At G = 300 K/m the 40 µm neck fills with ice and the air channel on its axis breaks "
  "into bubbles 3–8 µm across (a few ε); Π exceeds 1 in 57–65 % of the steps and every step is cut by "
  "the interface-CFL limit (dt 100–700 s against dtmax 1341 s), which is why these two runs did not "
  "finish. At G = 100 K/m, Π > 1 in 0.5–4 % of the steps.")

H("3.2 Interpretation", 2)
P("The observations agree with the diffusion-limited growth instability of a depositing interface "
  "(Mullins and Sekerka 1963, 1964). In a gradient the vapor arriving at the warm face has a relative "
  "gradient $g = s\\,G$. Parts of the face that stand ahead into the vapor collect more of it; curvature "
  "(the Gibbs–Thomson term d0χ) opposes this. For a flat face and a quasi-static vapor field, a "
  "perturbation of wavenumber k grows when $k^2 < g/d_0$. Across a pore of width W with no-flux walls, "
  "the symmetric mode (both contact lines ahead of the axis, the shape a θ < 90° meniscus already has) "
  "has $k = 2\\pi/W$, so the face is stable for")
P("$G < G_s = \\left(\\frac{2\\pi}{W}\\right)^{2} \\frac{d_0}{s}$", align=WD_ALIGN_PARAGRAPH.CENTER)
P("and the tilting mode ($k = \\pi/W$) lowers the bound by 4. A subliming face is stable, as the cold "
  "faces were. This adaptation of the Mullins–Sekerka result to the vapor field leaves out kinetics "
  "(which slow a mode but do not move its neutral wavenumber), the contact-angle condition and the "
  "finite lens length; it is an order-of-magnitude estimate, and the threshold must be measured (S7).")
T(["Pore width W (µm)", "40", "80", "125", "150", "300"], [
    ["$G_s$, symmetric mode (K/m)", "261", "65", "27", "18.5", "4.6"],
    ["tilting mode (K/m)", "65", "16", "6.7", "4.6", "1.2"],
], [2.1, 0.9, 0.9, 0.9, 0.9, 0.9], size=8.5)
P("A second scale is the gradient at which the gradient drive over a length l equals the Gibbs–Thomson "
  "drive of a feature of curvature χ, $G_x = d_0\\chi/(s\\,l)$: 0.5–0.9 K/m for a 150 µm feature over "
  "its own width, 1.7–3.3 K/m for 80 µm and 6.6–13 K/m for 40 µm (χ = 1/W to 2/W). Below $G_x$, curvature "
  "(ripening, the contact angle) controls the ice; above it, the gradient does. The type 2 ripening pair "
  "has its reversal near G* ≈ 1.6 K/m (SimulationPlan_v2).")
P("At G = 100 and 300 K/m the 150 µm pores are 5 and 16 times above $G_s$, and every capillary scale is far "
  "below the gradient drive. The runs therefore show gradient-driven, unstable deposition and not the "
  "capillary dynamics (meniscus-controlled migration, ripening, droplet drift) the study is meant to "
  "resolve. The summer 2026 runs used G ≤ 50 K/m (user), which is still 2.7 times $G_s$ for a 150 µm pore.")

H("3.3 Proposed gradients", 2)
box("Proposal (needs approval)", [
    "Production gradients **G = 1, 3 and 10 K/m** (with G = 0 controls), all below $G_s$ of the 150 µm pores. "
    "**G = 30 K/m** is kept as one point per type, above $G_s$ for 125–150 µm pores and below it for "
    "features of 80 µm and smaller, to show the transition. G = 100 and 300 K/m are dropped; the T6 runs "
    "at those values remain as the record of the unstable regime.",
    "S7 (section 11) measures the threshold in the straight 150 µm pore before the gradient runs of types "
    "1, 2, 6 and 7 are rebuilt.",
], AMBER, "E0A800", title_color=NAVY)
P("Lower gradients move ice more slowly. The estimate v ≈ D_v·dρvs/dT·G/ρ_ice gives 0.15, 0.44, 1.5 and "
  "4.4 µm/d at G = 1, 3, 10 and 30 K/m; the T6 lenses moved 2.3–2.7 times faster than this estimate "
  "(center-of-mass speed up to 39 µm/d at G = 100 K/m and 101 µm/d at G = 300 K/m in the straight pore). "
  "Run lengths below are set so that ice moves at least one feature size; they are estimates. Π stays "
  "below 0.04 at G ≤ 30 K/m (estimate), so resolution does not limit these runs.")
P("Whether lunar near-surface gradients in ice-bearing regolith fall in this range is not established "
  "here; a cited range is needed for the AGU framing (section 12).")

# ================================================================== 4
H("4. Type 1: one wall droplet", 1)
fig("type1", "Figure 2. Type 1 reference geometry: a 600 by 300 µm box with the regolith wall at y = 0 and "
    "one ice cap (arc radius 50 µm, x = 300 µm). Drawn from t1_box600x300um_capR50_x300um_th40_eps0.86um.",
    width=5.6)
P("**What the set tests.** A single ice cap on a wall in a gradient. With closed vapor the ice mass is "
  "fixed and the cap can only drift (1A, 1C); with fixed vapor at the ends it can also grow or shrink (1B); "
  "the dh twins (1D) impose the same vapor contrast without a temperature gradient, which separates the "
  "effect of the temperature field on the cap from that of the vapor supply. The measured quantities are "
  "the drift speed of the center of mass and of each contact line, the growth rate, and their dependence "
  "on G, cap size and position.")
runs_table([
    ["1A-0", "control, R 50 µm", "0", "closed", "10", "Ready"],
    ["1A-1…4", "R 50 µm; G ladder", "1, 3, 10, 30 (was 10, 30, 100, 300)", "closed", "60", "Revise"],
    ["1B-1…3", "x = 500, 300, 100 µm", "10 (was 100)", "res", "30", "Revise"],
    ["1B-4…6", "x = 500, 300, 100 µm", "10 (was 100)", "ls", "30", "Revise"],
    ["1B-7", "isothermal, open", "0", "res", "10", "Ready"],
    ["1B-8", "R 50 µm", "30", "ls", "10", "Ready"],
    ["1C-1, 2", "R 25, 100 µm", "10 (was 100)", "closed", "60", "Revise"],
    ["1C-3, 4", "R 25, 100 µm", "10 (was 100)", "ls", "30", "Revise"],
    ["1D-1", "vapor twin of 1B-5", "0 (Δh of G 10)", "dh", "30", "Revise"],
    ["1D-2", "vapor twin of 1B-8", "0 (Δh of G 30)", "dh", "10", "Ready"],
])

# ================================================================== 5
H("5. Type 2: several droplets", 1)
fig("type2", "Figure 3. Type 2 geometries: (a) the ripening pair (50 µm cap on the warm side, 20 µm cap), "
    "600 by 300 µm; (b) the five-cap chain (arc radii 20–60 µm), 900 by 300 µm. Other layouts: pair with the "
    "large cap on the cold side, three equal caps, three caps with a small center, two 40 µm caps 20, 80 "
    "and 200 µm apart.", width=5.6)
P("**What the set tests.** Competition between ripening and the gradient. In a closed isothermal box the "
  "small cap of a pair feeds the large one (2A); a gradient opposes or assists this depending on which cap "
  "is warmer, with a reversal expected near G* ≈ 1.6 K/m (2B). Chains (2C) test whether caps shield each "
  "other, and spacing (2D) how far that shielding reaches. The measured quantities are the area of each "
  "cap against time, the time for the smaller cap to vanish, and the gradient at which transfer reverses.")
runs_table([
    ["2A-1", "pair, closed", "0", "closed", "60", "Ready"],
    ["2A-2", "pair, open", "0", "res", "60", "Ready"],
    ["2B-1…5", "pair, large warm; G ladder", "0.8, 1.6, 3.2, 8, 32", "closed", "60", "Ready"],
    ["2B-6", "pair, large cold", "3.2", "closed", "60", "Ready"],
    ["2C-1…3", "chains: equal, small center, five random", "10 (was 100)", "closed", "60 (was 10)", "Revise"],
    ["2D-1…3", "gap 20, 80, 200 µm", "10 (was 100)", "closed", "60 (was 10)", "Revise"],
])

# ================================================================== 6
H("6. Type 3: lens in a flat pore", 1)
fig("type3", "Figure 4. Type 3 geometry: a 500 by 125 µm channel (also 60 µm wide) with a 150 µm lens, "
    "flat at the start. Drawn from channel_2D_H125um_L500um_eps0.86um_lens150.", width=6.0)
P("**What the set tested.** Whether a lens grows or shrinks at the rate set by its meniscus curvature, "
  "$v = 2 d_0\\cos\\theta / H$ in series with diffusion to the reservoirs, with nothing fitted. "
  "**Result:** the law holds to 1 % for θ = 30–150° and H = 125 and 60 µm; θ = 90° does not move; the closed "
  "control does not move. Report: REPORT_type3 in batch 2026-11-07__08.00.43_type3.")
runs_table([
    ["3A-1…5", "θ = 30, 60, 90, 120, 150°, H 125 µm", "0", "res (h = 1)", "120", "Done"],
    ["3B-1, 2", "θ = 60, 120°, H 60 µm", "0", "res (h = 1)", "120", "Done"],
    ["3C-1", "closed control", "0", "closed", "60", "Done"],
    ["3D-1, 2", "supply at x = 0 only, θ = 60, 120°", "0", "top", "120", "Done"],
], extra="Variable")
P("Open item: the circle-fit contact angle of a moving face is below Young's angle (25° at θ = 30°) while "
  "the speed follows Young's angle (section 12).")

H("7. Type 4: wedge", 1)
fig("type4", "Figure 5. Type 4 geometry: a wedge with wall tilt ω = 10° (also 2.5, 5 and 20°) and a lens "
    "where the gap is 100 µm. Drawn from wedge_2D_om10deg_gap100um_eps0.86um.", width=5.0)
P("**What the set tested.** Whether a taper moves a lens toward the narrow end, through the different "
  "curvature of its two faces (θ − ω and θ + ω). **Result:** yes in every wedge; the series law with "
  "$\\theta \\mp \\omega$ holds to 3–5 %. Every θ = 60° lens reached a vapor face, at day 8 for ω = 20°. Report: "
  "REPORT_type4 in batch 2026-11-07__08.00.47_type4.")
runs_table([
    ["4A-0", "flat channel reference", "0", "res (h = 1)", "180", "Done"],
    ["4A-1…4", "ω = 2.5, 5, 10, 20°", "0", "res (h = 1)", "180", "Done; 4A-4 margin short"],
    ["4B-1, 2", "ω = 10°, θ = 90, 85°", "0", "res (h = 1)", "180", "Done"],
], extra="Variable")

# ================================================================== 8
H("8. Type 5: lens at a throat", 1)
fig("type5", "Figure 6. Type 5 geometry: an 800 µm pore, 150 µm wide, with a 40 µm throat (also 80 µm) at "
    "x = 400 µm and a 120 µm lens at the throat. Ice from the phi_file IC "
    "(t5_pore800x150um_wt40um_lens120c_eps0.80um).")
P("**What the set tests.** A throat is a curvature trap: a lens that spans it has faces of different "
  "curvature unless it sits symmetrically. The set measures where the lens comes to rest with no gradient "
  "(5A-1, 5A-6), how far a gradient displaces it (5A), whether a lens started above or below the throat "
  "returns to it (5C), and whether a throat collects the ice of a nearby cap and grain (5B).")
runs_table([
    ["5A-1…5", "throat 40 µm; G ladder", "0, 1, 3, 10, 30", "res (G 0), ls", "120", "Ready"],
    ["5A-6…8", "throat 80 µm", "0, 3, 10", "res (G 0), ls", "120", "Ready"],
    ["5B-1, 2", "lens + 2 caps + grain", "0, 3", "closed", "120", "Ready"],
    ["5C-1, 2", "lens 100 µm above, below", "3", "ls", "120", "Ready"],
], extra="Variable")
P("All gradients are already at or below 30 K/m; 5A-5 (G = 30 K/m) is the above-threshold point for the "
  "150 µm pore. The set can be submitted as built.")

# ================================================================== 9
H("9. Type 6: lens through a constriction", 1)
fig("type6", "Figure 7. Type 6. (a) The pore that ran: 1100 µm, wide sections 150 µm, a 300 µm neck of "
    "40 µm (also 80 and 150 µm) and the lens 55 µm from the cold-side shoulder. (b) Proposed: 720 µm, a "
    "150 µm neck, the lens 20 µm from the shoulder; built in a scratch folder for this figure only.")
P("**What the set tests.** How much a constriction slows a lens driven by a gradient, and how the delay "
  "scales with the neck width. The measured quantities are the transit time through the neck against neck "
  "width, the mass flux through the neck against the face speed (section_flux), and the change of face "
  "shape in the neck.")
P("**Why it is redesigned.** At G = 100 and 300 K/m the lens split into wall films before reaching the "
  "neck in most runs, so the runs measured the unstable regime of section 3 rather than the delay of a "
  "lens. At G = 10 K/m the lens moves about 10–30 times more slowly; the proposed pore puts it next to "
  "the neck and shortens the neck to keep the run near 180 days (about 11 600 steps at dtmax, an "
  "estimated 35 h at 18 cores; submit with a 3-day limit).")
runs_table([
    ["T6 (7)", "neck 40, 80, 150 µm; θ 60°, one θ 90°", "100, 300", "ls", "60, 30",
     "Ran; superseded (record of the unstable regime)"],
    ["6A-1…3", "neck 40, 80, 150 µm, new pore", "10", "ls", "180", "Proposed"],
    ["6A-4…6", "neck 40, 80, 150 µm, new pore", "30", "ls", "60", "Proposed"],
    ["6B-1", "neck 40 µm, θ = 90° (flat start)", "30", "ls", "60", "Proposed"],
], extra="Variable")
P("The two unfinished G = 300 K/m runs are not extended.")

# ================================================================== 10
H("10. Type 7: rough regolith pore with mixed ice", 1)
fig("type7", "Figure 8. Type 7, layout A: a 2 mm pore with rough walls (35 bumps per wall, gap 50–147 µm, "
    "a 50 µm throat at x = 1 mm), three lenses at 60°, four caps and three grains at 40° (24.6 % ice). "
    "Ice from the phi_file IC (type7_mixed_A_2D_L2000um_eps0.86um). The collector is at x = 0.")
P("**What the set tests.** All the mechanisms of types 1–6 at once in one pore under harvesting: which body "
  "empties first, whether bodies near the collector shield those behind them, whether a throat holds ice "
  "back, and how the harvest rate falls as the dry layer grows. The thermal variants test self-cooling "
  "(7B, insulated dead end) and heating from the surface (7C) or from depth (7D). The measured quantities "
  "are the area of each body against time, the harvest rate at the collector, and the recapture "
  "fraction.")
runs_table([
    ["7A", "isothermal reference", "0", "collector 0.99", "60", "Ready (pre-check first)"],
    ["7B", "dead end insulated", "0", "collector 0.99", "60", "To build"],
    ["7C-1, 2", "warm top (collector side)", "10, 30 (was 100, 300)", "collector 0.99", "60", "To build"],
    ["7D-1, 2", "warm dead end", "10, 30 (was 100, 300)", "collector 0.99", "60", "To build"],
], extra="Variable")
P("Each run is 2.41 M DOF (about 48 cores at 50 000 DOF per core). Pre-check: 7A runs first and is "
  "inspected for Π and step size before 7B–7D are submitted. The ice fraction (24.6 % or 30 %) is "
  "open.")

# ================================================================== 11
H("11. Additional runs", 1)
fig("S", "Figure 9. (a) S7: straight 150 µm pore, 800 µm long, 150 µm lens with 60° faces, 75 µm from the "
    "cold end (proposed; scratch geometry). (b) S5: the 850 µm two-throat pore (throats 30 and 90 µm, both "
    "plugged), sealed at x = L_x (throat_2D_L850um_eps0.86um_narrowL_deadend).")
P("**S7, face-stability scan (proposed, first).** Measures the gradient below which the depositing face "
  "of a lens in a 150 µm pore keeps its meniscus shape, which fixes the production gradients (section 3). "
  "The measured quantity is the lead of the warm face's contact lines over its center (ice_bodies.py: "
  "face_lo minus x_min): stable if it relaxes to the 20 µm of the 60° meniscus, unstable if it grows. The "
  "T6 straight-pore runs (G = 100, 300 K/m) are the upper points.")
runs_table([
    ["S7-1…3", "θ 60°", "3, 10, 30", "ls", "120, 90, 40", "Proposed"],
    ["S7-4", "θ 90° (flat start)", "30", "ls", "40", "Proposed"],
], extra="Variable")
P("**S5, harvest cutoff (to build).** With a collector at h_open and a warm top, ice colder than the "
  "temperature at which ρvs equals the collector density cannot lose vapor to it; the cutoff lies "
  "$dz = \\ln(1/h_{\\mathrm{open}})/(s\\,G)$ beyond the T0 point. The plan's G = 300 and 1000 K/m are replaced "
  "by G = 10 and 30 K/m with h_open = 0.999 and 0.9995, which put dz at 1.04, 0.35, 0.52 and 0.17 mm, "
  "inside or near the 0.85 mm pore (plan_estimates.py, section 5).")
runs_table([
    ["S5-1…4", "h_open 0.999, 0.9995", "10, 30", "collector", "60", "To build"],
], extra="Variable")
P("**S2, throat as an equivalent length (to build).** The same plugged pore with and without one throat, "
  "collector at h_open = 0.99; the delay is expressed as the extra straight length that gives the same "
  "harvest time. Two runs, geometry from the type 5 builder.")
P("**Not scheduled:** S1 (ε pairs and a time-step halving for each new geometry; I1 covers the 0.86 µm "
  "mesh), S3 (insulated vs fixed dead end at three h_open), S4 (channel at θ = 30, 90, 150°), S6 "
  "(axisymmetric runs).")

# ================================================================== 12
H("12. Items to investigate", 1)
T(["Item", "Why it matters", "Next step"], [
    ["Gradient threshold for a stable depositing face", "sets G for types 1, 2, 6, 7 (section 3)", "S7"],
    ["Ice growth of 11–15 % in T6 before reaching x = 0", "the local-saturation ends add water; "
     "the source has not been checked", "water_budget.py on the T6 G = 100 K/m runs"],
    ["Moving-face contact angle below Young (T3)", "contact-angle reporting in moving cases",
     "compare with a dynamic-angle model; decide what to report"],
    ["contact_angle.py fails on runs with several bodies", "types 2, 5, 7 postprocessing", "fix the tool"],
    ["4A-4 (ω = 20°) narrow margin", "narrow-face speed unmeasurable after day 8",
     "optional rerun with a longer narrow end"],
    ["Droplet start angle 40°", "placeholder for liquid water on regolith", "literature (Zhang et al. 2025)"],
    ["Wall angle θ = 60° (D1)", "used in every run", "confirm"],
    ["Lunar gradient range", "framing of the chosen G", "find a cited range"],
], [2.2, 2.3, 2.1], size=8)

# ================================================================== 13
H("13. Order of work", 1)
for line in [
    "Submit type 5 (12 runs, as built) and the 7A pre-check.",
    "Build and run S7 (4 runs, about 1 M DOF each).",
    "From S7, fix the production gradients; rebuild the gradient runs of types 1 and 2 and build the "
    "type 6 redesign and 7B–7D; submit types 1, 2 and 6.",
    "After 7A is inspected, submit 7B–7D; then S5 and S2.",
]:
    B(line, style="List Number")

H("References", 2)
P("Kaempfer, T. U. and Plapp, M. (2009). Phase-field modeling of dry snow metamorphism. Phys. Rev. E 79, "
  "031502. doi:10.1103/PhysRevE.79.031502", size=8.5)
P("Mullins, W. W. and Sekerka, R. F. (1963). Morphological stability of a particle growing by diffusion or "
  "heat flow. J. Appl. Phys. 34, 323–329. doi:10.1063/1.1702607", size=8.5)
P("Mullins, W. W. and Sekerka, R. F. (1964). Stability of a planar interface during solidification of a "
  "dilute binary alloy. J. Appl. Phys. 35, 444–451. doi:10.1063/1.1713333", size=8.5)

doc.save(OUTPATH)
print("wrote", OUTPATH)

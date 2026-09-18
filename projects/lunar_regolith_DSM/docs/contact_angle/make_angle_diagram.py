#!/usr/bin/env python3
"""Diagram of how the contact angle is measured from a simulation.

Draws a real diffuse phase field -- the model's own tanh profile in the signed
distance to a circular cap -- rather than a cartoon, so the phi = 0.5 contour
shown is genuinely the level set the measurement uses.

Output: contact_angle_diagram.{pdf,png}, 5.4 x 3.25 in, all text >= 16 pt.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Arc, Polygon

# ---------------------------------------------------------------- style ----
# Aptos if it is installed, otherwise STIX Two Text -- the Times-based face the
# LaTeX math fonts are drawn from, so the two look consistent on a slide.
import matplotlib.font_manager as _fm
_have = {f.name for f in _fm.fontManager.ttflist}
FONT = next((f for f in ("Aptos", "Aptos Display", "STIX Two Text") if f in _have),
            "DejaVu Serif")
plt.rcParams.update({
    "font.family": FONT,
    "mathtext.fontset": "stix",
    "font.size": 16,
    "pdf.fonttype": 42,          # embed as TrueType so slides can edit/scale
    "ps.fonttype": 42,
})

THETA = 60.0                     # contact angle to draw, degrees
R     = 1.0                      # cap radius
EPS   = 0.030                    # interface decay length, in cap radii

ICE  = "#a9d8e0"      # pale glacial cyan -- reads as ice, not water
LINE = "#12283a"
SED  = "#6b5744"
TAN  = "#c2410c"

th  = np.radians(THETA)
yc  = -R * np.cos(th)            # centre below the wall => cap sits above it
xP  =  R * np.sin(th)            # right contact point

# ------------------------------------------------------------- the field ----
xr, yr = 1.255, 1.28
x = np.linspace(-xr, xr, 1400)
y = np.linspace(-0.30, -0.30 + yr, 760)
X, Y = np.meshgrid(x, y)
d   = R - np.hypot(X - 0.0, Y - yc)          # >0 inside the cap
phi = 0.5 * (1.0 + np.tanh(d / (2.0 * EPS)))
phi[Y < 0.0] = 0.0                            # the wall truncates the drop

fig, ax = plt.subplots(figsize=(5.4, 3.25))
cmap = LinearSegmentedColormap.from_list("ice", ["#ffffff", ICE])
ax.imshow(phi, extent=[x[0], x[-1], y[0], y[-1]], origin="lower",
          cmap=cmap, vmin=0, vmax=1, interpolation="bilinear", zorder=1)

# ------------------------------------------------------------- sediment ----
ax.add_patch(Polygon([[x[0], y[0]], [x[-1], y[0]], [x[-1], 0], [x[0], 0]],
                     closed=True, facecolor=SED, edgecolor="none", zorder=3))
ax.plot([x[0], x[-1]], [0, 0], color=LINE, lw=2.2, zorder=4)
ax.text(0.0, -0.155, "S E D I M E N T", color="white", zorder=5,
        va="center", ha="center", fontsize=16)

# ------------------------------------------------- the phi = 0.5 contour ----
ax.contour(X, Y, np.where(Y >= 0, phi, np.nan), levels=[0.5],
           colors=[LINE], linewidths=2.4, zorder=6)

# --------------------------------------------------------- tangent line ----
# tangent at P is perpendicular to the radius; it points up-left into the drop
tx, ty = -np.cos(th), np.sin(th)
L0, L1 = 0.13, 0.88
ax.plot([xP - L0 * tx, xP + L1 * tx], [ty * -L0, ty * L1],
        color=TAN, lw=2.2, ls=(0, (5, 3)), zorder=7)
ax.plot([xP], [0], "o", ms=7, color=LINE, zorder=8)

# ------------------------------------------------------------- the angle ----
ax.add_patch(Arc((xP, 0), 0.60, 0.60, angle=0, theta1=180.0 - THETA,
                 theta2=180.0, color=LINE, lw=1.8, zorder=8))
ax.text(xP + 0.400 * np.cos(np.radians(180 - THETA / 2)),
        0.400 * np.sin(np.radians(180 - THETA / 2)) + 0.010,
        r"$\theta$", fontsize=19, color=LINE, ha="center", va="center", zorder=9)

# ----------------------------------------------------------- phase labels ---
ax.text(-0.10, 0.17, r"$\phi_i$", fontsize=19, ha="center", va="center",
        color=LINE, zorder=7)
ax.text(-0.98, 0.72, r"$\phi_a$", fontsize=19, ha="center", va="center",
        color=LINE, zorder=7)

# ------------------------------------------------- contour + tangent notes --
ax.text(0.0, 0.615, r"$\phi_i = 0.5$", fontsize=16, color=LINE,
        ha="center", va="bottom", zorder=9)


ax.set_xlim(x[0], x[-1]); ax.set_ylim(y[0], y[-1])
ax.set_aspect("equal"); ax.axis("off")
fig.subplots_adjust(left=0.005, right=0.995, top=0.995, bottom=0.005)
for ext in ("pdf", "png"):
    fig.savefig(f"docs/contact_angle/contact_angle_diagram.{ext}",
                dpi=400, transparent=False, facecolor="white")
print("wrote docs/contact_angle/contact_angle_diagram.{pdf,png}")

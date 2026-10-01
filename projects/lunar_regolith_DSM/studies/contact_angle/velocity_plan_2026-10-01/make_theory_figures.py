#!/usr/bin/env python3
"""Sharp-interface predictions for the contact-angle / interface-velocity study.

Theory only -- no simulation output is read. These are the curves the four
planned batches (channel and wedge, each swept in theta and in sigma_inf) are
expected to land on, drawn before the runs so the comparison is a prediction
rather than a fit.

MODEL

Interface condition (the one wedge_gt_velocity.py verified on the solver):

    sigma_i = d0 * chi + beta * v_n                               (Gibbs-Thomson)

Quasi-steady vapour diffusion between the meniscus and a fixed-sigma wall,
with the mass balance rho_ice * v_n = D_v * rho_vs * d(sigma)/dn, adds a
diffusive resistance in series with the kinetic one:

    v_n = (sigma_inf - d0 * chi) / (beta + Z_D)

    channel:  Z_D = K * l                 l   = meniscus-to-wall distance
    wedge:    Z_D = K * r * ln(r_res/r)   radial diffusion about the apex
    K = rho_ice / (rho_vs * D_v)

CURVATURE  (chi > 0 where the ice is convex into the vapour; theta is measured
through the ice; alpha is the wedge half-angle; r is the meniscus position on
the centreline, measured from the virtual apex)

    channel:       chi = -2 cos(theta) / H
    wedge, inner:  chi = -(cos(theta) + sin(alpha)) / (r sin(alpha))
    wedge, outer:  chi = -(cos(theta) - sin(alpha)) / (r sin(alpha))

Both wedge forms reduce to the -1/r, +1/r of the apex-centred 90-degree band.

BETA. beta_sub0 = 7.9408e3 / alpha_c [s/m] at -20 C (Hertz-Knudsen). The
measured interface coefficient on the 2026-08-07 wedge_bc batch was
beta_sub0 + a1*a2*eps*(1/D_T + 1/D_v)*rho_ice/rho_vs, an ADDITIVE offset of
8.74e5 s/m at this eps: +22 % at alpha_c = 2e-3, +11 % at 1e-3, +110 % at 1e-2.
--beta-offset 0 draws the curves for a solver that removes it exactly.

Usage:
    venv_lunar/bin/python3 studies/contact_angle/velocity_plan_2026-10-01/make_theory_figures.py \\
        [--alpha 1e-3] [--beta-offset 8.74e5]

Figures 2-4 are written with an _ac<alpha> suffix; figure 1 is alpha-independent.
"""
import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Arc, Polygon
from scipy.integrate import solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "postprocess"))
from pplib import rho_vs                                    # noqa: E402

# --- parameters, as the -20 C contact-angle runs set them --------------------
T_C = -20.0
D0 = 1.0166e-9                     # -d0_sub0 [m]
BETA_HK1 = 7.9408e3                # beta_sub0 at alpha_c = 1 [s/m]
BETA_OFFSET = 0.22 * 3.9704e6      # measured thin-interface offset [s/m]
BETA = BETA_HK1 / 2e-3 + BETA_OFFSET   # reset from --alpha in __main__
TAG = ""
RHO_ICE = 919.0
D_V = 2.178e-5 * ((T_C + 273.15) / 273.15) ** 1.81         # VaporDiffus()
K = RHO_ICE / (float(rho_vs(T_C)) * D_V)                   # [s/m^2]

# channel_2D_H125um_eps0.86um_flat90
H = 125e-6
L_CH = 101e-6                      # meniscus-to-wall distance at t = 0
W_CH = 172.5e-6                    # bridge width at t = 0

# wedge_2D_L300um_eps0.86um_wall: apex 100 um left of the domain, slopes +/-0.25
ALPHA = np.arctan(0.25)
SA = np.sin(ALPHA)
R_L, R_R = 100e-6, 400e-6          # the two fixed-sigma walls
R1, R2 = 200e-6, 300e-6            # wedge_band IC

DAY = 86400.0
NM_DAY = 1e9 * DAY                 # m/s -> nm/day

# --- style -------------------------------------------------------------------
INK, MUTED, GRID = "#1c2430", "#5b6673", "#d9dee4"
BLUE, ORANGE, GREY = "#2a78d6", "#eb6834", "#7b8794"
ICE = "#cfe3f7"
THETA_RAMP = ["#9cc4f0", "#5a9be0", "#2a78d6", "#1c569c", "#0f3462"]
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 15, "axes.labelsize": 16,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
    "lines.linewidth": 2.2, "legend.frameon": False, "mathtext.default": "regular",
})


def save(fig, name):
    fig.savefig(os.path.join(HERE, name), dpi=200, transparent=True)
    plt.close(fig)
    print("wrote", name)


def panel(ax, letter):
    ax.text(-0.14, 1.04, f"({letter})", transform=ax.transAxes,
            fontsize=16, fontweight="bold", va="bottom")


# --- theory ------------------------------------------------------------------
def chi_channel(th):
    return -2.0 * np.cos(th) / H


def v_channel(th, sig, ell=L_CH):
    return (sig - D0 * chi_channel(th)) / (BETA + K * ell)


def chi_in(th, r):
    return -(np.cos(th) + SA) / (r * SA)


def chi_out(th, r):
    return -(np.cos(th) - SA) / (r * SA)


def v_in(th, sig, r):
    """Growth speed of the inner meniscus; growth moves it to SMALLER r."""
    return (sig - D0 * chi_in(th, r)) / (BETA + K * r * np.log(r / R_L))


def v_out(th, sig, r):
    return (sig - D0 * chi_out(th, r)) / (BETA + K * r * np.log(R_R / r))


def wedge_ode(th, sig, t_end, n=400):
    """Integrate the two meniscus positions; stops if either reaches a wall."""
    rhs = lambda t, y: [-v_in(th, sig, y[0]), v_out(th, sig, y[1])]
    hit_l = lambda t, y: y[0] - 1.02 * R_L
    hit_r = lambda t, y: 0.98 * R_R - y[1]
    gone = lambda t, y: y[1] - y[0] - 2e-6
    for ev in (hit_l, hit_r, gone):
        ev.terminal = True
    s = solve_ivp(rhs, (0, t_end), [R1, R2], events=(hit_l, hit_r, gone),
                  t_eval=np.linspace(0, t_end, n), rtol=1e-9, atol=1e-12)
    return s.t, s.y[0], s.y[1]


# --- figure 1: the two geometries -------------------------------------------
def fig_geometry():
    th = np.radians(60.0)
    fig, axs = plt.subplots(1, 2, figsize=(15, 4.6), gridspec_kw={"wspace": 0.08})
    um = 1e6

    # channel
    ax = axs[0]
    Lx, xl, xr = 375e-6, 101e-6, 274e-6
    R = H / (2 * np.cos(th))
    y = np.linspace(0, H, 80)
    bulge = np.sqrt(R ** 2 - (y - H / 2) ** 2) - R       # 0 on the centreline
    left = np.c_[xl + bulge, y]
    right = np.c_[xr - bulge, y][::-1]
    ax.add_patch(Polygon(np.r_[left, right] * um, fc=ICE, ec=BLUE, lw=2))
    ax.plot([0, Lx * um], [0, 0], color=INK, lw=3)
    ax.plot([0, Lx * um], [H * um, H * um], color=INK, lw=3)
    for x in (0, Lx * um):
        ax.plot([x, x], [0, H * um], color=ORANGE, lw=2.5, ls=(0, (4, 2)))
    cp = left[-1] * um
    ax.add_patch(Arc(cp, 44, 44, theta1=-60, theta2=0, color=INK, lw=1.5))
    ax.text(cp[0] + 27, cp[1] - 17, r"$\theta$", fontsize=17)
    ax.annotate("", (xl * um - 34, H * um / 2), (xl * um - 2, H * um / 2),
                arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.8))
    ax.annotate("", (xr * um + 34, H * um / 2), (xr * um + 2, H * um / 2),
                arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.8))
    ax.text(xl * um - 40, H * um / 2 + 8, r"$v_n$", ha="right", fontsize=16)
    ax.text(xr * um + 40, H * um / 2 + 8, r"$v_n$", ha="left", fontsize=16)
    ax.text((xl + xr) / 2 * um, H * um / 2, "ice", ha="center", va="center",
            fontsize=17, color=THETA_RAMP[4])
    ax.annotate("", (Lx * um - 28, 0), (Lx * um - 28, H * um),
                arrowprops=dict(arrowstyle="<->", color=MUTED, lw=1.2))
    ax.text(Lx * um - 34, H * um * 0.2, "H", ha="right", va="center", color=MUTED)
    ax.annotate("", (0, -14), (cp[0], -14),
                arrowprops=dict(arrowstyle="<->", color=MUTED, lw=1.2))
    ax.text(cp[0] / 2, -30, r"$\ell$", ha="center", color=MUTED, fontsize=17)
    ax.text(-8, H * um / 2, r"$\sigma_\infty$", ha="right", va="center",
            color=ORANGE, fontsize=18)
    ax.text(Lx * um + 8, H * um / 2, r"$\sigma_\infty$", ha="left", va="center",
            color=ORANGE, fontsize=18)
    ax.set_xlim(-45, Lx * um + 45)
    ax.set_ylim(-110, H * um + 100)

    # wedge
    ax = axs[1]
    ri, ro = R1, R2
    Ri = ri * SA / (np.cos(th) + SA)
    Ro = ro * SA / (np.cos(th) - SA)
    yi = np.linspace(-1, 1, 80) * Ri * np.cos(th - ALPHA)
    yo = np.linspace(-1, 1, 80) * Ro * np.cos(th + ALPHA)
    inner = np.c_[(ri - Ri) + np.sqrt(Ri ** 2 - yi ** 2), yi]
    outer = np.c_[(ro + Ro) - np.sqrt(Ro ** 2 - yo ** 2), yo][::-1]
    ax.add_patch(Polygon(np.r_[inner, outer] * um, fc=ICE, ec=BLUE, lw=2))
    for s in (1, -1):
        ax.plot([R_L * um, R_R * um], [s * R_L * um * 0.25, s * R_R * um * 0.25],
                color=INK, lw=3)
        ax.plot([0, R_L * um], [0, s * R_L * um * 0.25], color=MUTED, lw=1, ls=":")
    for r in (R_L, R_R):
        ax.plot([r * um] * 2, [-r * um * 0.25, r * um * 0.25],
                color=ORANGE, lw=2.5, ls=(0, (4, 2)))
    ax.plot([0, R_R * um], [0, 0], color=MUTED, lw=1, ls=":")
    ax.add_patch(Arc((0, 0), 110, 110, theta1=0, theta2=np.degrees(ALPHA),
                     color=MUTED, lw=1.2))
    ax.text(62, 4, r"$\alpha$", color=MUTED, fontsize=16)
    cp = outer[0] * um
    a0 = 180 + np.degrees(ALPHA)
    ax.add_patch(Arc(cp, 50, 50, theta1=a0, theta2=a0 + 60, color=INK, lw=1.5))
    ax.text(cp[0] - 40, cp[1] - 22, r"$\theta$", fontsize=17)
    ax.text((ri + ro) / 2 * um, 22, "ice", ha="center", va="center", fontsize=17,
            color=THETA_RAMP[4])
    for r, lab, dy in ((ri, r"$r_{in}$", -92), (ro, r"$r_{out}$", -116)):
        ax.annotate("", (0, dy), (r * um, dy),
                    arrowprops=dict(arrowstyle="<->", color=MUTED, lw=1.2))
        ax.plot([r * um] * 2, [0, dy], color=MUTED, lw=0.8, ls=":")
        ax.text(r * um / 2, dy + 5, lab, ha="center", color=MUTED, fontsize=16)
    xh = 362e-6                                  # local channel height H(x)
    ax.annotate("", (xh * um, -xh * um * 0.25), (xh * um, xh * um * 0.25),
                arrowprops=dict(arrowstyle="<->", color=MUTED, lw=1.2))
    ax.text(xh * um + 5, -42, "H(x)", ha="left", va="center", color=MUTED)
    ax.text(R_L * um - 8, 44, r"$\sigma_\infty$", ha="right", color=ORANGE, fontsize=18)
    ax.text(R_R * um + 8, 0, r"$\sigma_\infty$", ha="left", va="center",
            color=ORANGE, fontsize=18)
    ax.set_xlim(-20, R_R * um + 50)
    ax.set_ylim(-130, 120)

    for ax in axs:
        ax.set_aspect("equal")
        ax.axis("off")
    fig.text(0.25, 0.93, "Rectangular channel", ha="center", fontsize=18)
    fig.text(0.75, 0.93, "Wedge", ha="center", fontsize=18)
    fig.subplots_adjust(left=0.02, right=0.98, top=0.9, bottom=0.02)
    save(fig, "fig1_geometry.png")


# --- figure 2: channel -------------------------------------------------------
def fig_channel():
    fig, axs = plt.subplots(1, 2, figsize=(15, 5.6))
    th = np.radians(np.linspace(20, 160, 200))
    ax = axs[0]
    for sig, c, lab in ((2e-5, ORANGE, r"$\sigma_\infty=+2\times10^{-5}$"),
                        (0.0, GREY, r"$\sigma_\infty=0$"),
                        (-2e-5, BLUE, r"$\sigma_\infty=-2\times10^{-5}$")):
        v = v_channel(th, sig) * NM_DAY
        ax.plot(np.degrees(th), v, color=c)
        ax.text(162, v[-1], lab, color=INK, va="center", fontsize=14)
    ax.axhline(0, color=MUTED, lw=1)
    ax.set_xlim(20, 215)
    ax.set_xticks([30, 60, 90, 120, 150])
    ax.set_xlabel(r"contact angle $\theta$ [deg]")
    ax.set_ylabel(r"interface velocity $v_n$ [nm/day]")
    panel(ax, "a")

    ax = axs[1]
    sig = np.linspace(-3e-5, 3e-5, 100)
    for thd, c in zip((30, 60, 90, 120, 150), THETA_RAMP):
        v = v_channel(np.radians(thd), sig) * NM_DAY
        ax.plot(sig * 1e5, v, color=c)
        ax.text(3.1, v[-1], rf"$\theta={thd}°$", color=INK, va="center", fontsize=14)
    ax.axhline(0, color=MUTED, lw=1)
    ax.set_xlim(-3, 3.9)
    ax.set_xticks([-3, -2, -1, 0, 1, 2, 3])
    ax.set_xlabel(r"wall supersaturation $\sigma_\infty$  [$10^{-5}$]")
    ax.set_ylabel(r"interface velocity $v_n$ [nm/day]")
    panel(ax, "b")
    fig.tight_layout()
    save(fig, "fig2_channel%s.png" % TAG)


# --- figure 3: wedge curvature and the ODE ----------------------------------
def fig_wedge_ode():
    fig, axs = plt.subplots(1, 2, figsize=(15, 5.6))
    ax = axs[0]
    r = np.linspace(R_L, R_R, 200)
    for thd, c in zip((45, 90, 135), (THETA_RAMP[0], THETA_RAMP[2], THETA_RAMP[4])):
        th = np.radians(thd)
        ax.plot(r * 1e6, chi_in(th, r) * 1e-3, color=c)
        ax.plot(r * 1e6, chi_out(th, r) * 1e-3, color=c, ls=(0, (5, 2)))
        ax.text(404, 0.5e-3 * (chi_in(th, r[-1]) + chi_out(th, r[-1])),
                rf"$\theta={thd}°$", color=INK, va="center", fontsize=14)
    ax.plot([], [], color=MUTED, label="inner meniscus")
    ax.plot([], [], color=MUTED, ls=(0, (5, 2)), label="outer meniscus")
    ax.legend(loc="lower right", fontsize=14, ncol=2)
    ax.axhline(0, color=MUTED, lw=1)
    ax.set_xlim(100, 462)
    ax.set_ylim(-45, 45)
    ax.set_xlabel(r"meniscus position from apex $r$ [µm]")
    ax.set_ylabel(r"curvature $\chi$ [mm$^{-1}$]")
    panel(ax, "a")

    ax = axs[1]
    for thd, c in ((60, BLUE), (120, ORANGE)):
        t, ri, ro = wedge_ode(np.radians(thd), 0.0, 150 * DAY)
        ax.plot(t / DAY, (R1 - ri) * 1e6, color=c)
        ax.plot(t / DAY, (ro - R2) * 1e6, color=c, ls=(0, (5, 2)))
        ax.text(152, (R1 - ri[-1]) * 1e6, rf"$\theta={thd}°$ inner", va="center", fontsize=14)
        ax.text(152, (ro[-1] - R2) * 1e6, rf"$\theta={thd}°$ outer", va="center", fontsize=14)
        da = (ro[-1] ** 2 - ri[-1] ** 2) / (R2 ** 2 - R1 ** 2) - 1
        print(f"  wedge theta={thd}: r_in {ri[-1]*1e6:.1f}  r_out {ro[-1]*1e6:.1f} um"
              f"  ice area {100*da:+.1f} % at 150 d")
    ax.axhline(0, color=MUTED, lw=1)
    ax.set_xlim(0, 198)
    ax.set_xticks([0, 30, 60, 90, 120, 150])
    ax.set_xlabel("time [days]")
    ax.set_ylabel("meniscus advance into vapour [µm]")
    panel(ax, "b")
    fig.tight_layout()
    save(fig, "fig3_wedge_ode%s.png" % TAG)


# --- figure 4: wedge predictions --------------------------------------------
def fig_wedge():
    fig, axs = plt.subplots(1, 2, figsize=(15, 5.6))
    ax = axs[0]
    th = np.radians(np.linspace(20, 160, 200))
    vi, vo = v_in(th, 0.0, R1) * NM_DAY, v_out(th, 0.0, R2) * NM_DAY
    ax.plot(np.degrees(th), vi, color=BLUE, label=r"inner, $r_{in}=200$ µm")
    ax.plot(np.degrees(th), vo, color=ORANGE, label=r"outer, $r_{out}=300$ µm")
    ax.plot(np.degrees(th), v_channel(th, 0.0) * NM_DAY, color=GREY, lw=1.6, ls=":",
            label="channel, H = 125 µm")
    ax.legend(loc="upper right", fontsize=14)
    for x0, lab in ((90 + np.degrees(ALPHA), r"$90°+\alpha$"),
                    (90 - np.degrees(ALPHA), r"$90°-\alpha$")):
        ax.plot([x0], [0], "o", color=INK, ms=7, zorder=5)
    ax.annotate(r"$90°-\alpha$", (90 - np.degrees(ALPHA), 0), (38, -75), fontsize=14,
                arrowprops=dict(arrowstyle="-", color=MUTED, lw=1))
    ax.annotate(r"$90°+\alpha$", (90 + np.degrees(ALPHA), 0), (118, 40), fontsize=14,
                arrowprops=dict(arrowstyle="-", color=MUTED, lw=1))
    ax.axhline(0, color=MUTED, lw=1)
    ax.set_xticks([30, 60, 90, 120, 150])
    ax.set_xlabel(r"contact angle $\theta$ [deg]")
    ax.set_ylabel(r"growth velocity $v_n$ [nm/day]")
    panel(ax, "a")

    ax = axs[1]
    sig = np.linspace(-3e-5, 3e-5, 100)
    th = np.radians(60)
    vi, vo = v_in(th, sig, R1) * NM_DAY, v_out(th, sig, R2) * NM_DAY
    ax.plot(sig * 1e5, vi, color=BLUE, label=r"inner, $r_{in}=200$ µm")
    ax.plot(sig * 1e5, vo, color=ORANGE, label=r"outer, $r_{out}=300$ µm")
    ax.legend(loc="lower right", fontsize=14, title=r"$\theta=60°$", title_fontsize=14)
    for chi, c in ((chi_in(th, R1), BLUE), (chi_out(th, R2), ORANGE)):
        ax.plot([D0 * chi * 1e5], [0], "o", color=c, ms=8, mec="white", mew=1.5, zorder=5)
    ax.annotate(r"$v_n=0$ at $\sigma_\infty=d_0\chi$", (D0 * chi_in(th, R1) * 1e5, 0),
                (-2.9, 150), fontsize=14,
                arrowprops=dict(arrowstyle="-", color=MUTED, lw=1))
    ax.axhline(0, color=MUTED, lw=1)
    ax.set_xlim(-3, 3)
    ax.set_xticks([-3, -2, -1, 0, 1, 2, 3])
    ax.set_xlabel(r"wall supersaturation $\sigma_\infty$  [$10^{-5}$]")
    ax.set_ylabel(r"growth velocity $v_n$ [nm/day]")
    panel(ax, "b")
    fig.tight_layout()
    save(fig, "fig4_wedge%s.png" % TAG)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--alpha", type=float, default=1e-3, help="attachment coefficient alpha_c")
    ap.add_argument("--beta-offset", type=float, default=BETA_OFFSET,
                    help="additive thin-interface offset on beta [s/m]; 0 = fully corrected")
    args = ap.parse_args()
    BETA = BETA_HK1 / args.alpha + args.beta_offset
    TAG = "_ac%.0e" % args.alpha
    TAG = TAG.replace("e-0", "e-")
    print(f"alpha_c = {args.alpha:g}: beta_sub0 = {BETA_HK1/args.alpha:.3e}, "
          f"kinetic share of channel resistance = {BETA/(BETA + K*L_CH):.2f}")
    print(f"K = {K:.3e} s/m^2   beta_eff = {BETA:.3e} s/m   "
          f"Z_D(channel) = {K*L_CH:.3e} s/m   alpha = {np.degrees(ALPHA):.2f} deg")
    for thd in (60, 120):
        v = v_channel(np.radians(thd), 0.0)
        print(f"  channel theta={thd}: d0*chi = {D0*chi_channel(np.radians(thd)):+.2e}"
              f"  v = {v*NM_DAY:+.1f} nm/day"
              f"  ice {100*2*v*90*DAY/W_CH:+.1f} % at 90 d (fixed-l estimate)")
    fig_geometry()
    fig_channel()
    fig_wedge_ode()
    fig_wedge()

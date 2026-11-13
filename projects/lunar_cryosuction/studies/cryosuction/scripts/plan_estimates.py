"""Order-of-magnitude numbers for the paper-2 simulation plan (2026-11-06).

Every number quoted in the simulation proposal's "scales" boxes comes from
this script. Constants come from preprocess/mesh_params (which mirrors the C
solver), so rho_vs, D_v(T), d0 and beta0 are the solver's values at T0.
These are estimates for choosing parameters and run lengths, not results.

Sections
  1. B1 series-law slope: predicted with D_v at T0 and with D_v at 0 C, against
     the measured 5.73e12 s/m^2 (studies/cryosuction/RUNS.md, row 0c-B1).
  2. Ripening pair: Gibbs-Thomson drive, the gradient that cancels it, and a
     rough vanishing time (2D quasi-steady diffusion in series with kinetics).
  3. Gradient-driven migration speed of ice (and K&P Eq. 45 resolution number
     Pi = beta0 v eps / d0 at --eps; the design limit is 0.5), v ~ D_v * d(rho_vs)/dT * G / rho_ice
     (no conductivity-contrast factor; ice conducts ~100x better than air, so
     the real local gradient in the air differs by an O(1) factor).
  4. Harvest time of a filled channel, front moving from l0 to l1 under
     suction sigma: t = [beta0 (l1 - l0) + K (l1^2 - l0^2) / 2] * f / sigma,
     K = rho_ice / (D_v rho_vs), f = ice fraction of the pore cross-section.
  5. Warm-top (or warm-bottom) gradient runs: the collector holds
     rho_open = h_open * rho_vs(T0), with T0 the mid-domain temperature
     (src/cryosuction_main.c:1415). Ice colder than T_c, where
     rho_vs(T_c) = rho_open, cannot lose vapor to the collector. Cutoff
     distance on the cold side of the T0 point: dz = ln(1/h_open) / (s G).
  6. Gradient against capillarity (added 2026-11-11 after the type-6 runs).
     (a) Crossover: the gradient whose saturation change over a length l equals
     the Gibbs-Thomson change of a feature of curvature chi, G_x = d0 chi /
     (s l). Below G_x curvature (ripening, contact angle) controls the ice;
     above it the gradient does. (b) Stability of a depositing face
     (Mullins-Sekerka, quasi-static one-sided diffusion, kinetics left out:
     they slow the growth of a mode but do not move its neutral wavenumber).
     A face growing into a relative vapor gradient g = s G is unstable to
     modes with k^2 < g / d0. Across a pore of width W with no-flux walls the
     symmetric mode (both contact lines ahead of the axis, the shape a
     theta < 90 meniscus already has) has k = 2 pi / W, so the face is stable
     for G < G_s = (2 pi / W)^2 d0 / s; the antisymmetric (tilting) mode,
     k = pi / W, gives a 4x lower bound. A sublimating face is stable.
     These are linear, flat-face estimates; the threshold is to be measured.
  7. Choosing G (added 2026-11-12, after the theta-90 type-6 run).
     (a) The user's estimate: a 2D meniscus at angle theta in a pore of width W
     is an arc of radius R = W/(2 cos theta), curvature chi = 1/R, and its
     contact lines sit dx = R (1 - sin theta) ahead of its axis point. The
     Gibbs-Thomson shift d0 chi is the same everywhere on the arc; the gas
     saturation differs by s G dx between the two x positions. They are equal
     at G* = d0 chi / (s dx).
     (b) Gas gradient. Ice conducts 115x better than air (k 2.29 against
     0.02 W/m/K, src/cryosuction_main.c), so the temperature drop falls on the
     gas, and in a neck of width wc it is W/wc times steeper than in the wide
     sections (insulating walls, 1D series estimate; shoulders counted at the
     mean width). Measured in the 6B-1 run (wc 40, G 300): 164-174 K/m in the
     wide gas, 614-640 in the neck, 2 in the ice.
     (c) Growth rate of the symmetric face mode (k = 2 pi/W) with the ice
     isothermal: sigma = (D_v rho_vs/rho_ice) k (s G_gas - d0 k^2); measured
     1.18/day in 6B-1 against 0.90/day from this formula.
     (d) Cost: time for a face to move --move at v = D_v rho_vs s G_gas/rho_ice
     (the 6B-1 center moved at 24 um/day at G_gas 169 K/m, 0.14 um/day per
     K/m), steps at dt = min(dtmax, --cfl-step / v) (the interface-CFL cap held
     the fastest point to 0.10-0.16 um per step in T6) and --s-per-step.

    venv_cryosuction/bin/python studies/cryosuction/scripts/plan_estimates.py
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from preprocess.mesh_params import physics  # noqa: E402
from preprocess.mesh_params.constants import CONSTANTS as C  # noqa: E402

DAY = 86400.0


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--T0", type=float, default=-20.0, help="temperature [C]")
    p.add_argument("--alpha", type=float, default=1e-2, help="alpha_c")
    p.add_argument("--b1-slope", type=float, default=5.73e12, help="measured B1 slope of 1/v vs l [s/m^2]")
    p.add_argument("--b1-sigma", type=float, default=0.01, help="B1 suction 1 - h")
    p.add_argument("--pair", type=float, nargs=3, default=[20e-6, 50e-6, 200e-6],
                   help="ripening pair: R_small R_large center distance [m]")
    p.add_argument("--G", type=float, nargs="+", default=[1, 10, 100, 1000], help="gradients [K/m]")
    p.add_argument("--harvest", type=float, nargs="+", default=[0.99, 0.9, 0.5], help="open-end h values")
    p.add_argument("--l0", type=float, default=1e-4, help="harvest: start distance [m]")
    p.add_argument("--l1", type=float, nargs="+", default=[0.85e-3, 2e-3, 5e-3], help="harvest: end distances [m]")
    p.add_argument("--fill", type=float, default=1.0, help="harvest: ice fraction of the pore cross-section")
    p.add_argument("--cutoff-h", type=float, nargs="+", default=[0.99, 0.995, 0.999], help="h_open for section 5")
    p.add_argument("--widths", type=float, nargs="+", default=[40e-6, 80e-6, 125e-6, 150e-6, 300e-6],
                   help="section 6: pore widths or feature sizes [m]")
    p.add_argument("--eps", type=float, default=8.584e-7, help="eps for the resolution number Pi [m]")
    p.add_argument("--G7", type=float, nargs="+", default=[3, 5, 10, 15, 20, 30, 50, 100, 300],
                   help="section 7: nominal gradients [K/m]")
    p.add_argument("--theta", type=float, default=60.0, help="section 7a: contact angle [deg]")
    p.add_argument("--move", type=float, default=300e-6, help="section 7d: distance the face must move [m]")
    p.add_argument("--dtmax", type=float, default=1341.2, help="section 7d: dtmax = 2 tau_sub [s]")
    p.add_argument("--cfl-step", type=float, default=0.12e-6, help="section 7d: face travel per CFL-capped step [m]")
    p.add_argument("--s-per-step", type=float, default=17.0, help="section 7d: wall seconds per step (T6: 13-28)")
    p.add_argument("--wt", type=float, nargs="+", default=[30e-6, 40e-6, 60e-6, 80e-6, 100e-6, 120e-6],
                   help="section 8: throat widths [m]")
    p.add_argument("--W8", type=float, default=150e-6, help="section 8: pore width [m]")
    p.add_argument("--bump-R", type=float, default=250e-6, help="section 8: bump half-width (fixed) [m]")
    p.add_argument("--lens-len", type=float, default=120e-6, help="section 8: lens length on the axis [m]")
    p.add_argument("--Lx8", type=float, default=800e-6, help="section 8: pore length, throat at Lx/2 [m]")
    p.add_argument("--Lx8c", type=float, default=400e-6, help="section 8c: pore length (no flat ends) [m]")
    p.add_argument("--lens8c", type=float, default=100e-6, help="section 8c: initial lens length [m]")
    p.add_argument("--start8c", type=float, default=70e-6, help="section 8c: lens center start, cold side [m]")
    p.add_argument("--start8", type=float, default=100e-6,
                   help="section 8b: lens center start, this far on the COLD side of the throat [m]")
    a = p.parse_args()
    if not (-150.0 < a.T0 < -0.001):
        p.error(f"--T0 must be in (-150, -0.001) C, where the ice formulas apply (got {a.T0})")
    if not (0.0 < a.alpha <= 1.0):
        p.error(f"--alpha must be in (0, 1] (got {a.alpha})")
    if any(not (g > 0.0) for g in a.G):
        p.error("--G values must be > 0 (sections 2 and 5 divide by G)")
    Rs_, Rl_, L_ = a.pair
    if not (0.0 < Rs_ < Rl_):
        p.error(f"--pair needs 0 < R_small < R_large (got {Rs_:g}, {Rl_:g})")
    if not (L_ > Rs_ + Rl_):
        p.error(f"--pair: center distance {L_:g} m must exceed R_small + R_large = {Rs_ + Rl_:g} m "
                "(the grains would overlap)")
    if not (0.0 < a.fill <= 1.0):
        p.error(f"--fill must be in (0, 1] (got {a.fill})")
    if not (a.l0 >= 0.0):
        p.error(f"--l0 must be >= 0 (got {a.l0})")
    if any(not (0.0 < h <= 1.0) for h in a.cutoff_h):
        p.error("--cutoff-h values must be in (0, 1]")
    if any(not (0.0 < h < 1.0) for h in a.harvest):
        p.error("--harvest values must be in (0, 1): at h >= 1 the open end supplies vapor and "
                "nothing is harvested")
    if any(not (l1 > a.l0) for l1 in a.l1):
        p.error(f"--l1 values must exceed --l0 = {a.l0:g} m")

    T = a.T0
    rho_vs = physics.rho_vs_sat_kg_per_m3(T)
    D = physics.dv_m2_per_s(T)
    D0 = C.dv0_m2_per_s
    d0 = physics.capillary_length_m(T)
    beta0 = physics.beta_hk_s_per_m(T, a.alpha) * C.rho_ice_kg_per_m3 / rho_vs  # unscaled K&P beta0, as front_speed.py
    rho_i = C.rho_ice_kg_per_m3
    dT = 0.01
    s = (math.log(physics.rho_vs_sat_kg_per_m3(T + dT)) - math.log(physics.rho_vs_sat_kg_per_m3(T - dT))) / (2 * dT)
    K = rho_i / (D * rho_vs)
    print("# inputs")
    print(f"T0 = {T} C, alpha_c = {a.alpha}")
    print(f"rho_vs = {rho_vs:.5e} kg/m3, D_v(T0) = {D:.5e} m2/s, D_v(0 C) = {D0:.5e} m2/s, rho_ice = {rho_i}")
    print(f"d0 = {d0:.5e} m, beta0 = {beta0:.5e} s/m, s = dln(rho_vs)/dT = {s:.5e} 1/K, K = {K:.5e} s/m2")

    print("\n# 1. B1 series-law slope, 1/v = (beta0 + K l)/sigma")
    sl_T = K / a.b1_sigma
    sl_0 = rho_i / (D0 * rho_vs) / a.b1_sigma
    print(f"slope with D_v(T0) = {sl_T:.4e} s/m2 -> measured/predicted = {a.b1_slope / sl_T:.4f}")
    print(f"slope with D_v(0C) = {sl_0:.4e} s/m2 -> measured/predicted = {a.b1_slope / sl_0:.4f}")
    print(f"intercept beta0/sigma = {beta0 / a.b1_sigma:.4e} s/m;"
          f" a 4.2e7 intercept equals a shift in l of {(beta0 / a.b1_sigma - 4.2e7) / sl_T * 1e6:.2f} um")

    print("\n# 2. Ripening pair (2D discs)")
    Rs, Rl, L = a.pair
    dsig = d0 * (1 / Rs - 1 / Rl)
    gap = L - Rs - Rl
    print(f"R_small = {Rs:.2e}, R_large = {Rl:.2e}, center distance {L:.2e}, gap {gap:.2e} m")
    print(f"drive d0 (1/Rs - 1/Rl) = {dsig:.4e}; mean-curvature h for the pair = 1 + {d0 * 0.5 * (1 / Rs + 1 / Rl):.4e}")
    print(f"gradient that cancels it over the center distance: G* = {dsig / (s * L):.3f} K/m"
          f" (temperature difference {dsig / s * 1e3:.4f} mK)")
    v_dif = D * rho_vs * dsig / (rho_i * Rs * math.log(L / Rs))
    v_kin = dsig / beta0
    v = 1 / (1 / v_dif + 1 / v_kin)
    print(f"small-grain shrink speed: diffusion {v_dif:.3e}, kinetics {v_kin:.3e}, series {v:.3e} m/s"
          f" = {v * DAY * 1e6:.3f} um/day; Rs/v = {Rs / v / DAY:.1f} days (rough vanishing time)")

    print("\n# 3. Gradient-driven migration speed v ~ D rho_vs s G / rho_ice")
    for G in a.G:
        vm = D * rho_vs * s * G / rho_i
        print(f"G = {G:7.1f} K/m: v = {vm:.3e} m/s = {vm * DAY * 1e6:8.3f} um/day;"
              f" saturation change over 100 um = {s * G * 1e-4:.3e};"
              f" Pi = beta0 v eps/d0 = {beta0 * vm * a.eps / d0:.3f} at eps {a.eps:.3e} m")

    print(f"\n# 4. Harvest time, front from l0 = {a.l0:.2e} m, ice fill fraction {a.fill}")
    for h in a.harvest:
        sig = 1 - h
        for l1 in a.l1:
            t = (beta0 * (l1 - a.l0) + K * (l1 ** 2 - a.l0 ** 2) / 2) * a.fill / sig
            v_start = sig / (beta0 + K * a.l0)
            print(f"h = {h:5.3f}  l1 = {l1:.2e} m: t = {t / DAY:9.1f} days;"
                  f" starting front speed {v_start:.3e} m/s")
    print("\n# 5. Harvest cutoff on the cold side of the T0 point, dz = ln(1/h_open)/(s G)")
    for h in a.cutoff_h:
        for G in a.G:
            print(f"h_open = {h:6.4f}, G = {G:7.1f} K/m: dz = {math.log(1 / h) / (s * G) * 1e3:9.3f} mm")
    print("\n# 6a. Gradient-capillarity crossover G_x = d0 chi / (s l), chi = 2/W (2D meniscus at theta 60: cos60*2/W = 1/W)")
    for W in a.widths:
        for chi_name, chi in (("1/W", 1 / W), ("2/W", 2 / W)):
            print(f"W = {W * 1e6:6.1f} um, chi = {chi_name}: G_x over l = W: {d0 * chi / (s * W):8.3f} K/m;"
                  f" over l = 1 mm: {d0 * chi / (s * 1e-3):8.4f} K/m")
    print("\n# 6b. Depositing-face stability: G_s = (2 pi/W)^2 d0/s (symmetric), (pi/W)^2 d0/s (tilting)")
    for W in a.widths:
        Gs = (2 * math.pi / W) ** 2 * d0 / s
        Gt = (math.pi / W) ** 2 * d0 / s
        print(f"W = {W * 1e6:6.1f} um: G_s = {Gs:8.2f} K/m, tilting {Gt:8.2f} K/m")
    for G in a.G:
        lam = 2 * math.pi * math.sqrt(d0 / (s * G))
        print(f"G = {G:7.1f} K/m: marginal wavelength 2 pi sqrt(d0/(s G)) = {lam * 1e6:8.1f} um")

    th = math.radians(a.theta)
    print(f"\n# 7a. Meniscus at theta = {a.theta:g} deg: Gibbs-Thomson shift against the gas saturation"
          " difference over the face's x span")
    for W in a.widths:
        R = W / (2 * math.cos(th))
        dx = R * (1 - math.sin(th))
        gt = d0 / R
        print(f"W = {W * 1e6:6.1f} um: R = {R * 1e6:7.1f} um, dx = {dx * 1e6:6.2f} um, d0 chi = {gt:.3e};"
              f" G* = d0 chi/(s dx) = {gt / (s * dx):7.2f} K/m;  s G dx at G 10/30/100: "
              + ", ".join(f"{s * g * dx:.2e}" for g in (10, 30, 100)))

    Lx, Ln, Lsh, Lice, Wd = 1.1e-3, 300e-6, 200e-6, 150e-6, 150e-6
    print(f"\n# 7b. Gas gradient in the type-6 pore (Lx {Lx * 1e6:g}, neck {Ln * 1e6:g}, shoulders"
          f" {Lsh * 1e6:g}, ice {Lice * 1e6:g} um; ice resistance neglected)")
    fac = {}
    for wc in (40e-6, 80e-6, 150e-6):
        r_neck = Wd / wc
        r_sh = Wd / (0.5 * (Wd + wc))
        Lw = Lx - Ln - Lsh - Lice
        f_w = Lx / (Lw + r_sh * Lsh + r_neck * Ln)
        fac[wc] = f_w
        print(f"wc = {wc * 1e6:5.0f} um: G_gas(wide) = {f_w:.3f} G, G_gas(neck) = {f_w * r_neck:.3f} G")

    print("\n# 7c-d. Symmetric-mode growth (W = 150 um, wide gas) and cost to move the face"
          f" {a.move * 1e6:g} um (dtmax {a.dtmax:g} s, {a.cfl_step * 1e6:g} um/step CFL,"
          f" {a.s_per_step:g} s/step)")
    W = 150e-6
    k = 2 * math.pi / W
    pref = D * rho_vs / rho_i
    for G in a.G7:
        for wc in (40e-6, 150e-6):
            Gg = fac[wc] * G
            sig = pref * k * (s * Gg - d0 * k * k)
            v = pref * s * Gg
            t = a.move / v
            dt = min(a.dtmax, a.cfl_step / v)
            n = t / dt
            print(f"G = {G:6.1f}, wc {wc * 1e6:3.0f}: G_gas {Gg:6.1f} K/m, sigma = {sig * DAY:+7.3f} /day,"
                  f" v = {v * DAY * 1e6:6.2f} um/day, t_move = {t / DAY:7.1f} d,"
                  f" e-folds in t_move {sig * t:+6.1f}, dt {dt:6.0f} s, steps {n:8.0f},"
                  f" wall {n * a.s_per_step / 3600:6.1f} h")
    throat_arrest(a, d0, s, beta0, K)
    sloped_pore(a, d0, s, beta0, K)
    optimal_h(a, d0, s, beta0, K)
    return 0


def throat_arrest(a, d0, s, beta0, K):
    """8. Lens arrested in a throat, and the gradient that pulls it out (added 2026-11-13).

    Walls: y = b(x) and W - b(x), one C-infinity bump per wall, b = h exp(1 - 1/(1 - u^2)),
    u = (x - Lx/2)/R (preprocess/build_ice_layout.py), h = (W - wt)/2, R fixed. A face at x
    meets walls that open toward its gas side at the local slope angle omega, so its
    meniscus curvature is chi(x) = 2 cos(theta + omega)/g(x), g = W - 2b (as the type-4 law
    with theta -/+ omega). At G = 0 with vapor reservoirs at both ends the lens drifts
    toward the face of larger chi, so a lens straddling the throat is held at its center
    (a displacement lowers chi on the leading face and raises it on the trailing one).
    With a gradient and local-saturation ends, the ice is nearly isothermal (k_ice/k_air
    115), so the reservoirs differ from the ice by s dT on each side and dT_w + dT_c = G Lx.
    With equal gas-column resistances on both sides (lens centered in the pore; kinetic
    resistances equal), the drift is zero when s G Lx = d0 [chi(x_c) - chi(x_w)], x_w and
    x_c the warm and cold faces. The largest restoring term over all lens positions gives
    G_crit = d0 max[chi(x_c) - chi(x_w)] / (s Lx): below it the lens stays at a displaced
    rest position in the throat, above it the lens is pulled out toward the warm end.
    Approximations: lens length fixed (the ls lens grows), isothermal ice, 1D series
    resistances, Young angle at a moving contact line."""
    th = math.radians(a.theta)
    W, R, L, Lx = a.W8, a.bump_R, a.lens_len, a.Lx8
    print(f"\n# 8. Throat arrest: W {W * 1e6:g} um, bump half-width {R * 1e6:g} um (fixed), lens {L * 1e6:g} um,"
          f" Lx {Lx * 1e6:g} um, theta {a.theta:g} deg")
    print("   G_crit = d0 max[chi(x_c) - chi(x_w)] / (s Lx),  chi(x) = 2 cos(theta + omega)/g")
    x = [Lx / 2 + R * (k / 4000.0 - 1.0) * 0.999 for k in range(8001)]

    def wall(xx, h):
        u = (xx - Lx / 2) / R
        if abs(u) >= 1:
            return 0.0, 0.0
        b = h * math.exp(1 - 1 / (1 - u * u))
        db = b * (-2 * u / (1 - u * u) ** 2) / R
        return b, db

    for wt in a.wt:
        h = 0.5 * (W - wt)

        def chi(xx, h=h):
            b, db = wall(xx, h)
            return 2 * math.cos(th + math.atan(abs(db))) / (W - 2 * b)

        om = max(math.degrees(math.atan(abs(wall(xx, h)[1]))) for xx in x)
        # Lens displaced d toward x = 0 (the warm end). Gas columns l_w = Lx/2 - d - L/2 and
        # l_c = Lx/2 + d - L/2; the ice is isothermal, so the drop G Lx splits in proportion to
        # them. Reservoirs at h_eq rho_vs(T_end), h_eq = 1 - d0 chi of a centered face, so a
        # centered lens neither grows nor shrinks. Face speeds (type-3 series law):
        #   v_w = (s dT_w + d0 (chi_w - chi0)) / (beta0 + K l_w)   deposition, toward x = 0
        #   v_c = (s dT_c - d0 (chi_c - chi0)) / (beta0 + K l_c)   sublimation, toward x = 0
        # center speed v = (v_w + v_c)/2. G_crit = max over d of the G that makes v = 0.
        chi0 = chi(Lx / 2 + L / 2)

        def vel(d, G, chi=chi, chi0=chi0):
            lw, lc = Lx / 2 - d - L / 2, Lx / 2 + d - L / 2
            dTw, dTc = G * Lx * lw / (lw + lc), G * Lx * lc / (lw + lc)
            vw = (s * dTw + d0 * (chi(Lx / 2 - d - L / 2) - chi0)) / (beta0 + K * lw)
            vc = (s * dTc - d0 * (chi(Lx / 2 - d + L / 2) - chi0)) / (beta0 + K * lc)
            return 0.5 * (vw + vc)

        ds = [(R - L / 2) * k / 2000.0 for k in range(1, 2001)]   # warm face still on the bump
        g0 = [(-vel(d, 0.0) / (vel(d, 1.0) - vel(d, 0.0)), d) for d in ds]
        Gc, dpk = max(g0)
        print(f"wt = {wt * 1e6:5.1f} um: bump h {h * 1e6:5.1f} um, steepest wall {om:5.1f} deg"
              f"{' (theta + omega > 90: convex faces)' if a.theta + om >= 90 else ''};"
              f" chi(throat) {chi(Lx / 2):.3e}, chi(wide) {2 * math.cos(th) / W:.3e} 1/m;"
              f" G_crit {Gc:.3f} K/m (equal-column formula {d0 * max(chi(Lx / 2 - d + L / 2) - chi(Lx / 2 - d - L / 2) for d in ds) / (s * Lx):.3f}),"
              f" weakest hold at lens center {dpk * 1e6:5.1f} um off the throat")
        out = []
        for lab, G, d_start, stop in (("G 0: 40 -> 10 um off", 0.0, 40e-6, lambda d: d < 10e-6),
                                      ("4 G_crit: center -> past the hold", 4 * Gc, 0.0, lambda d, p=dpk: d > p),
                                      ("10 G_crit: center -> out of the bump", 10 * Gc, 0.0, lambda d: d > R - L / 2)):
            d, t, dt = d_start, 0.0, 6 * 3600.0
            while not stop(d) and t < 2e4 * DAY:
                d += vel(d, G) * dt
                t += dt
            out.append(f"{lab} {t / DAY:6.0f} d" + (" (>20000)" if t >= 2e4 * DAY else ""))
        v40 = [vel(d, 0.0) * DAY * 1e6 for d in (20e-6, 40e-6, dpk)]
        print(f"      G 0 restoring speed at 20, 40 um and the hold point: {v40[0]:+.4f}, {v40[1]:+.4f},"
              f" {v40[2]:+.4f} um/d; " + "; ".join(out))



def sloped_pore(a, d0, s, beta0, K):
    """8b. Type 5 redesign (user, 13 Nov): the walls slope toward the throat along the whole pore,
    gap g(xi) = sqrt(wt^2 + c^2 xi^2), xi = x - Lx/2, g(+-Lx/2) = W (smooth at the throat, slope
    never zero elsewhere; the throat width changes only the depth of the profile). The lens starts
    on the COLD side (x > Lx/2) and keeps its area (reservoirs at the humidity of the lens at its
    start, so it neither grows nor shrinks there); its length follows from the area. Center speed
    from the two face speeds as in section 8 (series law, isothermal ice, static angle). G pulls
    toward x = 0: at G = 0 the lens should come to rest at the throat; above G_crit it passes
    through and leaves on the warm side."""
    th = math.radians(a.theta)
    W, Lx, L0 = a.W8, a.Lx8, a.lens_len
    print(f"\n# 8b. Whole-pore slope: g = sqrt(wt^2 + c^2 xi^2), W {W * 1e6:g} um at the ends, Lx {Lx * 1e6:g} um,"
          f" lens {L0 * 1e6:g} um started {a.start8 * 1e6:g} um on the cold side, area conserved")
    for wt in a.wt:
        c = math.sqrt(W * W - wt * wt) / (Lx / 2)

        def g(xi, c=c, wt=wt):
            return math.sqrt(wt * wt + (c * xi) ** 2)

        def chi(xi, c=c, wt=wt):
            om = math.atan(0.5 * c * c * abs(xi) / math.sqrt(wt * wt + (c * xi) ** 2))
            return 2 * math.cos(th + om) / math.sqrt(wt * wt + (c * xi) ** 2)

        def area(xm, L, g=g):
            n = 200
            return sum(g(xm - L / 2 + (k + 0.5) * L / n) for k in range(n)) * L / n

        A = area(-a.start8, L0)                       # xm: lens center in xi (negative = warm side)

        def length(xm, A=A, area=area):
            lo, hi = 1e-6, Lx
            for _ in range(60):
                mid = 0.5 * (lo + hi)
                lo, hi = (mid, hi) if area(xm, mid) < A else (lo, mid)
            return 0.5 * (lo + hi)

        xs0 = -a.start8
        chi0 = 0.5 * (chi(xs0 - L0 / 2) + chi(xs0 + L0 / 2))

        def vel(xm, G, chi=chi, chi0=chi0, length=length):   # center speed toward x = 0 (xi decreasing)
            L = length(xm)
            xw, xc = xm - L / 2, xm + L / 2
            lw, lc = Lx / 2 + xw, Lx / 2 - xc
            if lw <= 5e-6 or lc <= 5e-6:
                return None
            dTw, dTc = G * Lx * lw / (lw + lc), G * Lx * lc / (lw + lc)
            vw = (s * dTw + d0 * (chi(xw) - chi0)) / (beta0 + K * lw)
            vc = (s * dTc - d0 * (chi(xc) - chi0)) / (beta0 + K * lc)
            return 0.5 * (vw + vc)

        # G_crit: largest G the warm half of the throat can hold (lens center on the warm side)
        gc = 0.0
        for k in range(1, 400):
            xm = -k * 0.5e-6
            v0, v1 = vel(xm, 0.0), vel(xm, 1.0)
            if v0 is None or v1 is None:
                break
            gc = max(gc, -v0 / (v1 - v0))
        line = []
        for f in (0.0, 0.5, 2.0, 5.0):
            G, xm, t, dt = f * gc, xs0, 0.0, 12 * 3600.0
            v = 0.0
            while t < 6000 * DAY:
                v = vel(xm, G)
                if v is None:
                    break
                xm -= v * dt
                t += dt
            state = (f"left the throat (warm face at the reservoir) at {t / DAY:5.0f} d" if v is None else
                     f"center at {-xm * 1e6:+6.1f} um after 6000 d, moving {v * DAY * 1e6:+.4f} um/d")
            line.append(f"{f:g} G_crit: {state}")
        Ls, Lt = length(xs0), length(0.0)
        print(f"wt {wt * 1e6:5.1f} um: wall slope at the ends {math.degrees(math.atan(0.5 * c * c * (Lx / 2) / W)):4.1f} deg;"
              f" lens {Ls * 1e6:5.1f} um at the start, {Lt * 1e6:5.1f} um centered on the throat;"
              f" G_crit {gc:.3f} K/m; start speed at G 0 {vel(xs0, 0.0) * DAY * 1e6:+.4f} um/d")
        print("      " + "; ".join(line))


def optimal_h(a, d0, s, beta0, K):
    """8c. End humidity that makes the lens grow or sublimate least (user, 13 Nov).

    Pore with no flat ends: gap g(xi) = sqrt(wt^2 + c^2 xi^2), g(+-Lx/2) = W. A lens with faces at
    xi_w < xi_c exchanges vapor with the two end reservoirs only. In area per unit depth, face i
    gains ice at Q_i = Delta_i / Rr_i, with Delta_i the supersaturation the face sees and
    Rr_i = beta0/g_i + K * int_column dx/g(x) (kinetic + diffusive resistance; for a uniform gap
    this is the type-3 law v = Delta/(beta0 + K l)). With ends at h * rho_vs(T_end):
        Delta_w = h - 1 + d0 chi_w + s dT_w,   Delta_c = h - 1 + d0 chi_c - s dT_c,
    so the lens area changes at dA/dt = sum_i (h - h_i)/Rr_i - s (dT_c/Rr_c - dT_w/Rr_w), with the
    face equilibrium h_i = 1 - d0 chi_i. dA/dt = a (h - h0(t)), where a = sum_i 1/Rr_i and
        h0(t) = [sum_i h_i/Rr_i + s (dT_c/Rr_c - dT_w/Rr_w)] / a        (instantaneous zero).
    A fixed h can cancel growth only on average. Candidates:
        h_LS  = int a^2 h0 dt / int a^2 dt   minimizes J(h) = int (dA/dt)^2 dt: J'(h) = 2 int a^2 (h - h0)
                dt = 0 there and J'' = 2 int a^2 dt > 0, so it is the unique minimum;
        h_net = int a h0 dt / int a dt      makes the net area change over the run zero;
        h_chi = 1 - d0 <chi>                the user's estimate: mean face curvature over the run.
    The trajectory depends weakly on h, so h is iterated (trajectory -> h_LS -> trajectory).
    Face speeds: v_w = Q_w/g_w toward x = 0, v_c = -Q_c/g_c toward x = 0; center = mean."""
    th = math.radians(a.theta)
    W, Lx = a.W8, a.Lx8c
    print(f"\n# 8c. Optimal end humidity: no flat ends, Lx {Lx * 1e6:g} um, W {W * 1e6:g} um at the ends,"
          f" lens {a.lens8c * 1e6:g} um started {a.start8c * 1e6:g} um on the cold side; 300 d per case")
    for wt in a.wt:
        c = math.sqrt(W * W - wt * wt) / (Lx / 2)

        def g(xi, c=c, wt=wt):
            return math.sqrt(wt * wt + (c * xi) ** 2)

        def chi_face(xi, side, c=c, g=g):
            """Face curvature 2 cos(theta + omega)/g. side = -1 for the warm face (its gas lies toward
            x = 0), +1 for the cold face. omega > 0 when the gap widens into the face's gas,
            omega < 0 when it narrows (the gas side runs back through the throat)."""
            dg = c * c * xi / g(xi)
            return 2 * math.cos(th + math.atan(0.5 * side * dg)) / g(xi)

        def col(x0, x1, g=g, n=200):                 # int_x0^x1 dx/g
            return sum(1 / g(x0 + (k + 0.5) * (x1 - x0) / n) for k in range(n)) * (x1 - x0) / n

        def area(xm, L, g=g, n=200):
            return sum(g(xm - L / 2 + (k + 0.5) * L / n) for k in range(n)) * L / n

        def length(xm, A, area=area):
            lo, hi = 1e-6, Lx
            for _ in range(50):
                mid = 0.5 * (lo + hi)
                lo, hi = (mid, hi) if area(xm, mid) < A else (lo, mid)
            return 0.5 * (lo + hi)

        def state(xm, A, h, G, g=g, chi_face=chi_face, col=col, length=length):
            L = length(xm, A)
            xw, xc = xm - L / 2, xm + L / 2
            if xw <= -Lx / 2 + 5e-6 or xc >= Lx / 2 - 5e-6:
                return None
            # heat and vapor cross the same gas columns: the drop G Lx splits as int dx/g (ice:
            # k_ice/k_air = 115, so its share is neglected)
            iw, ic = col(-Lx / 2, xw), col(xc, Lx / 2)
            dTw, dTc = G * Lx * iw / (iw + ic), G * Lx * ic / (iw + ic)
            Rw = beta0 / g(xw) + K * iw
            Rc = beta0 / g(xc) + K * ic
            hw, hc = 1 - d0 * chi_face(xw, -1), 1 - d0 * chi_face(xc, +1)
            aa = 1 / Rw + 1 / Rc
            h0 = (hw / Rw + hc / Rc + s * (dTc / Rc - dTw / Rw)) / aa
            Qw = (h - hw + s * dTw) / Rw
            Qc = (h - hc - s * dTc) / Rc
            v = 0.5 * (Qw / g(xw) - Qc / g(xc))       # center speed toward x = 0
            return v, Qw + Qc, aa, h0, 0.5 * (chi_face(xw, -1) + chi_face(xc, +1))

        def run(h, G, xs=a.start8c, days=300, state=state, area=area):
            xm, A = xs, area(xs, a.lens8c)
            A0, t, dt, hist = A, 0.0, 6 * 3600.0, []
            run.maxdev = 0.0
            while t < days * DAY:
                st = state(xm, A, h, G)
                if st is None:
                    break
                v, dA, aa, h0, cm = st
                hist.append((aa, h0, cm))
                xm -= v * dt
                A += dA * dt
                t += dt
                run.maxdev = max(run.maxdev, abs(A / A0 - 1))
            return xm, A / A0 - 1, hist, t

        # G_crit for this pore (lens held on the warm side of the throat), with h at the lens' start value
        hstart = state(a.start8c, area(a.start8c, a.lens8c), 1.0, 0.0)[3]
        gc = 0.0
        for k in range(1, 300):
            xm = -k * 0.5e-6
            A = area(a.start8c, a.lens8c)
            s0, s1 = state(xm, A, hstart, 0.0), state(xm, A, hstart, 1.0)
            if s0 is None or s1 is None:
                break
            gc = max(gc, -s0[0] / (s1[0] - s0[0]))
        print(f"wt {wt * 1e6:4.0f} um: end wall slope {math.degrees(math.atan(0.5 * c * c * (Lx / 2) / W)):4.1f} deg,"
              f" G_crit ~ {gc:.3f} K/m")
        def converge(rule, G, run=run, hstart=hstart):
            """h -> rule(path(h)) until |dh| < 1e-12 (at most 20 passes); returns h and passes."""
            h = hstart
            for k in range(20):
                hist = run(h, G)[2]
                hn = rule(hist)
                if abs(hn - h) < 1e-12:
                    return hn, k + 1
                h = hn
            return h, 20

        def rule_ls(hist):
            return sum(x[0] ** 2 * x[1] for x in hist) / sum(x[0] ** 2 for x in hist)

        def rule_net(hist):
            return sum(x[0] * x[1] for x in hist) / sum(x[0] for x in hist)

        for G in (0.0, 1.0, 2.0):
            h_ls, k_ls = converge(rule_ls, G)
            h_net, k_net = converge(rule_net, G)
            hist = run(h_net, G)[2]
            h_chi = 1 - d0 * sum(x[2] for x in hist) / len(hist)
            res = []
            for lab, hh in (("h_LS", h_ls), ("h_net", h_net), ("h_chi", h_chi), ("h 1", 1.0)):
                xm2, dA2, _, t2 = run(hh, G)
                res.append(f"{lab} {hh:.8f}: area {dA2 * 100:+6.2f} % (max {run.maxdev * 100:.2f} %), "
                           f"center {-xm2 * 1e6:+6.1f} um"
                           + (f" (left at {t2 / DAY:.0f} d)" if t2 < 300 * DAY - 1 else ""))
            print(f"   G {G:4.1f} K/m (passes LS {k_ls}, net {k_net}): " + "; ".join(res))


if __name__ == "__main__":
    raise SystemExit(main())

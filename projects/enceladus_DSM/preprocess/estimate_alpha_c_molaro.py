"""
estimate_alpha_c_molaro.py — derive the condensation coefficient alpha_c from
Molaro et al. (2019) Fig. 11 instead of assuming it.

WHY
---
alpha_c is the only genuinely free parameter in the sublimation model, and it is
uncertain by orders of magnitude in the literature. But Molaro measured a neck
growth rate, and neck growth is driven by a curvature difference we can compute
exactly. That turns alpha_c into something the data constrains directly.

THE MATH
--------
Gibbs-Thomson gives the equilibrium vapour density over a curved ice surface

    rho_eq(K) = rho_vs * (1 + d0*K),       d0 = gamma*V_m/(R_gas*T)

where K is the SUM of the principal curvatures, positive for convex ice
(Kaempfer & Plapp 2009 Eq. 11). For the two surfaces that exchange mass during
neck growth:

    grain surface (sphere)   K_g = 2/R
    neck (saddle)            K_n = 1/r - 1/rho

with r the neck radius and rho the fillet radius. rho follows from requiring the
fillet circle to be tangent to both spheres:

    sqrt(R^2 + (r+rho)^2) = R + rho   ->   rho = r^2 / (2*(R - r))

Because rho << r, K_n is large and negative: the neck is a strong vapour sink.
The driving force, in the dimensionless concentration u = drho_v/rho_i used by
the phase-field model, is

    du = d0' * dK,     d0' = d0 * rho_vs/rho_i,     dK = K_g - K_n

Vapour must cross two resistances in series to get from the grain surface to the
neck: attachment/detachment at the interface, and diffusion through the pore.

    du = (beta' + l/D_v) * v_n,     beta' = (1/alpha_c)*sqrt(2*pi*m/(k_B*T))

with l the diffusion path length and v_n = dr/dt the measured neck growth rate.

TWO LIMITS, BOTH COMPUTED
-------------------------
LIMIT A -- attachment-only (D_v -> infinity). Drops l entirely, so it needs no
geometric assumption, and gives a LOWER BOUND on alpha_c: any real transport
resistance means more attachment speed is needed, not less.

    alpha_c_min = sqrt(2*pi*m/(k_B*T)) * v_n / (d0' * dK)

LIMIT B -- the vapour-transport ceiling at alpha_c = 1 (beta' at its physical
minimum). If v_obs exceeds this, vapour transport alone cannot explain the
observed rate at ANY alpha_c.

    v_max = d0' * dK / (beta'_min + l/D_v),      l ~ rho

SATURATION -- v_n(alpha_c) over the literature range, showing how little headroom
is left above alpha_c = 0.1. Once Lambda = beta'*D_v/rho drops below 1 the
process is transport-limited and further attachment speed buys nothing.

HOW TO READ THE ANSWER
----------------------
alpha_c_min ~ 0.14 is a lower bound on what the MODEL needs, and therefore an
UPPER bound on the physical attachment coefficient: surface diffusion (which this
model does not have) would supply part of the flux and leave less for vapour.
The residual gap is the same ~50% vapour share that
docs/molaro_validation_synthesis.md section 4 arrives at independently.

Usage
-----
  python preprocess/estimate_alpha_c_molaro.py
  python preprocess/estimate_alpha_c_molaro.py --data inputs/validation/molaro2019_fig11_T-5.csv \
      --T0 -5 --R1 98.5e-6 --R2 76.5e-6 --outdir studies/molaro_2019
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))

# Single source of truth for the thermodynamics — do not re-derive here.
from comp_eps import (  # noqa: E402
    beta_HK, Dv_T, rho_vs_sat, capillary_length, _RHO_ICE, _M_H2O, _K_B,
)

_REPO = _HERE.parent

# Validated categorical pair (repo convention; passes the OKLab/CVD checks:
# protan 19.7, deutan 24.1, normal 28.5, both in the light lightness band).
C_MODEL, C_REF, C_INK, C_MUTED = "#3d74d9", "#d1495b", "#3f434a", "#8a8a8a"


# =========================================================================
# Geometry
# =========================================================================

def fillet_radius(r: float, R: float) -> float:
    """Fillet radius rho of the neck between two spheres of radius R.

    Exact tangency solution rho = r^2/(2(R-r)); the small-r form r^2/(2R) is
    what the resolution-floor convention sqrt(12*eps/R) is built on.
    """
    return r * r / (2.0 * (R - r))


def curvature_difference(r: float, R: float) -> float:
    """dK = K_grain - K_neck = 2/R - (1/r - 1/rho)  [1/m]. Positive: neck fills."""
    return 2.0 / R - (1.0 / r - 1.0 / fillet_radius(r, R))


# =========================================================================
# The two limits
# =========================================================================

def alpha_c_min(v_n: float, r: float, R: float, T_C: float) -> float:
    """LIMIT A: smallest alpha_c that can produce v_n with D_v -> infinity."""
    d0p = capillary_length(T_C) * rho_vs_sat(T_C) / _RHO_ICE
    c_th = math.sqrt(2.0 * math.pi * _M_H2O / (_K_B * (T_C + 273.15)))
    return c_th * v_n / (d0p * curvature_difference(r, R))


def v_n_of_alpha(alpha_c: float, r: float, R: float, T_C: float,
                 ell: float | None = None) -> float:
    """Neck growth rate with both resistances in series. ell defaults to rho."""
    d0p = capillary_length(T_C) * rho_vs_sat(T_C) / _RHO_ICE
    rho = fillet_radius(r, R)
    if ell is None:
        ell = rho
    return d0p * curvature_difference(r, R) / (beta_HK(T_C, alpha_c) + ell / Dv_T(T_C))


def lambda_kinetic(alpha_c: float, r: float, R: float, T_C: float) -> float:
    """Lambda = L*/rho with L* = beta_HK*D_v. >>1 attachment-limited, <<1 transport."""
    return beta_HK(T_C, alpha_c) * Dv_T(T_C) / fillet_radius(r, R)


# =========================================================================
# Data
# =========================================================================

def read_fig11(path: Path):
    """Return [(t_s, neck_width_m, lg_diam_m, sm_diam_m), ...] from the validation CSV."""
    rows = []
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            f = next(csv.reader([line]))
            rows.append((float(f[0]) * 60.0, float(f[1]) * 1e-6,
                         float(f[4]) * 1e-6, float(f[5]) * 1e-6))
    return rows


# =========================================================================
# Report
# =========================================================================

def analyse(rows, R: float, T_C: float):
    """Per-interval and integrated estimates. Returns (intervals, integrated)."""
    c_th = math.sqrt(2.0 * math.pi * _M_H2O / (_K_B * (T_C + 273.15)))
    out = []
    for (t0, w0, *_), (t1, w1, *_) in zip(rows, rows[1:]):
        r = 0.5 * (w0 + w1) / 2.0                 # mid-interval neck RADIUS
        v = (w1 - w0) / 2.0 / (t1 - t0)           # dr/dt
        rho = fillet_radius(r, R)
        vmax = v_n_of_alpha(1.0, r, R, T_C)
        out.append(dict(t0_s=t0, t1_s=t1, r_m=r, rho_m=rho,
                        dK_per_m=curvature_difference(r, R), v_obs=v,
                        alpha_c_min=alpha_c_min(v, r, R, T_C),
                        v_max_alpha1=vmax, v_over_vmax=v / vmax))
    # Integrated: first to last point.
    t0, w0 = rows[0][0], rows[0][1]
    t1, w1 = rows[-1][0], rows[-1][1]
    r = 0.5 * (w0 + w1) / 2.0
    v = (w1 - w0) / 2.0 / (t1 - t0)
    integ = dict(r_m=r, rho_m=fillet_radius(r, R),
                 dK_per_m=curvature_difference(r, R), v_obs=v,
                 alpha_c_min=alpha_c_min(v, r, R, T_C),
                 v_max_alpha1=v_n_of_alpha(1.0, r, R, T_C),
                 c_th=c_th)
    integ["v_over_vmax"] = v / integ["v_max_alpha1"]
    return out, integ


def saturation_sweep(r: float, R: float, T_C: float, alphas):
    ref = v_n_of_alpha(0.1, r, R, T_C)
    ceil = v_n_of_alpha(1.0, r, R, T_C)
    return [dict(alpha_c=a, beta_hk=beta_HK(T_C, a), v_n=v_n_of_alpha(a, r, R, T_C),
                 rel_to_0p1=v_n_of_alpha(a, r, R, T_C) / ref,
                 frac_of_ceiling=v_n_of_alpha(a, r, R, T_C) / ceil,
                 Lambda=lambda_kinetic(a, r, R, T_C)) for a in alphas]


def make_figure(intervals, integ, sweep, T_C, out_png):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11.0, 4.3))

    # --- left: alpha_c_min per interval, vs neck radius -------------------
    r_um = [d["r_m"] * 1e6 for d in intervals]
    a_min = [d["alpha_c_min"] for d in intervals]
    axL.axhspan(1.0, 10.0, color=C_MUTED, alpha=0.18, lw=0, zorder=0)
    axL.text(r_um[0], 1.35, "$\\alpha_c > 1$ is unphysical", fontsize=8,
             color=C_INK, va="bottom")
    axL.plot(r_um, a_min, "o-", color=C_MODEL, lw=2, ms=8, zorder=3,
             label="per interval")
    axL.axhline(integ["alpha_c_min"], color=C_REF, lw=2, ls="--", zorder=2,
                label=f"integrated = {integ['alpha_c_min']:.3f}")
    axL.set_yscale("log")
    axL.set_ylim(0.02, 10)
    axL.set_xlabel("neck radius r  [$\\mu$m]")
    axL.set_ylabel("$\\alpha_c$ lower bound  [–]")
    axL.set_title(f"Minimum $\\alpha_c$ from the observed rate  (T = {T_C:g} $\\degree$C)",
                  fontsize=10, color=C_INK)
    axL.grid(alpha=0.25, lw=0.6)
    axL.set_axisbelow(True)
    axL.legend(frameon=False, fontsize=9, loc="upper left")

    # --- right: v_n(alpha_c) saturation -----------------------------------
    al = [s["alpha_c"] for s in sweep]
    vn = [s["v_n"] for s in sweep]
    axR.plot(al, vn, "o-", color=C_MODEL, lw=2, ms=8, zorder=3, label="model $v_n(\\alpha_c)$")
    axR.axhline(integ["v_obs"], color=C_REF, lw=2, ls="--", zorder=2,
                label=f"observed = {integ['v_obs']:.2e} m/s")
    axR.set_xscale("log")
    axR.set_yscale("log")
    axR.set_xlabel("$\\alpha_c$  [–]")
    axR.set_ylabel("neck growth rate $v_n$  [m/s]")
    axR.set_title("$\\alpha_c$ saturates: transport takes over", fontsize=10, color=C_INK)
    axR.grid(alpha=0.25, lw=0.6)
    axR.set_axisbelow(True)
    axR.legend(frameon=False, fontsize=9, loc="lower right")
    # Direct-label the two endpoints only (never every point).
    for s, ha, dx in ((sweep[0], "left", 6), (sweep[-1], "right", -6)):
        axR.annotate(f"$\\Lambda$ = {s['Lambda']:.2f}", (s["alpha_c"], s["v_n"]),
                     textcoords="offset points", xytext=(dx, 10),
                     fontsize=8, color=C_INK, ha=ha)

    for ax in (axL, axR):
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.tick_params(colors=C_INK, labelsize=9)

    fig.tight_layout()
    fig.savefig(out_png, dpi=160)
    print(f"  wrote {out_png}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", type=Path,
                    default=_REPO / "inputs/validation/molaro2019_fig11_T-20.csv")
    ap.add_argument("--T0", type=float, default=-20.0, help="temperature [C]")
    ap.add_argument("--R1", type=float, default=101.0e-6, help="large grain radius [m]")
    ap.add_argument("--R2", type=float, default=72.5e-6, help="small grain radius [m]")
    ap.add_argument("--outdir", type=Path, default=_REPO / "studies/molaro_2019")
    args = ap.parse_args()

    R = 0.5 * (args.R1 + args.R2)
    T = args.T0
    rows = read_fig11(args.data)
    intervals, integ = analyse(rows, R, T)
    sweep = saturation_sweep(integ["r_m"], R, T,
                             [1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1, 1.0])

    print("=" * 78)
    print(f"  alpha_c from Molaro et al. (2019), T = {T:g} C")
    print("=" * 78)
    print(f"  data          : {args.data}")
    print(f"  R_ave         : {R*1e6:.2f} um   (R1 = {args.R1*1e6:.1f}, R2 = {args.R2*1e6:.1f})")
    print(f"  d0            : {capillary_length(T):.4e} m")
    print(f"  d0'           : {capillary_length(T)*rho_vs_sat(T)/_RHO_ICE:.4e} m")
    print(f"  D_v           : {Dv_T(T):.4e} m^2/s")
    print(f"  sqrt(2 pi m/kT): {integ['c_th']:.4e} s/m   (= beta_HK * alpha_c)")

    print("\n--- LIMIT A: alpha_c_min = sqrt(2 pi m/kT) * v_n / (d0' * dK) ---")
    print("  (attachment-only, D_v -> inf; no geometric assumption)")
    print("   interval[min]   r[um]  rho[um]     dK[1/m]   v_obs[m/s]  alpha_c_min  v_obs/v_max(a=1)")
    for d in intervals:
        print(f"   {d['t0_s']/60:5.0f}->{d['t1_s']/60:4.0f}  {d['r_m']*1e6:7.2f} {d['rho_m']*1e6:7.3f}"
              f" {d['dK_per_m']:11.3e}  {d['v_obs']:.3e}    {d['alpha_c_min']:8.4f}      {d['v_over_vmax']:7.2f}")
    print(f"\n   INTEGRATED  r={integ['r_m']*1e6:.2f} um  rho={integ['rho_m']*1e6:.3f} um"
          f"  v_obs={integ['v_obs']:.3e} m/s")
    print(f"   -> alpha_c_min = {integ['alpha_c_min']:.4f}")

    print("\n--- LIMIT B: vapour-transport ceiling at alpha_c = 1 (l = rho) ---")
    print(f"   v_max = {integ['v_max_alpha1']:.3e} m/s   ->   v_obs / v_max = {integ['v_over_vmax']:.2f}")
    if integ["v_over_vmax"] > 1.0:
        print("   v_obs EXCEEDS the ceiling: vapour transport alone cannot explain the")
        print("   observed rate at ANY alpha_c. The l = rho path length is crude (it")
        print("   overestimates the resistance by ~2.5x against measured runs), but the")
        print("   sign of the conclusion is robust.")

    print("\n--- SATURATION: v_n(alpha_c), Lambda = beta_HK*D_v/rho ---")
    print("   alpha_c   beta_HK[s/m]     v_n[m/s]   rel. to 0.1   frac of a=1   Lambda")
    for s in sweep:
        print(f"   {s['alpha_c']:7.3g} {s['beta_hk']:12.4e}  {s['v_n']:.4e}      {s['rel_to_0p1']:6.3f}"
              f"        {s['frac_of_ceiling']:6.3f}   {s['Lambda']:8.3f}")
    head = [s for s in sweep if s["alpha_c"] == 1.0][0]["rel_to_0p1"]
    print(f"\n   Headroom from alpha_c = 0.1 to the alpha_c = 1 ceiling: {head:.2f}x")

    print("\n--- HOW TO READ THIS ---")
    print(f"   alpha_c_min = {integ['alpha_c_min']:.3f} is a LOWER bound on what the model needs,")
    print("   and therefore an UPPER bound on the physical attachment coefficient:")
    print("   surface diffusion (absent from this model) would carry part of the flux")
    print("   and leave less for vapour. Combined with the saturation above, alpha_c is")
    print("   exhausted as a knob near 0.1 -- see docs/molaro_validation_synthesis.md sec. 4.")
    print("=" * 78)

    args.outdir.mkdir(parents=True, exist_ok=True)
    csv_path = args.outdir / "alpha_c_estimate.csv"
    with csv_path.open("w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["kind", "t0_s", "t1_s", "r_m", "rho_m", "dK_per_m", "v_obs_m_s",
                    "alpha_c_min", "v_max_alpha1_m_s", "v_obs_over_v_max"])
        for d in intervals:
            w.writerow(["interval", f"{d['t0_s']:.1f}", f"{d['t1_s']:.1f}",
                        f"{d['r_m']:.6e}", f"{d['rho_m']:.6e}", f"{d['dK_per_m']:.6e}",
                        f"{d['v_obs']:.6e}", f"{d['alpha_c_min']:.6f}",
                        f"{d['v_max_alpha1']:.6e}", f"{d['v_over_vmax']:.4f}"])
        w.writerow(["integrated", f"{rows[0][0]:.1f}", f"{rows[-1][0]:.1f}",
                    f"{integ['r_m']:.6e}", f"{integ['rho_m']:.6e}", f"{integ['dK_per_m']:.6e}",
                    f"{integ['v_obs']:.6e}", f"{integ['alpha_c_min']:.6f}",
                    f"{integ['v_max_alpha1']:.6e}", f"{integ['v_over_vmax']:.4f}"])
        w.writerow([])
        w.writerow(["saturation", "alpha_c", "beta_HK_s_m", "v_n_m_s", "rel_to_0p1",
                    "frac_of_ceiling", "Lambda"])
        for s in sweep:
            w.writerow(["saturation", f"{s['alpha_c']:.6g}", f"{s['beta_hk']:.6e}",
                        f"{s['v_n']:.6e}", f"{s['rel_to_0p1']:.4f}",
                        f"{s['frac_of_ceiling']:.4f}", f"{s['Lambda']:.4f}"])
    print(f"  wrote {csv_path}")
    make_figure(intervals, integ, sweep, T, args.outdir / "alpha_c_estimate.png")


if __name__ == "__main__":
    main()

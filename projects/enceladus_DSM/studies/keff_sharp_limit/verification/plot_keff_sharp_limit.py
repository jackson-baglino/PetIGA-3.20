#!/usr/bin/env python3
"""Gate and plot the k_eff finite-eps ladder against the closed-form prediction.

Reads keff_sharp_limit.csv (written by verify_keff_sharp_limit.sh), checks the
three predictions of
  projects/effective_thermal_cond/docs/calonne_to_phasefield_equivalence.tex
and writes the two-panel figure. Exits non-zero if any gate fails, so the driver
can be used as a regression check.

THE FIGURE IS IN RESISTIVITY ON PURPOSE
---------------------------------------
Panel (a) plots 1/k_11 against eps, not k_11. The prediction

    1/k_perp(eps) = <1/K>_sharp - (n_Gamma/L) * C * eps

is LINEAR in eps in resistivity and a hyperbola in conductivity. Plotted as a
line, the test has two independently predicted numbers -- intercept and slope,
neither fitted -- and a reader can see agreement or its absence. Plotted as
k_perp it is a curve whose agreement can only be eyeballed, which is how a 10%
slope error passes review.

Panel (b) plots k_00 against eps on an axis scaled to the tolerance, because the
prediction is that NOTHING happens: the tangential excess is exactly zero. On an
auto-scaled axis a flat line's numerical noise fills the frame and looks like
structure, so the band being compared against is drawn explicitly.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "preprocess"))
import figstyle as fs                                     # noqa: E402
import keff_laminate_analytic as ana                      # noqa: E402

# Gate tolerances. The parallel one is tight because the prediction is exact;
# it is set by the IC's quadrature error (~4e-4 relative at Nx=256, per the
# geometry file header), not by anything physical.
TOL_PARALLEL = 2.0e-3      # relative, k_00 vs phi*K_i + (1-phi)*K_a
TOL_PERP = 2.0e-2          # relative, k_11 vs the closed form, per rung
TOL_SLOPE = 5.0e-2         # relative, fitted vs predicted ladder slope
TOL_OFFDIAG = 1.0e-3       # |k_01| / k_00


def read_ladder(path: Path) -> list[dict]:
    with path.open() as fh:
        rows = [{k: v for k, v in r.items()} for r in csv.DictReader(fh)]
    if not rows:
        raise SystemExit(f"no rows in {path} -- did the runs produce k_eff.csv?")
    for r in rows:
        for k in ("k_00", "k_01", "k_10", "k_11", "phi_bar", "eps"):
            r[k] = float(r[k])
    return sorted(rows, key=lambda r: r["eps"])


def gate(name: str, ok: bool, detail: str, failures: list[str]) -> None:
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {detail}")
    if not ok:
        failures.append(name)


def check(rows: list[dict], L: float, n_gamma: int) -> list[str]:
    failures: list[str] = []
    eps = np.array([r["eps"] for r in rows])
    phi = np.array([r["phi_bar"] for r in rows])
    k00 = np.array([r["k_00"] for r in rows])
    k11 = np.array([r["k_11"] for r in rows])
    offd = np.maximum(np.abs([r["k_01"] for r in rows]),
                      np.abs([r["k_10"] for r in rows]))

    print("\nmeasured ladder (predictions evaluated at the MEASURED phi_bar):")
    print(f"  {'eps':>11} {'phi_bar':>9} {'k_00':>10} {'k_00 pred':>10} "
          f"{'k_11':>10} {'k_11 pred':>10} {'k11 err':>9}")
    for e, p, a, b in zip(eps, phi, k00, k11):
        kp = ana.k_parallel(p)
        kq = ana.k_perp(e, p, L, n_gamma)
        print(f"  {e:11.4e} {p:9.6f} {a:10.6f} {kp:10.6f} "
              f"{b:10.6f} {kq:10.6f} {abs(b - kq) / kq:8.2%}")

    # -- Gate 1: tangential excess is exactly zero, so k_00 must not move.
    print("\nGate 1  parallel component is eps-exact (Sigma_t = 0)")
    pred00 = np.array([ana.k_parallel(p) for p in phi])
    err00 = np.max(np.abs(k00 - pred00) / pred00)
    gate("k_00 vs arithmetic mean", err00 < TOL_PARALLEL,
         f"max relative error {err00:.2e} (tol {TOL_PARALLEL:.0e})", failures)
    spread = (k00.max() - k00.min()) / k00.mean()
    gate("k_00 flat across ladder", spread < TOL_PARALLEL,
         f"spread {spread:.2e} (tol {TOL_PARALLEL:.0e})", failures)

    # -- Gate 2: the normal excess, which is what the theory actually predicts.
    print("\nGate 2  perpendicular component follows the closed form (Sigma_n)")
    pred11 = np.array([ana.k_perp(e, p, L, n_gamma) for e, p in zip(eps, phi)])
    err11 = np.max(np.abs(k11 - pred11) / pred11)
    gate("k_11 vs closed form", err11 < TOL_PERP,
         f"max relative error {err11:.2e} (tol {TOL_PERP:.0e})", failures)

    icept_p, slope_p = ana.resistivity_line(float(phi.mean()), L, n_gamma)
    if len(eps) >= 2:
        slope_f, icept_f = np.polyfit(eps, 1.0 / k11, 1)
        slope_f = -slope_f          # fit is 1/k = a*eps + b with a = -slope
        rel_s = abs(slope_f - slope_p) / slope_p
        rel_i = abs(icept_f - icept_p) / icept_p
        gate("ladder slope", rel_s < TOL_SLOPE,
             f"fitted {slope_f:.4e} vs predicted {slope_p:.4e} "
             f"({rel_s:.2%}, tol {TOL_SLOPE:.0%})", failures)
        gate("ladder intercept", rel_i < TOL_PERP,
             f"fitted {icept_f:.4f} vs predicted {icept_p:.4f} "
             f"({rel_i:.2%}, tol {TOL_PERP:.0%})", failures)
        if rel_s > TOL_SLOPE and rel_i < TOL_PERP:
            print("    NOTE  a bad slope with a good intercept is the signature "
                  "of resolution\n          NOT scaling with eps -- check "
                  "EPS_PER_ELEM in the driver before\n          doubting the "
                  "theory.")

    # -- Gate 3: the cell solver's own isotropy check.
    print("\nGate 3  off-diagonal components vanish")
    rel_off = np.max(offd / k00)
    gate("|k_01|, |k_10| ~ 0", rel_off < TOL_OFFDIAG,
         f"max |k_0i|/k_00 = {rel_off:.2e} (tol {TOL_OFFDIAG:.0e})", failures)

    return failures


def figure(rows: list[dict], L: float, n_gamma: int, out: Path) -> None:
    eps = np.array([r["eps"] for r in rows])
    phi = np.array([r["phi_bar"] for r in rows])
    k00 = np.array([r["k_00"] for r in rows])
    k11 = np.array([r["k_11"] for r in rows])
    phi_m = float(phi.mean())
    icept, slope = ana.resistivity_line(phi_m, L, n_gamma)

    fig, (ax_r, ax_p) = plt.subplots(1, 2, figsize=(10.0, 4.2))

    # -- (a) resistivity: the predicted LINE, and the measured points on it ---
    e_fine = np.linspace(0.0, eps.max() * 1.08, 200)
    ax_r.plot(e_fine * 1e6, icept - slope * e_fine, color=fs.C[0], lw=2.0,
              zorder=2, label="predicted: no fitted parameters")
    ax_r.plot(eps * 1e6, 1.0 / k11, "o", color=fs.C[1], ms=7, zorder=3,
              markeredgecolor="white", markeredgewidth=1.0, label="measured")
    # clip_on=False: the star sits exactly on x = 0, and the axes would
    # otherwise slice it in half against the spine.
    ax_r.plot([0], [icept], "*", color=fs.INK, ms=13, zorder=4, clip_on=False,
              label=r"sharp limit $\langle 1/K\rangle$")
    fs.style(ax_r, r"interface decay length  $\varepsilon$  [$\mu$m]",
             r"resistivity  $1/k_{\perp}$  [m K W$^{-1}$]",
             "(a)  normal: the surface excess", logy=False)
    ax_r.set_xlim(left=0.0)
    ax_r.legend(fontsize=fs.FS_LEG, frameon=False, loc="lower left")
    ax_r.text(0.97, 0.94,
              f"slope  {slope:.3e}\nintercept  {icept:.4f}",
              transform=ax_r.transAxes, ha="right", va="top",
              fontsize=fs.FS_NOTE, color=fs.MUTED)

    # -- (b) the parallel component, on an axis scaled to the tolerance -------
    pred00 = ana.k_parallel(phi_m)
    ax_p.axhspan(pred00 * (1 - TOL_PARALLEL), pred00 * (1 + TOL_PARALLEL),
                 color=fs.C[2], alpha=0.16, zorder=1,
                 label=f"gate tolerance  $\\pm${TOL_PARALLEL:.0e}")
    ax_p.axhline(pred00, color=fs.C[2], lw=1.6, zorder=2,
                 label=r"predicted  $\bar\phi K_i + (1-\bar\phi)K_a$")
    ax_p.plot(eps * 1e6, k00, "o", color=fs.C[1], ms=7, zorder=3,
              markeredgecolor="white", markeredgewidth=1.0, label="measured")
    fs.style(ax_p, r"interface decay length  $\varepsilon$  [$\mu$m]",
             r"$k_{\parallel}$  [W m$^{-1}$K$^{-1}$]",
             "(b)  tangential: no excess at all", logy=False)
    ax_p.set_xlim(left=0.0)
    ax_p.set_ylim(pred00 * (1 - 4 * TOL_PARALLEL),
                  pred00 * (1 + 4 * TOL_PARALLEL))
    ax_p.legend(fontsize=fs.FS_LEG, frameon=False, loc="lower left")

    fig.tight_layout()
    fs.save(fig, out, "keff_sharp_limit", dpi=200)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--csv", type=Path,
                    default=Path(__file__).with_name("keff_sharp_limit.csv"))
    ap.add_argument("--L", type=float, default=1.0e-3, help="cell period [m]")
    ap.add_argument("--n-gamma", type=int, default=2,
                    help="interfaces per period (2 under -periodic 1)")
    args = ap.parse_args()

    rows = read_ladder(args.csv)
    print(ana.phi_bar_note())
    failures = check(rows, args.L, args.n_gamma)
    figure(rows, args.L, args.n_gamma, args.csv.parent)

    print()
    if failures:
        print(f"FAILED {len(failures)} gate(s): {', '.join(failures)}")
        return 1
    print("all gates passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

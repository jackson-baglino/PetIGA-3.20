#!/usr/bin/env python3
"""Gate and plot the disk-array eps ladder, arithmetic vs tensor law.

Reads keff_disk.csv (written by verify_keff_disk.sh). Exits non-zero if any gate
fails, so the driver doubles as a regression check.

The compared quantity is k = (k_00 + k_11)/2 against the sharp Rayleigh value
at the NOMINAL area fraction pi R^2/L^2. Unlike the laminate, the measured
phi_bar is the wrong argument here: on a curved interface the diffuse disk
genuinely holds O(eps^2) more ice than the sharp one, and that excess is part
of the bias being measured, not a discretisation error to be divided out.

GATES
  1  square symmetry, both laws: k_00 = k_11, k_01 = k_10 = 0
  2  arith: quadratic fit k(eps) = a + b eps + c eps^2 has a = k_sharp and
     b = the predicted first-order slope (dipole field, no fitted constants)
  3  tensor: at every rung its error is a small fraction of arith's, and at the
     finest rung it is small outright
  4  tensor: fitted a = k_sharp and fitted b ~ 0 (first-order bias removed)
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
import keff_disk_analytic as ana                          # noqa: E402

TOL_SYM = 1.0e-3           # |k00-k11|/k00 and |k01|/k00
TOL_ICEPT = 5.0e-3         # fitted intercept vs Rayleigh, relative
TOL_SLOPE_ARITH = 0.15     # fitted vs dipole slope; dipole field is O(f^4)-approximate
TOL_SLOPE_TENSOR = 0.10    # |tensor slope| as a fraction of the arith prediction
RATIO_TENSOR = 0.20        # |err_tensor| / |err_arith| per rung
TOL_FINEST = 5.0e-3        # |err_tensor| at the finest rung


def read(path: Path) -> dict[str, dict[str, np.ndarray]]:
    with path.open() as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        raise SystemExit(f"no rows in {path}")
    out: dict[str, dict[str, np.ndarray]] = {}
    for law in ("arith", "tensor"):
        rs = sorted((r for r in rows if r["interp"] == law),
                    key=lambda r: float(r["eps"]))
        if not rs:
            raise SystemExit(f"no {law} rows in {path}")
        out[law] = {k: np.array([float(r[k]) for r in rs])
                    for k in ("eps", "k_00", "k_01", "k_10", "k_11", "phi_bar")}
        out[law]["k"] = 0.5 * (out[law]["k_00"] + out[law]["k_11"])
    return out


def gate(name: str, ok: bool, detail: str, failures: list[str]) -> None:
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {detail}")
    if not ok:
        failures.append(name)


def fit(eps: np.ndarray, k: np.ndarray) -> tuple[float, float]:
    """Intercept and linear coefficient of k = a + b eps (+ c eps^2 if >= 3 rungs)."""
    deg = 2 if len(eps) >= 3 else 1
    coef = np.polyfit(eps, k, deg)
    return float(coef[-1]), float(coef[-2])


def check(d: dict, L: float, R: float) -> list[str]:
    failures: list[str] = []
    ks = ana.k_sharp(ana.area_fraction(R, L))
    s_pred = ana.arith_slope(R, L)
    a, t = d["arith"], d["tensor"]
    if not np.allclose(a["eps"], t["eps"]):
        raise SystemExit("arith and tensor ladders do not share eps rungs")

    err_a = a["k"] / ks - 1.0
    err_t = t["k"] / ks - 1.0
    print(f"\nsharp reference k = {ks:.7f}   predicted arith slope {s_pred:.2f}")
    print(f"  {'eps':>10} {'k arith':>10} {'err':>8} {'pred':>8} "
          f"{'k tensor':>10} {'err':>9}")
    for e, ka, ea, kt, et in zip(a["eps"], a["k"], err_a, t["k"], err_t):
        print(f"  {e:10.3e} {ka:10.7f} {ea:+8.2%} {s_pred * e / ks:+8.2%} "
              f"{kt:10.7f} {et:+9.3%}")

    print("\nGate 1  square symmetry")
    for law, x in d.items():
        asym = np.max(np.abs(x["k_00"] - x["k_11"]) / x["k_00"])
        offd = np.max(np.maximum(np.abs(x["k_01"]), np.abs(x["k_10"])) / x["k_00"])
        gate(f"{law}: k_00 = k_11", asym < TOL_SYM,
             f"max rel diff {asym:.1e} (tol {TOL_SYM:.0e})", failures)
        gate(f"{law}: k_01, k_10 ~ 0", offd < TOL_SYM,
             f"max |k_0i|/k_00 {offd:.1e} (tol {TOL_SYM:.0e})", failures)

    print("\nGate 2  arith law: first-order bias as predicted")
    icept_a, slope_a = fit(a["eps"], a["k"])
    gate("arith intercept", abs(icept_a / ks - 1) < TOL_ICEPT,
         f"{icept_a:.7f} vs {ks:.7f} ({icept_a / ks - 1:+.2%}, tol {TOL_ICEPT:.1%})",
         failures)
    gate("arith slope", abs(slope_a / s_pred - 1) < TOL_SLOPE_ARITH,
         f"{slope_a:.2f} vs predicted {s_pred:.2f} "
         f"({slope_a / s_pred - 1:+.1%}, tol {TOL_SLOPE_ARITH:.0%})", failures)

    print("\nGate 3  tensor law: error small against arith, rung by rung")
    ratio = np.max(np.abs(err_t) / np.abs(err_a))
    gate("|err tensor| / |err arith|", ratio < RATIO_TENSOR,
         f"worst {ratio:.3f} (tol {RATIO_TENSOR})", failures)
    gate("tensor error at finest eps", abs(err_t[0]) < TOL_FINEST,
         f"{err_t[0]:+.3%} at eps = {t['eps'][0]:.3e} (tol {TOL_FINEST:.1%})",
         failures)

    print("\nGate 4  tensor law: no first-order term")
    icept_t, slope_t = fit(t["eps"], t["k"])
    gate("tensor intercept", abs(icept_t / ks - 1) < TOL_ICEPT,
         f"{icept_t:.7f} vs {ks:.7f} ({icept_t / ks - 1:+.2%}, tol {TOL_ICEPT:.1%})",
         failures)
    gate("tensor slope", abs(slope_t) < TOL_SLOPE_TENSOR * s_pred,
         f"{slope_t:.2f} = {slope_t / s_pred:+.1%} of the arith slope "
         f"(tol {TOL_SLOPE_TENSOR:.0%})", failures)
    return failures


def figure(d: dict, L: float, R: float, out: Path) -> None:
    ks = ana.k_sharp(ana.area_fraction(R, L))
    s_pred = ana.arith_slope(R, L)
    a, t = d["arith"], d["tensor"]
    e_um = a["eps"] * 1e6
    e_fine = np.geomspace(a["eps"].min() * 0.8, a["eps"].max() * 1.25, 100)

    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    ax.plot(e_fine * 1e6, s_pred * e_fine / ks, color=fs.C[0], lw=1.6,
            label="arith, predicted first order")
    ax.plot(e_um, np.abs(a["k"] / ks - 1), "o", color=fs.C[0], ms=7,
            markeredgecolor="white", label="arith, measured")
    ax.plot(e_um, np.abs(t["k"] / ks - 1), "s", color=fs.C[1], ms=7,
            markeredgecolor="white", label="tensor, measured")
    ref = np.abs(t["k"][-1] / ks - 1)
    ax.plot(e_fine * 1e6, ref * (e_fine / t["eps"][-1]) ** 2, "--",
            color=fs.C[1], lw=1.0, label=r"$\propto\varepsilon^2$ guide")
    ax.set_xscale("log")
    ax.set_yscale("log")
    fs.style(ax, r"interface decay length  $\varepsilon$  [$\mu$m]",
             r"$|k/k_{\mathrm{sharp}} - 1|$",
             "disk array: arithmetic vs tensor law", logy=True)
    ax.legend(fontsize=fs.FS_LEG, frameon=False, loc="lower right")
    fig.tight_layout()
    fs.save(fig, out, "keff_disk", dpi=200)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", type=Path,
                    default=Path(__file__).with_name("keff_disk.csv"))
    ap.add_argument("--L", type=float, default=1.0e-3, help="cell period [m]")
    ap.add_argument("--R", type=float, default=2.5e-4, help="disk radius [m]")
    args = ap.parse_args()

    d = read(args.csv)
    failures = check(d, args.L, args.R)
    figure(d, args.L, args.R, args.csv.parent)
    print()
    if failures:
        print(f"FAILED {len(failures)} gate(s): {', '.join(failures)}")
        return 1
    print("all gates passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

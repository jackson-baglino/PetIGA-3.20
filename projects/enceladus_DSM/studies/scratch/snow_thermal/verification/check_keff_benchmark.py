#!/usr/bin/env python3
"""Analytic gate on the k_eff cell-problem solver.

For a LAYERED medium the effective conductivity tensor is known in closed form
without solving any PDE:

    k_parallel = (1/Ly) INT k(y) dy            (arithmetic mean of k)
    k_perp     = Ly / INT dy/k(y)              (harmonic mean of k)
    off-diagonals = 0

and that holds for the DIFFUSE profile, not just for a sharp interface. So the
exact expected answer can be computed by 1-D quadrature over the same tanh
profile the initial condition lays down, and compared with what the cell-problem
solver returns. Any discrepancy is the solver's (or the spline discretisation's)
-- there is no modelling error left in the comparison.

This is a stricter and more useful gate than comparing against the sharp
interface values 1.155000 / 0.0396538, which the diffuse profile is not supposed
to reproduce at finite eps.

Usage:
    python3 check_keff_benchmark.py [--csv PATH] [--tol 5e-3]
"""
import argparse
import csv
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
K_ICE, K_AIR, LY = 2.29, 0.02, 1.0e-3


def diffuse_profile(y, eps, Ly=LY):
    """The phase field FormInitialIceSlab2D writes: ice on [0, Ly/2], periodic,
    built from the signed distance to the nearest of the two interfaces."""
    sdf = np.where(y <= 0.5 * Ly,
                   -np.minimum(y, 0.5 * Ly - y),
                   np.minimum(y - 0.5 * Ly, Ly - y))
    return np.clip(0.5 - 0.5 * np.tanh(0.5 * sdf / eps), 0.0, 1.0)


def exact_layered(eps, n=4_000_001, Ly=LY):
    """Closed-form k_par, k_perp for the continuous diffuse profile."""
    y = np.linspace(0.0, Ly, n)
    k = diffuse_profile(y, eps) * K_ICE + (1.0 - diffuse_profile(y, eps)) * K_AIR
    phi_bar = np.trapezoid(diffuse_profile(y, eps), y) / Ly
    return phi_bar, np.trapezoid(k, y) / Ly, Ly / np.trapezoid(1.0 / k, y)


def sharp_harmonic(phi):
    return 1.0 / (phi / K_ICE + (1.0 - phi) / K_AIR)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default=os.path.join(HERE, "keff_benchmark.csv"))
    ap.add_argument("--tol", type=float, default=3.5e-3,
                    help="relative tolerance on k_perp at the finest mesh for each eps. "
                         "Set from observed mesh convergence (1.7e-3 at N=512), not "
                         "from a solver tolerance -- the residual here is the spline's "
                         "representation of the profile, not the linear solve.")
    args = ap.parse_args()

    rows = list(csv.DictReader(open(args.csv)))
    ok = True

    print("Gate 1 — k_parallel against the arithmetic mean at the MEASURED phi_bar")
    print("  (exact for any profile, since k is linear in phi)\n")
    print(f"  {'sweep':9} {'eps':>10} {'N':>5} {'k_00 measured':>16} {'predicted':>16} {'rel':>10}")
    for r in rows:
        phi = float(r["phi_bar"])
        pred = phi * K_ICE + (1 - phi) * K_AIR
        got = float(r["k00"])
        rel = abs(got - pred) / pred
        flag = "" if rel < 1e-5 else "   <-- FAIL"
        if rel >= 1e-5:
            ok = False
        print(f"  {r['sweep']:9} {float(r['eps']):10.4e} {int(r['Nx']):5d} "
              f"{got:16.9e} {pred:16.9e} {rel:10.2e}{flag}")

    print("\nGate 2 — off-diagonals vanish\n")
    worst_off = max(max(abs(float(r["k01"])), abs(float(r["k10"]))) for r in rows)
    print(f"  max |k_01|, |k_10| over all runs = {worst_off:.3e}   "
          f"{'PASS' if worst_off < 1e-12 else 'FAIL'}")
    if worst_off >= 1e-12:
        ok = False

    print("\nGate 3 — k_perp against the EXACT harmonic mean of the diffuse profile")
    print("  (the real solver test: no modelling error left in the comparison)")
    print("  Judged by CONVERGENCE, not by a fixed tolerance at every mesh. k_perp")
    print("  is a harmonic mean across a 114:1 contrast, so it is dominated by the")
    print("  low-k region and the spline's pointwise error there is amplified by")
    print("  1/k^2; a coarse mesh is expected to be off by ~0.5%.\n")
    print(f"  {'sweep':9} {'eps':>10} {'N':>5} {'k_11 measured':>16} {'exact diffuse':>16} {'rel':>10}")
    errs = {}
    for r in rows:
        eps = float(r["eps"])
        _, _, kperp_exact = exact_layered(eps)
        got = float(r["k11"])
        rel = abs(got - kperp_exact) / kperp_exact
        errs.setdefault(eps, []).append((int(r["Nx"]), rel))
        print(f"  {r['sweep']:9} {eps:10.4e} {int(r['Nx']):5d} "
              f"{got:16.9e} {kperp_exact:16.9e} {rel:10.2e}")

    # 3a: refining the mesh at fixed eps must reduce the error monotonically.
    for eps, seq in sorted(errs.items()):
        if len(seq) < 2:
            continue
        seq.sort()
        rels = [e for _, e in seq]
        mono = all(b <= a * 1.02 for a, b in zip(rels, rels[1:]))
        print(f"\n  3a  eps={eps:.4e}: error vs mesh "
              f"{' -> '.join(f'{e:.2e}' for e in rels)}   "
              f"{'monotone PASS' if mono else 'NOT MONOTONE -- FAIL'}")
        if not mono:
            ok = False

    # 3b: at the finest mesh available for each eps, the error must be small.
    print()
    for eps, seq in sorted(errs.items()):
        n, rel = max(seq)
        good = rel < args.tol
        print(f"  3b  eps={eps:.4e}  finest N={n:4d}  rel={rel:.2e}   "
              f"{'PASS' if good else 'FAIL'}  (tol {args.tol:.1e})")
        if not good:
            ok = False

    print("\nContext — how far the diffuse profile sits from the SHARP-interface limit")
    print("  (a model property, not an error: the band short-circuits the series"
          " resistance)\n")
    print(f"  {'eps':>10} {'band/layer':>11} {'k_perp diffuse':>16} {'k_perp sharp':>14} {'ratio':>7}")
    for eps in sorted({float(r["eps"]) for r in rows}):
        _, _, kperp_exact = exact_layered(eps)
        ks = sharp_harmonic(0.5)
        print(f"  {eps:10.4e} {9.2 * eps / (0.5 * LY):11.3f} "
              f"{kperp_exact:16.9e} {ks:14.9e} {kperp_exact / ks:7.3f}")

    print("\nRESULT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

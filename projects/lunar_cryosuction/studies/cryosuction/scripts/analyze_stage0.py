"""Stage-0 checks for the cryosuction study (studies/cryosuction/README.md).

Reads SSA_evo.dat (columns: ssa, tot_ice, t, step, dt, tot_air, tot_rhov,
tot_mass) and bflux.dat from each run folder of a downloaded batch and prints
every number it used. Nothing here is fitted except D_eff from test A3.

A  open face, pure curvature flow (dA/dt = -D * total turning):
   A3  full disc: D_eff = -(dA/dt)/(2 pi), least squares over the run.
   A1  disc centered at x = -d, d = 7.5 um: visible area x > 0 of the
       shrinking disc, R^2 = R0^2 - 2 D_eff t, against tot_ice. PASS if the
       relative error stays below --tol while R > d + 4 eps.
   A2  same prediction with the 90 deg (Neumann) face; reported, expected to
       miss (that shows A1 can tell the two apart).
   A4  sessile cap, theta = 60 deg on y = 0, circle center (-5 um, -R cos th):
       full area A = R^2 (th - sin th cos th) falls at 2 th D_eff; the visible
       part x > 0, y > 0 is integrated numerically.
B  xi invariance: ice change and mean outflow J_x0 of B2-B4 against B1.
   Also the IC transient: ice change over the first --t-transient seconds
   against the later steady loss rate.
C  pilot: ice left, outflow, and the mass balance J_x0 vs -rho_ice dA/dt.

    venv_cryosuction/bin/python studies/cryosuction/scripts/analyze_stage0.py \
        --batch <downloaded batch folder> [--tol 0.02]
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np

EPS = 8.5840e-7
RHO_ICE = 919.0  # src/cryosuction_main.c user.rho_ice


def find_run(batch: Path, geom: str, exp: str) -> Path | None:
    hits = [p for p in batch.rglob("SSA_evo.dat") if f"{geom}__{exp}" in str(p)]
    return hits[0].parent if hits else None


def ssa(run: Path):
    a = np.loadtxt(run / "SSA_evo.dat", ndmin=2)
    return a[:, 2], a[:, 1]  # t, tot_ice


def bflux(run: Path):
    f = run / "bflux.dat"
    if not f.exists():
        return None, None
    a = np.loadtxt(f, comments="#", ndmin=2)
    return a[:, 0], a[:, 2]  # t, J_vap_x0 (first pinned face is x0 in these runs)


def seg_area(R: float, d: float) -> float:
    """Area of a disc of radius R centered at x = -d that lies in x > 0."""
    if R <= d:
        return 0.0
    return R * R * math.acos(d / R) - d * math.sqrt(R * R - d * d)


def cap_visible(R: float, th: float, cx: float, n: int = 4000) -> float:
    """Area of the circle centered (cx, -R cos th) inside x > 0, y > 0."""
    cy = -R * math.cos(th)
    x = np.linspace(max(0.0, cx - R), cx + R, n)
    half = np.sqrt(np.clip(R * R - (x - cx) ** 2, 0.0, None))
    top = cy + half
    h = np.clip(top - np.maximum(cy - half, 0.0), 0.0, None)
    h[top <= 0.0] = 0.0
    return float(np.trapezoid(h, x))


def report(name, t, model, meas, valid, tol):
    rel = np.abs(meas[valid] - model[valid]) / np.abs(model[valid])
    worst = float(rel.max()) if rel.size else float("nan")
    verdict = "PASS" if worst < tol else "FAIL"
    print(f"  {name}: max rel. error {worst:.3e} over {valid.sum()} rows "
          f"(t <= {t[valid].max() if valid.any() else 0:.3e} s)  {verdict}")
    return verdict


def part_a(batch: Path, tol: float):
    print("\n== A: open face, curvature flow ==")
    exp, expc = "openface_T-20_curvflow_7d", "openface_T-20_curvflow_7d_closed"
    r3 = find_run(batch, "openface_2D_box60um_disc_interior", exp)
    if r3 is None:
        print("  A3 run not found; A skipped")
        return
    t3, a3 = ssa(r3)
    slope = np.polyfit(t3, a3, 1)[0]
    d_eff = -slope / (2 * math.pi)
    r0 = math.sqrt(a3[0] / math.pi)
    print(f"  A3 {r3.name}: R0 = {r0:.4e} m (from A(0)), D_eff = {d_eff:.4e} m2/s "
          f"(eps^2/tau_sub = {EPS**2 / 5.8604e3:.4e})")
    d = 7.5e-6
    for tag, e in (("A1 open", exp), ("A2 90deg", expc)):
        r = find_run(batch, "openface_2D_box60um_disc_offcenter", e)
        if r is None:
            print(f"  {tag}: run not found")
            continue
        t, a = ssa(r)
        R = np.sqrt(np.clip(r0 * r0 - 2 * d_eff * t, 0.0, None))
        model = np.array([seg_area(Ri, d) for Ri in R])
        report(tag, t, model, a, R > d + 4 * EPS, tol)
    r4 = find_run(batch, "openface_2D_box60x30um_cap", exp)
    if r4 is not None:
        t, a = ssa(r4)
        th, r0c, cx = math.pi / 3, 2.0e-5, -5.0e-6
        full0 = r0c * r0c * (th - math.sin(th) * math.cos(th))
        full = full0 - 2 * th * d_eff * t
        R = np.sqrt(np.clip(full, 0.0, None) / (th - math.sin(th) * math.cos(th)))
        model = np.array([cap_visible(Ri, th, cx) for Ri in R])
        report("A4 cap at the corner", t, model, a, model > 0.2 * model[0], tol)


def part_b(batch: Path, tol: float, t_tr: float):
    print("\n== B: time-scaling invariance ==")
    g = "channel_2D_H125um_eps0.86um_deadend"
    names = {"B1 xi_v 1e-3": "cryo_T-20_h0.990_iso_3d_xiv1e-3",
             "B2 xi_v 1e-4": "cryo_T-20_h0.990_iso_3d_xiv1e-4",
             "B3 xi_v 1e-2": "cryo_T-20_h0.990_iso_3d_xiv1e-2",
             "B4 xi_T = xi_v = 1e-3": "cryo_T-20_h0.990_iso_3d_xiv1e-3_xiT1e-3"}
    runs = {k: find_run(batch, g, v) for k, v in names.items()}
    if runs["B1 xi_v 1e-3"] is None:
        print("  B1 not found; B skipped")
        return
    t1, a1 = ssa(runs["B1 xi_v 1e-3"])
    da1 = a1[-1] - a1[0]
    tb, jb = bflux(runs["B1 xi_v 1e-3"])
    j1 = np.trapezoid(jb, tb) / (tb[-1] - tb[0]) if tb is not None else float("nan")
    print(f"  B1: ice change {da1:+.4e} m2 over {t1[-1]:.3e} s, mean J_x0 {j1:.4e} kg/(m s)")
    early = t1 <= t_tr
    if early.sum() >= 2 and t1[-1] > 2 * t_tr:
        rate_late = (a1[-1] - a1[np.searchsorted(t1, 2 * t_tr)]) / (t1[-1] - 2 * t_tr)
        tr = a1[early][-1] - a1[0]
        print(f"  IC transient: ice change in the first {t_tr:g} s = {tr:+.3e} m2; "
              f"the late rate over that span would give {rate_late * t_tr:+.3e} m2")
    for k, r in runs.items():
        if k.startswith("B1") or r is None:
            continue
        t, a = ssa(r)
        ai = np.interp(t1, t, a)
        rel = abs((ai[-1] - ai[0]) - da1) / abs(da1)
        tbk, jbk = bflux(r)
        jk = np.trapezoid(jbk, tbk) / (tbk[-1] - tbk[0]) if tbk is not None else float("nan")
        print(f"  {k}: ice change differs from B1 by {rel:.3e}; mean J_x0 {jk:.4e} "
              f"({abs(jk - j1) / abs(j1):.3e} from B1)  {'PASS' if rel < tol else 'FAIL'}")


def part_c(batch: Path):
    print("\n== C: pilot ==")
    for g, e in (("throat_2D_L850um_eps0.86um_narrowL_deadend", "cryo_T-20_h0.900_iso_5d"),
                 ("regolithpore_2D_L814um_eps0.86um_throat_deadend", "cryo_T-20_h0.900_iso_5d"),
                 ("throat_2D_L850um_eps0.86um_narrowL_deadend", "cryo_T-20_h0.900_G+1000_warmdead_5d"),
                 ("throat_2D_L850um_eps0.86um_narrowL_deadend", "cryo_T-20_h0.900_insdead_5d")):
        r = find_run(batch, g, e)
        if r is None:
            print(f"  {g}:{e} not found")
            continue
        t, a = ssa(r)
        tb, jb = bflux(r)
        out = np.trapezoid(jb, tb) if tb is not None else float("nan")
        lost = RHO_ICE * (a[0] - a[-1])
        print(f"  {g}:{e}\n    t_end {t[-1]:.3e} s, ice left {a[-1] / a[0]:.3f} of initial; "
              f"mass out through x0 {out:.4e} kg/m vs ice lost {lost:.4e} kg/m "
              f"(ratio {out / lost if lost else float('nan'):.4f})")


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--batch", type=Path, required=True)
    p.add_argument("--tol", type=float, default=0.02)
    p.add_argument("--t-transient", type=float, default=600.0)
    a = p.parse_args()
    print(f"# batch {a.batch}\n# tol {a.tol}, IC transient window {a.t_transient} s, eps {EPS}, rho_ice {RHO_ICE}")
    part_a(a.batch, a.tol)
    part_b(a.batch, a.tol, a.t_transient)
    part_c(a.batch)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

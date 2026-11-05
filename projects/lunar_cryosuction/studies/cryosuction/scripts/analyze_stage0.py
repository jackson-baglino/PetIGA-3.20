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
B  xi invariance: ice change and mean outflow J_x0 of B2-B4 against B1, and
   solver cost from solver_evo.dat (step t dt snes_its ksp_its ksp/snes
   rejections): steps, Newton and Krylov iterations, rejected steps.
   Also the IC transient: ice change over the first --t-transient seconds
   against the later steady loss rate.
C  pilot: ice left, outflow, and the mass balance J_x0 vs -rho_ice dA/dt.

    venv_cryosuction/bin/python studies/cryosuction/scripts/analyze_stage0.py \
        --batch <downloaded batch folder> [--tol 0.02]
"""
from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

import numpy as np

RHO_ICE = 919.0  # src/cryosuction_main.c user.rho_ice


def run_value(run: Path, key: str) -> float:
    """eps (from the staged opts) or tau_sub (from the solver header in outp.txt)."""
    if key == "eps":
        val = float("nan")
        for f in sorted(run.glob("*.opts")):
            for line in f.read_text().splitlines():
                if line.split()[:1] == ["-eps"]:
                    val = float(line.split()[1])
        return val
    for name in ("outp.txt", "out.txt"):
        f = run / name
        if f.exists():
            m = re.search(r"tau_sub\s*=\s*([0-9.eE+-]+)\s*s", f.read_text(errors="ignore"))
            if m:
                return float(m.group(1))
    return float("nan")


def find_run(batch: Path, geom: str, exp: str) -> Path | None:
    """HPC batch layout <geom>__<exp>/, or the local run-script layout
    <geom>/<YYYY-MM-DD__HH.MM.SS>_<exp>[_stage0_<test>]/ (newest wins)."""
    hits = []
    for p in sorted(batch.rglob("SSA_evo.dat")):
        run = p.parent
        if run.name == f"{geom}__{exp}":
            hits.append(run)
            continue
        name = re.sub(r"^\d{4}-\d\d-\d\d__\d\d\.\d\d\.\d\d_", "", run.name)
        name = re.sub(r"_stage0_[A-Z]\d$", "", name)
        if run.parent.name == geom and name == exp:
            hits.append(run)
    return hits[-1] if hits else None


def ssa(run: Path):
    a = np.loadtxt(run / "SSA_evo.dat", ndmin=2)
    return a[:, 2], a[:, 1]  # t, tot_ice


def bflux(run: Path):
    f = run / "bflux.dat"
    if not f.exists():
        return None, None
    a = np.loadtxt(f, comments="#", ndmin=2)
    return a[:, 0], a[:, 2]  # t, J_vap_x0 (first pinned face is x0 in these runs)


def resolution(run: Path) -> str:
    f = run / "resolution.dat"
    if not f.exists():
        return "resolution.dat missing"
    a = np.loadtxt(f, ndmin=2)
    flag = "  <-- CHECK" if (a[:, 4] > 1).any() or (a[:, 5] > 0.01).any() else ""
    return (f"max Pi {a[:, 4].max():.3f} (steps > 0.5: {(a[:, 4] > 0.5).sum()}, > 1: {(a[:, 4] > 1).sum()}); "
            f"max phi on open face {a[:, 5].max():.2e}{flag}")


def solver_cost(run: Path) -> str:
    f = run / "solver_evo.dat"
    if not f.exists():
        return "solver_evo.dat missing"
    a = np.loadtxt(f, ndmin=2)
    return (f"{a.shape[0]} steps, {int(a[:, 3].sum())} Newton its "
            f"({a[:, 3].mean():.2f}/step), {int(a[:, 4].sum())} Krylov its, "
            f"{int(a[:, 6].sum())} rejected")


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


def fig_a(batch: Path, out: Path, d_eff: float, r0: float):
    """Ice area against the exact curvature-flow law, A1-A4."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    exp = "openface_T-20_curvflow_22h_ac1e-2"
    th, r0c, cx, d = math.pi / 3, 2.0e-5, -5.0e-6, 7.5e-6
    full0 = r0c * r0c * (th - math.sin(th) * math.cos(th))
    cases = [("A3 full disc", "openface_2D_box60um_disc_interior", exp, "#52514e",
              lambda t: math.pi * max(r0 * r0 - 2 * d_eff * t, 0.0)),
             ("A1 cut disc, open face", "openface_2D_box60um_disc_offcenter", exp, "#eb6834",
              lambda t: seg_area(math.sqrt(max(r0 * r0 - 2 * d_eff * t, 0.0)), d)),
             ("A2 cut disc, 90° face", "openface_2D_box60um_disc_offcenter", exp + "_closed", "#6250d6",
              lambda t: seg_area(math.sqrt(max(r0 * r0 - 2 * d_eff * t, 0.0)), d)),
             ("A4 wall cap at corner, open face", "openface_2D_box60x30um_cap", exp, "#2a78d6",
              lambda t: cap_visible(math.sqrt(max(full0 - 2 * th * d_eff * t, 0.0)
                                              / (th - math.sin(th) * math.cos(th))), th, cx))]
    fig, ax = plt.subplots(figsize=(8, 5))
    for label, g, e, c, law in cases:
        r = find_run(batch, g, e)
        if r is None:
            continue
        t, a = ssa(r)
        tt = np.linspace(0, t[-1], 200)
        ax.plot(t / 86400, a * 1e12, color=c, lw=2, label=f"{label}: solver")
        if label.startswith("A1"):
            continue  # A1 and A2 share one exact law, drawn with A2
        ax.plot(tt / 86400, [law(x) * 1e12 for x in tt], color=c, lw=1.2, ls="--",
                label="A1 & A2 cut disc: exact law" if label.startswith("A2") else f"{label}: exact law")
    ax.set_xlabel("time [days]")
    ax.set_ylabel("ice area in the domain [µm²]")
    ax.set_title("Open-face test: pure curvature flow against the exact area law\n"
                 f"(D_eff = {d_eff:.3e} m²/s from A3; dashed = exact, solid = solver)", fontsize=10)
    ax.legend(fontsize=7.5, frameon=False, ncol=2, loc="upper right")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    fig.tight_layout()
    fig.savefig(out, dpi=180)
    print(f"  figure: {out}")


def part_a(batch: Path, tol: float, fig: Path | None = None):
    print("\n== A: open face, curvature flow ==")
    exp, expc = "openface_T-20_curvflow_22h_ac1e-2", "openface_T-20_curvflow_22h_ac1e-2_closed"
    r3 = find_run(batch, "openface_2D_box60um_disc_interior", exp)
    if r3 is None:
        print("  A3 run not found; A skipped")
        return
    t3, a3 = ssa(r3)
    slope = np.polyfit(t3, a3, 1)[0]
    d_eff = -slope / (2 * math.pi)
    r0 = math.sqrt(a3[0] / math.pi)
    eps, tau = run_value(r3, "eps"), run_value(r3, "tau_sub")
    print(f"  A3 {r3.name}: R0 = {r0:.4e} m (from A(0)), D_eff = {d_eff:.4e} m2/s "
          f"(eps^2/tau_sub = {eps**2 / tau:.4e}, eps {eps:.4e}, tau_sub {tau:.4e})")
    d = 7.5e-6
    for tag, e in (("A1 open", exp), ("A2 90deg", expc)):
        r = find_run(batch, "openface_2D_box60um_disc_offcenter", e)
        if r is None:
            print(f"  {tag}: run not found")
            continue
        t, a = ssa(r)
        R = np.sqrt(np.clip(r0 * r0 - 2 * d_eff * t, 0.0, None))
        model = np.array([seg_area(Ri, d) for Ri in R])
        report(tag, t, model, a, R > d + 4 * eps, tol)
    r4 = find_run(batch, "openface_2D_box60x30um_cap", exp)
    if r4 is not None:
        t, a = ssa(r4)
        th, r0c, cx = math.pi / 3, 2.0e-5, -5.0e-6
        full0 = r0c * r0c * (th - math.sin(th) * math.cos(th))
        full = full0 - 2 * th * d_eff * t
        R = np.sqrt(np.clip(full, 0.0, None) / (th - math.sin(th) * math.cos(th)))
        model = np.array([cap_visible(Ri, th, cx) for Ri in R])
        report("A4 cap at the corner", t, model, a, model > 0.2 * model[0], tol)
    if fig is not None:
        fig_a(batch, fig, d_eff, r0)


def part_b(batch: Path, tol: float, t_tr: float):
    print("\n== B: time-scaling invariance ==")
    g = "channel_2D_H125um_eps0.39um_deadend"
    names = {"B1 xi_v 1e-3": "cryo_T-20_h0.990_iso_3d_ac1e-2_xiv1e-3",
             "B2 xi_v 1e-4": "cryo_T-20_h0.990_iso_3d_ac1e-2_xiv1e-4",
             "B3 xi_v 1e-2": "cryo_T-20_h0.990_iso_3d_ac1e-2_xiv1e-2",
             "B4 xi_T = xi_v = 1e-3": "cryo_T-20_h0.990_iso_3d_ac1e-2_xiv1e-3_xiT1e-3"}
    runs = {k: find_run(batch, g, v) for k, v in names.items()}
    if runs["B1 xi_v 1e-3"] is None:
        print("  B1 not found; B skipped")
        return
    t1, a1 = ssa(runs["B1 xi_v 1e-3"])
    da1 = a1[-1] - a1[0]
    tb, jb = bflux(runs["B1 xi_v 1e-3"])
    j1 = np.trapezoid(jb, tb) / (tb[-1] - tb[0]) if tb is not None else float("nan")
    print(f"  B1: ice change {da1:+.4e} m2 over {t1[-1]:.3e} s, mean J_x0 {j1:.4e} kg/(m s)")
    for k, r in runs.items():
        if r is not None:
            print(f"  cost {k}: {solver_cost(r)}\n  resolution {k}: {resolution(r)}")
    early = t1 <= t_tr
    if early.sum() >= 2 and t1[-1] > 2 * t_tr:
        rate_late = (a1[-1] - a1[np.searchsorted(t1, 2 * t_tr)]) / (t1[-1] - 2 * t_tr)
        tr = a1[early][-1] - a1[0]
        print(f"  IC transient: ice change in the first {t_tr:g} s = {tr:+.3e} m2; "
              f"the late rate over that span would give {rate_late * t_tr:+.3e} m2")
    for k, r in runs.items():
        if k.startswith("B1") or r is None:
            continue
        # Compare over the window both runs reached (a crashed run ends early;
        # np.interp would otherwise hold its last value).
        t, a = ssa(r)
        tc = min(t[-1], t1[-1])
        ref = np.interp(tc, t1, a1) - a1[0]
        rel = abs((np.interp(tc, t, a) - a[0]) - ref) / abs(ref)
        tbk, jbk = bflux(r)

        def mean_j(tb_, jb_, tc=tc):
            m = tb_ <= tc
            return np.trapezoid(jb_[m], tb_[m]) / (tb_[m][-1] - tb_[m][0])
        jk, j1c = mean_j(tbk, jbk), mean_j(tb, jb)
        print(f"  {k}: over t <= {tc:.3e} s, ice change differs from B1 by {rel:.3e}; mean J_x0 "
              f"{jk:.4e} vs B1 {j1c:.4e} ({abs(jk - j1c) / abs(j1c):.3e})  {'PASS' if rel < tol else 'FAIL'}")


def part_c(batch: Path):
    print("\n== C: pilot ==")
    for g, e in (("throat_2D_L850um_eps0.86um_narrowL_deadend", "cryo_T-20_h0.990_iso_5d_ac1e-2"),
                 ("regolithpore_2D_L814um_eps0.86um_throat_deadend", "cryo_T-20_h0.998_iso_5d_ac1e-2"),
                 ("throat_2D_L850um_eps0.86um_narrowL_deadend", "cryo_T-20_h0.990_G+100_warmdead_5d_ac1e-2"),
                 ("throat_2D_L850um_eps0.86um_narrowL_deadend", "cryo_T-20_h0.990_insdead_5d_ac1e-2")):
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
              f"(ratio {out / lost if lost else float('nan'):.4f})\n    resolution: {resolution(r)}")


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--batch", type=Path, required=True)
    p.add_argument("--tol", type=float, default=0.02)
    p.add_argument("--t-transient", type=float, default=600.0)
    p.add_argument("--fig-a", type=Path, default=None, help="write the A area-law figure here")
    a = p.parse_args()
    print(f"# batch {a.batch}\n# tol {a.tol}, IC transient window {a.t_transient} s, rho_ice {RHO_ICE}")
    part_a(a.batch, a.tol, a.fig_a)
    part_b(a.batch, a.tol, a.t_transient)
    part_c(a.batch)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

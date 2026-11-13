#!/usr/bin/env python3
"""check_periodic_x.py: the deterministic checks of verify_periodic_x.sh.

Every number is computed here from the files the builder and the solver
wrote; nothing is taken from model output. Each check prints one row

    gate,case,quantity,measured,expected,tolerance,result

to stdout (the driver appends them to verify_periodic_x.csv).

INTEGRALS. Fields are read with pplib.read_sol_vec (which appends the p wrapped
control columns of a -periodic_x run) and evaluated with igakit at 3 x 3 Gauss
points per element, the rule PetIGA uses for p = 2 (p + 1 points), with the
geometric Jacobian from the mesh (igasol.dat):

    A   = int phi dA                                  ice area per unit depth
    X_c = circular mean of x weighted by phi           lens center (mod Lx)
    Q   = xi_v rho_ice int phi + int rho_v (1 - phi)   the water the scaled
          equations conserve: R_vap tested with 1 sums to
          d/dt int [(1 - phi) rho_v + xi_v rho_ice phi] = 0 when no vapor
          crosses the boundary (assembly.c, R[a][2]; the air_lim floor and the
          time discretization of the product rule are the only departures).

The circular mean, X_c = Lx/(2 pi) atan2(int phi sin(2 pi x/Lx), int phi cos(2 pi x/Lx)),
is the center of a body on a ring; it is continuous when the lens crosses the
seam. Displacements are wrapped to (-Lx/2, Lx/2].

Usage (from the project root):
    python check_periodic_x.py mesh   <geometry dir>
    python check_periodic_x.py runs   <runs.txt>      (lines: label run_dir)
"""
from __future__ import annotations

import math
import re
import struct
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "postprocess"))
sys.path.insert(0, str(ROOT / "preprocess"))
import pplib  # noqa: E402

ROWS: list[str] = []


def row(gate, case, quantity, measured, expected, tol, ok):
    r = f"{gate},{case},{quantity},{measured},{expected},{tol},{'PASS' if ok else 'FAIL'}"
    ROWS.append(r)
    print(r, flush=True)


# ------------------------------------------------------------------ mesh ----
def read_phi(path):
    with open(path, "rb") as f:
        assert f.read(8) == b"PHIIC001"
        mx, my = struct.unpack("<qq", f.read(16))
        return np.frombuffer(f.read(), dtype="<f8").reshape(my, mx, 3)


def check_mesh(geo_dir: Path):
    """Gate a: the periodic mesh and its phi file match across the seam."""
    from igakit.io import PetIGA
    for case in ("px_flat120x40_lens_seam_eps1.5um", "px_cd120_W40_wt16_lens_eps1.5um",
                 "px_cd120_W40_wt16_shift_lens_eps1.5um"):
        nrb = PetIGA().read(str(geo_dir / f"{case}.dat"))
        p, U = nrb.degree[0], np.asarray(nrb.knots[0])
        ncol = nrb.shape[0]
        N = ncol - p
        lx = float(re.search(r"^-Lx (\S+)", (geo_dir / f"{case}.opts").read_text(), re.M).group(1))
        # a1: knots are PetIGA's periodic uniform vector (petigaaxis.c IGAAxisInitUniform)
        want = (np.arange(N + 2 * p + 1) - p) / N
        e = float(np.abs(U - want).max())
        row("a1", case, "knot_vs_petiga_periodic", f"{e:.2e}", 0, "1e-14", e < 1e-14)
        # a2: last p control columns are the first p shifted by exactly Lx
        Cw = nrb.control
        dx = float(np.abs(Cw[N:, :, 0] - Cw[:p, :, 0] - lx).max())
        dy = float(np.abs(Cw[N:, :, 1] - Cw[:p, :, 1]).max())
        row("a2", case, "wrapped_ctrl_dx_minus_Lx_m", f"{dx:.2e}", 0, f"{1e-12 * lx:.0e}", dx < 1e-12 * lx)
        row("a2", case, "wrapped_ctrl_dy_m", f"{dy:.2e}", 0, "0 (copied)", dy == 0.0)
        # a3: both walls and their slope match across the seam (mapping at u = 0 and u = 1)
        v = np.array([0.0, 1.0])
        P0, P1 = nrb(0.0, v), nrb(1.0, v)
        G0 = nrb.gradient(u=np.array([0.0]), v=v)
        G1 = nrb.gradient(u=np.array([1.0]), v=v)
        jy = float(np.abs(P1[:, 1] - P0[:, 1]).max())
        jx = float(np.abs(P1[:, 0] - P0[:, 0] - lx).max())
        # slope dy/dx of each wall = (dy/du)/(dx/du)
        s0 = G0[0, :, 1, 0] / G0[0, :, 0, 0]
        s1 = G1[0, :, 1, 0] / G1[0, :, 0, 0]
        js = float(np.abs(s1 - s0).max())
        row("a3", case, "wall_y_jump_at_seam_m", f"{jy:.2e}", 0, f"{1e-12 * lx:.0e}", jy < 1e-12 * lx)
        row("a3", case, "seam_x_minus_Lx_m", f"{jx:.2e}", 0, f"{1e-12 * lx:.0e}", jx < 1e-12 * lx)
        row("a3", case, "wall_slope_jump_at_seam", f"{js:.2e}", 0, "1e-9", js < 1e-9)
        # a4: the phi file holds the N periodic nodes at the first N control points
        rec = read_phi(geo_dir / f"{case}_phi.bin")
        okshape = rec.shape[1] == N and rec.shape[0] == nrb.shape[1]
        row("a4", case, "phi_file_columns", rec.shape[1], N, "exact", okshape)
        e = float(np.abs(rec[..., 0].T - Cw[:N, :, 0]).max() + np.abs(rec[..., 1].T - Cw[:N, :, 1]).max())
        row("a4", case, "phi_file_xy_vs_ctrl_m", f"{e:.2e}", 0, "0", e == 0.0)
    # a5: the seam lens is the mid lens shifted by Lx/2 = N/2 columns, node for node
    for a, b in (("px_flat120x40_lens_mid_eps1.5um", "px_flat120x40_lens_seam_eps1.5um"),
                 ("px_cd120_W40_wt16_lens_eps1.5um", "px_cd120_W40_wt16_shift_lens_eps1.5um")):
        ra, rb = read_phi(geo_dir / f"{a}_phi.bin"), read_phi(geo_dir / f"{b}_phi.bin")
        N = ra.shape[1]
        e = float(np.abs(np.roll(ra[..., 2], -N // 2, axis=1) - rb[..., 2]).max())
        ey = float(np.abs(np.roll(ra[..., 1], -N // 2, axis=1) - rb[..., 1]).max())
        row("a5", b, "phi_vs_shifted_" + a.split("_")[1], f"{e:.2e}", 0, "1e-9", e < 1e-9)
        row("a5", b, "y_vs_shifted_" + a.split("_")[1] + "_m", f"{ey:.2e}", 0, "1e-15", ey < 1e-15)


# ------------------------------------------------------------------ runs ----
class Run:
    """Integrals of every snapshot of one run."""

    def __init__(self, label, rd):
        from igakit.io import PetIGA
        self.label, self.rd = label, Path(rd)
        out = (self.rd / "outp.txt").read_text(errors="replace")
        self.rho_ice = float(re.search(r"rho_ice\s+=\s+(\S+)", out).group(1))
        self.xi_v = float(re.search(r"xi_v\s+=\s+(\S+)", out).group(1))
        self.periodic = "periodic_x: x axis periodic" in out
        self.nprocs = int(re.search(r"Running on (\d+) process", out).group(1))
        nrb = PetIGA().read(str(self.rd / "igasol.dat"))
        self.nrb = nrb
        ux, wx = self._gauss(nrb, 0)
        uy, wy = self._gauss(nrb, 1)
        self.u, self.v = ux, uy
        C = nrb(ux, uy)
        G = nrb.gradient(u=ux, v=uy)
        det = G[..., 0, 0] * G[..., 1, 1] - G[..., 0, 1] * G[..., 1, 0]
        self.W = np.abs(det) * wx[:, None] * wy[None, :]
        self.X = C[..., 0]
        # parametric limits: [0, 1] for a builder mesh, [0, Lx] for the built-in rectangle
        U, p = np.asarray(nrb.knots[0]), nrb.degree[0]
        self.lx = float(nrb(U[-p - 1], 0.0)[0] - nrb(U[p], 0.0)[0])
        tmap = pplib.snapshot_times(str(self.rd))
        self.snaps = []
        for f in sorted(self.rd.glob("sol_*.dat")):
            st = int(re.search(r"sol_(\d+)", f.name).group(1))
            F = pplib.read_sol_vec(str(f), nrb)
            _, Fq = nrb(ux, uy, fields=F)
            phi, rv = Fq[..., 0], Fq[..., 2]
            A = float((self.W * phi).sum())
            k = 2 * math.pi / self.lx
            ang = math.atan2(float((self.W * phi * np.sin(k * self.X)).sum()),
                             float((self.W * phi * np.cos(k * self.X)).sum()))
            xc = (ang / k) % self.lx
            Q = self.xi_v * self.rho_ice * A + float((self.W * rv * (1.0 - phi)).sum())
            self.snaps.append((st, tmap.get(st), A, xc, Q))

    @staticmethod
    def _gauss(nrb, axis):
        p, U = nrb.degree[axis], np.asarray(nrb.knots[axis])
        br = np.unique(U[p:len(U) - p])
        g, w = np.polynomial.legendre.leggauss(p + 1)
        a, b = br[:-1, None], br[1:, None]
        u = (0.5 * (a + b) + 0.5 * (b - a) * g[None, :]).ravel()
        ww = (0.5 * (b - a) * w[None, :]).ravel()
        return u, ww

    def first(self):
        return self.snaps[0]

    def last(self):
        return self.snaps[-1]


def wrap(d, lx):
    return (d + 0.5 * lx) % lx - 0.5 * lx


def check_runs(listfile: Path):
    runs = {}
    for line in listfile.read_text().splitlines():
        if line.strip():
            lab, rd = line.split()
            runs[lab] = Run(lab, rd)
            r = runs[lab]
            print(f"# {lab}: {r.rd.name}, {r.nprocs} rank(s), periodic_x={r.periodic}, "
                  f"{len(r.snaps)} snapshots, t_end {r.last()[1]:.6e} s", flush=True)
            for st, t, A, xc, Q in r.snaps:
                print(f"#   step {st:5d} t {t:.6e}  A {A:.10e}  X_c {xc * 1e6:.6f} um  Q {Q:.10e}")
    eps = 1.5e-6
    P, O = runs["P_mid"], runs["O_mid"]

    # b: translation invariance and conservation in the flat periodic channel
    def drift(r):
        return wrap(r.last()[3] - r.first()[3], r.lx)

    def dq(r):
        return (r.last()[4] - r.first()[4]) / r.first()[4]
    dP, dO = drift(P), drift(O)
    qP, qO = dq(P), dq(O)
    # The lens sits mid-cell in a uniform channel: no direction is preferred, so the
    # drift is roundoff in the solve. 1e-4 eps (1.5e-10 m) is 4 orders below the
    # interface width and 2 orders below the d sign-check signal it must not mask.
    tol_drift = 1e-4 * eps
    row("b", "P_mid", "lens_center_drift_m", f"{dP:.3e}", 0, f"{tol_drift:.1e}", abs(dP) < tol_drift)
    row("b", "O_mid(ref)", "lens_center_drift_m", f"{dO:.3e}", 0, "info", True)
    # Conservation: the open sealed channel is the same scheme with no seam; the
    # periodic run must conserve Q as well as it does (2x, floor 1e-8: SNES rtol 1e-6
    # on a residual whose Q-part is ~1e-2 of its norm).
    tol_q = max(2 * abs(qO), 1e-8)
    row("b", "P_mid", "water_Q_rel_change", f"{qP:.3e}", 0, f"{tol_q:.1e}", abs(qP) < tol_q)
    row("b", "O_mid(ref)", "water_Q_rel_change", f"{qO:.3e}", 0, "info", True)
    # periodic vs open: the lens sees the same symmetric problem; areas should agree
    # closely but not exactly (the open channel's x faces are 90-degree Neumann ends)
    eA = abs(P.last()[2] - O.last()[2]) / O.last()[2]
    row("b", "P_mid_vs_O_mid", "final_area_rel_diff", f"{eA:.3e}", 0, "info", True)

    # c: seam crossing = the same configuration shifted by Lx/2
    for lab_a, lab_b, case in (("P_mid", "P_seam", "c1_flat"), ("CD", "CD_shift", "c2_cd")):
        a, b = runs[lab_a], runs[lab_b]
        sa = {s[0]: s for s in a.snaps}
        sb = {s[0]: s for s in b.snaps}
        common = sorted(set(sa) & set(sb))
        same_t = all(abs(sa[k][1] - sb[k][1]) <= 1e-9 * max(1.0, sa[k][1]) for k in common)
        row("c", case, "snapshot_steps_and_times_match", f"{len(common)}/{len(a.snaps)}",
            f"{len(a.snaps)}/{len(a.snaps)}", "exact", same_t and len(common) == len(a.snaps))
        eA = max(abs(sa[k][2] - sb[k][2]) / sa[k][2] for k in common)
        eX = max(abs(wrap(sb[k][3] - (sa[k][3] - 0.5 * a.lx), a.lx)) for k in common)
        # The two problems are one discrete problem with the nodes renumbered by N/2
        # columns; Newton (rtol 1e-6) and BiCGStab+ILU (rtol 1e-8) see a permuted
        # matrix, so results agree to the solver tolerance, not bitwise. 1e-6 in area
        # and 1e-4 eps in position bound that with margin; a seam that leaks or
        # distorts shows at the 1e-2 level.
        row("c", case, "max_area_rel_diff", f"{eA:.3e}", 0, "1e-6", eA < 1e-6)
        row("c", case, "max_center_diff_mod_Lx_m", f"{eX:.3e}", 0, f"{1e-4 * eps:.1e}", eX < 1e-4 * eps)

    # d: converging-diverging cell: the lens moves toward the throat (x = Lx/2)
    cd = runs["CD"]
    xt = 0.5 * cd.lx
    x0, x1 = cd.first()[3], cd.last()[3]
    toward = wrap(x1 - x0, cd.lx) * np.sign(wrap(xt - x0, cd.lx))
    row("d", "CD", "center_move_toward_throat_m", f"{toward:.3e}", "> 0", f"> {100 * tol_drift:.1e} (100x b tol)",
        toward > 100 * tol_drift)
    row("d", "CD", "distance_to_throat_start_end_um", f"{abs(wrap(xt - x0, cd.lx)) * 1e6:.3f} -> "
        f"{abs(wrap(xt - x1, cd.lx)) * 1e6:.3f}", "decreasing", "-", abs(wrap(xt - x1, cd.lx)) < abs(wrap(xt - x0, cd.lx)))
    # In the cd cell the lens moves, so every step carries a real vapor flux and
    # the per-step Newton tolerance (vapor residual accepted at -snes_atol 1e-6,
    # i.e. after one iteration) leaves a per-step error in Q. Bound: 1e-6 relative
    # over the run, the SNES rtol; gate h shows the change is that tolerance
    # (it shrinks when the solve is tightened), not a leak at the seam.
    qc = dq(cd)
    row("d", "CD", "water_Q_rel_change", f"{qc:.3e}", 0, "1e-6 (Newton tol; see h)", abs(qc) < 1e-6)
    swept = abs(wrap(x1 - x0, cd.lx)) * 1.6e-5        # lens shift x throat gap: area moved
    row("d", "CD", "Q_change_per_ice_moved", f"{qc * cd.first()[4] / (cd.xi_v * cd.rho_ice * swept):.3e}",
        0, "info", True)

    # h: tightened Newton (-snes_rtol 1e-10 -snes_atol 1e-30 -snes_stol 1e-14)
    if "CD_tight" in runs and "CD_shift_tight" in runs:
        a, b = runs["CD_tight"], runs["CD_shift_tight"]
        qa, qb = dq(a), dq(b)
        row("h", "CD_tight", "water_Q_rel_change", f"{qa:.3e}", f"< |CD| {abs(qc):.1e}/10", "10x smaller",
            abs(qa) < 0.1 * abs(qc))
        row("h", "CD_shift_tight", "water_Q_rel_change", f"{qb:.3e}", f"< |CD_shift| {abs(dq(runs['CD_shift'])):.1e}/10",
            "10x smaller", abs(qb) < 0.1 * abs(dq(runs["CD_shift"])))
        sa = {s[0]: s for s in a.snaps}
        common = [s for s in b.snaps if s[0] in sa]
        eA = max(abs(s[2] - sa[s[0]][2]) / sa[s[0]][2] for s in common)
        eX = max(abs(wrap(s[3] - (sa[s[0]][3] - 0.5 * a.lx), a.lx)) for s in common)
        row("h", "c2_tight", "snapshots_in_common", f"{len(common)}/{len(a.snaps)}", f"{len(a.snaps)}/{len(a.snaps)}",
            "exact", len(common) == len(a.snaps))
        row("h", "c2_tight", "max_area_rel_diff", f"{eA:.3e}", 0, "1e-8", eA < 1e-8)
        row("h", "c2_tight", "max_center_diff_mod_Lx_m", f"{eX:.3e}", 0, f"{1e-6 * eps:.1e}", eX < 1e-6 * eps)
        tow = wrap(a.last()[3] - a.first()[3], a.lx) * np.sign(wrap(xt - a.first()[3], a.lx))
        row("h", "CD_tight", "center_move_toward_throat_m", f"{tow:.3e}", "> 0; within 1 % of CD",
            "1 %", tow > 0 and abs(tow - toward) < 0.01 * abs(toward))

    # f: built-in rectangle (-Nx/-Ny, no -geom_file) = the flat periodic geom_file mesh
    if "RECT_mid" in runs:
        r = runs["RECT_mid"]
        eA = abs(r.last()[2] - P.last()[2]) / P.last()[2]
        eX = abs(wrap(r.last()[3] - P.last()[3], r.lx))
        same_t = abs(r.last()[1] - P.last()[1]) < 1e-9 * P.last()[1]
        row("f", "RECT_mid_vs_P_mid", "final_time_match", same_t, True, "exact", same_t)
        row("f", "RECT_mid_vs_P_mid", "final_area_rel_diff", f"{eA:.3e}", 0, "1e-6", eA < 1e-6)
        row("f", "RECT_mid_vs_P_mid", "final_center_diff_m", f"{eX:.3e}", 0, f"{1e-4 * eps:.1e}", eX < 1e-4 * eps)

    # g: MPI, ranks split along x: same result as the serial run
    if "CD_mpi" in runs:
        r = runs["CD_mpi"]
        sa = {s[0]: s for s in cd.snaps}
        common = [s for s in r.snaps if s[0] in sa]
        eA = max(abs(s[2] - sa[s[0]][2]) / sa[s[0]][2] for s in common)
        eX = max(abs(wrap(s[3] - sa[s[0]][3], r.lx)) for s in common)
        row("g", f"CD_mpi({r.nprocs} ranks)_vs_CD", "snapshots_in_common", f"{len(common)}/{len(cd.snaps)}",
            f"{len(cd.snaps)}/{len(cd.snaps)}", "exact", len(common) == len(cd.snaps))
        row("g", f"CD_mpi({r.nprocs} ranks)_vs_CD", "max_area_rel_diff", f"{eA:.3e}", 0, "1e-6", eA < 1e-6)
        row("g", f"CD_mpi({r.nprocs} ranks)_vs_CD", "max_center_diff_m", f"{eX:.3e}", 0, f"{1e-4 * eps:.1e}",
            eX < 1e-4 * eps)
        tow = wrap(r.last()[3] - r.first()[3], r.lx) * np.sign(wrap(xt - r.first()[3], r.lx))
        row("g", f"CD_mpi({r.nprocs} ranks)", "center_move_toward_throat_m", f"{tow:.3e}", "> 0", "-", tow > 0)


def main(argv):
    mode, arg = argv[1], Path(argv[2])
    if mode == "mesh":
        check_mesh(arg)
    else:
        check_runs(arg)
    return 0 if all(r.endswith("PASS") for r in ROWS) else 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))

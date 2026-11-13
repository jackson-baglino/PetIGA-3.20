#!/usr/bin/env python3
"""Deterministic checks of the thermal face conditions (-qflux_faces, -Twall_faces).

Called by verify_heatflux.sh with the run folders it produced. Each subcommand
compares a run against a closed-form solution and prints one CSV row per check:

    test,case,quantity,measured,expected,tolerance,result

Field values come from the spline itself: igasol.dat + sol_N.dat evaluated
with igakit (pplib's "nurbs" source), not the control net. Times come from
SSA_evo.dat via pplib.snapshot_times. Stored heat is integrated EXACTLY from
the control values: for a B-spline basis, int N_i dx = (U[i+p+1] - U[i])/(p+1),
so int T dA = sum_ij c_ij wx_i wy_j with no quadrature error.

Material constants are the solver defaults (src/cryosuction_main.c, user.* at
the top of main: k 2.29/0.02 W/m/K, rho 919/1.341 kg/m3, cp 1960/1044 J/kg/K,
latent heat 2.83e6 J/kg), overridden by any value in the run's staged opts.
"""
from __future__ import annotations

import argparse
import math
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..", "..", "postprocess")))
import pplib  # noqa: E402

DEFAULTS = {"-thcond_ice": 2.29, "-thcond_air": 0.02, "-rho_ice": 919.0,
            "-rho_air": 1.341, "-cp_ice": 1.96e3, "-cp_air": 1.044e3,
            "-xi_T": 1.0, "-temp": -20.0}
LAT_SUB = 2.83e6   # user.lat_sub, src/cryosuction_main.c


class Run:
    def __init__(self, rd: str):
        from igakit.io import PetIGA
        self.rd = rd
        self.opts = pplib.read_opts(rd)
        self.io = PetIGA()
        self.nrb = self.io.read(os.path.join(rd, "igasol.dat"))
        self.tmap = pplib.snapshot_times(rd)
        files = sorted(f for f in os.listdir(rd) if re.fullmatch(r"sol_\d+\.dat", f))
        self.snaps = sorted((self.tmap[int(f[4:9])], int(f[4:9]), f)
                            for f in files if int(f[4:9]) in self.tmap)
        self.Lx = self.opt("-Lx")
        self.Ly = self.opt("-Ly")

    def opt(self, key, default=None):
        d = DEFAULTS.get(key, default)
        return pplib.opt_float(self.opts, key, d)

    def phase(self, phi_mean):
        """k, rho*cp, alpha of the uniform phase present (phi 0 or 1)."""
        if phi_mean > 0.5:
            k, rc = self.opt("-thcond_ice"), self.opt("-rho_ice") * self.opt("-cp_ice")
        else:
            k, rc = self.opt("-thcond_air"), self.opt("-rho_air") * self.opt("-cp_air")
        return k, rc, k / rc

    def fields(self, fname):
        return self.io.read_vec(os.path.join(self.rd, fname), self.nrb)

    def evaluate(self, F, x, y):
        """T and phi at physical (x, y); knots are physical (no -geom_file)."""
        _, V = self.nrb(np.asarray(x, float), np.asarray(y, float), fields=F)
        return V[..., 1], V[..., 0]

    def weights(self):
        w = []
        for d in range(2):
            U, p = np.asarray(self.nrb.knots[d]), self.nrb.degree[d]
            n = len(U) - p - 1
            w.append(np.array([(U[i + p + 1] - U[i]) / (p + 1) for i in range(n)]))
        return w

    def integral(self, F, dof):
        wx, wy = self.weights()
        return float(np.einsum("ij,i,j->", F[..., dof], wx, wy))

    def nearest(self, t):
        return min(self.snaps, key=lambda s: abs(s[0] - t))

    def final(self):
        return self.snaps[-1]


def row(test, case, quantity, measured, expected, tol, ok):
    print(f"{test},{case},{quantity},{measured:.6e},{expected:.6e},{tol:.1e},{'PASS' if ok else 'FAIL'}")
    return ok


# ---------------------------------------------------------------- test a ----
def test_a(args):
    """Steady flux q on x0, T = T0 on x1: T = T0 + (q/k)(Lx - x)."""
    ok = True
    for rd in args.runs:
        r = Run(rd)
        t, st, f = r.final()
        F = r.fields(f)
        q, T0 = r.opt("-qflux"), r.opt("-temp")
        k, _, _ = r.phase(F[..., 0].mean())
        x = np.linspace(0.0, r.Lx, 401)
        y = np.linspace(0.0, r.Ly, 5)
        T, _ = r.evaluate(F, x, y)
        Tex = T0 + (q / k) * (r.Lx - x)[:, None]
        rise = q * r.Lx / k
        err = float(np.max(np.abs(T - Tex)) / rise)
        case = f"Nx={r.nrb.shape[0] - r.nrb.degree[0]}"
        ok &= row("a", case, "max|T-T_lin|/(qL/k)", err, 0.0, args.tol, err < args.tol)
        # Steady energy balance: the heat leaving through x1 equals q*Ly.
        bf = np.loadtxt(os.path.join(rd, "bflux.dat"))
        hdr = open(os.path.join(rd, "bflux.dat")).read().splitlines()[1].split()[1:]
        Qx1 = bf[-1, hdr.index("Q_heat_x1")]
        rel = abs(Qx1 - q * r.Ly) / (q * r.Ly)
        ok &= row("a", case, "|Q_out(x1)-q*Ly|/(q*Ly)", rel, 0.0, args.tol_flux, rel < args.tol_flux)
    return ok


# ---------------------------------------------------------------- test b ----
def slab_profile(x, t, q, k, a):
    """Constant flux q into a semi-infinite solid (Carslaw and Jaeger 1959, 2.9 (6))."""
    s = math.sqrt(a * t)
    return (2.0 * q / k) * (s / math.sqrt(math.pi) * np.exp(-x**2 / (4 * a * t))
                            - x / 2.0 * np.vectorize(math.erfc)(x / (2 * s)))


def test_b(args):
    ok = True
    for rd in args.runs:
        r = Run(rd)
        q, T0 = r.opt("-qflux"), r.opt("-temp")
        F0 = r.fields(r.snaps[0][2])
        k, rc, a = r.phase(F0[..., 0].mean())
        tag = "ice" if F0[..., 0].mean() > 0.5 else "air"
        tf = r.final()[0]
        for frac in (0.125, 0.25, 0.5, 1.0):
            t, st, f = r.nearest(frac * tf)
            F = r.fields(f)
            x = np.linspace(0.0, r.Lx, 1201)
            T, _ = r.evaluate(F, x, [r.Ly / 2])
            T = T[:, 0] - T0
            # Insulated far face x = Lx: method of images (Carslaw and Jaeger
            # 1959, 2.9 with images at 2nLx +- x); three pairs are exact to
            # far below the tolerance while sqrt(alpha t) < Lx.
            ex = slab_profile(x, t, q, k, a)
            for n in (1, 2, 3):
                ex = ex + slab_profile(2 * n * r.Lx - x, t, q, k, a) + slab_profile(2 * n * r.Lx + x, t, q, k, a)
            e_surf = abs(T[0] - ex[0]) / ex[0]
            e_prof = float(np.max(np.abs(T - ex)) / ex[0])
            case = f"{tag}_t={t:.4e}"
            ok &= row("b", case, "T(0,t)-T0 [K]", T[0], ex[0], args.tol, e_surf < args.tol)
            ok &= row("b", case, "max_x|T-T_ex|/(T_ex(0)-T0)", e_prof, 0.0, args.tol, e_prof < args.tol)
    return ok


# ---------------------------------------------------------------- test c ----
def test_c(args):
    """runs: COARSE FINE COARSE_DT. COARSE and FINE differ in h (x2), COARSE_DT
    is COARSE with a smaller time step. Every run is judged against the erfc
    (--tol). The convergence order h -> h/2 is judged only at the times where
    the error is spatial, i.e. where COARSE_DT changes it by < 5 %; at the
    other (early) times the time integration of the jump at t = 0 dominates
    and the order is printed as INFO."""
    ok = True
    errs = []
    for rd in args.runs:
        r = Run(rd)
        Ti = pplib.opt_float(r.opts, "-Twall_Tinit", r.opt("-temp"))
        dT, kap = r.opt("-Twall_dT"), r.opt("-Twall_kappa")
        F0 = r.fields(r.snaps[0][2])
        _, _, a = r.phase(F0[..., 0].mean())
        if abs(kap - a) / a > 1e-6:
            print(f"# {rd}: -Twall_kappa {kap} differs from the material alpha {a}", file=sys.stderr)
            return False
        tf = r.final()[0]
        nx = r.nrb.shape[0] - r.nrb.degree[0]
        # The largest step taken, from SSA_evo.dat. (Flags passed after `--`
        # are not in outp.txt, so read_opts cannot see a -dtmax override.)
        dtmax = float(np.nanmax(pplib.load_ssa(rd)[:, pplib.DT]))
        x = np.linspace(0.0, r.Lx, 601)
        y = np.linspace(0.0, r.Ly, 41)
        e_by_frac = []
        for frac in (1 / 64, 1 / 16, 0.25, 1.0):
            t, st, f = r.nearest(frac * tf)
            T, _ = r.evaluate(r.fields(f), x, y)
            ex = Ti + dT * np.vectorize(math.erfc)(x / (2 * math.sqrt(kap * t)))
            e = float(np.max(np.abs(T - ex[:, None])) / dT)
            e_by_frac.append(e)
            ok &= row("c", f"Nx={nx}_dtmax={dtmax:.1e}_t={t:.4e}", "max_xy|T-T_erfc|/dT", e, 0.0,
                      args.tol, e < args.tol)
        errs.append((nx, e_by_frac))
    if len(errs) == 3:
        (nc, ec), (nf, ef), (_, ed) = errs
        gated = 0
        for i, frac in enumerate((1 / 64, 1 / 16, 0.25, 1.0)):
            sens = abs(ed[i] - ec[i]) / ec[i]
            order = math.log(ec[i] / ef[i]) / math.log(nf / nc)
            case = f"t/t_final={frac:.4f}_dt_sensitivity={sens:.3f}"
            if sens < 0.05:
                gated += 1
                ok &= row("c", case, f"order Nx {nc}->{nf}", order, 2.0, 1.5, order > 1.5)
            else:
                print(f"c,{case},order Nx {nc}->{nf} (time error dominates),{order:.6e},2,,INFO")
        if gated == 0:
            print("c,order,no time with spatial error dominant,0,1,,FAIL")
            ok = False
    return ok


def test_c2(args):
    ok = True
    for rd in args.runs:
        r = Run(rd)
        t, st, f = r.final()
        Ti = pplib.opt_float(r.opts, "-Twall_Tinit", r.opt("-temp"))
        dT, L = r.opt("-Twall_dT"), r.opt("-Twall_length")
        x = np.linspace(0.0, r.Lx, 401)
        y = np.linspace(0.0, r.Ly, 41)
        T, _ = r.evaluate(r.fields(f), x, y)
        ex = Ti + dT * (1.0 - np.minimum(x, L) / L)
        e = float(np.max(np.abs(T - ex[:, None])) / dT)
        ok &= row("c2", f"t={t:.3e}", "max_xy|T-T_lin|/dT", e, 0.0, args.tol, e < args.tol)
    return ok


# ---------------------------------------------------------------- test d ----
def stored_heat(r, F):
    """rho*cp * int T dA - rho_i L int phi dA (per unit depth): the quantity
    whose rate equals the boundary heat input when rho*cp is uniform."""
    _, rc, _ = r.phase(F[..., 0].mean())
    return rc * r.integral(F, 1) - r.opt("-xi_T") * r.opt("-rho_ice") * LAT_SUB * r.integral(F, 0)


def test_d(args):
    """Energy budget. Flux runs: input xi_T q Ly (t - t0) is known exactly.
    Wall-temperature runs: input integrated in time from bflux.dat."""
    ok = True
    for rd in args.runs:
        r = Run(rd)
        F0 = r.fields(r.snaps[0][2])
        E0 = stored_heat(r, F0)
        if pplib.face_list(r.opts, "-qflux_faces"):
            q, t0 = r.opt("-qflux"), r.opt("-qflux_t0", 0.0)
            nface = len(pplib.face_list(r.opts, "-qflux_faces"))
            if nface != 1 or pplib.face_list(r.opts, "-qflux_faces")[0][0] != 0:
                print("# test d expects one x face with flux", file=sys.stderr)
                return False
            tag = "ice" if F0[..., 0].mean() > 0.5 else "air"
            xi = r.opt("-xi_T")
            # The integrator start leaves a fixed offset: TSALPHA starts with
            # U_t = 0 while the flux switches on at t = 0, so the first step
            # stores 20 % more than q Ly dt0 and the next few steps relax that
            # by a fraction. It does not grow, so it is reported (INFO) and the
            # budget is judged against the final input, where it is 2e-5 of it.
            t1, s1, f1 = r.snaps[1]
            off = stored_heat(r, r.fields(f1)) - E0 - xi * q * r.Ly * (t1 - t0)
            print(f"d,flux_{tag}_first_step,offset/(xi q Ly dt0),{off / (xi * q * r.Ly * (t1 - t0)):.6e},0,,INFO")
            tF = r.final()[0]
            QF = xi * q * r.Ly * (tF - t0)
            worst = 0.0
            for t, _st, f in r.snaps[1:]:
                dE = stored_heat(r, r.fields(f)) - E0
                worst = max(worst, abs(dE - xi * q * r.Ly * max(t - t0, 0.0)) / QF)
            ok &= row("d", f"flux_{tag}_all_snapshots", "max_t|dE-Q_in(t)|/Q_in(t_final)", worst, 0.0,
                      args.tol, worst < args.tol)
            dE = stored_heat(r, r.fields(r.final()[2])) - E0
            rel = abs(dE - QF) / QF
            ok &= row("d", f"flux_{tag}_final_t={tF:.4e}", "|dE-Q_in|/Q_in", rel, 0.0, args.tol, rel < args.tol)
        else:
            bf = np.loadtxt(os.path.join(rd, "bflux.dat"))
            hdr = open(os.path.join(rd, "bflux.dat")).read().splitlines()[1].split()[1:]
            cols = [i for i, h in enumerate(hdr) if h.startswith("Q_heat_")]
            tt = bf[:, 0]
            Pin = -bf[:, cols].sum(axis=1)          # bflux is outward; input = -outward
            Qcum = np.concatenate([[0.0], np.cumsum(0.5 * (Pin[1:] + Pin[:-1]) * np.diff(tt))])
            t, st, f = r.final()
            Qin = float(np.interp(t, tt, Qcum))
            dE = stored_heat(r, r.fields(f)) - E0
            ok &= row("d", f"twall_final_t={t:.4e}", "dE [J/m] vs trapezoid int of -k dT/dn", dE, Qin,
                      args.tol, abs(dE - Qin) / abs(Qin) < args.tol)
    return ok


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("test", choices=["a", "b", "c", "c2", "d"])
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--tol", type=float, required=True)
    ap.add_argument("--tol-flux", type=float, default=1e-6)
    args = ap.parse_args()
    fn = {"a": test_a, "b": test_b, "c": test_c, "c2": test_c2, "d": test_d}[args.test]
    sys.exit(0 if fn(args) else 1)


if __name__ == "__main__":
    main()

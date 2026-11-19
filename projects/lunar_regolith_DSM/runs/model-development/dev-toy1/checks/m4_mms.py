"""M4 check: manufactured solution and its source (spec.md M4-9, mms/mms_sources.py).

Claims tested:
  MM1 the closed form of f in M4-9 equals d_tau theta'_m - d_xi^2 theta'_m - Q for the theta'_m of M4-9.
  MM2 theta'_m with that f satisfies the strong form M4-1 (parsed) identically.
  MM3 theta'_m meets S29's end and initial data; the corners (0,0), (1,0) are compatible.
  MM4 A = B = 0 gives f = -Q.
  MM5 the sample table (theta'_m, f) at the stated defaults A, B, Q recomputes.
  MM6 mms_sources.py defines the same theta'_m and f = d_tau - d_xi^2 - Q (read as text, not executed).
  MM7 dimensional form: T = T0 theta'_m(x/L, kappa t/L^2) solves rho c_p T_t = k T_xx + q + q_m with
      q_m = (k T0/L^2) f(x/L, kappa t/L^2).
  MM8 an independent BE / second-difference solver on the MMS problem reproduces the stated errors
      (within 1 %) and observed orders (2 in h, 1 in dtau, within 0.05).

Usage: python m4_mms.py [--spec spec.md] [--ledger ledger.json] [--mms mms_sources.py]
"""
import os
import re
import sys

import numpy as np
import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import m3_common as C  # noqa: E402
import m4_common as M  # noqa: E402

xi, tau = C.S("xi"), C.S("tau")
A, B, Qs = C.S("A"), C.S("B"), C.S("Q")
THP = C.Fn("thetap")
THM = C.Fn("thetam")
F_ = C.Fn("f")


def parse_mms(spec):
    p0 = M.piece_eq(spec.block("M4-9_pre"), lambda p: p.rels == ["="] and isinstance(p.exprs[0], sp.Function)
                    and p.exprs[0].func == THM, "theta'_m = ...")
    th = p0.exprs[1]
    p1 = M.piece_eq(spec.block("M4-9"), lambda p: len(p.exprs) == 3 and p.rels == ["=", "="], "f = ... = ...")
    return th, p1.exprs[1], p1.exprs[2], p1.exprs[0]


def main():
    args = M.cli(__doc__, extra=[("--mms", M.DEFAULT_MMS)])
    spec = M.Spec(args.spec)
    L = M.Ledger(args.ledger)
    rep = M.Reporter("m4_mms.py", L, {"spec": str(spec.path), "spec_sha256": spec.sha256, "mms": args.mms})

    def mm1():
        th, fdef, fclosed, lhs = parse_mms(spec)
        fcalc = sp.diff(th, tau) - sp.diff(th, xi, 2) - Qs
        fdef_e = fdef.replace(THM, sp.Lambda((C.S("xi"), C.S("tau")), th)).doit()
        d1 = sp.simplify(fclosed - fcalc)
        d2 = sp.simplify(fdef_e - fcalc)
        lhs_ok = lhs == F_(xi, tau)
        return d1 == 0 and d2 == 0 and lhs_ok, {"theta'_m": str(th), "f closed form": str(fclosed),
                                                "closed form - computed": str(d1),
                                                "middle (definition) - computed": str(d2)}
    rep.check("MM1_source_closed_form", "the closed form of f in M4-9 equals d_tau th_m - d_xi^2 th_m - Q", mm1)

    def mm2():
        th, _, fclosed, _ = parse_mms(spec)
        blk = spec.block("M4-1")
        rhs = blk.split("\\equiv", 1)[1]
        p = M.piece_eq(rhs, lambda p: p.rels == ["="] and p.exprs[1] == 0 and p.exprs[0].has(THP), "PDE")
        r = p.exprs[0].replace(THP, sp.Lambda((C.S("xi"), C.S("tau")), th)).replace(
            F_, sp.Lambda((C.S("xi"), C.S("tau")), fclosed)).doit()
        d = sp.simplify(r)
        return d == 0, {"M4-1 residual on (theta'_m, f)": str(d)}
    rep.check("MM2_strong_form", "theta'_m with the M4-9 source satisfies the strong form M4-1", mm2)

    def mm3():
        th, _, fclosed, _ = parse_mms(spec)
        s29 = C.pieces(L.stmt("S29"))
        data = {}
        for pc in s29:
            if pc.ok and all(r == "=" for r in pc.rels):
                for e in pc.exprs[:-1]:
                    if isinstance(e, sp.Function) and e.func == THP:
                        data[e.args] = pc.exprs[-1]
        vals = {"left": sp.simplify(th.subs(xi, 0) - data[(0, tau)]),
                "right": sp.simplify(th.subs(xi, 1) - data[(1, tau)]),
                "init": sp.simplify(th.subs(tau, 0) - data[(xi, 0)])}
        corner = [sp.simplify((sp.diff(th, xi, 2) + Qs + fclosed).subs({xi: x0, tau: 0})) for x0 in (0, 1)]
        return all(v == 0 for v in vals.values()) and corner == [0, 0], {
            "theta'_m - S29 data": {k: str(v) for k, v in vals.items()},
            "corner d_xi^2 th + Q + f at (0,0),(1,0)": [str(c) for c in corner]}
    rep.check("MM3_data_and_corners", "theta'_m meets S29's end and initial data and is corner compatible", mm3)

    def mm4():
        _, _, fclosed, _ = parse_mms(spec)
        v = sp.simplify(fclosed.subs({A: 0, B: 0}))
        m = re.search(r"A = B = 0 gives f = (−?Q)", spec.text)
        claimed = M.uni(m.group(1)).subs(sp.Symbol("Q", real=True), Qs) if m else None
        return m is not None and sp.simplify(v - claimed) == 0 and v == -Qs, {"f(A=B=0)": str(v),
                                                                              "claimed": str(claimed)}
    rep.check("MM4_trivial_case", "A = B = 0 gives f = -Q", mm4)

    def mm5():
        th, _, fclosed, _ = parse_mms(spec)
        md = re.search(r"Defaults: A = ([\d.]+), B = ([\d.]+), and Q from the parameters \((\d+) in C1\)", spec.text)
        mq = re.search(r"Samples with A = ([\d.]+), B = ([\d.]+), Q = ([\d.]+):", spec.text)
        if not (md and mq) or md.groups()[:2] != mq.groups()[:2]:
            return False, "defaults / sample parameters not found or inconsistent"
        sub = {A: float(mq.group(1)), B: float(mq.group(2)), Qs: float(mq.group(3))}
        sec = spec.section(r"Manufactured solutions")
        rows = re.findall(r"^\| \(([\d.]+), ([\d.]+)\) \| ([\d.eE+-]+) \| ([\d.eE+-]+) \|$", sec, re.M)
        bad, det = [], []
        for x0, t0, tv, fv in rows:
            pt = {**sub, xi: float(x0), tau: float(t0)}
            tc = float(th.subs(pt))
            fc = float(fclosed.subs(pt))
            det.append({"(xi,tau)": (x0, t0), "theta'_m": [tc, float(tv)], "f": [fc, float(fv)]})
            if abs(tc - float(tv)) > 1e-12 or abs(fc - float(fv)) > 1e-12:
                bad.append((x0, t0))
        return len(rows) >= 2 and not bad, {"parameters": {str(k): v for k, v in sub.items()},
                                            "rows (computed, spec)": det, "failing": bad}
    rep.check("MM5_samples", "the sample table of theta'_m and f recomputes to 1e-12", mm5)

    def mm6():
        th, _, fclosed, _ = parse_mms(spec)
        src = open(args.mms).read()
        m = re.search(r"^theta_m\s*=\s*(.+)$", src, re.M)
        mf = re.search(r"^f\s*=\s*sp\.simplify\(sp\.diff\(theta_m, tau\) - sp\.diff\(theta_m, xi, 2\) - Q\)\s*$",
                       src, re.M)
        loc = {"sp": sp, "xi": xi, "tau": tau, "A": A, "B": B, "Q": Qs}
        expr = m.group(1).replace("sp.", "")
        ths = sp.sympify(expr, locals={"xi": xi, "tau": tau, "A": A, "B": B, "Q": Qs, "pi": sp.pi,
                                        "sin": sp.sin, "exp": sp.exp, "cos": sp.cos})
        d = sp.simplify(ths - th)
        _ = loc
        return d == 0 and mf is not None, {"mms_sources theta_m": str(ths), "minus spec theta'_m": str(d),
                                           "f defined as d_tau - d_xi^2 - Q": mf is not None}
    rep.check("MM6_script_matches_spec", "mms_sources.py builds the same theta'_m and the same f definition", mm6)

    def mm7():
        _, _, fclosed, _ = parse_mms(spec)
        m = re.search(r"q_m\(x, t\) = \((kT_0/L²)\)·f\(x/L, κt/L²\) W m⁻³", spec.text)
        if not m:
            return False, "dimensional MMS sentence not found"
        x, t = sp.symbols("x t", real=True)
        Ls, k, rho, cp, q, T0 = sp.symbols("L k rho c_p q T_0", positive=True)
        fac = M.uni(m.group(1)).subs({sp.Symbol("k", real=True): k, sp.Symbol("T0", real=True): T0,
                                      sp.Symbol("L", real=True): Ls})
        th, _, _, _ = parse_mms(spec)
        kap = k / (rho * cp)
        sub = {xi: x / Ls, tau: kap * t / Ls**2}
        Q_ = q * Ls**2 / (k * T0)
        T = T0 * th.subs(sub)
        qm = fac * fclosed.subs(Qs, Q_).subs(sub)
        r = sp.simplify(rho * cp * sp.diff(T, t) - k * sp.diff(T, x, 2) - q - qm)
        return r == 0, {"factor": str(fac), "rho c_p T_t - k T_xx - q - q_m": str(r)}
    rep.check("MM7_dimensional_source", "q_m = (k T0/L^2) f(x/L, kappa t/L^2) makes T0 theta'_m a solution of the "
              "dimensional equation", mm7)

    def mm8():
        th, _, fclosed, _ = parse_mms(spec)
        sub = {A: 0.5, B: 0.25, Qs: 1.0}
        thf = sp.lambdify((xi, tau), th.subs(sub), "numpy")
        ff = sp.lambdify((xi, tau), fclosed.subs(sub), "numpy")

        def run(N, dtau, tend):
            # independent implementation of (I - dtau D_h) U1 = U0 + dtau (Q + f(.,tau1)) + boundary
            h = 1.0 / N
            x = np.linspace(0, 1, N + 1)[1:-1]
            n = N - 1
            ab = np.zeros((3, n))
            ab[0, 1:] = -dtau / h**2
            ab[1, :] = 1 + 2 * dtau / h**2
            ab[2, :-1] = -dtau / h**2
            from scipy.linalg import solve_banded
            U = np.ones(n)
            nsteps = int(round(tend / dtau))
            for s in range(1, nsteps + 1):
                t1 = s * dtau
                rhs = U + dtau * (1.0 + ff(x, t1))
                rhs[0] += dtau / h**2
                rhs[-1] += dtau / h**2
                U = solve_banded((1, 1), ab, rhs)
            return np.max(np.abs(U - thf(x, nsteps * dtau)))

        sec = spec.section(r"Manufactured solutions")
        sec = sec.translate(str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹", "0123456789"))
        ms = re.search(r"Space, Δτ = 10⁻(\d+), N = (\d+) → (\d+) \| ([\d.e-]+) → ([\d.e-]+) \| ([\d.]+) \|", sec)
        mt = re.search(r"Time, N = (\d+), Δτ = 10⁻(\d+) → (\d+)×10⁻(\d+) \| ([\d.e-]+) → ([\d.e-]+) \| ([\d.]+) \|", sec)
        mr = re.search(r"At the reference N = (\d+), Δτ = 10⁻(\d+) the error is ([\d.e-]+)\.", sec)
        mtend = re.search(r"error = max over nodes at τ_end = ([\d.]+)", sec)
        if not (ms and mt and mr and mtend):
            return False, "order table not found"
        tend = float(mtend.group(1))
        dts = 10.0 ** -int(ms.group(1))
        e1, e2 = run(int(ms.group(2)), dts, tend), run(int(ms.group(3)), dts, tend)
        p_s = np.log(e1 / e2) / np.log(int(ms.group(3)) / int(ms.group(2)))
        Nt = int(mt.group(1))
        d1 = 10.0 ** -int(mt.group(2))
        d2 = int(mt.group(3)) * 10.0 ** -int(mt.group(4))
        f1, f2 = run(Nt, d1, tend), run(Nt, d2, tend)
        p_t = np.log(f1 / f2) / np.log(d1 / d2)
        er = run(int(mr.group(1)), 10.0 ** -int(mr.group(2)), tend)
        claims = {"space e1": (e1, float(ms.group(4))), "space e2": (e2, float(ms.group(5))),
                  "time e1": (f1, float(mt.group(5))), "time e2": (f2, float(mt.group(6))),
                  "reference": (er, float(mr.group(3)))}
        bad = [k for k, (c, s) in claims.items() if not M.rel_ok(c, s, 0.01)]
        if abs(p_s - 2) > 0.05 or abs(p_s - float(ms.group(6))) > 0.01:
            bad.append(f"space order {p_s}")
        if abs(p_t - 1) > 0.05 or abs(p_t - float(mt.group(7))) > 0.01:
            bad.append(f"time order {p_t}")
        return not bad, {"failures": bad, "computed vs spec": claims, "space order": p_s, "time order": p_t,
                         "tau_end": tend}
    rep.check("MM8_discrete_orders", "an independent BE/FD solve of the MMS problem gives the stated errors "
              "(1 %) and orders 2 in h and 1 in dtau", mm8)
    rep.finish()


if __name__ == "__main__":
    main()

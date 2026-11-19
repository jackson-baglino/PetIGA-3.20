"""M4 check: every Jacobian entry is the derivative of its residual entry (spec.md M4-6, M4-6b, M4-7).

Claims tested (all parsed from --spec; code_map.json from --code-map):
  J1  term by term: d(T_k)/dU_{j-1}, dU_j, dU_{j+1} of each residual-table term T1..T5 (with
      U = U^{n+1}, shift = d Udot / d U^{n+1} = 1/dtau) equals the Jacobian-table row of T_k; the
      'dF/dUdot' column is right.
  J2  the assembled-matrix sentence (diagonal shift + 2/h^2, off-diagonals -1/h^2) is the sum of the
      table rows and the symbolic derivative of the M4-6 row (interior and both boundary rows).
  J3  numerical full-Jacobian test on the residual parsed from M4-6 (N = 8 and 64): J v against the
      central difference of F along (v, shift v) for shift in {0, 10, 1e3, 1e5}, split dF/dU,
      dF/dUdot, boundary rows {1, N-1} and interior; and dR/dU^{n+1} of M4-6b for BE. With the MMS
      forcing on, J does not change. Sensitivity: a J without the shift term, and a J with one
      off-diagonal sign flipped, fail the same comparison.
  J4  matrix properties stated in M4-7 and the symbol table: symmetric, tridiagonal, eigenvalues of
      -D_h are (4/h^2) sin^2(l pi h/2), M-matrix (inverse >= 0), diagonal dominance margin = shift.
  J5  TSALPHA statements: shift_V = alpha_m/(alpha_f gamma dt) is d/dX1 of the (1/alpha_f)-scaled
      residual; rho_inf = 0.5 gives alpha_m = 5/6, alpha_f = 2/3, gamma = 2/3; BE is alpha = gamma = 1.
  J6  code_map.json Jacobian rows agree with the spec table.

Usage: python m4_jacobian.py [--spec spec.md] [--ledger ledger.json] [--code-map code_map.json]
"""
import json
import os
import re
import sys

import numpy as np
import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import m3_common as C  # noqa: E402
import m4_common as M  # noqa: E402

xi, tau = C.S("xi"), C.S("tau")
F_ = C.Fn("f")
Um, Uc, Up, Udc, Uoc = (M.U(n) for n in ("U_m", "U_c", "U_p", "Udot_c", "Uold_c"))
hs, dts, Qs, sh = (M.U(n) for n in ("h", "dtau", "Q", "shift"))
gL, gR, fj = M.U("g_L"), M.U("g_R"), M.U("fj")


def row46(spec):
    blk = spec.block("M4-6")
    rhs = re.split(r"F_j\(.*?\)\s*=", blk, maxsplit=1)[1]
    r6 = M.normalise_row(sp.expand(M.parse_latex_expr(rhs)))[0]
    r6 = r6.replace(F_, sp.Lambda((C.S("xi"), C.S("tau")), fj))
    return r6.subs({C.S("h"): hs, C.S("Q"): Qs}).xreplace({M.U("xi_c"): 0})


def res_terms(spec):
    sec = spec.section(r"Residual")

    def cell(tid):
        return spec.table_row(tid, within=sec)[1]

    def conv(s):
        s = s.replace("U_j^n", "Uold_c").replace("U_{j−1}", "U_m").replace("U_{j+1}", "U_p")
        s = s.replace("U_j", "U_c").replace("f(ξ_j, τ^{n+1})", "fj")
        return M.uni(s, extra_names=("Uold_c", "U_m", "U_p", "U_c", "fj"))
    m3 = re.search(r"(−?g_L/h²) in row 1, (−?g_R/h²) in row N−1", cell("T3"))
    return {
        "T1": conv(re.search(r"storage (.*)$", cell("T1")).group(1)),
        "T2": conv(re.search(r"interior neighbours: (.*?), using", cell("T2")).group(1)),
        "T3L": conv(m3.group(1)), "T3R": conv(m3.group(2)),
        "T4": conv(re.search(r"source (.*)$", cell("T4")).group(1)),
        "T5": conv(re.search(r"forcing (.*)$", cell("T5")).group(1)),
    }


def jac_cell(c):
    b = re.search(r"\*\*(.+?)\*\*", c)
    s = b.group(1) if b else c.split("(")[0].strip()
    if s.startswith("none") or s == "":
        return sp.Integer(0)
    return M.uni(s, extra_names=("shift",))


def jac_table(spec):
    sec = spec.section(r"Jacobian")
    out = {}
    for tid in ("T1", "T2", "T3", "T4", "T5"):
        cells = spec.table_row(tid + r"\b.*", within=sec)
        out[tid] = {"m": jac_cell(cells[1]), "c": jac_cell(cells[2]), "p": jac_cell(cells[3]),
                    "dFdUdot": cells[4]}
    return out


def matrix_sentence(spec):
    sec = spec.section(r"Jacobian")
    md = re.search(r"- diagonal: (.+?);", sec)
    mo = re.search(r"- sub-diagonal and super-diagonal: (.+?)\.\s*$", sec, re.M)
    return M.uni(md.group(1), extra_names=("shift",)), M.uni(mo.group(1), extra_names=("shift",))


def main():
    args = M.cli(__doc__, extra=[("--code-map", M.RUN + "/derivation/code_map.json")])
    spec = M.Spec(args.spec)
    L = M.Ledger(args.ledger)
    rep = M.Reporter("m4_jacobian.py", L, {"spec": str(spec.path), "spec_sha256": spec.sha256,
                                           "code_map": args.code_map})
    hu = sp.Symbol("h", real=True)
    dtu = sp.Symbol("dtau", real=True)
    Qu = sp.Symbol("Q", real=True)
    shu = sp.Symbol("shift", real=True)

    def to_m(e):  # unify symbols from uni() with the row symbols
        return e.subs({hu: hs, dtu: dts, Qu: Qs, shu: sh})

    # ------------------------------------------------------------------ J1
    def j1():
        T = res_terms(spec)
        J = jac_table(spec)
        bad, det = [], {}
        Unames = {"m": Um, "c": Uc, "p": Up}
        pairs = [("T1", T["T1"]), ("T2", T["T2"]), ("T3", T["T3L"]), ("T3", T["T3R"]), ("T4", T["T4"]),
                 ("T5", T["T5"])]
        for tid, term in pairs:
            term = to_m(term)
            for k, Uk in Unames.items():
                d = sp.diff(term, Uk)
                d = d.subs(dts, 1 / sh)  # shift = d Udot/dU^{n+1} = 1/dtau
                claimed = to_m(J[tid][k])
                if sp.simplify(d - claimed) != 0:
                    bad.append(f"{tid} d/dU_{k}: computed {d}, table {claimed}")
                det.setdefault(tid, {})[k] = {"computed": str(d), "table": str(claimed)}
        # dF/dUdot column: only T1 (dF_j/dUdot_j = 1)
        Fstor = T["T1"].subs(sp.Symbol("U_c", real=True) - sp.Symbol("Uold_c", real=True), 0)
        r6 = row46(spec)
        dFdUd = sp.diff(r6, Udc)
        col = J["T1"]["dFdUdot"]
        mcol = re.search(r"∂F_j/∂U̇_j = (\S+), times shift", col)
        ok_col = (col.startswith("yes") and mcol is not None and sp.sympify(mcol.group(1)) == dFdUd
                  and all(J[t]["dFdUdot"].startswith("no") for t in ("T2", "T3", "T4", "T5")))
        if not ok_col:
            bad.append(f"dF/dUdot column: spec says {col!r}; computed dF_j/dUdot_j = {dFdUd}")
        _ = Fstor
        return not bad, {"failures": bad, "entries": det, "dF_j/dUdot_j (from M4-6)": str(dFdUd)}
    rep.check("J1_term_by_term", "each Jacobian-table entry of T1..T5 is the derivative of the residual-table "
              "term, with shift = 1/dtau; only T1 has a dF/dUdot entry, equal to 1", j1)

    # ------------------------------------------------------------------ J2
    def j2():
        r6 = row46(spec)
        J = jac_table(spec)
        diag_s, off_s = (to_m(e) for e in matrix_sentence(spec))
        # derivative of R_j = F_j(U^{n+1}, (U^{n+1}-U^n)/dtau) w.r.t. U^{n+1}
        R = r6.subs(Udc, (Uc - Uoc) / dts)
        dR = {k: sp.simplify(sp.diff(R, u).subs(dts, 1 / sh)) for k, u in (("m", Um), ("c", Uc), ("p", Up))}
        # convention J = shift dF/dUdot + dF/dU
        conv = {k: sp.simplify(sh * sp.diff(r6, Udc) * (1 if k == "c" else 0) + sp.diff(r6, u))
                for k, u in (("m", Um), ("c", Uc), ("p", Up))}
        tab_sum = {k: sp.simplify(sum(to_m(J[t][k]) for t in J)) for k in ("m", "c", "p")}
        # boundary rows: U_0 = g_L, U_N = g_R are data -> no column; the remaining entries unchanged
        R1 = R.subs(Um, gL)
        RN = R.subs(Up, gR)
        b1 = [sp.simplify(sp.diff(R1, Uc).subs(dts, 1 / sh)), sp.simplify(sp.diff(R1, Up))]
        bN = [sp.simplify(sp.diff(RN, Um)), sp.simplify(sp.diff(RN, Uc).subs(dts, 1 / sh))]
        ok = (all(sp.simplify(dR[k] - conv[k]) == 0 for k in dR)
              and all(sp.simplify(dR[k] - tab_sum[k]) == 0 for k in dR)
              and sp.simplify(dR["c"] - diag_s) == 0 and sp.simplify(dR["m"] - off_s) == 0
              and sp.simplify(dR["p"] - off_s) == 0
              and sp.simplify(b1[0] - diag_s) == 0 and sp.simplify(b1[1] - off_s) == 0
              and sp.simplify(bN[0] - off_s) == 0 and sp.simplify(bN[1] - diag_s) == 0)
        return ok, {"dR/dU (m,c,p)": {k: str(v) for k, v in dR.items()},
                    "shift dF/dUdot + dF/dU": {k: str(v) for k, v in conv.items()},
                    "sum of table rows": {k: str(v) for k, v in tab_sum.items()},
                    "spec diagonal / off-diagonal": [str(diag_s), str(off_s)],
                    "row 1 (c, p)": [str(e) for e in b1], "row N-1 (m, c)": [str(e) for e in bN]}
    rep.check("J2_matrix", "dR_j/dU^{n+1} of M4-6b equals shift dF/dUdot + dF/dU, the table sum, and the stated "
              "diagonal shift + 2/h^2 and off-diagonals -1/h^2, in interior and boundary rows", j2)

    # ------------------------------------------------------------------ J3
    r6 = row46(spec)
    rowf = sp.lambdify((Um, Uc, Up, Udc, hs, Qs, fj), r6, "numpy")
    diag_e, off_e = (to_m(e) for e in matrix_sentence(spec))
    diagf = sp.lambdify((hs, sh), diag_e, "numpy")
    offf = sp.lambdify((hs, sh), off_e, "numpy")

    def Fvec(U, Ud, N, Q=1.0, fvals=None):
        h = 1.0 / N
        Ufull = np.concatenate(([1.0], U, [1.0]))
        f = np.zeros(N - 1) if fvals is None else fvals
        return np.asarray(rowf(Ufull[:-2], Ufull[1:-1], Ufull[2:], Ud, h, Q, f), dtype=float) * np.ones(N - 1)

    def Jmat(N, shift, mutate=None):
        h = 1.0 / N
        n = N - 1
        d = float(diagf(h, shift)) * np.ones(n)
        o = float(offf(h, shift)) * np.ones(n - 1)
        J = np.diag(d) + np.diag(o, 1) + np.diag(o, -1)
        if mutate == "drop_shift":
            J -= shift * np.eye(n)
        elif mutate == "flip_offdiag":
            J[n // 2, n // 2 + 1] *= -1
        return J

    def compare(N, shift, mutate=None, fvals=None, seed=0, eps=1e-3):
        rng = np.random.default_rng(seed)
        U = 1 + 0.1 * rng.standard_normal(N - 1)
        Ud = rng.standard_normal(N - 1)
        J = Jmat(N, shift, mutate)
        worst = {"full": 0.0, "dF/dU": 0.0, "dF/dUdot": 0.0, "boundary rows": 0.0, "interior rows": 0.0}
        for _ in range(5):
            v = rng.standard_normal(N - 1)
            fd = (Fvec(U + eps * v, Ud + shift * eps * v, N, fvals=fvals)
                  - Fvec(U - eps * v, Ud - shift * eps * v, N, fvals=fvals)) / (2 * eps)
            Jv = J @ v
            sc = np.linalg.norm(fd) + 1e-300
            worst["full"] = max(worst["full"], np.linalg.norm(Jv - fd) / sc)
            bidx = [0, N - 2]
            iidx = list(range(1, N - 2))
            worst["boundary rows"] = max(worst["boundary rows"],
                                         np.linalg.norm((Jv - fd)[bidx]) / (np.linalg.norm(fd[bidx]) + 1e-300))
            if iidx:
                worst["interior rows"] = max(worst["interior rows"],
                                             np.linalg.norm((Jv - fd)[iidx]) / (np.linalg.norm(fd[iidx]) + 1e-300))
            # parts
            fdU = (Fvec(U + eps * v, Ud, N, fvals=fvals) - Fvec(U - eps * v, Ud, N, fvals=fvals)) / (2 * eps)
            JU = Jmat(N, 0.0, None if mutate == "drop_shift" else mutate) @ v
            worst["dF/dU"] = max(worst["dF/dU"], np.linalg.norm(JU - fdU) / (np.linalg.norm(fdU) + 1e-300))
            fdUd = (Fvec(U, Ud + eps * v, N, fvals=fvals) - Fvec(U, Ud - eps * v, N, fvals=fvals)) / (2 * eps)
            # dF/dUdot implied by the spec's matrix: (J(shift) - J(0))/shift ; for shift=0 use shift=1
            s1 = shift if shift != 0 else 1.0
            JUd = ((Jmat(N, s1, mutate) - Jmat(N, 0.0, mutate if mutate != "drop_shift" else None)) / s1) @ v
            worst["dF/dUdot"] = max(worst["dF/dUdot"],
                                    np.linalg.norm(JUd - fdUd) / (np.linalg.norm(fdUd) + 1e-300))
        return worst

    TOL = 1e-8  # F is affine in (U, Udot): the central difference is exact up to roundoff

    def j3():
        res, bad = {}, []
        for N in (8, 64):
            for shift in (0.0, 1e1, 1e3, 1e5):
                w = compare(N, shift)
                res[f"N={N} shift={shift:g}"] = w
                if max(w.values()) > TOL:
                    bad.append(f"N={N} shift={shift:g}")
            # MMS forcing on: J unchanged
            fv = np.sin(np.arange(1, N)) * 3.0
            w = compare(N, 1e3, fvals=fv)
            res[f"N={N} shift=1e3 MMS on"] = w
            if max(w.values()) > TOL:
                bad.append(f"N={N} MMS on")
            # BE: R(U1) = F(U1, (U1-U0)/dtau); FD of R vs J(shift = 1/dtau)
            dtau = 1e-3
            rng = np.random.default_rng(3)
            U0 = 1 + 0.1 * rng.standard_normal(N - 1)
            U1 = U0 + 0.01 * rng.standard_normal(N - 1)
            v = rng.standard_normal(N - 1)
            R = lambda U: Fvec(U, (U - U0) / dtau, N)  # noqa: E731
            fd = (R(U1 + 1e-3 * v) - R(U1 - 1e-3 * v)) / 2e-3
            e = np.linalg.norm(Jmat(N, 1 / dtau) @ v - fd) / np.linalg.norm(fd)
            res[f"N={N} BE dR/dU1, dtau=1e-3"] = e
            if e > TOL:
                bad.append(f"N={N} BE")
        # sensitivity
        sens = {}
        for mut in ("drop_shift", "flip_offdiag"):
            w = compare(64, 1e1, mutate=mut)
            sens[mut] = w["full"]
            if w["full"] <= TOL * 100:
                bad.append(f"sensitivity {mut} not detected")
        return not bad, {"failures": bad, "tolerance": TOL, "worst relative errors": res,
                         "sensitivity (full rel. error of mutated J, must be large)": sens}
    rep.check("J3_full_jacobian_numeric", "J v equals the central difference of F along (v, shift v) for shifts "
              "0, 10, 1e3, 1e5 in every group, with MMS on, and for the BE residual; mutated J fail", j3)

    # ------------------------------------------------------------------ J4
    def j4():
        mlam = re.search(r"\| λ_ℓ \| eigenvalues of −D_h with zero end values: (\(4/h²\) sin²\(ℓπh/2\))", spec.text)
        if not mlam:
            return False, "eigenvalue formula not found"
        out, bad = {}, []
        for N in (8, 64):
            h = 1.0 / N
            A = Jmat(N, 0.0)
            ev = np.sort(np.linalg.eigvalsh(A))
            ell = np.arange(1, N)
            lam = np.sort((4 / h**2) * np.sin(ell * np.pi * h / 2) ** 2)
            e_ev = np.max(np.abs(ev - lam) / lam)
            for shift in (0.0, 10.0, 1e3):
                Js = Jmat(N, shift)
                sym = np.allclose(Js, Js.T)
                tri = np.allclose(np.triu(Js, 2), 0) and np.allclose(np.tril(Js, -2), 0)
                inv = np.linalg.inv(Js)
                mmat = inv.min() >= -1e-12 and np.all(np.diag(Js) > 0) and np.all(Js - np.diag(np.diag(Js)) <= 0)
                margin = np.min(np.abs(np.diag(Js)) - (np.abs(Js).sum(1) - np.abs(np.diag(Js))))
                out[f"N={N} shift={shift:g}"] = {"symmetric": bool(sym), "tridiagonal": bool(tri),
                                                 "M-matrix": bool(mmat), "min dominance margin": float(margin)}
                if not (sym and tri and mmat and margin >= shift - 1e-9 * (1 + shift)):
                    bad.append(f"N={N} shift={shift:g}")
            out[f"N={N} eig rel err"] = float(e_ev)
            if e_ev > 1e-10:
                bad.append(f"N={N} eigenvalues")
        return not bad, {"failures": bad, **out}
    rep.check("J4_matrix_properties", "J is symmetric tridiagonal, a nonsingular M-matrix with dominance margin "
              ">= shift, and -D_h has eigenvalues (4/h^2) sin^2(l pi h/2)", j4)

    # ------------------------------------------------------------------ J5
    def j5():
        sec = spec.section(r"Anchor solver")
        am, af, g, dt, rinf = sp.symbols("alpha_m alpha_f gamma Delta_t rho_inf", positive=True)
        X1, Xn, Vn = sp.symbols("X1 Xn Vn")
        V1 = (X1 - Xn) / (g * dt) - (1 - g) / g * Vn
        Xa = Xn + af * (X1 - Xn)
        Va = Vn + am * (V1 - Vn)
        a, b = sp.symbols("a b")
        Fsym = sp.Function("F")
        Fexpr = (1 / af) * Fsym(Xa, Va)
        dF = sp.diff(Fexpr, X1)
        FU = sp.Subs(sp.Derivative(Fsym(a, Va), a), a, Xa).doit()
        FV = sp.Subs(sp.Derivative(Fsym(Xa, b), b), b, Va).doit()
        ms = re.search(r"shift_V = (α_m/\(α_f γ Δt\))", spec.text)
        shiftV = am / (af * g * dt)
        chain_ok = sp.simplify(dF - (FU + shiftV * FV)) == 0
        dXa = sp.diff(Xa, X1)
        dVa = sp.diff(Va, X1)
        mdx = re.search(r"dX_a/dX_1 = α_f and dV_a/dX_1 = α_m/\(γΔt\)", sec)
        # rho_inf = 1/2 values
        mv = re.search(r"α_m = 0\.5\(3 − ρ∞\)/\(1 \+ ρ∞\) = (\d+/\d+), α_f = 1/\(1 \+ ρ∞\) = (\d+/\d+), and\s+γ = 1/2 \+ "
                       r"α_m − α_f = (\d+/\d+)", sec)
        r = sp.Rational(1, 2)
        amv = sp.Rational(1, 2) * (3 - r) / (1 + r)
        afv = 1 / (1 + r)
        gv = sp.Rational(1, 2) + amv - afv
        vals_ok = mv is not None and [sp.Rational(x) for x in mv.groups()] == [amv, afv, gv]
        be_shift = shiftV.subs({am: 1, af: 1, g: 1})
        return (chain_ok and ms is not None and mdx is not None and sp.simplify(dXa - af) == 0
                and sp.simplify(dVa - am / (g * dt)) == 0 and vals_ok and be_shift == 1 / dt), {
            "d/dX1 (1/af) F(Xa,Va) - (F_U + shift_V F_V)": str(sp.simplify(dF - (FU + shiftV * FV))),
            "dXa/dX1": str(dXa), "dVa/dX1": str(sp.simplify(dVa)),
            "rho_inf=0.5: (am, af, gamma) computed": [str(amv), str(afv), str(gv)],
            "spec values": list(mv.groups()) if mv else None, "BE shift": str(be_shift)}
    rep.check("J5_tsalpha", "shift_V = alpha_m/(alpha_f gamma dt) is the derivative of the 1/alpha_f-scaled "
              "residual; rho_inf = 0.5 gives (5/6, 2/3, 2/3); alpha = gamma = 1 is BE with shift 1/dt", j5)

    # ------------------------------------------------------------------ J6
    def j6():
        cm = json.load(open(args.code_map))
        terms = {e["term"]: e for e in cm}
        bad = []
        t1 = [t for t in terms if t.startswith("J from T1")]
        t2 = [t for t in terms if t.startswith("J from T2")]
        t345 = [t for t in terms if re.match(r"J from T[345]", t)]
        if not (t1 and "shift = 1/dtau" in t1[0] and "dF/dUdot = I" in t1[0]):
            bad.append(f"T1 entry: {t1}")
        m2 = re.search(r"diagonal ([+-]?\d+)/h\^2, sub- and super-diagonal ([+-]?\d+)/h\^2", t2[0]) if t2 else None
        if not (m2 and int(m2.group(1)) == 2 and int(m2.group(2)) == -1):
            bad.append(f"T2 entry: {t2}")
        if not (len(t345) == 3 and all(("no entry" in t or "zero" in t) for t in t345)):
            bad.append(f"T3-T5 entries: {t345}")
        resid = [t for t in terms if "(residual)" in t]
        ids = sorted(re.match(r"(T\d)", t).group(1) for t in resid)
        if ids != ["T1", "T2", "T3", "T4", "T5"]:
            bad.append(f"residual terms in code_map: {ids}")
        return not bad, {"failures": bad, "J rows": t1 + t2 + t345}
    rep.check("J6_code_map", "code_map.json lists T1..T5 and their Jacobian rows as the spec does", j6)
    rep.finish()


if __name__ == "__main__":
    main()

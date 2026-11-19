"""M4 check: strong form <-> weak form <-> FD rows <-> residual (spec.md M4-1 ... M4-6b).

Claims tested (each parsed from --spec, the strong form from --ledger S29):
  W1  M4-1 with f = 0 is S29 term by term; same Dirichlet and initial data.
  W2  M4-2 integrated by parts against w in H^1_0 gives back the M4-1 residual times w, with no
      leftover boundary term (so the ends are essential only, no natural condition).
  W3  the integration-by-parts identity written in the text has the right sign and bracket.
  W4  P1 Galerkin + lumped mass + nodal load, divided by h, gives the M4-3 row; stiffness and mass
      entries stated in the text.
  W5  D_h consistency: D_h u = u'' + (h^2/12) u'''' + O(h^4) as the text's (h^2/24)[..] says; exact on
      quadratics.
  W6  M4-4 is M4-3 under backward Euler.
  W7  M4-5 is M4-3 with Udot = 0, f = 0; the claimed nodal parabola solves M4-5 exactly.
  W8  M4-6 is M4-3 (index j), and M4-6b composed with M4-6 is M4-4.
  W9  the residual table T1..T5 adds up to the M4-4 row, for interior rows and for rows 1, N-1
      with the Dirichlet data.
  W10 units: (k T0/L^2) x M4-1 is the dimensional heat row rho c_p T_t - k T_xx - q - q_m in W m^-3.

Usage: python m4_weak_strong.py [--spec spec.md] [--ledger ledger.json]
"""
import os
import re
import sys

import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import m3_common as C  # noqa: E402
import m4_common as M  # noqa: E402

xi, tau = C.S("xi"), C.S("tau")
Qs = C.S("Q")
h = C.S("h")
THP = C.Fn("thetap")
F_ = C.Fn("f")
W_ = C.Fn("w")


def strong_m41(spec):
    blk = spec.block("M4-1")
    rhs = blk.split("\\equiv", 1)[1]
    p = M.piece_eq(rhs, lambda p: p.rels == ["="] and p.exprs[1] == 0 and p.exprs[0].has(THP),
                   "PDE residual = 0")
    return sp.expand(p.exprs[0]), rhs


def data_pieces(stmt_pre_text, parser):
    """Return {('end',0):val, ('end',1):val, ('init',):val} for theta' data in a statement."""
    out = {}
    for p in parser(stmt_pre_text):
        if not p.ok or any(r != "=" for r in p.rels):
            continue
        for e in p.exprs[:-1]:
            if isinstance(e, sp.core.function.AppliedUndef) and e.func == THP:
                a, b = e.args
                if a == 0 and b == tau:
                    out["left"] = p.exprs[-1]
                elif a == 1 and b == tau:
                    out["right"] = p.exprs[-1]
                elif a == xi and b == 0:
                    out["init"] = p.exprs[-1]
    return out


def main():
    args = M.cli(__doc__)
    spec = M.Spec(args.spec)
    L = M.Ledger(args.ledger)
    rep = M.Reporter("m4_weak_strong.py", L, {"spec": str(spec.path), "spec_sha256": spec.sha256})

    # ------------------------------------------------------------------ W1
    def w1():
        r41, rhs = strong_m41(spec)
        s29 = C.piece(L.stmt("S29"), 0).residual()
        r41_f0 = sp.expand(r41.subs(F_(xi, tau), 0))
        same = sp.simplify(r41_f0 - s29) == 0
        f_term = sp.expand(r41 - r41_f0)
        f_ok = sp.simplify(f_term + F_(xi, tau)) == 0
        d41 = data_pieces(rhs, M.latex_pieces)
        d29 = data_pieces(L.stmt("S29"), lambda s: C.pieces(s))
        data_ok = d41 == d29 and set(d29) == {"left", "right", "init"}
        return same and f_ok and data_ok, {
            "M4-1 residual": M.C.short(r41), "S29 residual": M.C.short(s29),
            "difference with f=0": M.C.short(sp.simplify(r41_f0 - s29)),
            "forcing term (must be -f)": M.C.short(f_term),
            "M4-1 data": {k: str(v) for k, v in d41.items()}, "S29 data": {k: str(v) for k, v in d29.items()}}
    rep.check("W1_m41_vs_S29", "M4-1 with f = 0 equals S29 term by term, with the same end and initial data", w1)

    # ------------------------------------------------------------------ W2
    a1, a2, a3, a4 = sp.symbols("a1:5", real=True)
    theta_test = sp.Lambda((C.S("xi"), C.S("tau")),
                           1 + a1 * sp.sin(sp.pi * xi) * tau + a2 * xi**2 * (1 - xi) * sp.exp(tau)
                           + a3 * xi**3 * tau**2 + a4 * sp.cos(3 * xi) * sp.sin(tau))
    f_test = sp.Lambda((C.S("xi"), C.S("tau")), sp.exp(xi) * tau + xi**2 * sp.cos(tau))
    w_h10 = [xi * (1 - xi), xi**2 * (1 - xi), sp.sin(sp.pi * xi), sp.sin(2 * sp.pi * xi) * (1 + xi),
             xi * (1 - xi) * sp.exp(xi)]
    w_free = [sp.Integer(1), xi, sp.exp(xi)]

    def plug_all(expr, wexpr):
        e = expr.replace(THP, theta_test).replace(F_, f_test)
        e = e.replace(W_, sp.Lambda(C.S("xi"), wexpr))
        # evaluate derivatives only; integrals are left for numerical quadrature
        return e.replace(lambda z: isinstance(z, sp.Derivative), lambda z: z.doit())

    def w2():
        blk = spec.block("M4-2")
        p = M.piece_eq(blk, lambda p: p.rels == ["="] and p.exprs[0].has(sp.Integral), "integral = 0")
        weak = p.exprs[0] - p.exprs[1]
        r41, _ = strong_m41(spec)
        bad, rows = [], []
        # numeric parameter values (a generic point) and quadrature with mpmath via evalf
        pts = [{a1: 0.7, a2: -1.3, a3: 0.45, a4: 2.1, Qs: 0.83, tau: 0.37},
               {a1: -0.2, a2: 0.9, a3: -1.7, a4: 0.3, Qs: -1.6, tau: 1.9}]
        for wexpr in w_h10 + w_free:
            for pt in pts:
                Iw = sp.N(plug_all(weak, wexpr).subs(pt), 30)
                integrand = (plug_all(r41, wexpr) * wexpr).subs(pt)
                strong_int = sp.N(sp.Integral(integrand, (xi, 0, 1)), 30)
                th = theta_test(xi, tau)
                bt = (sp.diff(th, xi) * wexpr)
                bterm = sp.N((bt.subs(xi, 1) - bt.subs(xi, 0)).subs(pt), 30)
                d = float(Iw - strong_int - bterm)
                scale = float(abs(Iw) + abs(strong_int) + abs(bterm) + 1)
                vanishes = wexpr.subs(xi, 0) == 0 and wexpr.subs(xi, 1) == 0
                rows.append({"w": str(wexpr), "point": {str(k): v for k, v in pt.items()},
                             "weak - int(strong*w) - [th_xi w]": d, "w in H10": vanishes,
                             "boundary term": float(bterm)})
                if abs(d) > 1e-12 * scale:  # quadrature roundoff; a wrong term gives O(1)
                    bad.append(str(wexpr))
        # sensitivity: the same comparison on a weak form with a flipped source sign must fail
        wm = w_h10[0]
        pt = pts[0]
        Imut = sp.N(plug_all(weak.subs(Qs, -Qs), wm).subs(pt), 30)
        Sint = sp.N(sp.Integral((plug_all(r41, wm) * wm).subs(pt), (xi, 0, 1)), 30)
        sens = abs(float(Imut - Sint))
        if sens < 1e-6:
            bad.append("sensitivity: flipped source sign not detected")
        # with the boundary term zero for H10 test functions, weak == integral(strong * w)
        return not bad, {"failing w": bad, "sensitivity (flipped Q) difference": sens, "per w": rows,
                         "note": "difference must be exactly the integration-by-parts bracket; for w in H10 "
                                 "the bracket is 0, so weak form == int(M4-1 residual * w) and, w being "
                                 "arbitrary, the strong form is recovered; no boundary integral is left"}
    rep.check("W2_weak_to_strong", "integrating M4-2 by parts gives the M4-1 residual times w, with no boundary "
              "term left for w in H^1_0", w2)

    def w2b():
        # Essential conditions: the lift value in "theta' - 1 in H^1_0" equals S29's end data
        m = re.search(r"find θ′\(·, τ\) with θ′ − (\S+) ∈ H¹₀\(0, 1\)", spec.text)
        d29 = data_pieces(L.stmt("S29"), lambda s: C.pieces(s))
        lift = sp.sympify(m.group(1)) if m else None
        ok = m is not None and lift == d29.get("left") == d29.get("right")
        no_flux = "no natural (flux) term" in spec.text
        return ok and no_flux, {"lift": str(lift), "S29 end data": {k: str(v) for k, v in d29.items()},
                                "text says no natural term": no_flux}
    rep.check("W2b_essential_bc", "the weak form's essential condition theta' - 1 in H^1_0 is S29's Dirichlet "
              "data and no natural condition is stated", w2b)

    # ------------------------------------------------------------------ W3
    def w3():
        m = re.search(r"Integration by parts gives (−?)∫₀¹ ∂_ξ²θ′ w dξ = ∫₀¹ ∂_ξθ′ ∂_ξw dξ ([−+]) \[∂_ξθ′ w\]₀¹",
                      spec.text)
        if not m:
            return False, "identity text not found in the expected form"
        sl = -1 if m.group(1) == "−" else 1
        sb = -1 if m.group(2) == "−" else 1
        pt = {a1: 0.7, a2: -1.3, a3: 0.45, a4: 2.1, tau: 0.37}
        th = theta_test(xi, tau).subs(pt)
        bad = []
        for wexpr in w_h10 + w_free:
            lhs = sl * sp.N(sp.Integral(sp.diff(th, xi, 2) * wexpr, (xi, 0, 1)), 30)
            rhs = sp.N(sp.Integral(sp.diff(th, xi) * sp.diff(wexpr, xi), (xi, 0, 1)), 30) + sb * sp.N(
                (sp.diff(th, xi) * wexpr).subs(xi, 1) - (sp.diff(th, xi) * wexpr).subs(xi, 0), 30)
            if abs(float(lhs - rhs)) > 1e-12 * (1 + abs(float(lhs))):
                bad.append(str(wexpr))
        return not bad, {"signs (lhs, bracket)": [sl, sb], "failing w": bad}
    rep.check("W3_ibp_identity", "the stated identity -int th'' w = int th' w' - [th' w] holds", w3)

    # ------------------------------------------------------------------ W4
    def w4():
        Um, Uc, Up, Udc = (M.U(n) for n in ("U_m", "U_c", "U_p", "Udot_c"))
        xc = M.U("xi_c")
        # hat functions on [xc-h, xc+h]; derivatives on each element
        s = sp.Symbol("s", real=True)
        phi_c_l, phi_c_r = (s - (xc - h)) / h, ((xc + h) - s) / h
        phi_m_l = ((xc) - s) / h
        phi_p_r = (s - xc) / h
        dl = lambda e: sp.diff(e, s)  # noqa: E731
        Kcc = sp.integrate(dl(phi_c_l)**2, (s, xc - h, xc)) + sp.integrate(dl(phi_c_r)**2, (s, xc, xc + h))
        Kcm = sp.integrate(dl(phi_c_l) * dl(phi_m_l), (s, xc - h, xc))
        Kcp = sp.integrate(dl(phi_c_r) * dl(phi_p_r), (s, xc, xc + h))
        mass_cons = sp.integrate(phi_c_l, (s, xc - h, xc)) + sp.integrate(phi_c_r, (s, xc, xc + h))
        load_Q = Qs * mass_cons  # exact for constant Q
        fval = F_(xc, tau)
        row = (mass_cons * Udc + Kcm * Um + Kcc * Uc + Kcp * Up - load_Q - mass_cons * fval) / h
        row = sp.simplify(sp.expand(row))
        parsed = M.piece_eq(spec.block("M4-3"), lambda p: p.rels == ["="] and p.exprs[1] == 0 and
                            p.exprs[0].has(Qs), "row = 0")
        prow, letters = M.normalise_row(sp.expand(parsed.exprs[0]))
        diff = sp.simplify(sp.expand(prow - row))
        # text claims
        mt = re.search(r"∫φ_i′φ_i′ = (\S+) and ∫φ_i′φ_\{i±1\}′ = (\S+?)\.", spec.text)
        mm = re.search(r"\*\*Lumped mass\.\*\* (\S+) on the diagonal", spec.text)
        tk = (M.uni(mt.group(1)), M.uni(mt.group(2))) if mt else (None, None)
        hsym = sp.Symbol("h", real=True)
        text_ok = (mt is not None and mm is not None and
                   sp.simplify(tk[0] - Kcc.subs(h, hsym)) == 0 and sp.simplify(tk[1] - Kcm.subs(h, hsym)) == 0
                   and sp.simplify(tk[1] - Kcp.subs(h, hsym)) == 0
                   and sp.simplify(M.uni(mm.group(1)) - mass_cons.subs(h, hsym)) == 0)
        return diff == 0 and text_ok, {"Galerkin row / h": M.C.short(row), "M4-3 row": M.C.short(prow),
                                       "difference": str(diff), "index letters": sorted(letters),
                                       "K_ii, K_i,i-1, K_i,i+1": [str(Kcc), str(Kcm), str(Kcp)],
                                       "row-sum mass": str(mass_cons),
                                       "text stiffness": [str(t) for t in tk],
                                       "text mass": mm.group(1) if mm else None}
    rep.check("W4_galerkin_to_fd", "P1 Galerkin with lumped mass and nodal load, divided by h, is the M4-3 row; "
              "stated stiffness 2/h, -1/h and lumped mass h are right", w4)

    # ------------------------------------------------------------------ W5
    def w5():
        m = re.search(r"D_h u = u″ \+ \(h²/(\d+)\)\[u⁗\(η₊\) \+ u⁗\(η₋\)\]", spec.text)
        if not m:
            return False, "consistency formula not found"
        c = int(m.group(1))
        x0 = sp.Symbol("x0", real=True)
        u = sp.Function("u")
        ser = sp.series((u(x0 - h) - 2 * u(x0) + u(x0 + h)) / h**2, h, 0, 4).removeO().doit()
        lead = sp.simplify(ser - sp.diff(u(x0), x0, 2))
        claimed = sp.Rational(2, c) * h**2 * sp.diff(u(x0), x0, 4)  # eta+ = eta- -> x0 as h -> 0
        ok_series = sp.simplify(lead - claimed) == 0
        p0, p1, p2 = sp.symbols("p0:3")
        quad = p0 + p1 * x0 + p2 * x0**2
        Dq = sp.simplify((quad.subs(x0, x0 - h) - 2 * quad + quad.subs(x0, x0 + h)) / h**2)
        return ok_series and sp.simplify(Dq - sp.diff(quad, x0, 2)) == 0, {
            "D_h u - u'' (to h^2)": str(lead), "claimed (h^2/%d)*2 u''''" % c: str(claimed),
            "D_h on quadratic minus u''": str(sp.simplify(Dq - sp.diff(quad, x0, 2)))}
    rep.check("W5_consistency", "D_h u = u'' + O(h^2) with the stated h^2/24 bracket; exact on quadratics", w5)

    # ------------------------------------------------------------------ helpers for rows
    def row43():
        p = M.piece_eq(spec.block("M4-3"), lambda p: p.rels == ["="] and p.exprs[1] == 0 and
                       p.exprs[0].has(Qs), "row = 0")
        return M.normalise_row(sp.expand(p.exprs[0]))[0]

    def row44():
        p = M.piece_eq(spec.block("M4-4"), lambda p: p.rels == ["="] and p.exprs[1] == 0 and
                       p.exprs[0].has(Qs), "row = 0")
        return M.normalise_row(sp.expand(p.exprs[0]))[0]

    dtau, taunew = C.S("dtau"), C.S("taunew")
    names = ("m", "c", "p")

    def be(expr):
        rep_ = {M.U("Udot_c"): (M.U("Unew_c") - M.U("Uold_c")) / dtau}
        for n in names:
            rep_[M.U("U_" + n)] = M.U("Unew_" + n)
        e = expr.xreplace(rep_)
        return e.replace(F_, sp.Lambda((C.S("xi"), C.S("tau")), F_(C.S("xi"), taunew)))

    # ------------------------------------------------------------------ W6
    def w6():
        r3, r4 = row43(), row44()
        d = sp.simplify(sp.expand(be(r3) - r4))
        return d == 0, {"BE(M4-3)": M.C.short(be(r3)), "M4-4": M.C.short(r4), "difference": str(d)}
    rep.check("W6_backward_euler", "M4-4 is M4-3 with Udot -> (U^{n+1}-U^n)/dtau and all else at tau^{n+1}", w6)

    # ------------------------------------------------------------------ W7
    def w7():
        p = M.piece_eq(spec.block("M4-5"), lambda p: p.rels == ["="] and p.exprs[1] == 0 and
                       p.exprs[0].has(Qs), "row = 0")
        r5 = M.normalise_row(sp.expand(p.exprs[0]))[0]
        expect = sp.expand(row43().xreplace({M.U("Udot_c"): 0}).replace(F_, sp.Lambda((C.S("xi"), C.S("tau")), 0)))
        d = sp.simplify(r5 - expect)
        # claimed nodal solution
        m = re.search(r"U_j = (1 \+ Qξ_j\(1 − ξ_j\)/2)", spec.text)
        if not m:
            return False, "steady nodal formula not found"
        Ujf = M.uni(m.group(1), extra_names=("xi_j",))
        xj = sp.Symbol("xi_j", real=True)
        Qr = sp.Symbol("Q", real=True)
        nod = {M.U("U_m"): Ujf.subs(xj, xi - h), M.U("U_c"): Ujf.subs(xj, xi),
               M.U("U_p"): Ujf.subs(xj, xi + h)}
        r5n = sp.simplify(r5.xreplace(nod).subs(Qr, Qs))
        ends = [sp.simplify(Ujf.subs(xj, 0) - 1), sp.simplify(Ujf.subs(xj, 1) - 1)]
        return d == 0 and r5n == 0 and ends == [0, 0], {
            "M4-5 minus M4-3(Udot=0,f=0)": str(d), "M4-5 on claimed nodal values": str(r5n),
            "claimed U_j": str(Ujf), "end values minus 1": [str(e) for e in ends]}
    rep.check("W7_steady_rows", "M4-5 is M4-3 with Udot = 0, f = 0, and U_j = 1 + Q xi_j(1-xi_j)/2 solves it "
              "exactly with U_0 = U_N = 1", w7)

    # ------------------------------------------------------------------ W8
    def w8():
        blk = spec.block("M4-6")
        rhs = re.split(r"F_j\(.*?\)\s*=", blk, maxsplit=1)[1]
        r6 = M.normalise_row(sp.expand(M.parse_latex_expr(rhs)))[0]
        d36 = sp.simplify(r6 - row43())
        blk_b = M.latex_pre(spec.block("M4-6b"))
        mb = re.search(r"R_j\(\\Unew\s*\)\s*=\s*F_j\\Big\((.*)\\Big\)", blk_b, re.S)
        if not mb:
            return False, {"M4-6 vs M4-3": str(d36), "error": "M4-6b form not found"}
        argstr = [a.strip() for a in re.split(r",\\,", mb.group(1))]
        a_tau = C._P(C.tokenize(argstr[0])).full_expr()
        a_U = C._P(C.tokenize(argstr[1])).full_expr()
        a_Ud = C._P(C.tokenize(argstr[2])).full_expr()
        comp = {}
        for n in names:
            comp[M.U("U_" + n)] = a_U.xreplace({C.S("Unew"): M.U("Unew_" + n), C.S("Uold"): M.U("Uold_" + n)})
        comp[M.U("Udot_c")] = a_Ud.xreplace({C.S("Unew"): M.U("Unew_c"), C.S("Uold"): M.U("Uold_c")})
        r6b = r6.xreplace(comp).replace(F_, sp.Lambda((C.S("xi"), C.S("tau")), F_(C.S("xi"), a_tau)))
        d6b = sp.simplify(sp.expand(r6b - row44()))
        return d36 == 0 and d6b == 0, {"M4-6 minus M4-3": str(d36), "M4-6b args": [str(a_tau), str(a_U), str(a_Ud)],
                                       "R_j minus M4-4": str(d6b)}
    rep.check("W8_residual_def", "M4-6 equals the M4-3 row and R_j of M4-6b equals the M4-4 row", w8)

    # ------------------------------------------------------------------ W9
    def w9():
        sec = spec.section(r"Residual")

        def cell(tid):
            return spec.table_row(tid, within=sec)[1]

        def conv(s):
            s = s.replace("U_j^n", "Uold_c").replace("U_{j−1}", "U_m").replace("U_{j+1}", "U_p")
            s = s.replace("U_j", "U_c").replace("f(ξ_j, τ^{n+1})", "fj")
            return M.uni(s, extra_names=("Uold_c", "U_m", "U_p", "U_c", "fj"))
        t1 = conv(re.search(r"storage (.*)$", cell("T1")).group(1))
        t2 = conv(re.search(r"interior neighbours: (.*?), using", cell("T2")).group(1))
        m3 = re.search(r"(−?g_L/h²) in row 1, (−?g_R/h²) in row N−1", cell("T3"))
        t3L, t3R = conv(m3.group(1)), conv(m3.group(2))
        t4 = conv(re.search(r"source (.*)$", cell("T4")).group(1))
        t5 = conv(re.search(r"forcing (.*)$", cell("T5")).group(1))
        Um, Uc, Up, fj = (M.U(n) for n in ("U_m", "U_c", "U_p", "fj"))
        gL, gR = M.U("g_L"), M.U("g_R")
        r4 = row44()
        r4 = r4.xreplace({M.U("Unew_" + n): M.U("U_" + n) for n in names})
        r4 = r4.replace(F_, sp.Lambda((C.S("xi"), C.S("tau")), fj))
        r4 = r4.xreplace({M.U("xi_c"): 0}).subs({dtau: M.U("dtau"), h: M.U("h"), Qs: M.U("Q")})
        tab = t1 + t2 + t4 + t5
        d_int = sp.simplify(sp.expand(tab - r4))
        d_row1 = sp.simplify(sp.expand((t1 + t2.subs(Um, 0) + t3L + t4 + t5) - r4.subs(Um, gL)))
        d_rowN = sp.simplify(sp.expand((t1 + t2.subs(Up, 0) + t3R + t4 + t5) - r4.subs(Up, gR)))
        return d_int == 0 and d_row1 == 0 and d_rowN == 0, {
            "T1..T5": [str(t) for t in (t1, t2, t3L, t3R, t4, t5)],
            "interior: table - M4-4": str(d_int), "row 1: table - M4-4(U_0=g_L)": str(d_row1),
            "row N-1: table - M4-4(U_N=g_R)": str(d_rowN)}
    rep.check("W9_term_table", "the residual table T1..T5 adds up to the M4-4 row in interior rows and in rows 1, "
              "N-1 with the Dirichlet data", w9)

    # ------------------------------------------------------------------ W10
    def w10():
        import pint
        ureg = pint.UnitRegistry()
        x, t = sp.symbols("x t", real=True)
        Ls, k, rho, cp, q, T0 = sp.symbols("L k rho c_p q T_0", positive=True)
        kap = k / (rho * cp)
        m = re.search(r"The dimensional row is \((k T_0/L²)\) × M4-1 in\s*W m⁻³", spec.text)
        if not m:
            return False, "dimensional factor sentence not found"
        fac = M.uni(m.group(1).replace(" ", ""))
        fac = fac.subs({sp.Symbol("k", real=True): k, sp.Symbol("T0", real=True): T0,
                        sp.Symbol("L", real=True): Ls})
        r41, _ = strong_m41(spec)
        th = sp.sin(2 * xi) * sp.exp(-tau) + xi**3 * tau + 1
        fx = sp.cos(xi) * tau**2
        r_nd = r41.replace(THP, sp.Lambda((C.S("xi"), C.S("tau")), th)).replace(
            F_, sp.Lambda((C.S("xi"), C.S("tau")), fx)).doit()
        r_nd = r_nd.subs(Qs, q * Ls**2 / (k * T0))
        Tdim = T0 * th.subs({xi: x / Ls, tau: kap * t / Ls**2})
        qm = (k * T0 / Ls**2) * fx.subs({xi: x / Ls, tau: kap * t / Ls**2})
        dim_row = rho * cp * sp.diff(Tdim, t) - k * sp.diff(Tdim, x, 2) - q - qm
        d = sp.simplify(dim_row - (fac * r_nd).subs({xi: x / Ls, tau: kap * t / Ls**2}))
        u = (1 * ureg("W/(m*K)") * ureg("K") / ureg("m")**2).to("W/m**3")
        u2 = (1 * ureg("kg/m**3") * ureg("J/(kg*K)") * ureg("m**2/s") * ureg("K") / ureg("m")**2).to("W/m**3")
        return d == 0, {"factor": str(fac), "dimensional row - factor*M4-1": str(d),
                        "units k T0/L^2": str(u.units), "units rho c_p kappa T0/L^2": str(u2.units)}
    rep.check("W10_units", "(k T0/L^2) times M4-1 is the dimensional row rho c_p T_t - k T_xx - q - q_m, in W m^-3",
              w10)
    rep.finish()


if __name__ == "__main__":
    main()

"""M4 check: dimensionless groups, scales, dt rule and end fluxes (spec.md grid table, scales table, M4-8,
end_fluxes), against the ledger and M1.

Claims tested:
  SC1 the grid-table formulas for kappa, Q, dtau agree with the ledger (S16, S27, S20) and give the M1
      reference values (m1_inputs.json: kappa, t_c, Q, t_slow, end flux) within 1 %.
  SC2 units (Pint): Q, dtau dimensionless; dt = dtau L^2/kappa in s; k T0/L^2 in W m^-3; (k T0/L) F' in
      W m^-2.
  SC3 every number of the a-priori scales table and the explicit limit recompute within 1 %.
  SC4 the dt rule M4-8: lambda_1 is the l = 1 eigenvalue; b_1 is S53's beta_1; the first-order expansion
      and its maximum at tau = 1/lambda_1 give dtau_max = 2 e eps/(|Q| b_1 lambda_1); n_steps, dtau, dt,
      tau_end formulas; the reference values (eps <-> dtau) recompute within 1 %. BE amplification is in (0,1).
  SC5 the measured ratios of BE error to the mode-1 prediction (0.982 ... 0.997) and the location of the
      maximum recompute with an independent eigen-decomposition (within 0.005 and 0.002).
  SC6 the fully discrete max |theta - theta_s| at tau = 1 (N = 64, dtau = 1e-3) is 7.0172e-6 (1 %), above
      the continuous S57 bound value 6.67e-6; the 1.4x margin.
  SC7 the half-cell end flux F'_0 is the integral energy balance over [0, h/2] (second order on a smooth
      solution), equals Q/2 on the steady nodal values, and (k T0/L)(Q/2) = qL/2 (S45).

Usage: python m4_scales.py [--spec spec.md] [--ledger ledger.json] [--m1 m1_inputs.json]
"""
import json
import math
import os
import re
import sys

import numpy as np
import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import m3_common as C  # noqa: E402
import m4_common as M  # noqa: E402

SUP = str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹⁻", "0123456789-")


def sci(s):
    s = s.strip().translate(SUP).replace("−", "-")
    s = re.sub(r"×\s*10(-?\d+)", r"e\1", s)
    s = re.sub(r"^10(-?\d+)$", r"1e\1", s)
    return float(s)


def main():
    args = M.cli(__doc__, extra=[("--m1", M.DEFAULT_M1)])
    spec = M.Spec(args.spec)
    L = M.Ledger(args.ledger)
    m1 = json.load(open(args.m1))
    rep = M.Reporter("m4_scales.py", L, {"spec": str(spec.path), "spec_sha256": spec.sha256, "m1": args.m1,
                                         "m1_reference_case": m1["reference_case"]})
    ref = {k: v[0] for k, v in m1["reference_case"].items()}
    R = sp.Symbol
    sym = {n: R(n, real=True) for n in ("k", "rho", "cp", "q", "L", "T0", "dt", "kappa")}
    vals = {sym["k"]: ref["k"], sym["rho"]: ref["rho"], sym["cp"]: ref["cp"], sym["q"]: ref["q"],
            sym["L"]: ref["L"], sym["T0"]: ref["T0"]}
    grid = spec.section(r"Discretization")

    def grid_formula(name):
        cells = spec.table_row(name, within=grid)
        return M.uni(cells[1], extra_names=("cp",))

    # ------------------------------------------------------------------ SC1
    def sc1():
        kap = grid_formula("κ")
        Qf = grid_formula("Q")
        dtf = grid_formula("Δτ")
        # ledger
        s16 = C.find_piece(L.stmt("S16"), C.lhs_is("kappa"), "kappa =").exprs[1]
        s27 = C.find_piece(L.stmt("S27"), C.lhs_is("Q"), "Q =").exprs[1]
        s20 = C.find_piece(L.stmt("S20"), C.lhs_is("tau"), "tau =").exprs[1]
        lsub = {C.S("k"): sym["k"], C.S("rho"): sym["rho"], C.S("c_p"): sym["cp"], C.S("q"): sym["q"],
                C.S("L"): sym["L"], C.S("T_0"): sym["T0"], C.S("kappa"): sym["kappa"]}
        d_kap = sp.simplify(kap - s16.xreplace(lsub))
        d_Q = sp.simplify(Qf - s27.xreplace(lsub))
        d_dt = sp.simplify(dtf - s20.xreplace(lsub).subs(C.S("t").xreplace({}), sym["dt"]).xreplace(
            {C.S("t"): sym["dt"]}))
        kv = float(kap.subs(vals))
        Qv = float(Qf.subs(vals))
        tc = ref["L"] ** 2 / kv
        cl = {k: v[0] for k, v in m1["claimed"].items()}
        scl = spec.section(r"A priori scales")
        tslow_spec = float(re.search(r"Slowest mode, continuous \| L²/\(π²κ\) \| ([\d.]+) s", scl).group(1))
        tc_spec = float(re.search(r"Domain diffusion time t_c \| L²/κ \(S20\) \| ([\d.]+) s", scl).group(1))
        checks = {"kappa": (kv, cl["kappa"]), "Q": (Qv, cl["Q"]), "t_c": (tc, cl["t_c"]),
                  "t_c spec vs M1": (tc_spec, cl["t_c"]), "t_slow spec vs M1": (tslow_spec, cl["t_slow"]),
                  "end flux (k T0/L)(Q/2) vs M1": (ref["k"] * ref["T0"] / ref["L"] * Qv / 2, cl["end_flux"]),
                  "peak Q/8 vs M1": (Qv / 8, cl["Q_over_8"])}
        bad = [k for k, (a, b) in checks.items() if not M.rel_ok(a, b, m1.get("rel_tol", 0.01))]
        sym_ok = d_kap == 0 and d_Q == 0 and d_dt == 0
        return sym_ok and not bad, {"kappa, Q, dtau (spec)": [str(kap), str(Qf), str(dtf)],
                                    "spec - ledger": [str(d_kap), str(d_Q), str(d_dt)],
                                    "values (computed, M1)": checks, "failures": bad}
    rep.check("SC1_groups_vs_ledger_and_M1", "kappa = k/(rho c_p), Q = q L^2/(k T0), dtau = kappa dt/L^2 match "
              "S16, S27, S20 and give M1's reference values", sc1)

    # ------------------------------------------------------------------ SC2
    def sc2():
        import pint
        u = pint.UnitRegistry()
        P = {"k": 1 * u("W/(m*K)"), "rho": 1 * u("kg/m**3"), "cp": 1 * u("J/(kg*K)"), "q": 1 * u("W/m**3"),
             "L": 1 * u("m"), "T0": 1 * u("K"), "dt": 1 * u("s")}
        kap = P["k"] / (P["rho"] * P["cp"])
        Q = P["q"] * P["L"] ** 2 / (P["k"] * P["T0"])
        dtau = kap * P["dt"] / P["L"] ** 2
        dt_back = dtau * P["L"] ** 2 / kap
        fac = P["k"] * P["T0"] / P["L"] ** 2
        flux = P["k"] * P["T0"] / P["L"]
        res = {"Q": Q.to_base_units().dimensionless, "dtau": dtau.to_base_units().dimensionless,
               "dt": dt_back.check("[time]"), "kT0/L^2": fac.check("[power]/[length]**3"),
               "kT0/L * F'": flux.check("[power]/[length]**2")}
        # the spec states the flux unit and residual unit
        t_ok = ("F_0 = (kT_0/L)·F′_0 W m⁻²" in spec.text
                and re.search(r"The dimensional row is \(k T_0/L²\) × M4-1 in\s+W m⁻³", spec.text) is not None)
        return all(res.values()) and t_ok, {k: bool(v) for k, v in res.items()} | {"spec unit text found": t_ok}
    rep.check("SC2_units", "Q and dtau are dimensionless, dt is a time, k T0/L^2 is W m^-3 and (k T0/L) F' is W m^-2",
              sc2)

    # ------------------------------------------------------------------ SC3
    def sc3():
        scl = spec.section(r"A priori scales")
        N = int(re.search(r"reference case C1, N = (\d+)", spec.text).group(1))
        kv, Lv = ref["k"] / (ref["rho"] * ref["cp"]), ref["L"]
        lam = lambda l: 4 * N**2 * math.sin(l * math.pi / (2 * N)) ** 2  # noqa: E731
        comp = {
            "t_c": Lv**2 / kv,
            "slow continuous": Lv**2 / (math.pi**2 * kv),
            "lambda_1": lam(1), "pi^2": math.pi**2,
            "mesh time": (Lv / N) ** 2 / kv,
            "lambda_max": lam(N - 1), "1/lambda_max": 1 / lam(N - 1),
            "stiffness ratio": lam(N - 1) / lam(1),
        }
        g = lambda pat: sci(re.search(pat, scl).group(1))  # noqa: E731
        claim = {
            "t_c": g(r"t_c \| L²/κ \(S20\) \| ([\d.]+) s"),
            "slow continuous": g(r"L²/\(π²κ\) \| ([\d.]+) s"),
            "lambda_1": g(r"λ_1 = ([\d.]+) \(π²"),
            "pi^2": g(r"\(π² = ([\d.]+)\)"),
            "mesh time": g(r"\(L/N\)²/κ \| ([\d.e-]+) s"),
            "lambda_max": g(r"λ_max = ([\d.]+)\)"),
            "1/lambda_max": g(r"\| ([\d.e-]+) \(λ_max ="),
            "stiffness ratio": g(r"Stiffness ratio λ_max/λ_1 \| \| (\d+) \|"),
        }
        mex = re.search(r"Δτ ≤ 2/λ_max = ([\d.e-]+)\.", spec.text)
        comp["explicit limit"] = 2 / lam(N - 1)
        claim["explicit limit"] = float(mex.group(1))
        bad = [k for k in comp if not M.rel_ok(comp[k], claim[k], 0.01)]
        return not bad, {"N": N, "computed vs spec": {k: [comp[k], claim[k]] for k in comp}, "failures": bad}
    rep.check("SC3_scale_table", "every entry of the a-priori scales table (and 2/lambda_max) recomputes within 1 %",
              sc3)

    # ------------------------------------------------------------------ SC4
    def sc4():
        blk = spec.block("M4-8_pre")
        ps = M.latex_pieces(blk)
        lam1 = next(p for p in ps if p.ok and p.exprs[0] == C.S("lambda_1")).exprs[1]
        dmax = next(p for p in ps if p.ok and p.exprs[0] == C.S("dtau_max") and len(p.exprs) == 2
                    and p.exprs[1] != C.S("tau_end")).exprs[1]
        b1 = next(p for p in ps if p.ok and p.exprs[0] == C.S("b_1")).exprs[1]
        # the Q = 0 branch sits in a parenthesised condition, which the piece splitter drops; read it as text
        q0 = re.search(r"Q\\neq0;\\ \\ \\Delta\\tau_\{\\max\}=\\tau_\{end\}\\ \\text\{if\}\\ Q=0\)", blk) is not None
        # the M3 parser reads a bare 'e' as the energy field e(x,t); here it is Euler's number (2e eps)
        dmax = dmax.replace(lambda z: isinstance(z, sp.core.function.AppliedUndef) and z.func.__name__ == "e",
                            lambda z: sp.E).subs(C.S("e"), sp.E)
        Ns, ell, hh = C.S("N"), sp.Symbol("ell", positive=True), sp.Symbol("h", positive=True)
        lam_l = (4 / hh**2) * sp.sin(ell * sp.pi * hh / 2) ** 2
        d_lam = sp.simplify(lam1 - lam_l.subs({ell: 1, hh: 1 / Ns}))
        # S53 beta_1
        s53 = L.stmt("S53")
        mb = re.search(r"= \\frac\{2\(1-\(-1\)\^\\ell\)\}\{\(\\ell\\pi\)\^3\}", s53)
        beta1 = sp.Integer(2) * (1 - (-1) ** 1) / sp.pi**3 if mb else None
        d_b1 = sp.simplify(b1 - beta1) if beta1 is not None else "S53 form not found"
        # independent derivation
        lam, dt, t_, Qa, eps = sp.symbols("lambda dtau tau Qa epsilon", positive=True)
        ratio_log = (t_ / dt) * (lam * dt - sp.log(1 + lam * dt))
        ser = sp.series(sp.exp(ratio_log), dt, 0, 2).removeO()
        lead = sp.simplify(ser - 1)
        e_t = b1 * lead * sp.exp(-lam * t_)
        tstar = sp.solve(sp.diff(e_t, t_), t_)
        emax = sp.simplify(e_t.subs(t_, tstar[0]))
        dsolve = sp.solve(sp.Eq(Qa * emax, eps), dt)[0]
        dclaim = dmax.subs({C.S("epsth"): eps, C.S("b_1"): b1, C.S("lambda_1"): lam}).subs(sp.Abs(C.S("Q")), Qa)
        d_rule = sp.simplify(dsolve - dclaim)
        # second block: n_steps, dtau, dt, tau_end
        ps2 = M.latex_pieces(spec.block("M4-8"))
        def get(lhs):
            return next(p for p in ps2 if p.ok and p.exprs[0] == C.S(lhs)).exprs[1]
        nst, dtv, dtstep, tend = get("n_steps"), get("dtau"), get("dtstep"), get("tau_end")
        ceil_ok = "\\lceil" in spec.block("M4-8") and sp.simplify(nst - C.S("tau_end") / C.S("dtau_max")) == 0
        f_ok = (sp.simplify(dtv - C.S("tau_end") / C.S("n_steps")) == 0
                and sp.simplify(dtstep - C.S("dtau") * C.S("L")**2 / C.S("kappa")) == 0
                and sp.simplify(tend - C.S("kappa") * C.S("t_end") / C.S("L")**2) == 0)
        # reference values
        N = 64
        try:
            l1 = float(lam1.subs(Ns, N))
            b1v = float(b1)
            eps_ref = b1v * l1 * 1e-3 / (2 * math.e)
            dmax_1e4 = float(dmax.subs({C.S("epsth"): 1e-4, sp.Abs(C.S("Q")): 1, C.S("b_1"): b1v,
                                        C.S("lambda_1"): l1}))
        except TypeError as exc:
            return False, {"error": str(exc), "lambda_1": str(lam1), "b_1": str(b1), "dtau_max": str(dmax),
                           "free symbols": [str(s) for s in dmax.free_symbols | lam1.free_symbols]}
        sec = spec.text
        m_e = re.search(r"Δτ = 10⁻³ corresponds to ε_θ′ = ([\d.e-]+) with Q = 1\. That is ([\d.]+) % of the peak rise 1/8",
                        sec)
        m_d = re.search(r"ε_θ′ = 10⁻⁴ gives Δτ_max = ([\d.e-]+)\.", sec)
        nums = {"eps for dtau=1e-3": (eps_ref, float(m_e.group(1))),
                "percent of 1/8": (100 * eps_ref / 0.125, float(m_e.group(2))),
                "dtau_max for eps=1e-4": (dmax_1e4, float(m_d.group(1)))}
        bad_n = [k for k, (a, b) in nums.items() if not M.rel_ok(a, b, 0.01 if "percent" not in k else 0.03)]
        # BE amplification
        amp = [1 / (1 + d * l) for d in (1e-5, 1e-3, 1.0) for l in (l1, 16374.13)]
        amp_ok = all(0 < a < 1 for a in amp)
        ok = (d_lam == 0 and d_b1 == 0 and d_rule == 0 and len(tstar) == 1 and sp.simplify(tstar[0] - 1 / lam) == 0
              and ceil_ok and f_ok and q0 and not bad_n and amp_ok)
        return ok, {"lambda_1 - (4/h^2)sin^2(pi h/2)|h=1/N": str(d_lam), "b_1 - S53 beta_1": str(d_b1),
                    "expansion ratio - 1": str(lead), "argmax tau": [str(t) for t in tstar],
                    "max error": str(emax), "dtau_max derived": str(dsolve), "spec": str(dclaim),
                    "derived - spec": str(d_rule), "n_steps/dtau/dt/tau_end ok": [ceil_ok, f_ok],
                    "Q = 0 branch stated": bool(q0), "numbers (computed, spec)": nums, "failures": bad_n,
                    "BE amplification in (0,1)": amp_ok}
    rep.check("SC4_dt_rule", "M4-8 follows from the mode-1 BE error expansion, its maximum at tau = 1/lambda_1, "
              "lambda_1 and b_1 = beta_1; the reference eps <-> dtau values recompute", sc4)

    # ------------------------------------------------------------------ discrete helpers
    def eig_setup(N):
        h = 1.0 / N
        n = N - 1
        A = (np.diag(2 * np.ones(n)) - np.diag(np.ones(n - 1), 1) - np.diag(np.ones(n - 1), -1)) / h**2
        w, V = np.linalg.eigh(A)
        x = np.arange(1, N) * h
        ths = np.linalg.solve(A, np.ones(n))  # discrete steady theta (Q = 1 scaling: theta)
        c = V.T @ (-ths)  # theta(0) - theta_s = -theta_s
        return w, V, x, ths, c

    # ------------------------------------------------------------------ SC5
    def sc5():
        N = 64
        w, V, x, ths, c = eig_setup(N)
        lam1 = w[0]
        b1 = 4 / math.pi**3
        m = re.search(r"has a maximum over τ that is ([\d.]+), ([\d.]+), ([\d.]+) and ([\d.]+) times the\s+mode-1 "
                      r"prediction, for Δτ = ([\d.]+), ([\d.]+), ([\d.]+) and ([\d.]+) × 10⁻³\. The maxima sit at "
                      r"τ = ([\d.]+) to ([\d.]+)", spec.text)
        if not m:
            return False, "measured-ratio sentence not found"
        g = [float(v) for v in m.groups()]
        ratios_s, dts_s, tlo, thi = g[0:4], [v * 1e-3 for v in g[4:8]], g[8], g[9]
        out, bad = [], []
        for rs, dt in zip(ratios_s, dts_s):
            nsteps = int(round(0.5 / dt))
            n_arr = np.arange(1, nsteps + 1)
            tau = n_arr * dt
            be = (1.0 / (1 + dt * w))[None, :] ** n_arr[:, None]
            ex = np.exp(-np.outer(tau, w))
            err = np.abs(((be - ex) * c[None, :]) @ V.T).max(axis=1)
            k = int(np.argmax(err))
            ratio = err[k] / (b1 * lam1 * dt / (2 * math.e))
            out.append({"dtau": dt, "ratio": ratio, "spec": rs, "tau at max": tau[k]})
            if abs(ratio - rs) > 0.005 or not (tlo - 0.002 <= tau[k] <= thi + 0.002):
                bad.append(dt)
        return not bad, {"per dtau": out, "failures": bad, "spec tau range": [tlo, thi]}
    rep.check("SC5_be_error_ratio", "the BE-only error maxima are 0.982..0.997 of the mode-1 prediction and sit near "
              "tau = 1/lambda_1", sc5)

    # ------------------------------------------------------------------ SC6
    def sc6():
        N, dt, nsteps = 64, 1e-3, 1000
        w, V, x, ths, c = eig_setup(N)
        theta = ths + V @ (c * (1.0 / (1 + dt * w)) ** nsteps)
        exact_s = x * (1 - x) / 2
        e = np.max(np.abs(theta - exact_s))
        m = re.search(r"maximum of \|θ − θ_s\| is ([\d.e-]+)\s*\(`m4_checks.py` item 9\), slightly above the continuous "
                      r"([\d.e-]+) \(S57\)", spec.text)
        mm = re.search(r"would have a margin of only ([\d.]+)×", spec.text)
        cont = 4 / math.pi**3 * math.exp(-math.pi**2)
        ok = (m is not None and M.rel_ok(e, float(m.group(1)), 0.01) and M.rel_ok(cont, float(m.group(2)), 0.01)
              and e > cont and mm is not None and abs(1e-5 / e - float(mm.group(1))) < 0.05)
        return ok, {"computed max|theta-theta_s|": e, "spec": m.group(1) if m else None,
                    "continuous mode-1 value": cont, "margin 1e-5/e": 1e-5 / e,
                    "spec margin": mm.group(1) if mm else None}
    rep.check("SC6_discrete_floor", "fully discrete max|theta - theta_s| at tau = 1 is 7.0172e-6 (1 %), above the "
              "continuous 6.67e-6, margin 1.4x to 1e-5", sc6)

    # ------------------------------------------------------------------ SC7
    def sc7():
        m = re.search(r"At ξ = 0: F′_0 = (.+?)\.\n", spec.text)
        if not m:
            return False, "end flux formula not found"
        s = m.group(1)
        s = s.replace("U_0^prev", "U0p").replace("f(0, τ)", "f0").replace("U_0", "U0").replace("U_1", "U1")
        Ff = M.uni(s, extra_names=("U0p", "f0", "U0", "U1"))
        h_, Q_, dt_ = (sp.Symbol(n, real=True) for n in ("h", "Q", "dtau"))
        U0, U1, U0p, f0 = (sp.Symbol(n, real=True) for n in ("U0", "U1", "U0p", "f0"))
        # (a) steady nodal values: U0 = 1, U1 = 1 + Q h (1-h)/2, f = 0
        st = sp.simplify(Ff.subs({U0: 1, U0p: 1, U1: 1 + Q_ * h_ * (1 - h_) / 2, f0: 0}))
        # (b) half-cell balance: on a smooth solution the formula approximates d_xi theta'(0) to O(h^2)
        xi_, t_ = sp.symbols("xi tau", real=True)
        th = 1 + sp.sin(3 * xi_) * sp.exp(-t_) + xi_**3 * t_  # any smooth function, with f from M4-1
        Qv = 0.7
        f_ = sp.diff(th, t_) - sp.diff(th, xi_, 2) - Qv
        t0 = 0.3
        errs = []
        for hv in (0.1, 0.05, 0.025):
            # U_0 is not fixed in this test; use the exact time derivative in place of the BE difference
            Fv = (th.subs({xi_: hv, t_: t0}) - th.subs({xi_: 0, t_: t0})) / hv + (hv / 2) * (Qv + f_.subs({xi_: 0, t_: t0})) \
                - (hv / 2) * sp.diff(th, t_).subs({xi_: 0, t_: t0})
            Fv_spec = Ff.subs({U1: th.subs({xi_: hv, t_: t0}), U0: th.subs({xi_: 0, t_: t0}), Q_: Qv,
                               f0: f_.subs({xi_: 0, t_: t0}), h_: hv,
                               U0p: th.subs({xi_: 0, t_: t0}) - sp.diff(th, t_).subs({xi_: 0, t_: t0}) * dt_})
            Fv_spec = sp.simplify(Fv_spec)
            ex = sp.diff(th, xi_).subs({xi_: 0, t_: t0})
            errs.append((float(abs(Fv_spec - ex)), float(abs(Fv - ex))))
        orders = [math.log(errs[i][0] / errs[i + 1][0], 2) for i in range(2)]
        # (c) dimensional: (k T0 / L) * Q/2 = q L / 2 with Q = q L^2/(k T0) (S45)
        k, T0, Ls, q = sp.symbols("k T_0 L q", positive=True)
        dim = sp.simplify((k * T0 / Ls) * (q * Ls**2 / (k * T0)) / 2 - q * Ls / 2)
        s45 = C.find_piece(L.stmt("S45"), lambda p: len(p.exprs) == 3, "F_0 = F_L = ...").exprs[-1]
        d45 = sp.simplify(s45.xreplace({C.S("q"): q, C.S("L"): Ls}) - q * Ls / 2)
        ok = sp.simplify(st - Q_ / 2) == 0 and all(o > 1.8 for o in orders) and dim == 0 and d45 == 0
        return ok, {"F'_0 (spec)": str(Ff), "on steady nodal values": str(st),
                    "|F'_0 - d_xi theta'(0)| at h = .1,.05,.025": [e[0] for e in errs], "observed orders": orders,
                    "(kT0/L)(Q/2) - qL/2": str(dim), "S45 - qL/2": str(d45)}
    rep.check("SC7_end_flux", "the half-cell flux F'_0 approximates d_xi theta'(0) to O(h^2), equals Q/2 on the steady "
              "nodal values, and (k T0/L)(Q/2) = qL/2 = S45", sc7)
    rep.finish()


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""interface_width_estimate.py — would a narrower interface raise the exponent?

    venv_enceladus/bin/python studies/molaro_2019/demmenie/interface_width_estimate.py

An ESTIMATE from the one saturated run, not a test. The run's neck velocity is
divided by the Gibbs–Thomson prediction (predict_neck_growth.v_neck); the
excess is fitted as a function of (diffuse band) / (neck fillet radius),

    run / prediction = 1 + k (band / fillet)^n ,   band = 9.19 eps (1-99 %),
    fillet = x^2 / (2 (R - x)),

and the prediction is then integrated with that excess for eps, eps/2, eps/4
and 0, and fitted with the same C (t + t0)^a as the run.

What it cannot tell apart: in a single run band/fillet falls monotonically as
the neck grows, so ANY early excess that decays with neck size (for instance
the initial union-of-spheres neck relaxing towards the fillet shape the
prediction assumes) would fit as well. Only a run at eps/2 separates the two.

Writes prediction/interface_width_estimate.txt.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import predict_neck_growth as P                                    # noqa: E402
import analyze_demmenie as ad                                       # noqa: E402

RUN = next((Path.home() / "SimulationResults/HPC_results/enceladus_DSM/GrainPairSintering"
            / "batch_2026-10-06__09.21.04_demmenie_mirror").glob("molaro_2D_*"))
R, T_C, ALPHA, EPS, BAND = 84.4e-6, -20.0, 1e-3, 0.12e-6, 9.19


def main():
    t, w, _ = ad.load(RUN)
    keep = t >= 11 * ad.read_tau_sub(RUN)
    t, w = t[keep], w[keep]
    th = P.thermo(T_C)
    x = w * 1e-6 / 2
    ratio = (2 / np.pi) * np.gradient(x, t) / np.array([P.v_neck(xi, R, th, ALPHA) for xi in x])
    rho = x ** 2 / (2 * (R - x))
    f = lambda br, k, n: 1 + k * br ** n
    s = slice(1, len(t) - 1)
    (k, n), _ = curve_fit(f, BAND * EPS / rho[s], ratio[s], p0=(1.0, 1.0))
    L = ["run velocity / Gibbs-Thomson prediction along the run:"]
    for i in (1, 5, 15, 30, 50, 70, 95):
        L.append(f"  t {t[i] / 3600:5.1f} h   w {w[i]:5.1f} um   fillet {rho[i] * 1e6:5.2f} um = {rho[i] / (BAND * EPS):4.1f} bands   ratio {ratio[i]:.2f}")
    rms = float(np.sqrt(np.mean((f(BAND * EPS / rho[s], k, n) - ratio[s]) ** 2)))
    L += ["", f"fit: ratio = 1 + {k:.2f} (band/fillet)^{n:.2f}   (rms {rms:.3f})", "",
          "prediction integrated with that excess, fitted as C (t + t0)^a over the same samples:"]

    def integrate(eps):
        def rhs(tt, y):
            r = y[0] ** 2 / (2 * (R - y[0]))
            return [(np.pi / 2) * P.v_neck(y[0], R, th, ALPHA) * (1 + (k * (BAND * eps / r) ** n if eps > 0 else 0.0))]
        return 2e6 * solve_ivp(rhs, (t[0], t[-1]), [x[0]], t_eval=t, rtol=1e-9, atol=1e-14).y[0]

    for lab, e in (("eps as run (0.12 um)", EPS), ("eps / 2", EPS / 2), ("eps / 4", EPS / 4), ("sharp limit", 0.0)):
        ww = integrate(e)
        L.append(f"  {lab:22s} final w {ww[-1]:.2f} um   a = {ad.fit_free(t, ww)['a']:.3f}")
    L.append(f"  {'the run itself':22s} final w {w[-1]:.2f} um   a = {ad.fit_free(t, w)['a']:.3f}")
    txt = "\n".join(L)
    print(txt)
    (HERE / "prediction").mkdir(exist_ok=True)
    (HERE / "prediction/interface_width_estimate.txt").write_text(txt + "\n")


if __name__ == "__main__":
    main()

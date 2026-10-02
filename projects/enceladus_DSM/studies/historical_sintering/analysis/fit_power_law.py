#!/usr/bin/env python3
"""fit_power_law.py -- fit the Thomas et al. (1994) neck data to x/a = C t^a.

Run from enceladus_DSM/:
    python studies/historical_sintering/analysis/fit_power_law.py

Two forms, each over two windows:
  d_fixed  x/a = C t^a        least squares in log-log (clock as reported)
  d_free   x/a = C (t+t0)^a   t0 free, Demmenie's protocol -- absorbs the unknown
                              contact time, at the price of a wider CI on a
Windows: all 11 points, and points 3-11 (x/a >= 0.181), the window the
relaxed simulation resolves. 95 % CIs. Also prints point-to-point local slopes.
Writes data/thomas_powerlaw_fits.csv.
"""

import csv
from pathlib import Path

import numpy as np
from scipy import stats
from scipy.optimize import curve_fit

DATA = Path(__file__).resolve().parent.parent / "data"


def load():
    lines = [l for l in open(DATA / "thomas1994_T-20_r120um.csv")
             if not l.startswith("#")]
    return np.loadtxt(lines[1:], delimiter=",", unpack=True)


def d_fixed(t, u):
    r = stats.linregress(np.log(t), np.log(u))
    ci = r.stderr * stats.t.ppf(0.975, len(t) - 2)
    return dict(a=r.slope, a_ci=ci, C=np.exp(r.intercept), t0_h=0.0,
                t0_ci=0.0, R2=r.rvalue**2)


def d_free(t, u):
    f = lambda t, C, t0, a: C * (t + t0)**a
    p, c = curve_fit(f, t, u, p0=[0.07, 1.0, 0.35],
                     bounds=([0, -t.min() + 1e-3, 0.01], [10, 1e3, 2]),
                     maxfev=20000)
    ci = np.sqrt(np.diag(c)) * stats.t.ppf(0.975, len(t) - 3)
    res = u - f(t, *p)
    return dict(a=p[2], a_ci=ci[2], C=p[0], t0_h=p[1], t0_ci=ci[1],
                R2=1 - res.var() / u.var())


def main():
    t, u = load()
    rows = []
    for win, s in [("all (1-11)", slice(None)), ("pts 3-11", slice(2, None))]:
        for form, fn in [("d_fixed", d_fixed), ("d_free", d_free)]:
            r = fn(t[s], u[s])
            rows.append(dict(window=win, form=form, n=len(t[s]),
                             **{k: round(float(v), 4) for k, v in r.items()}))
    for r in rows:
        print(r)
    slopes = np.diff(np.log(u)) / np.diff(np.log(t))
    print("local slopes:", np.round(slopes, 3))
    with open(DATA / "thomas_powerlaw_fits.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()

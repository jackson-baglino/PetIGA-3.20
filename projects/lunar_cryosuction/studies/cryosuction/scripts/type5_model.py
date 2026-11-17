#!/usr/bin/env python3
"""type5_model.py: reduced model of an ice lens in a type-5 pore (throat), used to design the
type-5 runs and as the prediction their results are compared with.

The model is the one in the SI draft "Supplementary Methods: vapor boundary conditions for the
throat simulations" (iCloud lunar_regolith_DSM/SupplementalInfo_paper2/SI_throat_vapor_BC.docx,
draft 3.2, Steps 1-5 and Eq. S7). It was first coded in plan_estimates.py section 8c
(optimal_h); that section now calls this module.

Coordinates. xi = x - L/2: the throat is at xi = 0, the WARM end (x = 0) at xi = -L/2, the cold
end at xi = +L/2. A lens center xm > 0 is on the COLD side, xm < 0 on the WARM side. Every print
says "cold side" or "warm side" in words (Table S1 of draft 3.2 had the sign of the G = 0 rows
reversed).

Pores (gap g(xi), widest W, throat wt, length L):
  hyp     g = sqrt(wt^2 + c^2 xi^2), c = sqrt(W^2 - wt^2)/(L/2): g = W at both ends (non-periodic;
          the gradient runs, ends at a fixed vapor density). wt = W gives the straight pore.
  sinhyp  g = sqrt(wt^2 + (W^2 - wt^2) sin^2(pi xi/L)): one cell of a periodic chain (G = 0 runs).

Lens faces. Face i spans the local gap g_i and meets the walls at theta; with the walls inclined
by omega_i it is an arc of curvature chi_i = 2 cos(theta + omega_i)/g_i (S1), with
tan omega_i = side g'(xi_i)/2, side = -1 for the face whose gas lies toward -xi (the warm face
w), +1 for the other (cold face c). Equilibrium humidity h_i = 1 - d0 chi_i (S2).

Ends model (hyp pores, S3-S6). Face i gains ice area at Q_i = Delta_i / R_i with
R_i = beta0/g_i + K int_column dx/g; Delta_w = h - h_w + s dT_w, Delta_c = h - h_c - s dT_c, the
temperature drop G L split between the two gas columns in proportion to int dx/g (isothermal
ice). dA/dt = Q_w + Q_c. The center moves at the mean face speed. h_net (S6) is the constant end
humidity that returns the lens to its initial area at t_final, iterated with the path.
G_crit is the largest G the warm half of the throat can hold the lens against.

Periodic model (sinhyp, G = 0, S7). The two faces exchange vapor through one gas column that runs
through the neighboring cell: Q_1 = -Q_2 = (h_2 - h_1)/(beta0/g_1 + beta0/g_2 + K int dx/g).
Area is conserved; the lens moves toward the narrower gap.

Approximations (SI): isothermal ice, quasi-steady 1D gas columns, face length = gap, static
contact angle, linearized Gibbs-Thomson, lens area = int g dx between the two face apexes.

    venv_cryosuction/bin/python studies/cryosuction/scripts/type5_model.py ends --wt 40e-6 80e-6 --G 0 1 2
    venv_cryosuction/bin/python studies/cryosuction/scripts/type5_model.py ends --wt 150e-6 --G 1   # straight pore
    venv_cryosuction/bin/python studies/cryosuction/scripts/type5_model.py periodic --wt 40e-6 80e-6 --days 300
"""
from __future__ import annotations

import argparse
import math
import sys
from dataclasses import dataclass
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from preprocess.mesh_params import physics  # noqa: E402
from preprocess.mesh_params.constants import CONSTANTS as C  # noqa: E402

DAY = 86400.0
W_DEF, L_DEF, LENS_DEF, START_DEF, THETA_DEF = 150e-6, 400e-6, 100e-6, 70e-6, 60.0


@dataclass(frozen=True)
class Consts:
    d0: float      # capillary length [m]
    s: float       # dln(rho_vs)/dT [1/K]
    beta0: float   # K&P beta0 [s/m]
    K: float       # rho_ice/(D_v rho_vs) [s/m^2]


def consts(T0: float = -20.0, alpha: float = 1e-2) -> Consts:
    """The solver's constants at T0 (same expressions as plan_estimates.py)."""
    rho_vs = physics.rho_vs_sat_kg_per_m3(T0)
    D = physics.dv_m2_per_s(T0)
    d0 = physics.capillary_length_m(T0)
    beta0 = physics.beta_hk_s_per_m(T0, alpha) * C.rho_ice_kg_per_m3 / rho_vs
    dT = 0.01
    s = (math.log(physics.rho_vs_sat_kg_per_m3(T0 + dT))
         - math.log(physics.rho_vs_sat_kg_per_m3(T0 - dT))) / (2 * dT)
    return Consts(d0, s, beta0, C.rho_ice_kg_per_m3 / (D * rho_vs))


def side_words(xm: float) -> str:
    """'34.0 um on the cold side' etc. for a lens center xm (xi, m)."""
    if abs(xm) < 0.05e-6:
        return "at the throat"
    return f"{abs(xm) * 1e6:.1f} um on the {'cold' if xm > 0 else 'warm'} side"


class EndsPore:
    """hyp pore (or straight, wt = W) with fixed vapor at both ends: plan_estimates section 8c."""

    def __init__(self, wt, k: Consts, W=W_DEF, L=L_DEF, theta=THETA_DEF, lens=LENS_DEF, start=START_DEF):
        self.wt, self.k, self.W, self.L, self.lens, self.start = wt, k, W, L, lens, start
        self.th = math.radians(theta)
        self.c = math.sqrt(max(W * W - wt * wt, 0.0)) / (L / 2)

    def g(self, xi):
        return math.sqrt(self.wt * self.wt + (self.c * xi) ** 2)

    def chi_face(self, xi, side):
        """2 cos(theta + omega)/g; side -1 warm face (gas toward x = 0), +1 cold face."""
        dg = self.c * self.c * xi / self.g(xi)
        return 2 * math.cos(self.th + math.atan(0.5 * side * dg)) / self.g(xi)

    def col(self, x0, x1, n=200):                  # int_x0^x1 dx/g
        return sum(1 / self.g(x0 + (k + 0.5) * (x1 - x0) / n) for k in range(n)) * (x1 - x0) / n

    def area(self, xm, Ln, n=200):
        return sum(self.g(xm - Ln / 2 + (k + 0.5) * Ln / n) for k in range(n)) * Ln / n

    def length(self, xm, A):
        lo, hi = 1e-6, self.L
        for _ in range(50):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if self.area(xm, mid) < A else (lo, mid)
        return 0.5 * (lo + hi)

    def state(self, xm, A, h, G):
        """(center speed toward x = 0, dA/dt, a, h0, mean face curvature) or None at an end."""
        k, L = self.k, self.L
        Ln = self.length(xm, A)
        xw, xc = xm - Ln / 2, xm + Ln / 2
        if xw <= -L / 2 + 5e-6 or xc >= L / 2 - 5e-6:
            return None
        iw, ic = self.col(-L / 2, xw), self.col(xc, L / 2)
        dTw, dTc = G * L * iw / (iw + ic), G * L * ic / (iw + ic)
        Rw = k.beta0 / self.g(xw) + k.K * iw
        Rc = k.beta0 / self.g(xc) + k.K * ic
        hw, hc = 1 - k.d0 * self.chi_face(xw, -1), 1 - k.d0 * self.chi_face(xc, +1)
        aa = 1 / Rw + 1 / Rc
        h0 = (hw / Rw + hc / Rc + k.s * (dTc / Rc - dTw / Rw)) / aa
        Qw = (h - hw + k.s * dTw) / Rw
        Qc = (h - hc - k.s * dTc) / Rc
        v = 0.5 * (Qw / self.g(xw) - Qc / self.g(xc))
        return v, Qw + Qc, aa, h0, 0.5 * (self.chi_face(xw, -1) + self.chi_face(xc, +1))

    def run(self, h, G, days=300.0, dt=6 * 3600.0, path=False):
        """Integrate the lens. Returns dict: xm (final center, xi), dA (final area/initial - 1),
        maxdev, t (s; < days if a face reached an end), hist [(a, h0, chi_mean)], and with path=True
        the trajectory [(t, xm, A/A0)]."""
        xm, A = self.start, self.area(self.start, self.lens)
        A0, t, hist, maxdev, traj = A, 0.0, [], 0.0, [(0.0, xm, 1.0)]
        while t < days * DAY:
            st = self.state(xm, A, h, G)
            if st is None:
                break
            v, dA, aa, h0, cm = st
            hist.append((aa, h0, cm))
            xm -= v * dt
            A += dA * dt
            t += dt
            maxdev = max(maxdev, abs(A / A0 - 1))
            if path:
                traj.append((t, xm, A / A0))
        return {"xm": xm, "dA": A / A0 - 1, "maxdev": maxdev, "t": t, "hist": hist, "traj": traj}

    def h_start(self):
        return self.state(self.start, self.area(self.start, self.lens), 1.0, 0.0)[3]

    def g_crit(self):
        """Largest G the warm half of the throat holds the lens against (h at the start value)."""
        hs, A, gc = self.h_start(), self.area(self.start, self.lens), 0.0
        for kk in range(1, 300):
            xm = -kk * 0.5e-6
            s0, s1 = self.state(xm, A, hs, 0.0), self.state(xm, A, hs, 1.0)
            if s0 is None or s1 is None:
                break
            gc = max(gc, -s0[0] / (s1[0] - s0[0]))
        return gc

    def converge(self, rule, G, days=300.0):
        """h -> rule(path(h)) until |dh| < 1e-12 (at most 20 passes); returns (h, passes)."""
        h = self.h_start()
        for kk in range(20):
            hn = rule(self.run(h, G, days)["hist"])
            if abs(hn - h) < 1e-12:
                return hn, kk + 1
            h = hn
        return h, 20

    def h_net(self, G, days=300.0):
        return self.converge(rule_net, G, days)


def rule_ls(hist):
    return sum(x[0] ** 2 * x[1] for x in hist) / sum(x[0] ** 2 for x in hist)


def rule_net(hist):
    return sum(x[0] * x[1] for x in hist) / sum(x[0] for x in hist)


class PeriodicCell:
    """sinhyp periodic cell, G = 0, faces coupled through the neighboring cell (S7)."""

    def __init__(self, wt, k: Consts, W=W_DEF, L=L_DEF, theta=THETA_DEF, lens=LENS_DEF, start=START_DEF):
        self.wt, self.k, self.W, self.L, self.lens, self.start = wt, k, W, L, lens, start
        self.th = math.radians(theta)
        self.a2 = W * W - wt * wt

    def g(self, xi):
        return math.sqrt(self.wt * self.wt + self.a2 * math.sin(math.pi * xi / self.L) ** 2)

    def dg(self, xi):
        return self.a2 * math.sin(2 * math.pi * xi / self.L) * (math.pi / self.L) / (2 * self.g(xi))

    def chi_face(self, xi, side):
        return 2 * math.cos(self.th + math.atan(0.5 * side * self.dg(xi))) / self.g(xi)

    def _int(self, f, a, b, n=400):
        return sum(f(a + (kk + 0.5) * (b - a) / n) for kk in range(n)) * (b - a) / n

    def area(self, xm, Ln):
        return self._int(self.g, xm - Ln / 2, xm + Ln / 2)

    def length(self, xm, A):
        lo, hi = 1e-6, 0.95 * self.L
        for _ in range(50):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if self.area(xm, mid) < A else (lo, mid)
        return 0.5 * (lo + hi)

    def velocity(self, xm, A):
        """center velocity d xm/dt [m/s] (negative = toward the throat for xm > 0)."""
        k = self.k
        Ln = self.length(xm, A)
        x1, x2 = xm - Ln / 2, xm + Ln / 2               # face 1: gas toward -xi; face 2: toward +xi
        h1, h2 = 1 - k.d0 * self.chi_face(x1, -1), 1 - k.d0 * self.chi_face(x2, +1)
        col = self._int(lambda x: 1 / self.g(x), x2, x1 + self.L)
        Q1 = (h2 - h1) / (k.beta0 / self.g(x1) + k.beta0 / self.g(x2) + k.K * col)
        # face 1 gains Q1 and advances toward -xi by Q1/g1; face 2 loses Q1 and recedes toward -xi
        return -0.5 * (Q1 / self.g(x1) + Q1 / self.g(x2))

    def run(self, days=300.0, dt=6 * 3600.0):
        xm, A = self.start, self.area(self.start, self.lens)
        t, traj = 0.0, [(0.0, xm)]
        while t < days * DAY - 1:
            xm += self.velocity(xm, A) * dt
            t += dt
            traj.append((t, xm))
        return traj


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("mode", choices=["ends", "periodic"])
    p.add_argument("--wt", type=float, nargs="+", default=[40e-6, 80e-6], help="throat widths [m] (= W: straight)")
    p.add_argument("--G", type=float, nargs="+", default=[0.0, 1.0, 2.0], help="gradients [K/m] (ends mode)")
    p.add_argument("--Gfrac", type=float, nargs="*", default=[],
                   help="also run these multiples of each pore's G_crit (ends mode)")
    p.add_argument("--days", type=float, default=300.0)
    p.add_argument("--W", type=float, default=W_DEF)
    p.add_argument("--L", type=float, default=L_DEF)
    p.add_argument("--lens", type=float, default=LENS_DEF, help="initial lens length on the axis [m]")
    p.add_argument("--start", type=float, default=START_DEF, help="initial center, cold side of the throat [m]")
    p.add_argument("--theta", type=float, default=THETA_DEF)
    p.add_argument("--T0", type=float, default=-20.0)
    p.add_argument("--alpha", type=float, default=1e-2)
    a = p.parse_args()
    k = consts(a.T0, a.alpha)
    print(f"# d0 {k.d0:.5e} m, s {k.s:.5e} 1/K, beta0 {k.beta0:.5e} s/m, K {k.K:.5e} s/m2; W {a.W * 1e6:g} um,"
          f" L {a.L * 1e6:g} um, lens {a.lens * 1e6:g} um started {a.start * 1e6:g} um on the cold side,"
          f" theta {a.theta:g}, {a.days:g} d")
    for wt in a.wt:
        if a.mode == "ends":
            pore = EndsPore(wt, k, a.W, a.L, a.theta, a.lens, a.start)
            gc = pore.g_crit() if wt < a.W else 0.0
            print(f"wt {wt * 1e6:5.1f} um: G_crit {gc:.4f} K/m")
            for G in list(a.G) + [f * gc for f in a.Gfrac]:
                h, n = pore.h_net(G, a.days)
                r = pore.run(h, G, a.days)
                left = f", a face reached an end at {r['t'] / DAY:.0f} d" if r["t"] < a.days * DAY - 1 else ""
                print(f"   G {G:7.4f} K/m: h_net {h:.10f} ({n} passes); area at t_final {r['dA'] * 100:+.2f} %"
                      f" (largest {r['maxdev'] * 100:.2f} %); lens center {side_words(r['xm'])}{left}")
        else:
            cell = PeriodicCell(wt, k, a.W, a.L, a.theta, a.lens, a.start)
            traj = cell.run(a.days)
            marks = [d for d in (10, 30, 100, 300, 600, 1000) if d <= a.days]
            out = []
            for d in marks:
                t, xm = min(traj, key=lambda q: abs(q[0] - d * DAY))
                out.append(f"{d} d: {side_words(xm)}")
            print(f"wt {wt * 1e6:5.1f} um, periodic: " + "; ".join(out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

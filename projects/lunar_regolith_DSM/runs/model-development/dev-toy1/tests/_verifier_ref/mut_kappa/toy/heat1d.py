"""Planted error: kappa = k/rho (c_p missing)."""
from pathlib import Path

_src = (Path(__file__).resolve().parents[2] / "ok" / "toy" / "heat1d.py").read_text()
exec(compile(_src, "ok_heat1d", "exec"))


def _scales(params):  # noqa: F811
    L, k, rho, cp, q, T0 = (float(params[s]) for s in ("L", "k", "rho", "cp", "q", "T0"))
    kappa = k / rho
    Q = q * L * L / (k * T0)
    return kappa, Q, L * L / kappa, k * T0 / (L * L)

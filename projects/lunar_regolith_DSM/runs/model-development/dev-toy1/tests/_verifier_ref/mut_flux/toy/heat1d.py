"""Planted error: end fluxes without the half-cell source term (one-sided difference only)."""
from pathlib import Path

_src = (Path(__file__).resolve().parents[2] / "ok" / "toy" / "heat1d.py").read_text()
exec(compile(_src, "ok_heat1d", "exec"))


def end_fluxes(U, U_prev, dtau, h, Q, f_end):  # noqa: F811
    F0 = (U[1] - U[0]) / h - h / 2 * (U[0] - U_prev[0]) / dtau
    F1 = (U[-2] - U[-1]) / h - h / 2 * (U[-1] - U_prev[-1]) / dtau
    return F0, F1

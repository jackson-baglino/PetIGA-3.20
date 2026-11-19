"""Planted error: right Dirichlet datum enters row N-1 as 0.9 instead of 1."""
from pathlib import Path

_src = (Path(__file__).resolve().parents[2] / "ok" / "toy" / "heat1d.py").read_text()
exec(compile(_src, "ok_heat1d", "exec"))


def residual(tau, U, Udot, n, params):  # noqa: F811
    _, Q, _, _ = _scales(params)
    h = 1.0 / n
    Ue = np.concatenate(([1.0], np.asarray(U, float), [0.9]))
    F = np.asarray(Udot, float) - (Ue[:-2] - 2 * Ue[1:-1] + Ue[2:]) / h ** 2 - Q
    m = params.get("mms")
    if m is not None:
        xi = np.arange(1, n) * h
        F = F - mms_forcing(xi, tau, m["A"], m["B"], Q)
    return F

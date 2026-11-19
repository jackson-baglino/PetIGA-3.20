"""Planted error: MMS forcing evaluated at tau^n instead of tau^{n+1} (still first order in time)."""
from pathlib import Path

_src = (Path(__file__).resolve().parents[2] / "ok" / "toy" / "heat1d.py").read_text()
exec(compile(_src, "ok_heat1d", "exec"))


def step_be(U, tau_new, dtau, n, params):  # noqa: F811
    shift = 1.0 / dtau
    tau_used = tau_new - dtau
    R = residual(tau_used, U, np.zeros_like(U), n, params)
    sub, diag, sup = jacobian_tridiag(shift, n, params)
    return U - thomas_solve(sub, diag, sup, R)

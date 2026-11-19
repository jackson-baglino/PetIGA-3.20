"""Planted error: the Jacobian omits the shift term (dF/dUdot)."""
from pathlib import Path

_src = (Path(__file__).resolve().parents[2] / "ok" / "toy" / "heat1d.py").read_text()
exec(compile(_src, "ok_heat1d", "exec"))


def jacobian_tridiag(shift, n, params):  # noqa: F811
    h = 1.0 / n
    m = n - 1
    return np.full(m, -1 / h ** 2), np.full(m, 2 / h ** 2), np.full(m, -1 / h ** 2)

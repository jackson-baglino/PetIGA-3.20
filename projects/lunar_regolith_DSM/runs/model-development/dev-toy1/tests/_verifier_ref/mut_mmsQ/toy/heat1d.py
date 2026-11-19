"""Planted error: MMS forcing without its -Q term."""
from pathlib import Path

_src = (Path(__file__).resolve().parents[2] / "ok" / "toy" / "heat1d.py").read_text()
exec(compile(_src, "ok_heat1d", "exec"))
_ok_mms_forcing = mms_forcing


def mms_forcing(xi, tau, A, B, Q):  # noqa: F811
    return _ok_mms_forcing(xi, tau, A, B, Q) + Q

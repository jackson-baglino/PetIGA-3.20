"""Shared postprocessing helpers.

Every function here previously existed as two to five copy-pasted variants
across the scripts in this directory, and every one of them had drifted:

    load_ssa        5 copies, 5 different contracts for the SAME file
    rho_vs          3 copies (and all three were wrong -- see below)
    read_vts        3 copies
    step_times      3 copies
    auto_time_unit  2 copies

Import from here instead of re-deriving. Running a script directly puts this
directory on sys.path, so ``from pplib import rho_vs`` just works.
"""

from __future__ import annotations

import base64
import os
import re
import struct
import xml.etree.ElementTree as ET

import numpy as np


# ---------------------------------------------------------------------------
# SSA_evo.dat
# ---------------------------------------------------------------------------
# Column layout written by monitoring.c. The last three are only present in
# runs produced after those quadrature integrals were added; older files stop
# at DT, and files older still stop at STEP.
SSA, ICE, TIME, STEP, DT, AIR, RHOV, MASS = range(8)

SSA_COLUMNS = """\
  0 sub_interf/eps   ice-air interface density proxy
  1 tot_ice          integrated ice volume
  2 t                [s]
  3 step
  4 dt               [s]   (NaN for legacy 4-column files)
  5 tot_air          only in runs newer than the monitoring.c quadrature change
  6 tot_rhov               "
  7 tot_mass               "
"""


def load_ssa(path: str, min_cols: int = 4):
    """Load SSA_evo.dat, returning an (N, >=5) array or None.

    `path` may be the file itself or the run directory containing it.

    A legacy 4-column file is padded with a NaN dt column so callers can index
    DT unconditionally. Rows are dropped only when one of the first four
    columns is NaN -- never on a NaN in dt or in the newer integral columns,
    which is what an earlier variant of this function got wrong: it filtered on
    `isnan(row).any()` and so discarded *every* row of a padded legacy file.

    Rows are then deduplicated by step number, keeping the last occurrence, to
    guard against repeated monitor calls after a timestep retry.
    """
    if os.path.isdir(path):
        path = os.path.join(path, "SSA_evo.dat")
    if not os.path.isfile(path):
        return None

    try:
        data = np.genfromtxt(path, dtype=float, comments="#", invalid_raise=False)
    except Exception:
        return None

    if data.ndim == 1:
        data = data[np.newaxis, :]
    if data.size == 0 or data.shape[1] < min_cols:
        return None

    if data.shape[1] == 4:
        data = np.hstack([data, np.full((len(data), 1), np.nan)])

    data = data[~np.isnan(data[:, :4]).any(axis=1)]
    if len(data) == 0:
        return None

    steps = data[:, STEP].astype(int)
    _, last = np.unique(steps[::-1], return_index=True)
    return data[np.sort(len(steps) - 1 - last)]


# ---------------------------------------------------------------------------
# Material properties
# ---------------------------------------------------------------------------
# Empirical saturation-pressure coefficients, copied from RhoVS_I in
# src/material_properties.c. Keep them in step with that function: postprocess
# reporting a different rho_vs than the solver used is worse than useless.
_K = (-0.5865e4, 0.2224e2, 0.1375e-1, -0.3403e-4, 0.2697e-7, 0.6918)
_PATM = 101325.0
_BB = 0.62
RHO_AIR_DEFAULT = 1.341          # -rho_air default in lunar_main.c


def rho_vs(T_C, rho_air: float = RHO_AIR_DEFAULT):
    """Saturation vapour density over ice [kg/m^3] at T_C [deg C].

    This mirrors the solver's RhoVS_I exactly.

    The three previous copies of this function all used

        3.25e-3 * exp(-6150 / T_K)

    which is roughly 1e10 times too small -- the prefactor should have been
    ~3e7 for that Clausius-Clapeyron form, so the exponent looks like a typo.
    Because supersaturation is computed as (rhov - rho_vs) / rho_vs, the error
    did not cancel: a saturated run, which should report 0, reported ~9e9.
    Any supersaturation figure produced before 2026-08-07 is wrong.
    """
    T_K = np.asarray(T_C, dtype=float) + 273.15
    Pvs = np.exp(_K[0] / T_K + _K[1] + _K[2] * T_K + _K[3] * T_K ** 2
                 + _K[4] * T_K ** 3 + _K[5] * np.log(T_K))
    return rho_air * _BB * Pvs / (_PATM - Pvs)


def supersaturation(rhov, T_C, rho_air: float = RHO_AIR_DEFAULT):
    """(rhov - rho_vs) / rho_vs, zero where rho_vs is non-positive."""
    rvs = rho_vs(T_C, rho_air)
    return np.where(rvs > 0, (np.asarray(rhov) - rvs) / rvs, 0.0)


# ---------------------------------------------------------------------------
# VTS snapshots
# ---------------------------------------------------------------------------
def _decode(da) -> np.ndarray:
    """Decode a base64 appended-format DataArray element."""
    raw = base64.b64decode("".join(da.text.split()))
    n = struct.unpack("<Q", raw[:8])[0]
    return np.frombuffer(raw[8:8 + n], dtype=np.float64)


def read_vts(fn, want=None):
    """Read a solV_*.vts snapshot.

    Returns (fields, X, Y) where X and Y are (ny, nx) meshgrids. Pass `want`
    to keep only named point-data arrays; the default reads all of them.

    Callers that want 1-D coordinate vectors should slice: X[0, :], Y[:, 0].
    """
    root = ET.parse(fn).getroot()
    grid = root.find(".//StructuredGrid")
    ext = [int(v) for v in grid.get("WholeExtent").split()]
    nx, ny = ext[1] - ext[0] + 1, ext[3] - ext[2] + 1

    pts = _decode(root.find(".//Points/DataArray")).reshape(ny, nx, 3)

    fields = {}
    for da in root.findall(".//PointData/DataArray"):
        name = da.get("Name")
        if want is None or name in want:
            fields[name] = _decode(da).reshape(ny, nx)

    return fields, pts[:, :, 0], pts[:, :, 1]


def step_of(fn) -> int:
    """Step number from a solV_NNNNN.vts filename."""
    return int(re.search(r"solV_(\d+)\.vts", str(fn)).group(1))


def step_times(path) -> dict:
    """Map step -> time [s] from the monitor tables in outp.txt.

    `path` may be the run directory or outp.txt itself. Returns {} when the
    file is absent, so callers can fall back to SSA_evo.dat.
    """
    if os.path.isdir(path):
        path = os.path.join(path, "outp.txt")
    tmap: dict[int, float] = {}
    if not os.path.isfile(path):
        return tmap

    pat = re.compile(r"^\s+(\d+)\s+\|\s+([0-9.eE+-]+)\s+\|")
    with open(path, errors="replace") as fh:
        for line in fh:
            if line.count("|") != 8:
                continue
            m = pat.match(line)
            if m:
                tmap[int(m.group(1))] = float(m.group(2))
    return tmap


# ---------------------------------------------------------------------------
# Axis formatting
# ---------------------------------------------------------------------------
_TIME_UNITS = {"s": 1.0, "min": 60.0, "h": 3600.0, "d": 86400.0}


def auto_time_unit(t_max_sec: float) -> str:
    """Pick a sensible x-axis time unit for a run of this length."""
    if t_max_sec <= 600:
        return "s"
    if t_max_sec <= 7200:
        return "min"
    if t_max_sec <= 3 * 86400:
        return "h"
    return "d"


def in_time_unit(t_sec, unit: str):
    """Convert seconds to `unit` (one of s, min, h, d)."""
    return np.asarray(t_sec, dtype=float) / _TIME_UNITS[unit]

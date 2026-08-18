"""
verify_grain_shrinkage.py — unit gate for postprocess/grain_shrinkage.py.

There is no analytic check available on a real run (the grains are deforming),
so this synthesises an axisymmetric snapshot whose answer is known exactly: two
tangent spheres of radii R0 and R1 on the symmetry axis, with the model's own
logistic interface profile phi = 1/(1 + exp(-(R - d)/eps)).

The volume the script should recover for each sphere is (4/3)*pi*R^3, so the
sphere-equivalent radius it reports must come back as R itself. That exercises
all three things that can silently go wrong:

  * the 2*pi*y axisymmetric measure (a planar integral would be off by orders
    of magnitude, and a 1/(2*pi) slip by 6.28x)
  * the split at the neck plane (mass must not be double-counted or dropped)
  * the volume -> equivalent-radius conversion

The first case supplies no neck_width.csv, so it exercises the DEFAULT split
plane -- which is how the midpoint-between-centres bug was caught: for unequal
radii that plane sits 14 um from the contact and moved 1.3 % of the small
grain's radius, half the shrinkage signal we are trying to fit. Two further
cases shift the plane +-5 um to confirm that once it is near the contact the
answer is insensitive to it.

Run:  python studies/molaro_2019/verification/verify_grain_shrinkage.py
Writes verification.csv next to this file and exits non-zero on failure.
"""

from __future__ import annotations

import base64
import csv
import math
import struct
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
_REPO = _HERE.parent.parent.parent

TOL_PCT = 0.5          # equivalent radius must land within 0.5 % of analytic


def b64_block(arr: np.ndarray) -> str:
    """Same encoding plot_fields.py uses: UInt64 byte count + Float64 payload."""
    data = np.asarray(arr, dtype="<f8").tobytes()
    return base64.b64encode(struct.pack("<Q", len(data)) + data).decode("ascii")


def write_vts(path: Path, X, Y, phi):
    ny, nx = phi.shape
    coords = np.zeros((ny, nx, 3))
    coords[:, :, 0], coords[:, :, 1] = X, Y
    ext = f"0 {nx-1} 0 {ny-1} 0 0"
    path.write_text("\n".join([
        '<?xml version="1.0"?>',
        '<VTKFile type="StructuredGrid" version="0.1" '
        'byte_order="LittleEndian" header_type="UInt64">',
        f'  <StructuredGrid WholeExtent="{ext}">',
        f'    <Piece Extent="{ext}">',
        '      <Points>',
        '        <DataArray type="Float64" NumberOfComponents="3" format="binary">',
        f'          {b64_block(coords.ravel())}',
        '        </DataArray>',
        '      </Points>',
        '      <PointData>',
        f'        <DataArray type="Float64" Name="IcePhase" format="binary">'
        f'{b64_block(phi.ravel())}</DataArray>',
        '      </PointData>',
        '    </Piece>',
        '  </StructuredGrid>',
        '</VTKFile>',
    ]) + "\n")


def build_run(tmp: Path, R0, R1, Lz, Lr, eps, n_per_eps=4):
    """Synthesise a run folder holding one tangent-pair snapshot."""
    d = R0 + R1
    c0 = Lz / 2 - d + R0
    c1 = c0 + d
    h = eps / n_per_eps
    nx, ny = int(Lz / h) + 1, int(Lr / h) + 1
    x = np.linspace(0.0, Lz, nx)
    y = np.linspace(0.0, Lr, ny)
    X, Y = np.meshgrid(x, y)

    # Union signed distance, then the model's equilibrium logistic profile.
    s0 = np.hypot(X - c0, Y) - R0
    s1 = np.hypot(X - c1, Y) - R1
    phi = 1.0 / (1.0 + np.exp(np.minimum(np.minimum(s0, s1) / eps, 700.0)))

    (tmp / "vtkOut").mkdir(parents=True, exist_ok=True)
    write_vts(tmp / "vtkOut" / "solV_00000.vts", X, Y, phi)
    (tmp / "geom.opts").write_text(
        f"-axisym 1\n-ice_grain_cx {c0:.8e},{c1:.8e}\n"
        f"-ice_grain_R  {R0:.8e},{R1:.8e}\n-eps {eps:.8e}\n")
    # step_times() reads this table; one row is enough.
    (tmp / "outp.txt").write_text(
        "  step |        t |       dt\n"
        "     0 | 0.000e+00 | 1.000e-04\n")
    return c0, c1


def run_case(label, R0, R1, Lz, Lr, eps, split_shift=0.0):
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        c0, c1 = build_run(tmp, R0, R1, Lz, Lr, eps)
        if split_shift:
            # Force an off-centre split by planting a neck_width.csv.
            (tmp / "neck_width.csv").write_text(
                "t_s,neck_width_m,x_neck_m\n"
                f"0.0,0.0,{c0 + R0 + split_shift:.8e}\n")
        out = subprocess.run(
            [sys.executable, str(_REPO / "postprocess/grain_shrinkage.py"),
             str(tmp), "--no-plot"],
            capture_output=True, text=True)
        if out.returncode != 0:
            print(out.stdout + out.stderr)
            raise SystemExit(f"{label}: grain_shrinkage.py failed")
        with (tmp / "grain_shrinkage.csv").open() as fh:
            row = next(csv.DictReader(fh))

    R_lg, R_sm = float(row["R_large_m"]), float(row["R_small_m"])
    exp_lg, exp_sm = max(R0, R1), min(R0, R1)
    e_lg = 100.0 * (R_lg / exp_lg - 1.0)
    e_sm = 100.0 * (R_sm / exp_sm - 1.0)
    ok = abs(e_lg) < TOL_PCT and abs(e_sm) < TOL_PCT
    print(f"  {label:34s} large {R_lg*1e6:7.3f} vs {exp_lg*1e6:7.3f} um ({e_lg:+.3f} %)"
          f"   small {R_sm*1e6:7.3f} vs {exp_sm*1e6:7.3f} um ({e_sm:+.3f} %)"
          f"   [{'PASS' if ok else 'FAIL'}]")
    return dict(case=label, R_large_m=R_lg, R_small_m=R_sm,
                expect_large_m=exp_lg, expect_small_m=exp_sm,
                err_large_pct=e_lg, err_small_pct=e_sm, passed=ok)


def main():
    R0, R1, eps = 7.25e-5, 1.01e-4, 2.5852e-7
    Lz, Lr = 6.80e-4, 3.40e-4
    print("verify_grain_shrinkage.py — analytic tangent-pair gate")
    print(f"  R0 = {R0*1e6:.1f} um   R1 = {R1*1e6:.1f} um   eps = {eps:.4e} m")
    print(f"  domain {Lz*1e6:.0f} x {Lr*1e6:.0f} um, axisymmetric,"
          f" 4 points per eps\n")
    results = [
        # No neck_width.csv -> exercises the DEFAULT split plane. This is the
        # case that caught the midpoint-between-centres bug: for these radii
        # that plane sits 14 um from the contact and biased the small grain by
        # +1.3 %. The radical plane puts it exactly on the contact.
        run_case("default split (radical plane)", R0, R1, Lz, Lr, eps),
        run_case("split shifted +5 um", R0, R1, Lz, Lr, eps, split_shift=+5e-6),
        run_case("split shifted -5 um", R0, R1, Lz, Lr, eps, split_shift=-5e-6),
    ]
    with (_HERE / "verification.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(results[0]), lineterminator="\n")
        w.writeheader()
        for r in results:
            w.writerow(r)
    print(f"\n  wrote {_HERE / 'verification.csv'}")
    if not all(r["passed"] for r in results):
        raise SystemExit("FAILED")
    print("  all cases within "
          f"{TOL_PCT} % of the analytic sphere radii")


if __name__ == "__main__":
    main()

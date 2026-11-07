"""Write the geometry, experiment and tests files for plan types 3 and 4.

Type 3  lens in a flat pore (contact angle): pore 500 um long, H = 125 um
        (60 um for 3B), lens 150 um long, centered, flat faces at t = 0
        (the -bridge_* slab trick of channel_2D_H125um_eps0.86um_flat90.opts).
Type 4  lens in a wedge (wall slope): half-angle omega = 2.5, 5, 10, 20 deg,
        gap at the lens center 100 um at t = 0, every lens 1.5e-8 m^2 (the
        area of a 100 x 150 um lens); 4A-0 is a flat 100 um pore with the
        same lens. The wedge lens is the apex-centered band (90 deg at t = 0),
        built by preprocess/build_geometry_wedge.py.

eps = 8.5840e-7 m for every run: the velocity-study value, so the 3A reruns
double as a regression check against it. comp_eps at the type-3 speed allows
up to 2.96e-6 m (B-HEAT bound) but then the realised attachment coefficient
is only 0.61 of the target; at 8.584e-7 it is within the comp_eps 0.9 check.
tau_sub = 670.6 s (mesh_params.derived_pf_params), dtmax = 2 tau_sub.

Contact angle from Young: cos(theta) = (gamma_as - gamma_is)/gamma_ia with
gamma_ia = 0.109, gamma_as = 0.3 (the stage-0 values), so
gamma_is = 0.3 - 0.109 cos(theta).

Every number is printed. Run from the project root:

    venv_cryosuction/bin/python studies/cryosuction/scripts/make_type34_inputs.py [--dry-run]
"""
from __future__ import annotations

import argparse
import math
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from preprocess.mesh_params import physics  # noqa: E402
from preprocess.mesh_params.constants import CONSTANTS as C  # noqa: E402

T0, ALPHA, EPS = -20.0, 1e-2, 8.5840e-7
BETA0 = 7.9408e5                       # comp_eps beta_sub at -20 C, alpha_c 1e-2
G_IA, G_AS = 0.109, 0.3
DAY = 86400.0
GEO = ROOT / "inputs/geometry/cryosuction"
EXP = ROOT / "inputs/experiment/cryosuction"
MESH = ROOT / "inputs/geometry/meshes"
STUDY = ROOT / "studies/cryosuction"

LENS_L, LENS_AREA, GAP0 = 1.5e-4, 1.5e-8, 1.0e-4
END_GAP = 1.75e-4                      # type 3: pore length (500 um) minus lens, split evenly
WIDE_MARGIN, NARROW_MARGIN, MIN_THROAT = 6.0e-5, 6.0e-5, 2.5e-5


def tau_sub(eps):
    d = physics.derived_pf_params(eps, T0, physics.capillary_length_m(T0),
                                  physics.beta_hk_s_per_m(T0, ALPHA),
                                  C.k_ice_w_per_m_k / C.c_ice_j_per_m3_k,
                                  C.k_air_w_per_m_k / C.c_air_j_per_m3_k, physics.dv_m2_per_s(T0))
    return d["tau_sub"]


def gamma_is(theta_deg):
    return G_AS - G_IA * math.cos(math.radians(theta_deg))


def write(path: Path, text: str, dry: bool):
    print(f"  {'(dry) ' if dry else ''}wrote {path.relative_to(ROOT)}")
    if not dry:
        path.write_text(text)


# ---------------------------------------------------------------- geometry
def channel_geometry(name, H, L, x_a, x_b, dtmax, purpose, dry):
    """Flat pore with a flat-faced slab x_a < x < x_b (crossed far discs)."""
    R = 50.0 * max(H, 1.25e-4)             # face sag (H/2)^2/(2R) < 0.04 um
    cxL, cxR = x_a + R, x_b - R
    Nx = math.ceil(L * math.sqrt(2) / EPS)
    Ny = math.ceil(H * math.sqrt(2) / EPS)
    sag = (0.5 * H) ** 2 / (2 * R)
    text = f"""# =============================================================================
# geometry/cryosuction/{name}.opts   (make_type34_inputs.py, 2026-11-07)
#
# {purpose}
# Flat pore {L*1e6:.0f} x {H*1e6:.0f} um, walls y0 and y1. Lens = slab {x_a*1e6:.1f}-{x_b*1e6:.1f} um
# ({(x_b-x_a)*1e6:.0f} um long, area {(x_b-x_a)*H:.4e} m^2) with flat faces (90 deg at t = 0): the
# intersection of two far, crossed discs (R = {R:.3e} m, face sag {sag*1e6:.3f} um), as in
# contactangle/channel_2D_H125um_eps0.86um_flat90.opts.
# eps {EPS:.4e} m (velocity-study value; see make_type34_inputs.py), h = eps/sqrt(2).
# dtmax = 2 tau_sub = {dtmax:.4e} s.
# =============================================================================
-ic_type multi_grains
-dim 2
-Nx {Nx}
-Ny {Ny}
-Nz 1
-Lx {L:.6e}
-Ly {H:.6e}
-Lz 0
-bridge_cxL {cxL:.6e}
-bridge_rL  {R:.6e}
-bridge_cxR {cxR:.6e}
-bridge_rR  {R:.6e}
-bridge_cy  {0.5*H:.6e}
-eps {EPS:.4e}
-eps_valid_temp {T0:g}
-delt_t 1.0e-4
-periodic 0
-wall_faces y0,y1
-dtmax {dtmax:.4e}
"""
    write(GEO / f"{name}.opts", text, dry)
    print(f"    Nx {Nx} Ny {Ny} ({(Nx+2)*(Ny+2)*3/1e6:.2f} M DOF)")


def wedge_geometry(omega_deg, dtmax, dry):
    m = math.tan(math.radians(omega_deg))
    w = math.radians(omega_deg)
    r_c = 0.5 * GAP0 / m                      # apex distance of the lens center
    dr = LENS_AREA / (2.0 * w * r_c)          # band thickness for the target area
    narrow = min(NARROW_MARGIN, (GAP0 - MIN_THROAT) / (2 * m) - 0.5 * dr)
    x_c = 0.5 * dr + narrow
    Lx = dr + narrow + WIDE_MARGIN
    w_t = GAP0 - 2 * m * x_c
    w_m = GAP0 + 2 * m * (Lx - x_c)
    tag = f"om{omega_deg:g}"
    name = f"wedge_2D_om{omega_deg:g}deg_gap100um_eps0.86um"
    print(f"  omega {omega_deg:g} deg: r_c {r_c*1e6:.1f} um, band {dr*1e6:.1f} um "
          f"(r {(r_c-dr/2)*1e6:.1f}-{(r_c+dr/2)*1e6:.1f}), Lx {Lx*1e6:.1f} um, width {w_t*1e6:.1f} -> "
          f"{w_m*1e6:.1f} um, narrow margin {narrow*1e6:.1f} um")
    dat = MESH / f"wedge_{tag}_gap100um.dat"
    opts = GEO / f"{name}.opts"
    cmd = [sys.executable, str(ROOT / "preprocess/build_geometry_wedge.py"),
           "--Lx", f"{Lx:.6e}", "--Ly", f"{w_m:.6e}", "--w-throat", f"{w_t:.6e}",
           "--w-mouth", f"{w_m:.6e}", "--ice-shape", "band", "--ice-x", f"{x_c/Lx:.8f}",
           "--band-dr", f"{dr:.6e}", "--eps", f"{EPS:.4e}", "--tag", tag,
           "--out", str(dat), "--opts", str(opts)]
    print("    " + " ".join(cmd[1:]).replace(str(ROOT) + "/", ""))
    if dry:
        return name
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=ROOT)
    if r.returncode != 0:
        raise SystemExit(f"wedge builder failed for omega {omega_deg}:\n{r.stdout[-1500:]}\n{r.stderr[-1500:]}")
    head = (f"# Plan type 4 (make_type34_inputs.py, 2026-11-07): wedge half-angle {omega_deg:g} deg,\n"
            f"# gap {GAP0*1e6:.0f} um at the lens center at t = 0, lens area {LENS_AREA:.2e} m^2 (90 deg band).\n"
            "# Walls are the mesh's y0 and y1 faces; x0 (narrow) and x1 (wide) are vapor faces.\n")
    body = opts.read_text()
    opts.write_text(head + body + f"-wall_faces y0,y1\n-dtmax {dtmax:.4e}\n")
    return name


# ---------------------------------------------------------------- experiments
def experiment(name, theta, vapor, days, purpose, dry):
    if vapor == "res":
        vap = ("-flag_BC_rhovfix 1\n-rhovfix_axis 0                         # both x faces, h = 1\n"
               "-rhovfix_lo 1.000000\n-rhovfix_hi 1.000000\n")
    elif vapor == "top":
        vap = ("-flag_BC_rhovfix 1\n-rhovfix_faces x0                       # supply at x0 only; x1 zero flux\n"
               "-rhovfix_lo 1.000000\n")
    elif vapor == "closed":
        vap = "-flag_BC_rhovfix 0                      # no vapor in or out: water conserved\n"
    else:
        raise ValueError(vapor)
    gis = gamma_is(theta)
    text = f"""# experiment/cryosuction/{name}.opts   (make_type34_inputs.py, 2026-11-07)
#
# {purpose}
# Isothermal at -20 C (T fixed on the boundary), theta = {theta:g} deg from Young
# (gamma_is = 0.3 - 0.109 cos theta), vapor set "{vapor}", {days:g} days.
# -dtmax is in the GEOMETRY file (2 tau_sub depends on its eps).
#
-temp -20.0
-humidity 1.0
-grad_temp0 0.0,0.0,0.0
-flag_BC_Tfix 1
{vap}-phi_open_bc 0                          # 90-degree phi face on every non-wall face
-t_final {days*DAY:.6e}
-n_out {int(days) + 1}                              # one snapshot a day
-beta_sub0 {BETA0:.4e}                     # K&P beta0 at -20 C, alpha_c = 1e-2
-gamma_ia {G_IA:g}
-gamma_as {G_AS:g}
-gamma_is {gis:.7f}                     # cos(theta) = {math.cos(math.radians(theta)):+.6f}
-dtCFL 1
-outp 0
"""
    write(EXP / f"{name}.opts", text, dry)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true", help="print what would be written")
    a = ap.parse_args()
    dry = a.dry_run
    tau = tau_sub(EPS)
    dtmax = 2.0 * tau
    print(f"eps {EPS:.4e} m, tau_sub {tau:.4e} s, dtmax {dtmax:.4e} s; d0 {physics.capillary_length_m(T0):.4e} m")

    print("\nType 3 geometry")
    L3 = LENS_L + 2 * END_GAP
    channel_geometry("channel_2D_H125um_L500um_eps0.86um_lens150", 1.25e-4, L3, END_GAP, END_GAP + LENS_L,
                     dtmax, "Plan type 3 (3A, 3C, 3D): lens in a flat 125 um pore.", dry)
    channel_geometry("channel_2D_H60um_L500um_eps0.86um_lens150", 6.0e-5, L3, END_GAP, END_GAP + LENS_L,
                     dtmax, "Plan type 3B: the same lens length in a 60 um pore (1/H scaling).", dry)
    print("\nType 4 geometry")
    L4 = LENS_L + 2 * WIDE_MARGIN
    channel_geometry("channel_2D_H100um_L270um_eps0.86um_lens150", GAP0, L4, WIDE_MARGIN, WIDE_MARGIN + LENS_L,
                     dtmax, "Plan type 4A-0: flat 100 um pore, the omega = 0 reference for the wedges.", dry)
    wedges = {om: wedge_geometry(om, dtmax, dry) for om in (2.5, 5.0, 10.0, 20.0)}

    print("\nExperiments")
    for th in (30, 60, 90, 120, 150):
        experiment(f"t3res_T-20_theta{th}_120d_ac1e-2", th, "res", 120,
                   "Plan type 3A/3B: lens grows or shrinks by contact angle alone; h = 1 at both ends.", dry)
    experiment("t3closed_T-20_theta60_60d_ac1e-2", 60, "closed", 60,
               "Plan type 3C: closed pore; nothing should move (control).", dry)
    for th in (60, 120):
        experiment(f"t3top_T-20_theta{th}_120d_ac1e-2", th, "top", 120,
                   "Plan type 3D: supply at x0 only, x1 closed; the lens should migrate.", dry)
    for th in (60, 85, 90):
        experiment(f"t4res_T-20_theta{th}_180d_ac1e-2", th, "res", 180,
                   "Plan type 4: lens in a wedge or the flat reference; h = 1 at both ends.", dry)

    H125, H60, H100 = ("channel_2D_H125um_L500um_eps0.86um_lens150", "channel_2D_H60um_L500um_eps0.86um_lens150",
                       "channel_2D_H100um_L270um_eps0.86um_lens150")
    lines3 = ["# Plan type 3 (make_type34_inputs.py). RUNS.md IDs in the comments."]
    lines3 += [f"{H125}:t3res_T-20_theta{th}_120d_ac1e-2   # 3A-{i}" for i, th in enumerate((30, 60, 90, 120, 150), 1)]
    lines3 += [f"{H60}:t3res_T-20_theta{th}_120d_ac1e-2   # 3B-{i}" for i, th in enumerate((60, 120), 1)]
    lines3 += [f"{H125}:t3closed_T-20_theta60_60d_ac1e-2   # 3C-1"]
    lines3 += [f"{H125}:t3top_T-20_theta{th}_120d_ac1e-2   # 3D-{i}" for i, th in enumerate((60, 120), 1)]
    lines4 = ["# Plan type 4 (make_type34_inputs.py). RUNS.md IDs in the comments."]
    lines4 += [f"{H100}:t4res_T-20_theta60_180d_ac1e-2   # 4A-0"]
    lines4 += [f"{wedges[om]}:t4res_T-20_theta60_180d_ac1e-2   # 4A-{i}" for i, om in enumerate((2.5, 5.0, 10.0, 20.0), 1)]
    lines4 += [f"{wedges[10.0]}:t4res_T-20_theta90_180d_ac1e-2   # 4B-1",
               f"{wedges[10.0]}:t4res_T-20_theta85_180d_ac1e-2   # 4B-2"]
    print("\nTests files")
    write(STUDY / "type3_tests.txt", "\n".join(lines3) + "\n", dry)
    write(STUDY / "type4_tests.txt", "\n".join(lines4) + "\n", dry)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

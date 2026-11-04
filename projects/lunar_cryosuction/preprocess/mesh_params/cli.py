"""Command line for the mesh_params tool.

  python -m preprocess.mesh_params --Lx 6e-4 --Ly 6e-4 --T0 -10 --alpha 1e-3          # readable summary
  python -m preprocess.mesh_params --json '{"domain_x_m": 6e-4, "alpha_c": 1e-3}'     # JSON in, JSON out
  python -m preprocess.mesh_params --Lx 6e-4 --T0 -20 --alpha_range 1e-4 1e-1 12      # alpha_c sweep
  python -m preprocess.mesh_params --Lx 6e-4 --T0 -5 --alpha 1e-3 --patch inputs/geometry/x.opts

preprocess/comp_eps.py forwards here, so the old flags keep working.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

from .tool import MeshParamsInput, ToolInputError, compute_mesh_params, sweep_alpha_c, TOOL_SPEC


def format_text(r: dict) -> str:
    i, m, rs, b, d, t = r["inputs"], r["mesh"], r["resolution"], r["eps_bounds_m"], r["derived"], r["thermo"]
    L = []
    L.append(f"mesh_params v{r['version']}   T = {i['temperature_c']:.1f} C   alpha_c = {i['alpha_c_used']:.3e}   "
             f"safety = {i['safety']:.2f}   channel = {i['thermal_channel']}")
    L.append("")
    L.append("eps bounds [m]  (KP eqs. 43-45 and geometric)")
    for k in ("B-HEAT", "B-VAPOR", "B-KINETIC", "B-CURV"):
        L.append(f"  {k:<10s} {b[k]:.4e}" + ("   <- binding" if k == b["binding"] else ""))
    L.append("")
    L.append(f"eps = {r['eps_m']:.4e} m     h = eps/sqrt(2) = {m['h_m']:.4e} m")
    n = f"Nx = {m['nx']}" + (f", Ny = {m['ny']}" if m["ny"] else "") + (f", Nz = {m['nz']}" if m["nz"] else "")
    L.append(f"{n}     elements total {m['elements_total']:,}     DOFs (p=2, C=1, 3 DOF) ~ {m['dofs_two_phase_3dof']:,}")
    L.append(f"visible band: {rs['elements_across_5_95']:.1f} elements across 5-95 %, {rs['elements_across_1_99']:.1f} across 1-99 %  "
             f"(eps is a decay length, not the band width)")
    L.append("")
    L.append("derived (MF SI Eq. 9):  tau_sub = {:.4e} s   lam_sub = {:.4e}   M_sub = {:.4e}   alpha_src = {:.4e}".format(
        d["tau_sub_s"], d["lam_sub"], d["M_sub"], d["alpha_src_per_s"]))
    L.append(f"  kinetic fraction of tau_sub {d['kinetic_fraction']:.1%}   min beta_eff/beta_target {d['min_beta_ratio']:.3f}   "
             f"Pe = {t['peclet']:.2e}")
    L.append(f"thermo: rho_vs = {t['rho_vs_kg_per_m3']:.4e} kg/m3   Dv(air) = {t['dv_m2_per_s']:.4e} m2/s   "
             f"beta_sub = {t['beta_sub_s_per_m']:.4e} s/m   d0 = {t['d0_m']:.4e} m")
    L.append("")
    L.append("checks")
    for c in r["checks"]:
        L.append(f"  [{c['status'].upper():<4s}] {c['name']}: {c['message']}")
    L.append("")
    L.append("solver options: " + " ".join(f"{k} {v}" for k, v in r["solver_options"].items()))
    return "\n".join(L)


def patch_opts(path: Path, options: dict) -> None:
    """Set (or append) exactly the options in `options` in an opts file. Only -eps and -N* are ever passed."""
    text = path.read_text()
    for key, val in options.items():
        pat = rf"^({re.escape(key)}\s+)\S+(.*)$"
        text = re.sub(pat, rf"{key} {val}\2", text, flags=re.M) if re.search(pat, text, flags=re.M) else text.rstrip() + f"\n{key} {val}\n"
    path.write_text(text)
    print(f"patched {path}: " + " ".join(f"{k} {v}" for k, v in options.items()), file=sys.stderr)


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description="Compute eps, mesh and derived phase-field parameters (KP 2009, MF 2024).")
    ap.add_argument("--json", help="JSON object with MeshParamsInput fields (domain_x_m, ...); prints JSON")
    ap.add_argument("--json-file", type=Path, help="Same as --json, read from a file")
    ap.add_argument("--spec", action="store_true", help="Print the tool spec (JSON schema) and exit")
    ap.add_argument("--Lx", type=float, help="Domain length x [m]")
    ap.add_argument("--Ly", type=float, default=0.0)
    ap.add_argument("--Lz", type=float, default=0.0)
    ap.add_argument("--Rave", type=float, default=3.0e-5, help="Smallest curvature radius [m]")
    ap.add_argument("--T0", type=float, default=-20.0, help="Temperature [C]")
    ap.add_argument("--alpha", type=float, default=None, help="alpha_c")
    ap.add_argument("--alpha_range", nargs=3, metavar=("LO", "HI", "N"), default=None)
    ap.add_argument("--sigma_surf", type=float, default=None)
    ap.add_argument("--vn", type=float, default=1.0e-9)
    ap.add_argument("--vn_feature", type=float, default=None, metavar="R_FEAT")
    ap.add_argument("--Dchannel", choices=("mean", "ice"), default="mean")
    ap.add_argument("--xiv", type=float, default=None)
    ap.add_argument("--safety", type=float, default=0.5)
    ap.add_argument("--dim", type=int, default=None)
    ap.add_argument("--patch", type=Path, default=None, help="Opts file: set -eps and -Nx/-Ny/-Nz only")
    ap.add_argument("--quiet", action="store_true", help="Print only eps and the element counts, one per line")
    ap.add_argument("--as-json", action="store_true", help="Print JSON instead of text for the flag-style call")
    return ap


def main(argv=None) -> int:
    ap = build_parser()
    a = ap.parse_args(argv)
    try:
        if a.spec:
            print(json.dumps(TOOL_SPEC, indent=2))
            return 0
        if a.json or a.json_file:
            raw = json.loads(a.json if a.json else a.json_file.read_text())
            print(json.dumps(compute_mesh_params(MeshParamsInput(**raw)), indent=2))
            return 0
        if a.Lx is None:
            ap.error("give --Lx (flag style) or --json / --json-file")
        n = sum(x is not None for x in (a.alpha, a.alpha_range, a.sigma_surf))
        if n != 1:
            ap.error("give exactly one of --alpha, --alpha_range, --sigma_surf")
        base = dict(domain_x_m=a.Lx, domain_y_m=a.Ly, domain_z_m=a.Lz, grain_radius_m=a.Rave, temperature_c=a.T0,
                    front_velocity_m_per_s=a.vn, feature_radius_m=a.vn_feature, safety=a.safety,
                    thermal_channel=a.Dchannel, xi_v=a.xiv, dimension=a.dim)
        if a.alpha_range:
            lo, hi, k = float(a.alpha_range[0]), float(a.alpha_range[1]), int(a.alpha_range[2])
            rows = sweep_alpha_c(MeshParamsInput(**base), lo, hi, k)
            print(f"{'alpha_c':>10s} {'eps [m]':>11s} {'Nx':>6s} {'binding':>10s} {'tau_sub [s]':>12s} {'min beta ratio':>15s}")
            for al, r in rows:
                print(f"{al:10.3e} {r['eps_m']:11.4e} {r['mesh']['nx']:6d} {r['eps_bounds_m']['binding']:>10s} "
                      f"{r['derived']['tau_sub_s']:12.4e} {r['derived']['min_beta_ratio']:15.3f}")
            return 0
        r = compute_mesh_params(MeshParamsInput(**base, alpha_c=a.alpha, sigma_surf=a.sigma_surf))
    except (ToolInputError, json.JSONDecodeError, TypeError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 2
    if a.quiet:
        print(f"{r['eps_m']:.4e}")
        for k in ("nx", "ny", "nz"):
            if r["mesh"][k]:
                print(r["mesh"][k])
    elif a.as_json:
        print(json.dumps(r, indent=2))
    else:
        print(format_text(r))
    if a.patch:
        patch_opts(a.patch, r["solver_options"])
    return 0

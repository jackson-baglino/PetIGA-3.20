#!/usr/bin/env python3
"""
PETSc Resource Planner — env-sourced
-----------------------------------
Estimate memory, ranks, nodes, and mem-per-cpu for a PETSc job.

This version **sources** one of your `./configs/*.env` files (without editing it)
so its key=value pairs become available as environment variables inside Python.
You can optionally override grid dimensions from the command line.

Primary inputs (CLI):
  --env-name <name>     (loads ./configs/<name>.env by default)
  --env <path>          (alternative: load env from an explicit path)
  --Nx --Ny [--Nz]      (optional overrides for grid; if omitted, read from env)
  --dim {2,3}           (defaults to 2 if Nz not provided; 3 means Nz is required)

Examples
--------
# Use configs/grainReadFile-35_s1-10.env and read Nx,Ny,Nz from it
python petsc_resource_planner.py --env-name grainReadFile-35_s1-10

# Same, but override Nz on the CLI (e.g., try a thinner domain)
python petsc_resource_planner.py --env-name grainReadFile-35_s1-10 --dim 3 --Nz 64

# Provide a full path to an env file
python petsc_resource_planner.py --env ./configs/grainReadFile-35_s1-10.env

Notes
-----
- We assume 3 fields by default (your case) and 1 scalar dof per element per field.
  Adjust with --fields and --dofs-per-elem-per-field if needed.
- We assume Nx*Ny*Nz counts **elements**. If your discretization differs, tune
  dofs-per-elem-per-field accordingly.
- We load the env file into os.environ, but we do NOT modify the file.
"""

from __future__ import annotations
import argparse
import math
import os
from pathlib import Path
from typing import Dict, Tuple

# ------------------------
# Utilities
# ------------------------

def humanize_bytes(num_bytes: float) -> str:
    units = ["B","KB","MB","GB","TB","PB"]
    val = float(num_bytes)
    for u in units:
        if val < 1024 or u == units[-1]:
            return f"{val:.2f} {u}"
        val /= 1024.0


def source_env_file(env_path: Path) -> Dict[str, str]:
    """"Source" a simple KEY=VALUE .env into the current process environment.
    Does not execute shell; just parses lines and injects into os.environ.
    Returns a dict of parsed key->value (original casing preserved).
    """
    if not env_path.exists():
        raise FileNotFoundError(f"Env file not found: {env_path}")
    parsed: Dict[str, str] = {}
    for raw in env_path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith('#'):
            continue
        if '=' not in line:
            continue
        k, v = line.split('=', 1)
        k = k.strip()
        v = v.strip()
        # Strip optional quotes
        if (v.startswith('"') and v.endswith('"')) or (v.startswith("'") and v.endswith("'")):
            v = v[1:-1]
        # Save into env (also keep both original and upper/lower variants for convenience)
        os.environ[k] = v
        os.environ.setdefault(k.upper(), v)
        os.environ.setdefault(k.lower(), v)
        parsed[k] = v
    return parsed


def get_int_env(name: str, default: int | None = None) -> int | None:
    s = os.environ.get(name)
    if s is None:
        s = os.environ.get(name.upper())
    if s is None:
        s = os.environ.get(name.lower())
    if s is None:
        return default
    try:
        return int(float(s))
    except Exception:
        return default


# ------------------------
# Memory / ranks planning
# ------------------------

def estimate_memory_bytes(total_dofs: int,
                          nnz_per_row: float,
                          scalar_bytes: int,
                          index_bytes: int,
                          ksp_pc_overhead: float,
                          num_vectors: int,
                          safety_factor: float) -> dict:
    rows = total_dofs
    nnz = rows * nnz_per_row

    mat_values = nnz * scalar_bytes
    mat_col_idx = nnz * index_bytes
    mat_row_ptr = (rows + 1) * index_bytes
    mat_bytes = mat_values + mat_col_idx + mat_row_ptr

    mat_total = mat_bytes * ksp_pc_overhead
    vec_total = num_vectors * rows * scalar_bytes
    total = (mat_total + vec_total) * safety_factor

    return {
        "rows": rows,
        "nnz": nnz,
        "mat_bytes": mat_bytes,
        "mat_total_bytes": mat_total,
        "vec_total_bytes": vec_total,
        "total_bytes": total,
    }


def plan_ranks(total_bytes: float,
               total_rows: int,
               target_rows_per_rank: int,
               min_rows_per_rank: int,
               cores_per_node: int | None,
               node_mem_gb: float | None,
               max_ranks: int):
    total_gb = total_bytes / (1024**3)

    # Step 1: choose ranks by DOFs per rank
    ntasks = max(1, min(max_ranks, math.ceil(total_rows / max(target_rows_per_rank, 1))))

    # Compute implied mem per rank
    mem_per_rank_gb = total_gb / ntasks

    # Optional sanity warning (not fatal): rows per rank too small
    rows_per_rank = total_rows / ntasks
    warn_over_parallel = rows_per_rank < max(1, min_rows_per_rank)

    # If no cluster info, return cluster-agnostic plan
    if not cores_per_node or not node_mem_gb:
        return {
            "ntasks": ntasks,
            "nodes": None,
            "ntasks_per_node_list": None,
            "mem_per_cpu_gb": mem_per_rank_gb,
            "rows_per_rank": rows_per_rank,
            "warn_over_parallel": warn_over_parallel,
            "ntasks_per_node": None,
        }

    # Step 2: pack onto nodes with both core & memory constraints
    rpn_core = cores_per_node
    rpn_mem = max(1, int(node_mem_gb // max(mem_per_rank_gb, 1e-9)))
    rpn = max(1, min(rpn_core, rpn_mem))

    nodes = math.ceil(ntasks / rpn)
    base = ntasks // nodes
    rem = ntasks % nodes
    ntasks_per_node_list = [base + 1] * rem + [base] * (nodes - rem)

    # Ensure per-node memory respected (repartition if needed)
    max_rpn_actual = max(ntasks_per_node_list) if ntasks_per_node_list else rpn
    if mem_per_rank_gb * max_rpn_actual > node_mem_gb + 1e-9:
        rpn_limit = max(1, int(node_mem_gb // max(mem_per_rank_gb, 1e-9)))
        rpn = min(rpn, rpn_limit)
        nodes = math.ceil(ntasks / rpn)
        base = ntasks // nodes
        rem = ntasks % nodes
        ntasks_per_node_list = [base + 1] * rem + [base] * (nodes - rem)

    return {
        "ntasks": ntasks,
        "nodes": nodes,
        "ntasks_per_node_list": ntasks_per_node_list,
        "mem_per_cpu_gb": mem_per_rank_gb,
        "rows_per_rank": rows_per_rank,
        "warn_over_parallel": warn_over_parallel,
        "ntasks_per_node": rpn if rem == 0 else None,
    }


# ------------------------
# CLI
# ------------------------

def parse_args():
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # Primary interface the user asked for
    ap.add_argument("--env", type=str, default=None, help="Path to a .env file (e.g., ./configs/foo.env)")
    ap.add_argument("--env-name", type=str, default=None, help="Name (without .env) to load from --configs-dir")
    ap.add_argument("--configs-dir", type=str, default="./configs", help="Directory containing .env files")
    ap.add_argument("--dim", type=int, choices=[2,3], default=None, help="Problem dimension; if omitted, infer from presence of Nz")
    ap.add_argument("--Nx", type=int, default=None, help="Override Nx (elements in x)")
    ap.add_argument("--Ny", type=int, default=None, help="Override Ny (elements in y)")
    ap.add_argument("--Nz", type=int, default=None, help="Override Nz (elements in z; required if --dim 3 and env lacks Nz)")

    # Expert knobs (retain from earlier version)
    ap.add_argument("--fields", type=int, default=3, help="Number of coupled fields (e.g., 3)")
    ap.add_argument("--dofs-per-elem-per-field", type=float, default=1.0,
                    help="Scalar dofs per element per field (depends on discretization)")
    ap.add_argument("--nnz-per-row", type=float, default=60.0,
                    help="Average nonzeros per scalar row")
    ap.add_argument("--scalar-bytes", type=int, default=8, help="Bytes per scalar (8 for double)")
    ap.add_argument("--index-bytes", type=int, default=8, help="Bytes per matrix index (8 for 64-bit)")
    ap.add_argument("--ksp-pc-overhead", type=float, default=1.20, help="Multiplier for PC/Krylov overhead")
    ap.add_argument("--num-vectors", type=int, default=6, help="Work vectors budget")
    ap.add_argument("--safety-factor", type=float, default=2.0, help="Global safety multiplier")

    # Cluster knobs
    ap.add_argument("--node-mem", type=float, default=192.0, help="Node memory in GB")
    ap.add_argument("--cores-per-node", type=int, default=40, help="MPI ranks per node")
    ap.add_argument("--target-mem-per-cpu", type=float, default=None, help="[DEPRECATED] Prefer --desired-mem-per-rank")
    ap.add_argument("--max-ranks", type=int, default=2000, help="Upper bound for ranks to consider")

    # New heuristic knobs
    ap.add_argument("--target-rows-per-rank", type=int, default=150_000, help="Heuristic target scalar rows (DOFs) per MPI rank. Drives ntasks selection.")
    ap.add_argument("--min-rows-per-rank", type=int, default=80_000, help="Lower bound before we warn about over-parallelization.")

    return ap.parse_args()


def resolve_env_path(args) -> Path | None:
    if args.env:
        return Path(args.env)
    if args.env_name:
        return Path(args.configs_dir) / f"{args.env_name}.env"
    return None


# ------------------------
# Main
# ------------------------

def main():
    args = parse_args()

    env_path = resolve_env_path(args)
    if env_path is not None:
        source_env_file(env_path)

    # Determine dimensions from CLI overrides first, then env
    nx = args.Nx if args.Nx is not None else get_int_env('Nx')
    ny = args.Ny if args.Ny is not None else get_int_env('Ny')
    nz_env = get_int_env('Nz')
    nz = args.Nz if args.Nz is not None else nz_env

    # Infer dim if not provided
    dim = args.dim
    if dim is None:
        dim = 3 if nz is not None and nz > 1 else 2

    if nx is None or ny is None:
        raise SystemExit("Nx and Ny must be provided (via CLI or env file).")
    if dim == 3:
        if nz is None:
            raise SystemExit("Nz must be provided for dim=3 (via CLI or env file).")
    else:
        nz = 1

    elements = nx * ny * nz
    total_dofs = int(elements * args.fields * args.dofs_per_elem_per_field)

    est = estimate_memory_bytes(
        total_dofs=total_dofs,
        nnz_per_row=args.nnz_per_row,
        scalar_bytes=args.scalar_bytes,
        index_bytes=args.index_bytes,
        ksp_pc_overhead=args.ksp_pc_overhead,
        num_vectors=args.num_vectors,
        safety_factor=args.safety_factor,
    )

    plan = plan_ranks(
        total_bytes=est["total_bytes"],
        total_rows=est["rows"],
        target_rows_per_rank=args.target_rows_per_rank,
        min_rows_per_rank=args.min_rows_per_rank,
        cores_per_node=args.cores_per_node if args.cores_per_node else None,
        node_mem_gb=args.node_mem if args.node_mem else None,
        max_ranks=args.max_ranks,
    )

    # Report
    print("=== PETSc Resource Plan ===")
    if env_path is not None:
        print(f"Env file                 : {env_path}")
    print(f"Grid (elements)          : Nx={nx}, Ny={ny}, Nz={nz}  (dim={dim})")
    print(f"Fields × dofs/elem/field : {args.fields} × {args.dofs_per_elem_per_field}")
    print(f"Total scalar DOFs        : {est['rows']:,}")
    print(f"Avg nonzeros per row     : {args.nnz_per_row:.1f}")
    print(f"Estimated matrix nnz     : {est['nnz']:,}")
    print("")
    print(f"Matrix bytes (CSR approx): {humanize_bytes(est['mat_bytes'])}")
    print(f"Matrix+PC overhead bytes : {humanize_bytes(est['mat_total_bytes'])}  (x{args.ksp_pc_overhead:.2f})")
    print(f"Vectors total bytes      : {humanize_bytes(est['vec_total_bytes'])}  ({args.num_vectors} vectors)")
    print(f"-------------------------------------------")
    print(f"TOTAL (with safety x{args.safety_factor:.2f}): {humanize_bytes(est['total_bytes'])}")
    print("")
    print("=== Independent Plan (cluster-agnostic) ===")
    print(f"Target rows/rank         : {args.target_rows_per_rank:,}")
    print(f"Planned MPI ranks (ntasks): {plan['ntasks']}")
    print(f"Rows per rank (actual)   : {int(plan['rows_per_rank']):,}")
    print(f"Suggested mem-per-cpu (GB): {plan['mem_per_cpu_gb']:.2f}")
    if plan['warn_over_parallel']:
        print("[WARN] Rows per rank below min threshold; consider fewer ranks to avoid overhead.")

    if plan['nodes'] is not None:
        print("")
        print("=== Node Partition (with cluster constraints) ===")
        print(f"Cores per node           : {args.cores_per_node}")
        print(f"Node memory (GB)         : {args.node_mem:.1f}")
        print(f"Nodes                    : {plan['nodes']}")
        if plan['ntasks_per_node_list'] and len(plan['ntasks_per_node_list']) > 1:
            csv = ",".join(str(x) for x in plan['ntasks_per_node_list'])
            print(f"--ntasks-per-node list   : {csv}")
        else:
            # uniform
            print(f"--ntasks-per-node        : {plan['ntasks_per_node'] or (plan['ntasks_per_node_list'][0] if plan['ntasks_per_node_list'] else args.cores_per_node)}")
        print(f"--mem-per-cpu            : {plan['mem_per_cpu_gb']:.2f}G")

        # Example SLURM line
        print("")
        print("# Example:")
        if plan['ntasks_per_node_list'] and len(plan['ntasks_per_node_list']) > 1:
            csv = ",".join(str(x) for x in plan['ntasks_per_node_list'])
            print("sbatch -J myPETScJob \\")
            print(f"  --nodes={plan['nodes']} --ntasks-per-node={csv} \\")
            print(f"  --mem-per-cpu={plan['mem_per_cpu_gb']:.2f}G --time=1-00:00:00 my_run.sh")
        else:
            npp = plan['ntasks_per_node'] or (plan['ntasks_per_node_list'][0] if plan['ntasks_per_node_list'] else args.cores_per_node)
            print("sbatch -J myPETScJob \\")
            print(f"  --nodes={plan['nodes']} --ntasks-per-node={npp} \\")
            print(f"  --mem-per-cpu={plan['mem_per_cpu_gb']:.2f}G --time=1-00:00:00 my_run.sh")
    else:
        print("")
        print("# Example (cluster-agnostic):")
        print("sbatch -J myPETScJob \\")
        print(f"  --ntasks={plan['ntasks']} --mem-per-cpu={plan['mem_per_cpu_gb']:.2f}G --time=1-00:00:00 my_run.sh")

if __name__ == "__main__":
    main()

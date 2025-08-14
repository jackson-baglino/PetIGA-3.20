#!/usr/bin/env python3
"""
PETSc Resource Planner
----------------------
Estimate (a) total memory, (b) number of MPI ranks, (c) nodes, and
(d) mem-per-cpu for a PETSc job given problem size & sparsity assumptions.

USAGE (examples):
  # Quick: elements only (assume 3 fields, 1 dof/elem/field, nnz/row=60)
  python petsc_resource_planner.py --elements 2_000_000

  # Provide direct total DoFs (recommended when you know it)
  python petsc_resource_planner.py --dofs 6_000_000

  # Customize sparsity and cluster defaults
  python petsc_resource_planner.py --dofs 6e6 --nnz-per-row 60 --node-mem 192 --cores-per-node 40

  # Target a specific mem-per-cpu (GB)
  python petsc_resource_planner.py --dofs 6e6 --target-mem-per-cpu 3.5

NOTES
-----
- By default, we assume 64-bit indices and 64-bit reals (8 bytes each).
- MPIAIJ roughly stores diag/offdiag CSR; here we estimate with a simple model:
    matrix_bytes ≈ nnz * (scalar_bytes + index_bytes) + (rows+1)*index_bytes
  Row-pointer storage is typically small versus nnz terms; we still include it.
- Total memory includes: matrix + (krylov/PC overhead) + multiple vectors.
  Vector budget defaults to 6 vectors (x, b, r, z, p, w); adjust if needed.
- A global safety factor (>1) is applied to cover extra workspace, assembly,
  halo/exchange buffers, monitoring, etc.
- Ranks are chosen so that estimated memory per rank ≤ target mem-per-cpu (if given)
  or ≤ node_mem/cores_per_node otherwise.
"""

import math
import argparse

def humanize_bytes(num_bytes: float) -> str:
    units = ["B","KB","MB","GB","TB","PB"]
    val = float(num_bytes)
    for u in units:
        if val < 1024 or u == units[-1]:
            return f"{val:.2f} {u}"
        val /= 1024.0

def parse_args():
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # Problem size
    ap.add_argument("--elements", type=float, default=None,
                    help="Number of elements (if DOFs unknown). Underscores ok, e.g., 2_000_000.")
    ap.add_argument("--dofs", type=float, default=None,
                    help="Total scalar DOFs (preferred if known). Overrides --elements.")
    ap.add_argument("--fields", type=int, default=3, help="Number of coupled fields (e.g., 3).")
    ap.add_argument("--dofs-per-elem-per-field", type=float, default=1.0,
                    help="Estimate of scalar dofs contributed per element per field (depends on discretization).")
    ap.add_argument("--nnz-per-row", type=float, default=60.0,
                    help="Average nonzeros per matrix row (per scalar dof). Increase for stronger coupling / higher order.")
    # Storage model
    ap.add_argument("--scalar-bytes", type=int, default=8, help="Bytes per scalar (8 for double).")
    ap.add_argument("--index-bytes", type=int, default=8, help="Bytes per matrix index (8 for 64-bit PETSc).")
    ap.add_argument("--ksp-pc-overhead", type=float, default=1.20,
                    help="Multiplier for PC/Krylov & assembly overhead applied to matrix bytes.")
    ap.add_argument("--num-vectors", type=int, default=6, help="How many work vectors to budget (solution, rhs, residual, search dirs, etc.).")
    ap.add_argument("--safety-factor", type=float, default=2.0, help="Global safety factor on total bytes.")
    # Cluster defaults
    ap.add_argument("--node-mem", type=float, default=192.0, help="Node memory in GB.")
    ap.add_argument("--cores-per-node", type=int, default=40, help="MPI ranks per node (if 1 rank per core).")
    ap.add_argument("--target-mem-per-cpu", type=float, default=None,
                    help="If set, target mem-per-cpu (GB) that you want to stay under.")
    # Rank hint
    ap.add_argument("--max-ranks", type=int, default=2000, help="Upper bound for ranks to consider.")
    return ap.parse_args()

def estimate_memory_bytes(total_dofs: int,
                          nnz_per_row: float,
                          scalar_bytes: int,
                          index_bytes: int,
                          ksp_pc_overhead: float,
                          num_vectors: int,
                          safety_factor: float) -> dict:
    rows = total_dofs
    nnz = rows * nnz_per_row

    # Matrix storage (simplified CSR + values; we include row-pointer bytes).
    mat_values = nnz * scalar_bytes
    mat_col_idx = nnz * index_bytes
    mat_row_ptr = (rows + 1) * index_bytes
    mat_bytes = mat_values + mat_col_idx + mat_row_ptr

    # Krylov/PC overhead (ILU/Jacobi/ASM/etc.). This is a coarse multiplier.
    mat_total = mat_bytes * ksp_pc_overhead

    # Vector storage (num_vectors * rows * scalar_bytes)
    vec_total = num_vectors * rows * scalar_bytes

    # Total with safety
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
               node_mem_gb: float,
               cores_per_node: int,
               target_mem_per_cpu_gb: float | None,
               max_ranks: int):
    # Available per-core memory if not specified
    default_mem_per_cpu_gb = node_mem_gb / cores_per_node
    mem_per_cpu_gb = target_mem_per_cpu_gb or default_mem_per_cpu_gb

    # Minimum ranks to satisfy memory
    required_ranks = math.ceil((total_bytes / (1024**3)) / mem_per_cpu_gb)
    required_ranks = max(1, min(required_ranks, max_ranks))

    # Fit into nodes
    nodes = math.ceil(required_ranks / cores_per_node)
    ntasks = nodes * cores_per_node
    mem_per_cpu_gb = max(mem_per_cpu_gb, (total_bytes / (1024**3)) / ntasks)  # ensure enough per-cpu

    return {
        "required_ranks": required_ranks,
        "nodes": nodes,
        "ntasks": ntasks,
        "ntasks_per_node": cores_per_node,
        "mem_per_cpu_gb": mem_per_cpu_gb,
    }

def main():
    args = parse_args()

    # Compute total DOFs
    if args.dofs is not None:
        total_dofs = int(args.dofs)
    elif args.elements is not None:
        total_dofs = int(args.elements * args.fields * args.dofs_per_elem_per_field)
    else:
        raise SystemExit("Provide either --dofs or --elements. See --help.")

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
        node_mem_gb=args.node_mem,
        cores_per_node=args.cores_per_node,
        target_mem_per_cpu_gb=args.target_mem_per_cpu,
        max_ranks=args.max_ranks,
    )

    # Report
    print("=== PETSc Resource Plan ===")
    print(f"Total scalar DOFs         : {est['rows']:,}")
    print(f"Avg nonzeros per row      : {args.nnz_per_row:.1f}")
    print(f"Estimated matrix nnz      : {est['nnz']:,}")
    print("")
    print(f"Matrix bytes (CSR approx) : {humanize_bytes(est['mat_bytes'])}")
    print(f"Matrix+PC overhead bytes  : {humanize_bytes(est['mat_total_bytes'])}  (x{args.ksp_pc_overhead:.2f})")
    print(f"Vectors total bytes       : {humanize_bytes(est['vec_total_bytes'])}  ({args.num_vectors} vectors)")
    print(f"-------------------------------------------")
    print(f"TOTAL (with safety x{args.safety_factor:.2f}): {humanize_bytes(est['total_bytes'])}")
    print("")
    print("=== Cluster Assumptions ===")
    print(f"Node memory (GB)          : {args.node_mem:.1f}")
    print(f"Cores per node            : {args.cores_per_node}")
    print(f"Target mem-per-cpu (GB)   : {args.target_mem_per_cpu if args.target_mem_per_cpu else args.node_mem/args.cores_per_node:.2f}")
    print("")
    print("=== Suggested SLURM Settings ===")
    print(f"--nodes                   : {plan['nodes']}")
    print(f"--ntasks-per-node        : {plan['ntasks_per_node']}")
    print(f"--ntasks                 : {plan['ntasks']}")
    print(f"--mem-per-cpu            : {plan['mem_per_cpu_gb']:.2f}G")
    print("")
    print("# Example:")
    print("sbatch -J myPETScJob \\")
    print(f"  --nodes={plan['nodes']} --ntasks-per-node={plan['ntasks_per_node']} \\")
    print(f"  --mem-per-cpu={plan['mem_per_cpu_gb']:.2f}G --time=1-00:00:00 my_run.sh")
    print("")
    print("# Tips:")
    print("# - If you know a tighter nnz-per-row, reduce it to get a less conservative plan.")
    print("# - Increase num-vectors if you use advanced preconditioners that need extra work vectors.")
    print("# - For 32-bit indices (rare on your setup), set --index-bytes=4 to reduce matrix footprint.")
    print("# - If your cluster has different cores-per-node/memory, set --cores-per-node / --node-mem.")
    print("# - To hard-cap mem-per-cpu, pass --target-mem-per-cpu GB; ranks will scale to satisfy it.")

if __name__ == "__main__":
    main()

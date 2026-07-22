#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# scripts/lib/alloc.sh — single source of truth for job-allocation constants.
#
# Sourced by the Studio and HPC run/submit scripts. Previously each of six
# scripts hard-coded these values and they drifted out of sync twice
# (2026-07-12 and 2026-07-15, when submit_regression.sh and run_batch_tests.sh
# were left at the old 10000 while the rest moved to 40000). Change them HERE
# and nowhere else.
#
# Each is set with := so an environment override wins, e.g.
#   TARGET_DOFS_PER_CORE=80000 ./scripts/HPC/submit_permafrost.sh ...
# lets large-domain runs tune the target without editing any script.
# ---------------------------------------------------------------------------

# Target unknowns (DOF) per MPI rank for implicit solves. PETSc's healthy band
# is ~20k-100k unknowns/rank; below ~20k, reductions and halo exchange
# dominate, and ASM+ILU weakens as the subdomain count grows. These runs are
# step-limited, so wall time is nearly flat in rank count and the allocation
# SIZE is what costs — so we target the upper part of the band. 50k keeps rank
# counts (hence core-hours and queue wait) modest; --half-cores doubles the
# per-rank load to ~100k, at the top of the demonstrated-good range
# (108k/rank ran ~7 s/step at 1.3M DoFs).
: "${TARGET_DOFS_PER_CORE:=50000}"

# MPI ranks per node on the Caltech Resnick cluster. 32 is the safe count
# across the icelake|skylake|cascadelake constraint. MAX_TASKS_PER_NODE is the
# same value under the name the batch/regression planners expect.
: "${NTASKS_PER_NODE:=32}"
: "${MAX_TASKS_PER_NODE:=${NTASKS_PER_NODE}}"

# Cap for local (Studio) runs — physical cores on the dev Mac.
: "${MAX_LOCAL_CORES:=12}"

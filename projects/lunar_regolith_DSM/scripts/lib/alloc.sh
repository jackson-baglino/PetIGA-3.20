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
#   TARGET_DOFS_PER_CORE=80000 ./scripts/HPC/submit_lunar.sh ...
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
#
# Raised 50k -> 200k on 2026-10-01, from the scaling test run in enceladus_DSM
# (same src/; ../enceladus_DSM/studies/keff_sintering/scaling/README.md). The
# phase-field step does not strong-scale: t ~ 12.7 s + 1765/P per step at 24M
# DoF, so core-hours per run are flat to +/-10 % from 150k to 400k DoF/core and
# 3-5x lower than at 60k. 200k is the bottom of that valley.
#
# What this does NOT buy is wall time. Cost is flat because the step slows as
# ranks are removed, so a run takes longer at 200k than it did at 50k. For the
# small contact-angle meshes (382k-495k DoF: 2-3 ranks now, 8-10 before) expect
# roughly 2-4x the wall time of the September batches -- size --time for that,
# or override the target for one submission:
#   TARGET_DOFS_PER_CORE=50000 ./scripts/HPC/submit_batch.sh ...
: "${TARGET_DOFS_PER_CORE:=200000}"

# MPI ranks per node on the Caltech Resnick cluster. 32 is the safe count
# across the icelake|skylake|cascadelake constraint. MAX_TASKS_PER_NODE is the
# same value under the name the batch/regression planners expect.
: "${NTASKS_PER_NODE:=32}"
: "${MAX_TASKS_PER_NODE:=${NTASKS_PER_NODE}}"

# Lower edge of the acceptable tasks-per-node band. Used by plan_alloc below to
# rebalance instead of leaving a nearly-empty last node.
: "${MIN_TASKS_PER_NODE:=28}"

# Cap for local (Studio) runs — physical cores on the dev Mac.
: "${MAX_LOCAL_CORES:=12}"

# ---------------------------------------------------------------------------
# mem_per_cpu <total_dofs> <nprocs>  ->  echoes the --mem-per-cpu to request, e.g. "2G"
#
# The flat 1G in run_lunar.sh's #SBATCH header was sized for ~50k DoF per core.
# Peak RSS per rank measured in the same scaling test (sacct MaxRSS, always on
# rank 0): 630 MB at 100k DoF/core, 860 MB at 198k, and OOM-killed at 394k
# under 1G. Fit, to ~10 %:
#     peak = 0.20 GB + 2.35 GB per 1M DoF on the rank + 8 B x total DoF
# Requested = 1.5 x peak, rounded UP to whole GB, never below 1G. Memory is not
# billed on Resnick (per-core-hour only), so the margin is free.
# MEM_PER_CPU=<N>G overrides.
# ---------------------------------------------------------------------------
mem_per_cpu() {
    local total="$1" np="$2"
    if [[ -n "${MEM_PER_CPU:-}" ]]; then echo "$MEM_PER_CPU"; return; fi
    awk -v N="$total" -v P="$np" 'BEGIN{
        D = N / (P > 0 ? P : 1)
        peak = 0.20 + 2.35e-6 * D + 8e-9 * N
        g = int(1.5 * peak); if (g < 1.5 * peak) g++
        if (g < 1) g = 1
        printf "%dG\n", g }'
}

# ---------------------------------------------------------------------------
# plan_alloc <total_ranks>  ->  echoes "<nodes> <tasks_per_node>"
#
# The naive plan (nodes = ceil(ranks/NTASKS_PER_NODE), tasks-per-node fixed at
# the max) wastes a whole node whenever ranks is just over a multiple: 33 ranks
# becomes 32+1, and the second node sits ~97% idle while still being billed.
# This spreads the ranks evenly instead — 33 ranks becomes 2x17 — keeping
# tasks-per-node inside [MIN_TASKS_PER_NODE, NTASKS_PER_NODE] when it can.
#
# Ported from the band-clamping planner in the legacy dry_snow_metamorphism
# batch submitter, which is the only place this logic existed.
# ---------------------------------------------------------------------------
plan_alloc() {
    local ranks="$1"
    (( ranks < 1 )) && ranks=1

    local nodes tpn
    nodes=$(( (ranks + NTASKS_PER_NODE - 1) / NTASKS_PER_NODE ))
    (( nodes < 1 )) && nodes=1
    tpn=$(( (ranks + nodes - 1) / nodes ))      # even spread over those nodes

    # If the even spread drops below the band, drop nodes until it climbs back
    # in (a single node holding the whole job is always acceptable).
    if (( tpn < MIN_TASKS_PER_NODE )); then
        local n p
        for (( n = nodes; n >= 1; n-- )); do
            p=$(( (ranks + n - 1) / n ))
            if (( p >= MIN_TASKS_PER_NODE && p <= NTASKS_PER_NODE )); then
                nodes=$n; tpn=$p; break
            fi
        done
    fi

    (( tpn > NTASKS_PER_NODE )) && tpn=$NTASKS_PER_NODE
    (( tpn < 1 )) && tpn=1
    echo "$nodes $tpn"
}

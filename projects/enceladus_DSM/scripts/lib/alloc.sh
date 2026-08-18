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
#   TARGET_DOFS_PER_CORE=80000 ./scripts/HPC/submit_enceladus.sh ...
# lets large-domain runs tune the target without editing any script.
# ---------------------------------------------------------------------------

# Target unknowns (DOF) per MPI rank for implicit solves. PETSc's healthy band
# is ~20k-100k unknowns/rank; below ~20k, reductions and halo exchange
# dominate, and ASM+ILU weakens as the subdomain count grows. These runs are
# step-limited, so wall time is nearly flat in rank count and the allocation
# SIZE is what costs — so we target the upper part of the band.
#
# Raised 50k -> 80k on 2026-07-31: the k_eff study's 2mm packings at safety=0.5
# are 24M DoF, which at 50k/rank demanded 480 ranks and a correspondingly long
# queue wait. 80k puts the same run at 300 ranks while staying inside the band
# (--half-cores would take it to ~160k, above the demonstrated-good 108k/rank
# that ran ~7 s/step at 1.3M DoFs, so do not combine the two on large jobs).
#
# Raised 80k -> 100k on 2026-08-18 for the Molaro 2019 campaign: 100k is the top
# of PETSc's healthy band and these runs are step-limited, so shrinking the
# allocation costs almost no wall time and cuts both the core-hour bill and the
# queue wait. 100k is still BELOW the 108k/rank that demonstrably ran ~7 s/step,
# so it is inside measured territory rather than extrapolated.
#
# CAUTION: per-rank memory scales with this, and 100k leaves less headroom than
# 80k did. Runs using -keff carry a SECOND large operator (the scalar corrector
# matrix plus its GAMG hierarchy) on top of the phase-field Jacobian and its
# ASM/ILU(3) factor. Check peak RSS on a small run before submitting a large one
# at this target. Do NOT combine with --half-cores on a large job: that would put
# it at ~200k/rank, well outside the band.
: "${TARGET_DOFS_PER_CORE:=100000}"

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

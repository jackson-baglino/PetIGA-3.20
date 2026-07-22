#!/usr/bin/env bash
# =============================================================================
# hpc_cost.sh — Caltech Resnick HPC cost utilities
#
# Source this file; do not execute directly.
#
# Rates: https://www.hpc.caltech.edu/rates
# Fiscal year resets October 1.
#
# Tiered rates (aggregate FY spend):
#   Tier 1: $0.012/core-hour  (spend ≤ $6,400)
#   Tier 2: $0.007/core-hour  (spend $6,401–$24,000)
#   Tier 3: $0.006/core-hour  (spend > $24,000)
#
# 1 CPU core-hour = 1 compute unit.
# =============================================================================

HPC_ACCOUNT="${HPC_ACCOUNT:-rubyfu}"
HPC_USER="${HPC_USER:-jbaglino}"

# ---------------------------------------------------------------------------
# Internal: Caltech fiscal year starts October 1.
# ---------------------------------------------------------------------------
_hpc_fy_start() {
    local m y
    m=$(date +%m | sed 's/^0//') # strip leading zero for arithmetic
    y=$(date +%Y)
    if (( m >= 10 )); then
        printf '%s-10-01T00:00:00\n' "$y"
    else
        printf '%s-10-01T00:00:00\n' "$(( y - 1 ))"
    fi
}

# ---------------------------------------------------------------------------
# Internal: compute_cost <core_hours>
#   Total spend in USD for <core_hours>, accounting for tier transitions.
# ---------------------------------------------------------------------------
_hpc_compute_cost() {
    awk -v ch="$1" 'BEGIN {
        t1 = 6400.0  / 0.012        # 533333.33 ch — end of tier 1
        t2 = t1 + 17600.0 / 0.007  # 3047619.04 ch — end of tier 2
        if      (ch <= t1) c = ch * 0.012
        else if (ch <= t2) c = t1*0.012 + (ch-t1)*0.007
        else               c = t1*0.012 + (t2-t1)*0.007 + (ch-t2)*0.006
        printf "%.2f", c
    }'
}

# ---------------------------------------------------------------------------
# Internal: marginal_rate <core_hours_fy>
#   Returns the per-core-hour rate string at a given FY core-hour total.
# ---------------------------------------------------------------------------
_hpc_marginal_rate() {
    awk -v ch="$1" 'BEGIN {
        t1 = 6400.0  / 0.012
        t2 = t1 + 17600.0 / 0.007
        if      (ch < t1) print "0.0120"
        else if (ch < t2) print "0.0070"
        else              print "0.0060"
    }'
}

# ---------------------------------------------------------------------------
# Internal: format_duration <seconds>
# ---------------------------------------------------------------------------
_hpc_fmt_dur() {
    local s="$1"
    printf '%dh %02dm %02ds' "$(( s/3600 ))" "$(( (s%3600)/60 ))" "$(( s%60 ))"
}

# ---------------------------------------------------------------------------
# fy_core_hours [<exclude_job_id>]
#   Query sacct for total CPU core-hours charged to $HPC_ACCOUNT since the
#   start of the current fiscal year.  Pass a SLURM job ID to exclude it
#   (avoids double-counting the currently-running job if it appears in sacct).
#   Returns "0" if sacct is unavailable.
# ---------------------------------------------------------------------------
fy_core_hours() {
    local exclude="${1:-}"
    if ! command -v sacct &>/dev/null; then
        echo "0"
        return
    fi
    local fy_start
    fy_start=$(_hpc_fy_start)
    # CPUTimeRAW = ncpus * elapsed_seconds; divide by 3600 for core-hours.
    sacct -A "$HPC_ACCOUNT" \
          --starttime="$fy_start" \
          -X --format=JobID,CPUTimeRAW --noheader --parsable2 \
          --state=COMPLETED,FAILED,CANCELLED,TIMEOUT,NODE_FAIL \
          2>/dev/null \
        | awk -F'|' -v excl="$exclude" '
            BEGIN { s = 0 }
            { if (excl != "" && $1 == excl) next
              if ($2 + 0 > 0) s += $2 }
            END { printf "%.0f", s / 3600 }'
}

# ---------------------------------------------------------------------------
# hpc_cost_pre_submit <nprocs>
#   Print estimated cost-per-hour before job submission (login node).
# ---------------------------------------------------------------------------
hpc_cost_pre_submit() {
    local nprocs="$1"
    local fy_ch rate cph fy_cost

    fy_ch=$( fy_core_hours )
    rate=$( _hpc_marginal_rate "$fy_ch" )
    cph=$( awk "BEGIN{printf \"%.2f\", $nprocs * $rate}" )
    fy_cost=$( _hpc_compute_cost "$fy_ch" )

    echo "------------------------------------------------------------"
    echo "  HPC Cost Estimate  (account: ${HPC_ACCOUNT} / user: ${HPC_USER})"
    printf "  %-26s %.0f core-hours (\$%s)\n" \
           "FY usage to date:"   "$fy_ch"   "$fy_cost"
    printf "  %-26s \$%s / core-hour\n" \
           "Marginal rate:"      "$rate"
    printf "  %-26s \$%s / hour of wall time\n" \
           "Rate for ${nprocs} cores:"  "$cph"
    echo "------------------------------------------------------------"
}

# ---------------------------------------------------------------------------
# hpc_cost_post_job <nprocs> <wall_secs> [<slurm_job_id>]
#   Print job cost and FY total at the end of a SLURM job.
#   <slurm_job_id> is excluded from the sacct FY query to avoid double-counting
#   (the job may or may not have appeared as COMPLETED yet).
# ---------------------------------------------------------------------------
hpc_cost_post_job() {
    local nprocs="$1"
    local wall_secs="$2"
    local jid="${3:-}"

    # Core-hours consumed by this job
    local job_ch
    job_ch=$( awk "BEGIN{printf \"%.2f\", $nprocs * $wall_secs / 3600}" )

    # FY usage from all *other* completed jobs (exclude this one)
    local fy_ch_prev
    fy_ch_prev=$( fy_core_hours "$jid" )

    local fy_ch_total fy_cost_prev fy_cost_total job_cost rate_now
    fy_ch_total=$( awk "BEGIN{printf \"%.2f\", $fy_ch_prev + $job_ch}" )
    fy_cost_prev=$( _hpc_compute_cost "$fy_ch_prev" )
    fy_cost_total=$( _hpc_compute_cost "$fy_ch_total" )
    job_cost=$( awk "BEGIN{printf \"%.2f\", $fy_cost_total - $fy_cost_prev}" )
    rate_now=$( _hpc_marginal_rate "$fy_ch_total" )

    echo ""
    echo "========================================================================="
    echo "  Resnick HPC Cost Summary"
    echo "  Account: ${HPC_ACCOUNT}   User: ${HPC_USER}"
    echo "  Rates: https://www.hpc.caltech.edu/rates  |  FY resets October 1"
    echo "  -----------------------------------------------------------------------"
    printf "  %-24s %s\n"          "Wall time:"            "$( _hpc_fmt_dur "$wall_secs" )"
    printf "  %-24s %d\n"          "MPI ranks:"            "$nprocs"
    printf "  %-24s %.2f ch\n"     "Job core-hours:"       "$job_ch"
    printf "  %-24s \$%s\n"        "Job cost:"             "$job_cost"
    echo "  -----------------------------------------------------------------------"
    printf "  %-24s %.0f ch\n"     "FY core-hours (total):" "$fy_ch_total"
    printf "  %-24s \$%s\n"        "FY total cost:"        "$fy_cost_total"
    printf "  %-24s \$%s/core-hour\n" "Marginal rate now:"  "$rate_now"
    echo "========================================================================="
}

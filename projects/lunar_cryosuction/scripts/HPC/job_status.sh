#!/usr/bin/env bash
# job_status.sh: one line per SLURM job log, with how its latest leg ended and the last simulated
# time. Run from the project root on the HPC (squeue marks running jobs) or locally on a
# downloaded output_files/ folder.
#
#   ./scripts/HPC/job_status.sh [dir] [name pattern, e.g. '7S*' or '*t7s_*']
#
# dir defaults to output_files/ and is searched recursively, so it can be one batch's log folder
# (output_files/<batch>/) or the whole tree. Two layouts are read:
#   per-run logs  <dir>/<batch>/<run ID>.out|.err   (submit_batch.sh since 2026-11-20): every leg
#                 of a run appends to the same files after a "===== LEG job <id> ..." line; the
#                 job ID and the state come from the last leg only.
#   per-job logs  <dir>/<geom>__<exp>.o<id>|.e<id>   (older jobs): one file pair per job.
#
# States (first match wins):
#   RUNNING      the job ID is in squeue
#   TIME_LIMIT   SLURM stopped it at the --time limit (resume with resume_run.sh)
#   STALLED      the solver's stall guard stopped it: the field stopped changing (usually all ice gone)
#   CANCELLED    cancelled by the user (scancel)
#   FINISHED_OK  reached t_final ("Simulation completed successfully")
#   CRASHED      ended with an error (PETSc error, signal, node fault)
#   NO_END_YET   no end marker: still running when the logs were copied, or the log was cut off
d=${1:-output_files}; pat=${2:-*}
q=$(command -v squeue >/dev/null && squeue -h -u "$USER" -o %i)

from_line() {   # line number of the last "===== LEG job" header in a file (1 if none)
  local n; n=$(grep -n '^===== LEG job ' "$1" 2>/dev/null | tail -1 | cut -d: -f1); echo "${n:-1}"
}

status() {     # status <job id> <stdout file> <stderr file> <first line of the leg in each>
  local id=$1 o=$2 e=$3 lo=${4:-1} le=${5:-1}
  local O E; O=$(mktemp); E=$(mktemp)
  tail -n +"$lo" "$o" > "$O"; { [ -f "$e" ] && tail -n +"$le" "$e"; } > "$E"
  if [ -n "$id" ] && echo "$q" | grep -qx "$id"; then echo RUNNING
  elif grep -q 'DUE TO TIME LIMIT' "$E"; then echo TIME_LIMIT
  elif grep -q 'has not changed for' "$E"; then echo STALLED
  elif grep 'JOB [0-9]* ON .* CANCELLED' "$E" | grep -vq 'DUE'; then echo CANCELLED
  elif grep -q 'Simulation completed successfully' "$O"; then echo FINISHED_OK
  elif grep -q 'simulation errors\|PETSC ERROR\|Terminate' "$O" "$E"; then echo CRASHED
  else echo NO_END_YET; fi
  rm -f "$O" "$E"
}

simtime() {    # last simulated time in days from the diagnostics table
  grep -E '^ *[0-9]+ \| +[0-9.e+-]+ \| [0-9.e+-]+ \|' "$1" | tail -1 | awk -F'|' '{printf "%.1f", $2/86400}'
}

{
  # per-run logs
  find "$d" -type f -name "$pat.out" 2>/dev/null | while read -r o; do
    e=${o%.out}.err; name=$(basename "${o%.out}")
    id=$(grep -o '^===== LEG job [0-9]*' "$o" | tail -1 | awk '{print $4}')
    legs=$(grep -c '^===== LEG job ' "$o")
    st=$(status "$id" "$o" "$e" "$(from_line "$o")" "$(from_line "$e")")
    printf '%-12s %8s  t=%6s d  leg %s  %s/%s\n' "$st" "${id:-?}" "$(simtime "$o")" "${legs:-?}" \
      "$(basename "$(dirname "$o")")" "$name"
  done
  # per-job logs (older layout)
  find "$d" -type f -name "$pat.o[0-9]*" 2>/dev/null | while read -r o; do
    id=${o##*.o}; e=${o%.o*}.e$id; name=$(basename "${o%.o*}")
    t=$(simtime "$o")
    printf '%-12s %8s  t=%6s d  %s\n' "$(status "$id" "$o" "$e")" "$id" "${t:-?}" "$name"
  done
} | sort

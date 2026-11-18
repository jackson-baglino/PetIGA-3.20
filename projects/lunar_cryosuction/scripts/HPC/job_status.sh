#!/usr/bin/env bash
# job_status.sh: one line per SLURM job log in output_files/, with how it ended and the last
# simulated time. Run from the project root on the HPC (squeue marks running jobs) or locally on a
# downloaded output_files/ folder.
#
#   ./scripts/HPC/job_status.sh [output_files dir] [name pattern, e.g. '*t7s_*']
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
for o in "$d"/$pat.o[0-9]*; do
  [ -e "$o" ] || continue
  id=${o##*.o}; e=${o%.o*}.e$id; name=$(basename "${o%.o*}")
  t=$(grep -E '^ *[0-9]+ \| +[0-9.e+-]+ \| [0-9.e+-]+ \|' "$o" | tail -1 | awk -F'|' '{printf "%.1f", $2/86400}')
  if echo "$q" | grep -qx "$id"; then s=RUNNING
  elif grep -q 'DUE TO TIME LIMIT' "$e" 2>/dev/null; then s=TIME_LIMIT
  elif grep -q 'has not changed for' "$e" 2>/dev/null; then s=STALLED
  elif grep 'JOB [0-9]* ON .* CANCELLED' "$e" 2>/dev/null | grep -vq 'DUE'; then s=CANCELLED
  elif grep -q 'Simulation completed successfully' "$o"; then s=FINISHED_OK
  elif grep -q 'simulation errors\|PETSC ERROR\|Terminate' "$o" "$e" 2>/dev/null; then s=CRASHED
  else s=NO_END_YET; fi
  printf '%-12s %8s  t=%6s d  %s\n' "$s" "$id" "${t:-?}" "$name"
done | sort

#!/usr/bin/env bash
# =============================================================================
# run_batch_measure.sh — measure a two-grain sintering batch ON THE CLUSTER and
# emit one summary.csv, so only the CSVs need downloading.
#
# WHY THIS IS SEPARATE FROM run_batch_postprocess.sh
# --------------------------------------------------
# run_batch_postprocess.sh is the generic figure sweep and runs LOCALLY on a
# downloaded batch. It does not run neck_width.py or grain_shrinkage.py, and
# for these runs the snapshots are far too large to move: a dom3 arm is 6.9M
# elements x 80 log-spaced snapshots. So: convert and measure in place, bring
# back kilobytes.
#
#   1. HPC:   ./scripts/HPC/submit_batch.sh --tag molaro2019_humidity \
#                 --tests-file studies/molaro_2019/batches/humidity_fit.txt
#   2. HPC:   bash postprocess/run_batch_measure.sh $SCRATCH/enceladus_DSM/batch_<...>
#   3. Read   <batch>/summary.csv, pick the 2-3 arms worth looking at
#   4. Local: rsync only those run folders down
#
# WHAT summary.csv HOLDS (one row per arm)
#   run, geometry, experiment, alpha_c, humidity, Lz_m, Lr_m, eps_m,
#   neck_w_final_um, neck_w_at_78min_um, t_star_s,
#   R_large_init_um, R_large_end_um, dR_large_pct, dR_small_pct, n_snapshots
#
#   t_star_s is the CLOCK SHIFT: the time at which the model neck first reaches
#   Molaro's first measured width (32.81 um). Our run starts from tangency and
#   theirs starts at an unknown time after contact, so comparing raw clocks is
#   meaningless -- everything is anchored on t_star.
#
# Usage:
#   bash run_batch_measure.sh /path/to/run        # ONE run (the usual case)
#   bash run_batch_measure.sh /path/to/batch      # a batch parent, fans out
#   bash run_batch_measure.sh                     # from inside either
#   ANCHOR_NECK_UM=32.81 bash run_batch_measure.sh /path/to/run
# =============================================================================
set -uo pipefail

if [[ $# -ge 1 ]]; then
    BATCH_DIR="$(cd "$1" && pwd)"
else
    BATCH_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi

if   [[ -d "$BATCH_DIR/postprocess" ]]; then
    POSTPROCESS="$BATCH_DIR/postprocess"
elif [[ -d "$(dirname "$BATCH_DIR")/postprocess" ]]; then
    POSTPROCESS="$(dirname "$BATCH_DIR")/postprocess"
else
    echo "Could not find postprocess/ next to or inside $BATCH_DIR" >&2
    exit 1
fi

PYTHON="$(command -v python3 || command -v python || true)"
[[ -n "$PYTHON" ]] || { echo "no python on PATH" >&2; exit 1; }

ANCHOR_NECK_UM="${ANCHOR_NECK_UM:-32.81}"   # Molaro's first measured width
SUMMARY="$BATCH_DIR/summary.csv"

echo "============================================================"
echo "  Batch measurement"
echo "  Batch dir   : $BATCH_DIR"
echo "  postprocess : $POSTPROCESS"
echo "  anchor neck : ${ANCHOR_NECK_UM} um"
echo "============================================================"

n_ok=0; n_fail=0; n_skip=0

# Accept a SINGLE run directory as well as a batch parent. The Molaro campaign
# runs one arm at a time (studies/molaro_2019/RUNBOOK.md), so the common case is
# pointing this at one run; a batch parent still fans out as before.
if [[ -f "$BATCH_DIR/igasol.dat" ]]; then
    RUNS=("$BATCH_DIR")
    SUMMARY="$BATCH_DIR/summary.csv"
    echo "  (single run directory)"
else
    RUNS=("$BATCH_DIR"/*/)
fi

printf 'run,geometry,experiment,alpha_c,humidity,Lz_m,Lr_m,eps_m,' > "$SUMMARY"
printf 'neck_w_final_um,neck_w_at_78min_um,t_star_s,' >> "$SUMMARY"
printf 'R_large_init_um,R_large_end_um,dR_large_pct,dR_small_pct,n_snapshots\n' >> "$SUMMARY"

for run in "${RUNS[@]}"; do
    run="${run%/}"; name="$(basename "$run")"
    case "$name" in inputs_snapshot|src_snapshot|postprocess) continue ;; esac
    [[ -f "$run/igasol.dat" ]] || { echo "  skip $name (no igasol.dat)"; ((n_skip++)); continue; }

    echo ""
    echo "------------------------------------------------------------"
    echo "> $name"

    axisym=$(awk '$1=="-axisym"{print $2; exit}' "$run"/*.opts 2>/dev/null | head -n1)
    ax_flag=""; [[ "${axisym:-0}" == "1" ]] && ax_flag="--axisym"

    "$PYTHON" "$POSTPROCESS/plot_fields.py"     --dir "$run"            2>&1 | sed 's/^/    /'
    "$PYTHON" "$POSTPROCESS/neck_width.py"      "$run" $ax_flag         2>&1 | sed 's/^/    /'
    "$PYTHON" "$POSTPROCESS/grain_shrinkage.py" "$run" --no-plot        2>&1 | sed 's/^/    /'

    if [[ ! -f "$run/neck_width.csv" ]]; then
        echo "    no neck_width.csv — arm not summarised"; ((n_fail++)); continue
    fi

    # Reduce this arm to one row. Kept in python so the anchor interpolation
    # and the opts parsing are not reimplemented in awk.
    ANCHOR_NECK_UM="$ANCHOR_NECK_UM" "$PYTHON" - "$run" "$SUMMARY" <<'PY' 2>&1 | sed 's/^/    /'
import csv, os, sys, glob
from pathlib import Path
run, summary = Path(sys.argv[1]), Path(sys.argv[2])
anchor = float(os.environ["ANCHOR_NECK_UM"]) * 1e-6

def opt(key, default=""):
    for f in sorted(run.glob("*.opts")):
        for line in f.read_text().splitlines():
            p = line.split("#", 1)[0].split()
            if len(p) >= 2 and p[0] == key:
                return p[1]
    return default

neck = list(csv.DictReader((run / "neck_width.csv").open()))
t  = [float(r["t_s"]) for r in neck]
w  = [float(r["neck_width_m"]) for r in neck]

def interp_time_at(target):
    """First time the neck reaches `target` (linear between samples)."""
    for (t0, w0), (t1, w1) in zip(zip(t, w), zip(t[1:], w[1:])):
        if w0 <= target <= w1 and w1 > w0:
            return t0 + (target - w0) * (t1 - t0) / (w1 - w0)
    return float("nan")

def interp_width_at(tt):
    for (t0, w0), (t1, w1) in zip(zip(t, w), zip(t[1:], w[1:])):
        if t0 <= tt <= t1 and t1 > t0:
            return w0 + (tt - t0) * (w1 - w0) / (t1 - t0)
    return float("nan")

t_star = interp_time_at(anchor)
# Their 78-min span, measured from OUR clock zero-point t_star.
w_78 = interp_width_at(t_star + 78 * 60.0) if t_star == t_star else float("nan")

gs = run / "grain_shrinkage.csv"
dR_lg = dR_sm = R_lg0 = R_lg1 = float("nan")
if gs.is_file():
    g = list(csv.DictReader(gs.open()))
    if len(g) >= 2:
        R_lg0, R_lg1 = float(g[0]["R_large_m"]), float(g[-1]["R_large_m"])
        R_sm0, R_sm1 = float(g[0]["R_small_m"]), float(g[-1]["R_small_m"])
        dR_lg = 100.0 * (R_lg1 / R_lg0 - 1.0) if R_lg0 else float("nan")
        dR_sm = 100.0 * (R_sm1 / R_sm0 - 1.0) if R_sm0 else float("nan")

# Identify the staged opts files by CONTENT, not by name: the geometry file
# is the one that sets the mesh, the experiment file the one that sets the
# run length. Name matching breaks the moment a file is renamed.
def which(*keys):
    for f in sorted(run.glob("*.opts")):
        heads = {ln.split("#", 1)[0].split()[0]
                 for ln in f.read_text().splitlines() if ln.split("#", 1)[0].split()}
        if heads & set(keys):
            return f.stem
    return ""

geom = which("-ice_grain_cx", "-Nx")
exp_ = which("-t_final")

row = [run.name, geom, exp_, opt("-alpha_c0"), opt("-humidity"),
       opt("-Lx"), opt("-Ly"), opt("-eps"),
       f"{w[-1]*1e6:.4f}", f"{w_78*1e6:.4f}", f"{t_star:.2f}",
       f"{R_lg0*1e6:.4f}", f"{R_lg1*1e6:.4f}", f"{dR_lg:.4f}", f"{dR_sm:.4f}", len(neck)]
with summary.open("a", newline="") as fh:
    csv.writer(fh, lineterminator="\n").writerow(row)
print(f"t* = {t_star:.0f} s   neck at t*+78min = {w_78*1e6:.2f} um"
      f"   dR_large = {dR_lg:+.2f} %   dR_small = {dR_sm:+.2f} %")
PY
    ((n_ok++))
done

echo ""
echo "============================================================"
echo "  measured $n_ok, failed $n_fail, skipped $n_skip"
echo "  summary: $SUMMARY"
echo ""
echo "  Molaro -20 C targets:  neck 32.81 -> 64.78 um over 78 min"
echo "                         large grain -2.93 % (fit; the caption says -3 %)"
echo "                         small grain: inside its own error bar, do not fit"
echo "============================================================"
column -s, -t "$SUMMARY" 2>/dev/null || cat "$SUMMARY"

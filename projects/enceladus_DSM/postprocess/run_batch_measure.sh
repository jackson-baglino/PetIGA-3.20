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
#   R_large_at_tstar_um, R_large_at_78min_um, dR_large_pct, dR_small_pct,
#   dR_large_fullrun_pct, dR_small_fullrun_pct, n_snapshots
#
#   dR_large_pct / dR_small_pct are over [t*, t*+78 min] -- the SAME window as
#   the neck, and the one Molaro's -2.93 % refers to. The *_fullrun_pct columns
#   are the whole run and are NOT comparable to their targets; they are carried
#   only so the difference is visible. Fitting humidity against the full-run
#   number lands it 1.54x too saturated on a 120-min run (see the note by the
#   computation below).
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
printf 'R_large_at_tstar_um,R_large_at_78min_um,dR_large_pct,dR_small_pct,' >> "$SUMMARY"
printf 'dR_large_fullrun_pct,dR_small_fullrun_pct,n_snapshots\n' >> "$SUMMARY"

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
dR_lg_full = dR_sm_full = float("nan")
if gs.is_file():
    g = list(csv.DictReader(gs.open()))
    if len(g) >= 2:
        gt  = [float(r["t_s"]) for r in g]
        gRl = [float(r["R_large_m"]) for r in g]
        gRs = [float(r["R_small_m"]) for r in g]

        def lerp(xs, ys, x):
            if x <= xs[0]:  return ys[0]
            if x >= xs[-1]: return ys[-1]
            for (x0, y0), (x1, y1) in zip(zip(xs, ys), zip(xs[1:], ys[1:])):
                if x0 <= x <= x1 and x1 > x0:
                    return y0 + (x - x0) * (y1 - y0) / (x1 - x0)
            return ys[-1]

        # THE SHRINKAGE MUST BE READ OVER THE SAME 78-MIN WINDOW AS THE NECK.
        # Molaro's -2.93 % is their large grain's least-squares slope over 78
        # minutes. Reading it over the whole run instead compares a 120-min
        # model change against a 78-min target, and because R_large(t) is
        # linear to R^2 = 1.0000 at fixed humidity that is a clean 120/78 =
        # 1.54x overstatement -- which propagates straight into the Newton
        # step and lands the fitted humidity 1.54x too saturated. Measured on
        # the 2026-09-03 batch: arm 1 reads -2.92 % full-run (looks like a
        # bullseye) against -1.90 % over the anchored window (35 % short).
        if t_star == t_star:
            t_end = t_star + 78 * 60.0
            R_lg0, R_lg1 = lerp(gt, gRl, t_star), lerp(gt, gRl, t_end)
            R_sm0, R_sm1 = lerp(gt, gRs, t_star), lerp(gt, gRs, t_end)
            dR_lg = 100.0 * (R_lg1 / R_lg0 - 1.0) if R_lg0 else float("nan")
            dR_sm = 100.0 * (R_sm1 / R_sm0 - 1.0) if R_sm0 else float("nan")
            if gt[-1] < t_end:
                print(f"    ! run ends {(t_end-gt[-1])/60:.1f} min before t*+78 min;"
                      f" dR_* are EXTRAPOLATED")
        # Kept beside it so the two windows can never be mistaken for each other.
        dR_lg_full = 100.0 * (gRl[-1] / gRl[0] - 1.0) if gRl[0] else float("nan")
        dR_sm_full = 100.0 * (gRs[-1] / gRs[0] - 1.0) if gRs[0] else float("nan")

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
       f"{R_lg0*1e6:.4f}", f"{R_lg1*1e6:.4f}", f"{dR_lg:.4f}", f"{dR_sm:.4f}",
       f"{dR_lg_full:.4f}", f"{dR_sm_full:.4f}", len(neck)]
with summary.open("a", newline="") as fh:
    csv.writer(fh, lineterminator="\n").writerow(row)
print(f"t* = {t_star:.0f} s   neck at t*+78min = {w_78*1e6:.2f} um"
      f"   dR_large = {dR_lg:+.2f} %   dR_small = {dR_sm:+.2f} %"
      f"   (full run: {dR_lg_full:+.2f} % / {dR_sm_full:+.2f} %)")
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

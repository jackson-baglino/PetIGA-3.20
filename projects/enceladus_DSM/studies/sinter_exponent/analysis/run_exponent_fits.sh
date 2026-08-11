#!/usr/bin/env bash
# =============================================================================
# run_exponent_fits.sh — regenerate every fit table and figure in this folder.
#
# This script IS the reproducible record of the numbers quoted in
# studies/sinter_exponent/README.md. Run it from the project root:
#
#     bash studies/sinter_exponent/analysis/run_exponent_fits.sh
#
# Outputs land under results/<topic>/. Stage 1 (the committed Molaro comparison) needs no simulation — it reads the
# neck CSVs already under $RESULTS_BASE/_neck_csv and the validation data in
# inputs/validation. Stage 2 is skipped unless the pilot runs exist.
# =============================================================================
set -euo pipefail


HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
cd "$ROOT"

PY="${PY:-$ROOT/venv_enceladus/bin/python}"
RESULTS_BASE="${RESULTS_BASE:-$HOME/SimulationResults/enceladus_DSM/scratch}"
STUDY="$(cd "$HERE/.." && pwd)"
RES="$STUDY/results"

mkdir -p "$RES/molaro_prenecked" "$RES/demmenie"

FIT="$PY postprocess/fit_neck_growth.py"
DATA20="inputs/validation/molaro2019_fig11_T-20.csv"

# --- 1. Molaro: three committed eps arms vs the -20 C data -------------------
# Expect two "no fittable window" notes: at eps = 0.60 and 0.86 um the
# resolution floor sqrt(12*eps/R) sits above the entire model curve. That is
# the finding, not a failure.
echo "== 1/3  Molaro eps arms vs data =="
$FIT "$RESULTS_BASE/_neck_csv/"*.csv "$DATA20" \
    --results-base "$RESULTS_BASE" \
    --labels "model eps0.86um (loose),model eps0.60um (mid),model eps0.35um (strict),Molaro 2019 data (-20C)" \
    --out "$RES/molaro_prenecked" --prefix molaro --demmenie

# --- 2. The same model arm, sampled where the experiment sampled -------------
# The protocol control. Nothing about the model changes here; only the times
# at which it is read. If this and the run above disagree, that difference is
# the sampling, not the physics.
echo "== 2/3  strict arm resampled at the data's times =="
$FIT "$RESULTS_BASE/_neck_csv/"*epsstrict*.csv "$DATA20" \
    --results-base "$RESULTS_BASE" --resample-at "$DATA20" \
    --labels "model eps0.35um (strict),Molaro 2019 data (-20C)" \
    --out "$RES/molaro_prenecked" --prefix molaro_resampled --demmenie

# --- 3. Demmenie replication, once the runs exist ----------------------------
echo "== 3/3  Demmenie 1 mm spheres =="
demmenie_csvs=()
for geom in sinter_2D_L2200um_eps0.64um_axisym_D1mm_tangent \
            sinter_2D_L2200um_eps0.32um_axisym_D1mm_tangent; do
    # Newest run of each arm; each arm may have several timestamped folders.
    newest="$(ls -1dt "$RESULTS_BASE/$geom"/*/ 2>/dev/null | head -1 || true)"
    [ -n "$newest" ] || continue
    [ -f "$newest/neck_width.csv" ] || \
        $PY postprocess/neck_width.py "$newest" --axisym || continue
    demmenie_csvs+=("$newest/neck_width.csv")
done

if [ "${#demmenie_csvs[@]}" -eq 0 ]; then
    echo "   (skipped — no sinter runs under $RESULTS_BASE yet;"
    echo "    submit studies/sinter_exponent/batches/gate.txt first)"
else
    $FIT "${demmenie_csvs[@]}" --out "$RES/demmenie" --prefix demmenie --demmenie
fi

echo
echo "outputs -> $RES"

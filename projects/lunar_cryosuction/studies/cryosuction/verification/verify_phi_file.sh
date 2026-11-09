#!/usr/bin/env bash
# =============================================================================
# studies/cryosuction/verification/verify_phi_file.sh
#
# Gates for -ic_type phi_file (ice field per control point, written by
# preprocess/build_ice_layout.py). Serial, IC only; calls the binary directly
# like the other gates. The multi-rank read is checked by a normal run through
# the run script (compare its SSA_evo.dat row 0 with P1 here).
#
#   P1  Area     -- a small rough pore with a plug, a cap and a grain: the
#                   solver's initial ice area (SSA_evo.dat row 0, column 2)
#                   matches the builder's sharp area within 0.5 %.
#   P2  Fields   -- phi spans [0, 1] as loaded (printed by the solver).
#   P3  Guards   -- a file for a different mesh, a truncated file, a missing
#                   file and a missing -phi_file all exit non-zero.
#
# Writes verify_phi_file.csv and .log next to this script.
# =============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
EXEC="$ROOT/lunar_cryosuction"
PY="$ROOT/venv_cryosuction/bin/python"
CSV="$HERE/verify_phi_file.csv"; LOG="$HERE/verify_phi_file.log"
WORK="${TMPDIR:-/tmp}/verify_phi_file.$$"; mkdir -p "$WORK"; trap 'rm -rf "$WORK"' EXIT
: > "$LOG"; echo "gate,case,value,result" > "$CSV"; fail=0
say() { printf '%s\n' "$*" | tee -a "$LOG"; }
rec() { echo "$1,$2,$3,$4" >> "$CSV"; [[ $4 == FAIL ]] && fail=1; say "  $1 $2 ($3): $4"; }
[[ -x "$EXEC" ]] || { echo "missing $EXEC -- run 'make' first" >&2; exit 2; }

cat > "$WORK/layout.json" <<'JSON'
{"name": "pf_test", "Lx": 3.0e-4, "mean_gap": 8.0e-5,
 "walls": {"seed_bottom": 3, "seed_top": 4, "n_bumps": 7, "R_range": [2.0e-5, 3.5e-5],
           "height_ratio_range": [0.3, 0.5], "overlap_range": [0.15, 0.3]},
 "theta_ice_init": 40,
 "plugs": [{"x0": 2.2e-4, "x1": 3.0e-4}],
 "caps": [{"x": 6.0e-5, "wall": "bottom", "R": 3.0e-5}],
 "grains": [{"x": 1.4e-4, "y": "axis", "R": 1.2e-5}]}
JSON
"$PY" "$ROOT/preprocess/build_ice_layout.py" "$WORK/layout.json" --eps 1.5e-6 --outdir "$WORK" --no-pdf \
    > "$WORK/build.txt" 2>&1
want=$(awk '/sharp area/{for(i=1;i<=NF;i++) if($i=="area") print $(i+1)}' "$WORK/pf_test.opts" | tr -d ',')
say "builder sharp area: $want"
EXP=(-temp -20 -humidity 1.0 -t_final 0 -pf_output 0)
runic() {  # name extra-args
  local d="$WORK/$1"; shift; mkdir -p "$d"
  folder="$d" "$EXEC" -options_file "$ROOT/inputs/solver.opts" -options_file "$WORK/pf_test.opts" \
      "${EXP[@]}" "$@" > "$d/out.txt" 2>&1; echo $? > "$d/rc"
}

say ""; say "=== P1/P2: IC area and phi range ==="
runic ok
got=$(awk 'NR==1{print $2}' "$WORK/ok/SSA_evo.dat" 2>/dev/null)
ok=$(awk -v g="$got" -v w="$want" 'BEGIN{r=(g-w)/w; if(r<0)r=-r; print (g!="" && r<0.005)?"PASS":"FAIL"}')
rec P1 "ice_area_vs_sharp_${want}" "${got:-none}" "$ok"
grep -q "phi in \[0.0000, 1.0000\]" "$WORK/ok/out.txt" && rec P2 phi_range "[0,1]" PASS || rec P2 phi_range other FAIL

say ""; say "=== P3: guards ==="
head -c 1000 "$WORK/pf_test_phi.bin" > "$WORK/short.bin"
runic badmesh -Lx 3.0e-4 -geom_file "" -Nx 50 -Ny 20
runic short   -phi_file "$WORK/short.bin"
runic missing -phi_file "$WORK/nope.bin"
runic noopt   -phi_file ""
for c in badmesh short missing noopt; do
  say "    $c: $(grep -m1 -E 'PETSC ERROR: .*(phi_file|control points|PHIIC)' "$WORK/$c/out.txt" | head -c 160)"
  rc=$(cat "$WORK/$c/rc"); [[ $rc != 0 ]] && rec P3 "$c" "rc=$rc" PASS || rec P3 "$c" "rc=$rc" FAIL
done
say ""; [[ $fail == 0 ]] && say "verify_phi_file: ALL PASS" || say "verify_phi_file: FAIL"
exit $fail

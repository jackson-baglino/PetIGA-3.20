#!/usr/bin/env bash
# =============================================================================
# studies/cryosuction/verification/verify_concave_bridge.sh [BASE_REF]
#
# Gates for the opt-in IC options -bridge_concaveL / -bridge_concaveR, which
# turn a -bridge_* plug face from convex (ice inside the disc) to concave (ice
# outside the disc, beyond its centre). Like verify_wall_bc.sh these call the
# binary directly: they give pass/fail numbers, not results.
#
#   C1  Regression  -- a convex bridge run (options unset) is bit-for-bit the
#                      BASE_REF binary (default HEAD), and setting both options
#                      to 0 is bit-for-bit the unset run.
#   C2  Geometry    -- in a flat channel of half-gap a, a lens with two concave
#                      faces of radius R = a/cos(60 deg) has the sharp area
#                      g*L + 2*E, E = 2aR - a*sqrt(R^2-a^2) - R^2*asin(a/R);
#                      its convex twin (same apexes, same R) has g*L - 2*E.
#                      The IC ice area (SSA_evo.dat row 0, column 2) must match
#                      each within 1 % (diffuse profile and spline projection).
#   C3  Guards      -- misuse exits non-zero: an entry other than 0/1, a list
#                      of the wrong length, and apexes that leave no ice.
#
# Writes verify_concave_bridge.csv and .log next to this script.
# =============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
EXEC="$ROOT/lunar_cryosuction"
CSV="$HERE/verify_concave_bridge.csv"
LOG="$HERE/verify_concave_bridge.log"
BASE_REF="${1:-HEAD}"
WORK="${TMPDIR:-/tmp}/verify_concave_bridge.$$"
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT
: > "$LOG"
echo "gate,case,value,result" > "$CSV"
fail=0
say() { printf '%s\n' "$*" | tee -a "$LOG"; }
rec() { echo "$1,$2,$3,$4" >> "$CSV"; [[ $4 == FAIL ]] && fail=1; say "  $1 $2 ($3): $4"; }

[[ -x "$EXEC" ]] || { echo "missing $EXEC -- run 'make' first" >&2; exit 2; }

# Flat channel 150 x 50 um (a = 25 um), walls at 60 deg. h = eps/sqrt(2).
A=2.5e-5; R=5.0e-5             # R = a / cos(60 deg)
COMMON=(-options_file "$ROOT/inputs/solver.opts" -periodic 0 -pf_monitor 1 -pf_output 0
        -dim 2 -Nx 212 -Ny 71 -Lx 1.5e-4 -Ly 5.0e-5 -eps 1.0e-6
        -ic_type multi_grains -temp -20 -humidity 1.0
        -wall_faces y0,y1 -contact_angle_deg 60 -bridge_cy 2.5e-5
        -bridge_rL $R -bridge_rR $R)
# Both lenses: apexes on the axis at x = 50 and 100 um (axis length L = 50 um).
CONCAVE=(-bridge_cxL 0.0    -bridge_cxR 1.5e-4 -bridge_concaveL 1 -bridge_concaveR 1)
CONVEX=( -bridge_cxL 1.0e-4 -bridge_cxR 5.0e-5)

runcase() {  # name binary t_final extra-args...
  local name=$1 bin=$2 tf=$3; shift 3
  local dir="$WORK/$name"; mkdir -p "$dir"
  folder="$dir" "$bin" "${COMMON[@]}" -t_final "$tf" -delt_t 1.0e-2 "$@" > "$dir/out.txt" 2>&1
  echo $? > "$dir/rc"
}
same() { cmp -s "$WORK/$1/SSA_evo.dat" "$WORK/$2/SSA_evo.dat" && [[ -s "$WORK/$1/SSA_evo.dat" ]]; }
area0() { awk 'NR==1{print $2; exit}' "$WORK/$1/SSA_evo.dat"; }

say "BASE_REF = $BASE_REF ($(git -C "$ROOT" rev-parse --short=8 "$BASE_REF"))"

# ---------------------------------------------------------------- C1 ----
say ""; say "=== C1: convex bridge, base vs new binary; explicit 0 vs unset ==="
SNAP="$WORK/base"; mkdir -p "$SNAP"
if git -C "$ROOT/.." archive "$BASE_REF" lunar_cryosuction | tar -x -C "$SNAP" \
   && ( cd "$SNAP/lunar_cryosuction" && make >>"$LOG" 2>&1 ); then
  runcase c1_base "$SNAP/lunar_cryosuction/lunar_cryosuction" 2.0e2 "${CONVEX[@]}"
  runcase c1_new  "$EXEC" 2.0e2 "${CONVEX[@]}"
  if same c1_base c1_new; then rec C1 convex_bitwise_vs_base identical PASS
  else rec C1 convex_bitwise_vs_base differs FAIL; fi
else
  rec C1 convex_bitwise_vs_base build_failed SKIP
fi
runcase c1_zero "$EXEC" 2.0e2 "${CONVEX[@]}" -bridge_concaveL 0 -bridge_concaveR 0
if same c1_new c1_zero; then rec C1 explicit_zero_bitwise identical PASS
else rec C1 explicit_zero_bitwise differs FAIL; fi

# ---------------------------------------------------------------- C2 ----
say ""; say "=== C2: IC ice area vs the sharp shape ==="
runcase c2_concave "$EXEC" 1.0e-2 "${CONCAVE[@]}"
runcase c2_convex  "$EXEC" 1.0e-2 "${CONVEX[@]}"
grep -q "faces: left CONCAVE, right CONCAVE" "$WORK/c2_concave/out.txt" \
  && rec C2 concave_banner printed PASS || rec C2 concave_banner missing FAIL
for kind in concave convex; do
  got=$(area0 c2_$kind)
  read -r want ok < <(awk -v a=$A -v r=$R -v got="$got" -v k=$kind 'BEGIN{
      L = 5.0e-5; g = 2*a
      E = 2*a*r - a*sqrt(r*r - a*a) - r*r*atan2(a, sqrt(r*r - a*a))
      w = (k == "concave") ? g*L + 2*E : g*L - 2*E
      rel = (got - w) / w; if (rel < 0) rel = -rel
      printf "%.6e %s\n", w, (got != "" && rel < 0.01) ? "PASS" : "FAIL"}')
  rec C2 "${kind}_area_sharp=${want}" "${got:-none}" "$ok"
done

# ---------------------------------------------------------------- C3 ----
say ""; say "=== C3: guards ==="
runcase c3_bad   "$EXEC" 1.0e-2 "${CONCAVE[@]:0:4}" -bridge_concaveL 2 -bridge_concaveR 1
runcase c3_len   "$EXEC" 1.0e-2 "${CONCAVE[@]:0:4}" -bridge_concaveL 1,1 -bridge_concaveR 1
runcase c3_empty "$EXEC" 1.0e-2 -bridge_cxL 1.0e-4 -bridge_cxR 5.0e-5 -bridge_concaveL 1 -bridge_concaveR 1
for c in c3_bad c3_len c3_empty; do
  rc=$(cat "$WORK/$c/rc")
  [[ $rc != 0 ]] && rec C3 "$c" "rc=$rc" PASS || rec C3 "$c" "rc=$rc" FAIL
done

say ""
if [[ $fail == 0 ]]; then say "verify_concave_bridge: ALL PASS"; else say "verify_concave_bridge: FAIL"; fi
exit $fail

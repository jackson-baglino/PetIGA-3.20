#!/usr/bin/env bash
# studies/cryosuction/scripts/check_open_bc.sh <tests-file>          (plan tool T5)
# studies/cryosuction/scripts/check_open_bc.sh <geometry> <experiment> [-- extra flags]
#
# The plan drops the open phase boundary: every run is built so ice never
# reaches a vapor-reservoir face, and those faces carry -phi_open_bc 0 (a
# 90-degree face for phi). The solver's default is -phi_open_bc 1
# (src/cryosuction_main.c:833), so a file that forgets the flag silently runs
# with ice allowed to flow out. This check fails such a run before it is
# submitted.
#
# For each geom:exp pair it merges the options the way the solver reads them
# (inputs/solver.opts < geometry < experiment < extra flags; last setting wins)
# and works out the reservoir faces the solver would pin:
#   none           when -flag_BC_rhovfix is off, or -periodic is on (the solver
#                  turns rhovfix off for periodic runs, cryosuction_main.c:977);
#   -rhovfix_faces when set; else both faces of -rhovfix_axis when >= 0; else
#                  every face; minus the y0 face in -axisym runs (the axis).
# A reservoir face that is also listed in -wall_faces is a wall, not open.
# If any open reservoir face remains and -phi_open_bc is not off: FAIL.
#
# Exit status 0 when every pair passes, 1 otherwise, 2 on usage errors.
set -uo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
source "$ROOT/scripts/lib/opts.sh"

# merged <file>... -- <extra flags>: "flag value" lines, last setting wins.
merged() {
    local files=() extra=()
    while [[ $# -gt 0 && $1 != "--" ]]; do files+=("$1"); shift; done
    [[ $# -gt 0 ]] && shift
    extra=("$@")
    {
        for f in "${files[@]}"; do sed 's/#.*//' "$f"; done
        # Extra flags: one flag per line, its value (if any) beside it.
        local i=0 n=${#extra[@]}
        while (( i < n )); do
            local k="${extra[$i]}" v=""
            if (( i + 1 < n )) && { [[ ${extra[$((i+1))]} != -* ]] || [[ ${extra[$((i+1))]} =~ ^-?[0-9.]+([eE][-+]?[0-9]+)?$ ]]; }; then
                v="${extra[$((i+1))]}"; i=$((i+2))
            else
                i=$((i+1))
            fi
            echo "$k $v"
        done
    } | awk '$1 ~ /^-[A-Za-z]/ { v=""; for (i=2;i<=NF;i++) v = v (i>2?" ":"") $i; val[$1]=v; seen[$1]=1 }
             END { for (k in seen) print k "\t" val[k] }'
}

get() { awk -F'\t' -v k="$1" '$1==k {print $2; f=1} END {if(!f) print "__unset__"}' <<<"$2" | tr -d '"'"'"; }

is_true() {   # PETSc boolean: bare flag, 1, true, yes, on
    case "$(echo "$1" | tr '[:upper:]' '[:lower:]')" in
        ""|1|true|yes|on|t|y) return 0 ;;
        *) return 1 ;;
    esac
}

check_pair() {
    local g="$1" e="$2"; shift 2
    local gf ef
    gf=$(resolve_opts "$ROOT/inputs/geometry" "$g") || { echo "FAIL $g:$e  geometry not found"; return 1; }
    ef=$(resolve_opts "$ROOT/inputs/experiment" "$e") || { echo "FAIL $g:$e  experiment not found"; return 1; }
    local m; m=$(merged "$ROOT/inputs/solver.opts" "$gf" "$ef" -- "$@")
    local rv per ax faces walls axisym pob dim
    rv=$(get -flag_BC_rhovfix "$m"); per=$(get -periodic "$m"); ax=$(get -rhovfix_axis "$m")
    faces=$(get -rhovfix_faces "$m"); walls=$(get -wall_faces "$m"); axisym=$(get -axisym "$m")
    pob=$(get -phi_open_bc "$m"); dim=$(get -dim "$m")
    [[ $dim == __unset__ ]] && dim=2

    local pinned=()
    if [[ $rv != __unset__ ]] && is_true "$rv" && ! { [[ $per != __unset__ ]] && is_true "$per" && [[ $per != "" ]] && [[ $per != 0 ]]; }; then
        local axes=(x y z)
        if [[ $faces != __unset__ && -n $faces ]]; then
            IFS=', ' read -r -a pinned <<<"$faces"
        elif [[ $ax != __unset__ && $ax -ge 0 ]]; then
            pinned=("${axes[$ax]}0" "${axes[$ax]}1")
        else
            local l
            for (( l=0; l<dim; l++ )); do pinned+=("${axes[$l]}0" "${axes[$l]}1"); done
        fi
    fi
    local open=() f
    for f in ${pinned[@]+"${pinned[@]}"}; do
        f=$(echo "$f" | tr '[:upper:]' '[:lower:]')
        [[ -z $f ]] && continue
        if [[ $axisym != __unset__ ]] && is_true "$axisym" && [[ $f == y0 ]]; then continue; fi
        [[ ",${walls//[[:space:]]/}," == *",$f,"* ]] && continue
        open+=("$f")
    done
    local shown_pob="$pob"; [[ $pob == __unset__ ]] && shown_pob="unset (solver default 1)"
    if (( ${#open[@]} == 0 )); then
        echo "OK   $g:$e  no open reservoir face"
        return 0
    fi
    if [[ $pob != __unset__ ]] && ! is_true "$pob"; then
        echo "OK   $g:$e  open reservoir faces ${open[*]}, -phi_open_bc $pob"
        return 0
    fi
    echo "FAIL $g:$e  open reservoir faces ${open[*]} with -phi_open_bc $shown_pob; add -phi_open_bc 0"
    return 1
}

[[ $# -ge 1 ]] || { sed -n '2,4p' "$0" >&2; exit 2; }
fail=0
if [[ $# -eq 1 ]]; then
    [[ -f $1 ]] || { echo "no such tests file: $1" >&2; exit 2; }
    while IFS= read -r line || [[ -n $line ]]; do
        line="${line%%#*}"; line="$(echo "$line" | tr -d '[:space:]')"
        [[ -z $line ]] && continue
        [[ $line == *:* ]] || { echo "FAIL $line  not a geometry:experiment line"; fail=1; continue; }
        check_pair "${line%%:*}" "${line#*:}" || fail=1
    done < "$1"
else
    g="$1"; e="$2"; shift 2
    [[ ${1:-} == "--" ]] && shift
    check_pair "$g" "$e" "$@" || fail=1
fi
exit $fail

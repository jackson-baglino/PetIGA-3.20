#!/usr/bin/env bash
# lint_c.sh - run clang-tidy (checks in .clang-tidy) on the C solver.
#
#   scripts/lint_c.sh                 # every src/*.c
#   scripts/lint_c.sh src/assembly.c  # chosen files
#
# Report only: nothing is changed. clang-tidy comes from the pip wheel in
# venv_lunar (requirements.txt), so Homebrew and the compilers are untouched.
# Include paths are the makefile's plus the MPI/PETSc arch directory that
# mpicc adds, so PETSC_DIR, PETSC_ARCH and PETIGA_DIR must be set, as for make.
# Exit status is clang-tidy's: nonzero if a file fails to parse.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
: "${PETSC_DIR:?set PETSC_DIR}" "${PETSC_ARCH:?set PETSC_ARCH}" "${PETIGA_DIR:?set PETIGA_DIR}"

TIDY="$ROOT/venv_lunar/bin/clang-tidy"
[[ -x "$TIDY" ]] || { echo "clang-tidy missing: venv_lunar/bin/pip install -r requirements.txt" >&2; exit 1; }

FLAGS=(
    -I"$PETSC_DIR/include" -I"$PETSC_DIR/include/petsc/private"
    -I"$PETSC_DIR/$PETSC_ARCH/include"
    -I"$PETIGA_DIR/include" -I"$PETIGA_DIR/$PETSC_ARCH/include"
    -I"$ROOT/include"
)
[[ "$(uname)" == Darwin ]] && FLAGS+=(-isysroot "$(xcrun --show-sdk-path)")

if [[ $# -gt 0 ]]; then FILES=("$@"); else FILES=(src/*.c); fi
echo "clang-tidy $("$TIDY" --version | sed -n 's/.*LLVM version //p'), config $ROOT/.clang-tidy"
echo "files: ${FILES[*]}"
"$TIDY" --quiet "${FILES[@]}" -- "${FLAGS[@]}"

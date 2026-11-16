#!/usr/bin/env bash
# The only way a model-line agent runs the solver. Caps and lock: run_test.py (CLAUDE.md, Hard limits).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec "$HERE/../../venv_workflows/bin/python" "$HERE/run_test.py" "$@"

#!/usr/bin/env bash
# PostToolUse hook (Write|Edit): run ruff on a Python file Claude just wrote
# or edited in this project, and hand any findings back to Claude as context.
# Report only: never changes the file, never blocks. Config: ruff.toml
# (--force-exclude keeps its excluded folders, e.g. scratch/, out of it).
ROOT="${CLAUDE_PROJECT_DIR:?}"
f=$(jq -r '.tool_input.file_path // .tool_response.filePath // empty')
case "$f" in "$ROOT"/*.py) ;; *) exit 0 ;; esac
[ -f "$f" ] || exit 0
out=$("$ROOT/venv_lunar/bin/ruff" check --quiet --force-exclude --output-format concise "$f" 2>&1) && exit 0
jq -n --arg o "$out" '{hookSpecificOutput: {hookEventName: "PostToolUse",
  additionalContext: ("ruff (report only) on the file just edited:\n" + $o)}}'
exit 0

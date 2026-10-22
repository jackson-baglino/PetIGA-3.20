#!/bin/bash
# SessionStart hook: append this session's ID to .claude/sessions.tsv and print it
# so Claude has it in context (SessionStart stdout is added to the session).
# Columns: local time, session_id, source (startup|resume|clear|compact), git branch, cwd.
# Resume any listed session with:  claude --resume <session_id>
set -u
dir="$(cd "$(dirname "$0")/.." && pwd)"
log="$dir/sessions.tsv"
input="$(cat)"
field() { printf '%s' "$input" | /usr/bin/python3 -c "import json,sys; print(json.load(sys.stdin).get('$1',''))" 2>/dev/null; }
sid="$(field session_id)"
src="$(field source)"
cwd="$(field cwd)"
branch="$(git -C "${cwd:-$dir}" rev-parse --abbrev-ref HEAD 2>/dev/null)"
[ -f "$log" ] || printf 'time\tsession_id\tsource\tbranch\tcwd\n' > "$log"
printf '%s\t%s\t%s\t%s\t%s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$sid" "$src" "$branch" "$cwd" >> "$log"
echo "Session ID: $sid (source: $src). Logged in .claude/sessions.tsv."
exit 0

# Permission files: literature review

These files are the literature-review workflow's permissions. Claude Code or
macOS enforces each one; nothing here is a description only. Every edit to
`.claude/` asks, so only the user changes them. The reasons for each rule are
in `workflows/literature-review/WORKFLOW.md` §4 and `PERMISSIONS.md` §5.

| File | What it is | Enforced by |
|---|---|---|
| `agents/<agent>.settings.json` | one per agent: the tools it may use, deny rules, sandbox | Claude Code: the driver passes it with `--settings`, together with `--setting-sources ""` so no other settings apply |
| `sessions/review-mode.json` | the main session whenever it works with library or run files: no web, no reads of untrusted stores | Claude Code, started as `claude --settings .claude/permissions/literature-review/sessions/review-mode.json` |
| `hooks/review_read_allowlist.py` | review mode's allow-list for Read, Grep and Glob (permission rules cannot express an allow-list) | Claude Code: a PreToolUse hook; exit 2 blocks the call |
| `launcher.allow` | the only scripts `wf_launch` may run, each with its profile | `wf_launch` (to be built) |
| `profiles/*.sb` | macOS sandbox profiles: `offline` (no network), `api`, `anthropic`, `zotero-web` (HTTPS out only), `zotero-local` (Zotero on this machine only); writes only to the run folder | macOS `sandbox-exec` |
| `inputs.json` | each agent's and session's input classes: `untrusted`, `user_private`, `public` | the driver checks every file it hands an agent; the transcript audit checks the reads |
| `accepted_risks.json` | entries the user allows to hold all three legs, with the reason | `check_permissions.py` |

The main session's **normal profile** is the project's `.claude/settings.json`
(and the local `settings.local.json`): it denies reads of the pipeline's
stores, so library text is reachable only in review mode.

**Check:** `venv_workflows/bin/python scripts/literature_review/check_permissions.py`
works out each entry's legs of the lethal trifecta from these files and fails
if one holds all three without an accepted reason, or if a structural rule is
broken (an agent may write, a sandbox may fail open, a store deny or the hook
is missing). `--write` refreshes the tables in `PERMISSIONS.md` and
`WORKFLOW.md`. Tests: `scripts/literature_review/tests/`.

**Verified on 2026-10-28:**
- the profiles under `sandbox-exec` (offline blocks all network and writes
  outside the run folder; api reaches HTTPS only; zotero-local reaches only
  `localhost:23119`);
- review mode in a live headless session (an allowed file read; a library
  store refused to Read; a shell read of it blocked by the sandbox).

**Known limits:**
- macOS profiles cannot filter by host name, so the API scripts enforce their
  host lists themselves.
- Deny rules do not cover a Python script that opens files itself; the sandbox
  does.
- Path-scoped `Write(...)` rules are ignored by Claude Code (it warns); use
  `Edit(...)`, which covers every file-editing tool.

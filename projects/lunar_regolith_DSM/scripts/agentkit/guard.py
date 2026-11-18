#!/usr/bin/env python3
"""agentkit guard: the PreToolUse hook for every script-run agent (doc 29 §4).

Claude Code runs this before each tool call and passes the call as JSON on stdin. Exit 0 lets the call run;
exit 2 blocks it, and what is printed to stderr goes back to the agent as the reason. The launcher sets the
policy through environment variables the agent cannot change:

    AGENTKIT_WORKSPACE    the agent's working folder: every file path must resolve inside it
    AGENTKIT_BASH_ALLOW   JSON list of regular expressions; a Bash command must match one in full
    AGENTKIT_GUARD_LOG    a JSON-lines file outside the workspace; every decision is appended

Fails safe: any error, including a missing variable or unreadable input, blocks the call (exit 2).
Writes to PROVENANCE.md are always blocked, wherever the file is.
"""
import json
import os
import re
import sys
import time
from pathlib import Path

FILE_TOOLS = {"Read": "file_path", "Write": "file_path", "Edit": "file_path", "NotebookEdit": "notebook_path"}
SEARCH_TOOLS = {"Glob", "Grep"}
WRITE_TOOLS = {"Write", "Edit", "NotebookEdit"}
ALWAYS_ALLOWED = {"StructuredOutput"}           # how the agent returns its typed answer
SHELL_META = re.compile(r"[;&|`$<>\n\\]")       # no chaining, substitution, redirection or escapes
PROTECTED_NAMES = {"PROVENANCE.md"}


def inside(path, workspace, cwd):
    p = Path(path)
    if not p.is_absolute():
        p = Path(cwd) / p
    real, ws = os.path.realpath(p), os.path.realpath(workspace)
    return real == ws or real.startswith(ws + os.sep)


def decide(event, workspace, bash_allow):
    """(allowed, reason) for one tool call."""
    tool = event.get("tool_name")
    inp = event.get("tool_input") or {}
    cwd = event.get("cwd") or workspace
    if tool in ALWAYS_ALLOWED:
        return True, "typed answer"
    if tool in FILE_TOOLS:
        path = inp.get(FILE_TOOLS[tool])
        if not path:
            return False, f"{tool} without a path"
        if tool in WRITE_TOOLS and Path(path).name in PROTECTED_NAMES:
            return False, f"{Path(path).name} is written only by the launcher"
        if not inside(path, workspace, cwd):
            return False, f"{path} is outside your working folder {workspace}"
        return True, "inside the working folder"
    if tool in SEARCH_TOOLS:
        pattern = str(inp.get("pattern", ""))
        if tool == "Glob" and (pattern.startswith(("/", "~")) or ".." in Path(pattern).parts):
            return False, f"Glob pattern {pattern!r} reaches outside your working folder"
        if not inside(inp.get("path") or cwd, workspace, cwd):
            return False, f"{tool} path {inp.get('path')} is outside your working folder {workspace}"
        return True, "search inside the working folder"
    if tool == "Bash":
        cmd = str(inp.get("command", "")).strip()
        if SHELL_META.search(cmd):
            return False, "one plain command per call: no ; & | ` $ < > or backslashes"
        if any(re.fullmatch(rx, cmd) for rx in bash_allow):
            return True, "command on the allow list"
        return False, f"command not allowed; the only commands you may run match: {bash_allow}"
    return False, f"tool {tool!r} is not allowed for this agent"


def log(entry):
    path = os.environ.get("AGENTKIT_GUARD_LOG")
    if path:
        with open(path, "a") as f:
            f.write(json.dumps(entry) + "\n")


def main():
    event = {}
    try:
        event = json.load(sys.stdin)
        workspace = os.environ["AGENTKIT_WORKSPACE"]
        bash_allow = json.loads(os.environ["AGENTKIT_BASH_ALLOW"])
        allowed, reason = decide(event, workspace, bash_allow)
    except Exception as e:                       # fail safe: a broken guard blocks
        allowed, reason = False, f"guard error ({type(e).__name__}: {e}); call blocked"
    try:
        log({"time": time.strftime("%Y-%m-%dT%H:%M:%S"), "tool": event.get("tool_name"),
             "input": event.get("tool_input"), "allowed": allowed, "reason": reason})
    except Exception as e:
        allowed, reason = False, f"guard could not log ({e}); call blocked"
    if allowed:
        return 0
    sys.stderr.write(f"Blocked by the agentkit guard: {reason}\n")
    return 2


if __name__ == "__main__":
    sys.exit(main())

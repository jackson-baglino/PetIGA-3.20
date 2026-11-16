#!/usr/bin/env python3
"""Model-development line, PreToolUse guard for headless agents (WORKFLOW.md §3.2, "Path guards").

Loaded only through the agents' --settings files in .claude/permissions/model-development/agents/; it does nothing
in a normal session (no MD_AGENT in the environment). md_agent.py passes the roots for each call:
  MD_READ_ROOTS, MD_WRITE_ROOTS, MD_DENY_ROOTS  JSON lists of absolute folders
  MD_BASH_ALLOW                                 JSON list of regexes; a Bash command must match one
Exit 2 blocks the call and tells the agent why.
"""
import json
import os
import re
import sys


def under(path, roots):
    p = os.path.realpath(path)
    return any(p == r or p.startswith(r.rstrip("/") + "/") for r in roots)


def block(msg):
    sys.stderr.write("Blocked by the model-line guard: " + msg + "\n")
    sys.exit(2)


def main():
    if not os.environ.get("MD_AGENT"):
        return 0
    data = json.load(sys.stdin)
    tool = data.get("tool_name")
    inp = data.get("tool_input") or {}
    cwd = os.environ.get("MD_CWD") or data.get("cwd") or os.getcwd()
    reads = json.loads(os.environ.get("MD_READ_ROOTS", "[]"))
    writes = json.loads(os.environ.get("MD_WRITE_ROOTS", "[]"))
    denies = json.loads(os.environ.get("MD_DENY_ROOTS", "[]"))
    if tool == "Bash":
        cmd = str(inp.get("command", ""))
        allow = json.loads(os.environ.get("MD_BASH_ALLOW", "[]"))
        if re.search(r"[;&|`]|\$\(", cmd.replace("2>&1", "")):
            block("one command per call; no ;, &&, |, backticks or $( ).")
        if not any(re.fullmatch(a, cmd.strip()) for a in allow):
            block(f"command not in this stage's list: {cmd[:120]}")
        return 0
    p = inp.get("file_path") or inp.get("path") or inp.get("notebook_path")
    if p is None:
        if tool in ("Grep", "Glob"):
            p = cwd
        else:
            return 0
    p = p if os.path.isabs(p) else os.path.join(cwd, p)
    if under(p, denies):
        block(f"{p} is not open to {os.environ['MD_AGENT']} in this stage.")
    if tool in ("Write", "Edit", "NotebookEdit") and not under(p, writes):
        block(f"{os.environ['MD_AGENT']} may write only under {writes}.")
    if tool in ("Read", "Grep", "Glob") and not under(p, reads):
        block(f"{p} is outside this stage's read roots.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

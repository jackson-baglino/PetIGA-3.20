#!/usr/bin/env python3
"""Start the B1 foreman session: `claude --agent beth-foreman` with its definition taken from
workflows/literature-review/agents/beth-foreman.md (WORKFLOW.md §0).

    venv_workflows/bin/python scripts/literature_review/foreman.py [extra claude arguments]

Every literature-review agent is defined in that one folder. Claude Code finds session agents only in
.claude/agents/, so this passes the definition with --agents (a JSON object: the frontmatter fields plus the body
as "prompt"; checked with Claude Code 2.1.316 on 2026-11-01 that the definition's PreToolUse hook applies) and
selects it with --agent. --print-json prints that object instead of starting a session.
"""
import json
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import agent_call  # noqa: E402

NAME = "beth-foreman"


def definition(name=NAME):
    a = agent_call.load(name)
    d = {"description": a["description"], "prompt": a["system_prompt"], "tools": a["tools"], "model": a["model"],
         "permissionMode": a["permission_mode"]}
    if a["disallowed"]:
        d["disallowedTools"] = a["disallowed"]
    if a["skills"]:
        d["skills"] = a["skills"]
    if a["hooks"]:
        d["hooks"] = a["hooks"]
    return {name: d}


def main(argv):
    agents = json.dumps(definition())
    if "--print-json" in argv:
        print(agents)
        return 0
    os.chdir(agent_call.PROJECT)
    os.execvp("claude", ["claude", "--agents", agents, "--agent", NAME] + argv)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

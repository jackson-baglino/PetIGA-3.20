#!/usr/bin/env python3
"""One headless Claude Code call in a clean folder (WORKFLOW.md §4).

The agent sees only the files copied into the folder, may use only Read, writes nothing
and returns JSON checked against a schema. Its settings file denies every other tool.
"""
import json
import subprocess
import time
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[2]
AGENTS = PROJECT / ".claude" / "permissions" / "literature-review" / "agents"
SCHEMAS = PROJECT / "workflows" / "literature-review" / "schemas"
PROMPTS = PROJECT / "workflows" / "literature-review" / "prompts"
MODELS = {"morty": "claude-sonnet-5-5", "beth": "claude-opus-5-5", "jerry": "claude-haiku-5-5"}


def prompt(name, **fields):
    text = (PROMPTS / f"{name}.md").read_text()
    for k, v in fields.items():
        text = text.replace("{" + k + "}", v)
    return text


def call(folder, text, schema_name, agent, timeout=3600):
    """(structured output or None, meta) for one call run in `folder`."""
    schema = json.loads((SCHEMAS / f"{schema_name}.schema.json").read_text())
    model = MODELS[agent.split("-")[0]]
    cmd = ["claude", "-p", text, "--model", model, "--output-format", "json",
           "--tools", "Read", "--allowedTools", "Read", "--permission-mode", "dontAsk",
           "--setting-sources", "", "--strict-mcp-config", "--no-session-persistence",
           "--disable-slash-commands", "--settings", str(AGENTS / f"{agent}.settings.json"),
           "--json-schema", json.dumps(schema)]
    t0 = time.time()
    p = subprocess.run(cmd, cwd=folder, capture_output=True, text=True, timeout=timeout)
    meta = {"agent": agent, "model": model, "seconds": round(time.time() - t0, 1), "exit": p.returncode}
    if p.returncode != 0:
        return None, meta
    out = json.loads(p.stdout)
    meta.update({"usage": out.get("usage"), "num_turns": out.get("num_turns"),
                 "is_error": out.get("is_error")})
    return out.get("structured_output"), meta

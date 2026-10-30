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


def call(folder, text, schema_name, agent, timeout=3600, ledger=None, **tag):
    """(structured output or None, meta) for one call run in `folder`.

    With `ledger`, meta (tagged with `tag`, e.g. step=...) is appended to that JSONL file as
    soon as the call returns, failed calls included, so token usage survives a crash."""
    schema = json.loads((SCHEMAS / f"{schema_name}.schema.json").read_text())
    model = MODELS[agent.split("-")[0]]
    cmd = ["claude", "-p", text, "--model", model, "--output-format", "json",
           "--tools", "Read", "--allowedTools", "Read", "--permission-mode", "dontAsk",
           "--setting-sources", "", "--strict-mcp-config", "--no-session-persistence",
           "--disable-slash-commands", "--settings", str(AGENTS / f"{agent}.settings.json"),
           "--json-schema", json.dumps(schema)]
    t0 = time.time()
    p = subprocess.run(cmd, cwd=folder, capture_output=True, text=True, timeout=timeout)
    meta = {**tag, "agent": agent, "model": model, "seconds": round(time.time() - t0, 1), "exit": p.returncode}
    try:                                    # a failed call may still print its JSON, with usage
        out = json.loads(p.stdout)
    except ValueError:
        out = {}
    meta.update({"usage": out.get("usage"), "num_turns": out.get("num_turns"),
                 "is_error": out.get("is_error")})
    if p.returncode != 0:
        meta["error"] = (str(out.get("result") or p.stderr or "")).strip()[:300]
    if ledger is not None:
        with open(ledger, "a") as f:
            f.write(json.dumps({"time": time.strftime("%Y-%m-%dT%H:%M:%S"), **meta}) + "\n")
    return (out.get("structured_output") if p.returncode == 0 else None), meta


TOKENS = ("input_tokens", "cache_creation_input_tokens", "cache_read_input_tokens", "output_tokens")


def usage_total(metas):
    """Token counts summed over calls; calls that returned no usage are counted separately."""
    tot = {k: 0 for k in TOKENS}
    for m in metas:
        for k in TOKENS:
            tot[k] += (m.get("usage") or {}).get(k) or 0
    tot["calls"] = len(metas)
    tot["calls_without_usage"] = sum(1 for m in metas if not m.get("usage"))
    return tot


def ledger_calls(ledger):
    """Every call recorded in a ledger file, oldest first."""
    ledger = Path(ledger)
    return [json.loads(l) for l in ledger.read_text().splitlines() if l.strip()] if ledger.exists() else []

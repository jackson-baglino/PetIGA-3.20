#!/usr/bin/env python3
"""One headless Claude Code call in a clean folder (WORKFLOW.md §0 and §4).

The agent is defined in workflows/literature-review/agents/<agent>.md (subagent format: frontmatter with the pinned
model and tools; the body is the system prompt). The work order is a template in prompts/<name>.md with {fields}.
The call runs through scripts/agentlib.py: the agent sees only the files in the folder, may use only its own tools,
writes nothing and returns JSON checked against a schema; its settings file in
.claude/permissions/literature-review/agents/ adds the deny rules and the sandbox. The transcript is audited, a
failed attempt is retried once, and every attempt is a row in runs.db next to the ledger.
"""
import json
import sys
import uuid
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT / "scripts"))
import agentlib  # noqa: E402

DEFS = PROJECT / "workflows" / "literature-review" / "agents"
AGENTS = PROJECT / ".claude" / "permissions" / "literature-review" / "agents"
SCHEMAS = PROJECT / "workflows" / "literature-review" / "schemas"
PROMPTS = PROJECT / "workflows" / "literature-review" / "prompts"
MODELS = {"morty": "claude-sonnet-5-5", "beth": "claude-opus-5-5", "jerry": "claude-haiku-5-5"}


def prompt(name, **fields):
    text = (PROMPTS / f"{name}.md").read_text()
    for k, v in fields.items():
        text = text.replace("{" + k + "}", v)
    return text


def load(agent):
    """The agent's definition; its model must be the one MODELS pins for its name prefix."""
    a = agentlib.load_agent(DEFS / f"{agent}.md", agent)
    if a["model"] != MODELS[agent.split("-")[0]]:
        raise agentlib.AgentError(f"{agent}: model {a['model']} differs from MODELS")
    return a


def call(folder, text, schema_name, agent, timeout=3600, ledger=None, **tag):
    """(structured output or None, meta) for one call run in `folder`, retried once on failure.

    With `ledger`, meta (tagged with `tag`, e.g. step=...) is appended to that JSONL file as soon as the call
    returns, failed calls included, so token usage survives a crash; runs.db and the transcripts go in the
    ledger's folder."""
    schema = json.loads((SCHEMAS / f"{schema_name}.schema.json").read_text())
    a = load(agent)
    stage = str(tag.get("step", "call"))
    label = "_".join(f"{k}{v}" for k, v in tag.items() if k != "step").replace("/", "-").replace(" ", "")[:80] \
        or "call"
    label = f"{label}_{uuid.uuid4().hex[:8]}"
    rdir = Path(ledger).parent if ledger is not None else None
    try:
        out, info = agentlib.call_agent(rdir, stage, a, label, text, schema, folder=folder,
                                        settings=AGENTS / f"{agent}.settings.json", timeout=timeout)
    except agentlib.AgentError as e:
        out, info = None, e.info
    meta = {**tag, "agent": agent, "model": a["model"], "seconds": info.get("seconds"),
            "exit": 0 if out is not None else 1, "usage": info.get("usage"), "num_turns": info.get("num_turns"),
            "is_error": info.get("is_error"), "attempts": len(info.get("attempts", []))}
    if out is None:
        meta["error"] = info.get("error")
    if ledger is not None:
        with open(ledger, "a") as f:
            f.write(json.dumps({"time": agentlib._now()[:19], **meta}) + "\n")
    return out, meta


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

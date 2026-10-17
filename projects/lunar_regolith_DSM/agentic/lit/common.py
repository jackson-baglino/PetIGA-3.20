"""Shared helpers for the literature pipeline: agent definitions, headless Claude calls, token ledger.

Agent definitions live in .claude/agents/<name>.md (frontmatter: name, description, model, tools). The same file serves
interactive subagents and the scripted workers below, so there is one source of truth for every prompt.
Workers are run with `claude -p` from a clean directory with no tools, no MCP and no settings, which cuts the
per-call overhead from about 27k tokens to about 1k.
"""
from __future__ import annotations

import json
import os
import re
import subprocess
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]                         # project root
AGENT_DIR = ROOT / ".claude" / "agents"
STATE = Path(os.environ.get("AGENT_STATE_DIR", ROOT / ".agent_state"))
REVIEW = Path(os.environ.get("LIT_REVIEW_DIR", ROOT / "research" / "Review"))   # records, memos (iCloud, never in git)
CLEAN_CWD = Path.home() / ".cache" / "lit_workers" / "clean"
LEDGER = STATE / "ledger.jsonl"


def parse_agent(name: str) -> dict:
    text = (AGENT_DIR / f"{name}.md").read_text()
    m = re.match(r"^---\n(.*?)\n---\n(.*)$", text, re.S)
    meta = dict(line.split(":", 1) for line in m.group(1).splitlines() if ":" in line)
    return {"name": name, "model": meta.get("model", "sonnet").strip(), "system": m.group(2).strip(),
            "description": meta.get("description", "").strip()}


def call_claude(agent: str, user_text: str, *, model: str | None = None, schema: dict | None = None,
                label: str = "", timeout: int = 1800) -> tuple[object, dict]:
    """Run one headless worker. Returns (parsed_output, usage_dict). Raises RuntimeError on failure."""
    a = parse_agent(agent)
    m = model or a["model"]
    CLEAN_CWD.mkdir(parents=True, exist_ok=True)
    cmd = ["claude", "-p", "--model", m, "--output-format", "json", "--no-session-persistence", "--tools", "",
           "--strict-mcp-config", "--setting-sources", "", "--disable-slash-commands", "--system-prompt", a["system"]]
    if schema is not None:
        cmd += ["--json-schema", json.dumps(schema)]
    t0 = time.time()
    p = subprocess.run(cmd, input=user_text, capture_output=True, text=True, cwd=CLEAN_CWD, timeout=timeout)
    if p.returncode != 0:
        raise RuntimeError(f"claude -p failed ({p.returncode}): {p.stderr[-400:]}")
    # stdout may start with warning lines; the JSON object is the last line that parses
    out = None
    for line in reversed(p.stdout.strip().splitlines()):
        try:
            out = json.loads(line)
            break
        except json.JSONDecodeError:
            continue
    if out is None:
        out = json.loads(p.stdout[p.stdout.index("{"):])
    if out.get("is_error"):
        raise RuntimeError(f"claude reported an error: {str(out.get('result'))[:300]}")
    u = out.get("usage", {})
    usage = {"agent": agent, "model": m, "label": label, "input": u.get("input_tokens", 0),
             "cache_read": u.get("cache_read_input_tokens", 0), "cache_write": u.get("cache_creation_input_tokens", 0),
             "output": u.get("output_tokens", 0), "cost_usd_list": out.get("total_cost_usd"),
             "seconds": round(time.time() - t0, 1), "ts": time.strftime("%Y-%m-%dT%H:%M:%S")}
    STATE.mkdir(parents=True, exist_ok=True)
    with open(LEDGER, "a") as f:
        f.write(json.dumps(usage) + "\n")
    parsed = out.get("structured_output")
    if parsed is None:
        txt = (out.get("result") or "").strip()
        txt = re.sub(r"^```(?:json)?\s*|\s*```$", "", txt)
        parsed = json.loads(txt) if txt.startswith(("{", "[")) else txt
    return parsed, usage

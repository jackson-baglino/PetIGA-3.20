"""Run a project agent headless for the units-verification workflow, and audit what it did.

call_agent() takes the agent's definition from .claude/agents/<name>.md (model, tools, body as system prompt),
copies the given images into a fresh empty folder outside the project (so no CLAUDE.md, memory or units table is
in reach) and runs

    claude -p <task> --model <model> --tools <tools> --permission-mode dontAsk --setting-sources ""
           --strict-mcp-config --no-session-persistence --disable-slash-commands --system-prompt <body>
           --output-format stream-json --verbose --json-schema <schema>

from that folder. Reading outside the folder is refused by the harness in dontAsk mode. The transcript is then
audited: the session model must be the pinned one, and every tool call must be a Read of one of the given images.
The structured output is validated against the schema. A failed call is retried once, then reported.
Each attempt is a row in <run>/runs.db; each transcript is <run>/transcripts/<stage>/<label>_try<k>.jsonl.
"""
import json
import shutil
import sqlite3
import subprocess
import tempfile
import time
from pathlib import Path

import uvlib

AGENTS_DIR = uvlib.ROOT / ".claude" / "agents"
TIMEOUT_S = 1200


class AgentError(RuntimeError):
    pass


def load_agent(name):
    path = AGENTS_DIR / f"{name}.md"
    text = path.read_text()
    if not text.startswith("---\n"):
        raise AgentError(f"{path}: no frontmatter")
    head, body = text[4:].split("\n---\n", 1)
    meta = {}
    for line in head.splitlines():
        k, _, v = line.partition(":")
        meta[k.strip()] = v.strip()
    if meta.get("name") != name:
        raise AgentError(f"{path}: name {meta.get('name')!r} does not match the file")
    tools = [t.strip() for t in meta.get("tools", "").split(",") if t.strip()]
    return {"name": name, "model": meta["model"], "tools": tools, "system_prompt": body.strip(),
            "path": uvlib.rel(path), "sha256": uvlib.sha256(path)}


def validate(obj, schema, where="$"):
    """Minimal JSON-schema check (type, enum, required, properties, additionalProperties, items). Returns errors."""
    errs = []
    types = schema.get("type")
    if types is not None:
        types = types if isinstance(types, list) else [types]
        pytypes = {"object": dict, "array": list, "string": str, "boolean": bool, "null": type(None),
                   "integer": int, "number": (int, float)}
        ok = any(isinstance(obj, pytypes[t]) and not (t in ("integer", "number") and isinstance(obj, bool))
                 for t in types)
        if not ok:
            return [f"{where}: expected {types}, got {type(obj).__name__}"]
    if "enum" in schema and obj not in schema["enum"]:
        errs.append(f"{where}: {obj!r} not in {schema['enum']}")
    if isinstance(obj, dict):
        props = schema.get("properties", {})
        for k in schema.get("required", []):
            if k not in obj:
                errs.append(f"{where}: missing {k!r}")
        if schema.get("additionalProperties") is False:
            errs += [f"{where}: unexpected {k!r}" for k in obj if k not in props]
        for k, v in obj.items():
            if k in props:
                errs += validate(v, props[k], f"{where}.{k}")
    if isinstance(obj, list) and "items" in schema:
        for i, v in enumerate(obj):
            errs += validate(v, schema["items"], f"{where}[{i}]")
    return errs


def _db(rdir):
    con = sqlite3.connect(Path(rdir) / "runs.db")
    con.execute("""CREATE TABLE IF NOT EXISTS calls (
        id INTEGER PRIMARY KEY, stage TEXT, agent TEXT, model TEXT, label TEXT, attempt INTEGER,
        input_sha256 TEXT, transcript TEXT, status TEXT, problems TEXT, started TEXT, seconds REAL,
        input_tokens INTEGER, output_tokens INTEGER, cache_read_tokens INTEGER, cache_creation_tokens INTEGER)""")
    return con


def _abs(path, clean):
    """A Read path as absolute: a relative path is relative to the clean folder (the session's working folder)."""
    p = Path(str(path))
    return str(p if p.is_absolute() else Path(clean) / p)


def _audit(events, agent, clean, image_names):
    """Problems with what the session did: wrong model, a tool other than Read, a Read of another file."""
    problems, result, model = [], None, None
    allowed = {str(Path(clean) / n) for n in image_names} | {str(Path(clean).resolve() / n) for n in image_names}
    for e in events:
        if e.get("type") == "system" and e.get("subtype") == "init":
            model = e.get("model")
        if e.get("type") == "assistant":
            for c in e.get("message", {}).get("content", []):
                if c.get("type") != "tool_use" or c.get("name") == "StructuredOutput":
                    continue
                if c.get("name") not in agent["tools"]:
                    problems.append(f"used tool {c.get('name')}")
                elif c["name"] == "Read" and _abs(c.get("input", {}).get("file_path"), clean) not in allowed:
                    problems.append(f"read {c.get('input', {}).get('file_path')}")
        if e.get("type") == "result":
            result = e
    if model != agent["model"]:
        problems.append(f"session model {model!r}, pinned {agent['model']!r}")
    if result is None:
        problems.append("no result event")
    elif result.get("subtype") != "success" or result.get("is_error"):
        problems.append(f"result {result.get('subtype')}: {str(result.get('result'))[:200]}")
    return problems, result


def call_agent(rdir, stage, agent_name, label, task, images, schema, check=None, retries=1):
    """Run one agent task. images: {name in the clean folder: source path}. check(output) -> extra problems.

    Returns (output, info). Raises AgentError when every attempt failed; info lists the problems.
    """
    agent = load_agent(agent_name)
    rdir = Path(rdir)
    tdir = rdir / "transcripts" / stage
    tdir.mkdir(parents=True, exist_ok=True)
    input_sha = uvlib.sha256_bytes(json.dumps(
        {"agent": agent["sha256"], "task": task, "schema": schema,
         "images": {n: uvlib.sha256(p) for n, p in sorted(images.items())}}, sort_keys=True).encode())
    history = []
    for attempt in range(1, retries + 2):
        started, t0 = uvlib.now(), time.time()
        with tempfile.TemporaryDirectory(prefix="uv_clean_") as clean:
            for n, p in images.items():
                shutil.copyfile(p, Path(clean) / n)
            cmd = ["claude", "-p", task, "--model", agent["model"], "--tools", ",".join(agent["tools"]),
                   "--permission-mode", "dontAsk", "--setting-sources", "", "--strict-mcp-config",
                   "--no-session-persistence", "--disable-slash-commands", "--system-prompt", agent["system_prompt"],
                   "--output-format", "stream-json", "--verbose", "--json-schema", json.dumps(schema)]
            try:
                p = subprocess.run(cmd, cwd=clean, capture_output=True, text=True, timeout=TIMEOUT_S)
                out_text, err_text = p.stdout, p.stderr
            except subprocess.TimeoutExpired as e:
                out_text, err_text = (e.stdout or b"").decode() if isinstance(e.stdout, bytes) else (e.stdout or ""), "timeout"
            transcript = tdir / f"{label}_try{attempt}.jsonl"
            transcript.write_text(out_text)
            events = []
            for line in out_text.splitlines():
                try:
                    events.append(json.loads(line))
                except json.JSONDecodeError:
                    pass
            problems, result = _audit(events, agent, clean, list(images))
        if err_text and not events:
            problems.append(f"stderr: {err_text.strip()[:300]}")
        output = (result or {}).get("structured_output")
        if output is None and not problems:
            problems.append("no structured output")
        if output is not None:
            problems += validate(output, schema)
            if not problems and check:
                problems += check(output)
        usage = (result or {}).get("usage", {}) or {}
        con = _db(rdir)
        con.execute("INSERT INTO calls (stage, agent, model, label, attempt, input_sha256, transcript, status, problems,"
                    " started, seconds, input_tokens, output_tokens, cache_read_tokens, cache_creation_tokens)"
                    " VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                    (stage, agent_name, agent["model"], label, attempt, input_sha, uvlib.rel(transcript),
                     "ok" if not problems else "failed", json.dumps(problems), started, round(time.time() - t0, 1),
                     usage.get("input_tokens"), usage.get("output_tokens"), usage.get("cache_read_input_tokens"),
                     usage.get("cache_creation_input_tokens")))
        con.commit()
        con.close()
        history.append({"attempt": attempt, "problems": problems, "transcript": uvlib.rel(transcript)})
        if not problems:
            return output, {"agent": agent_name, "model": agent["model"], "agent_sha256": agent["sha256"],
                            "input_sha256": input_sha, "attempts": history}
    raise AgentError(json.dumps({"label": label, "attempts": history}))

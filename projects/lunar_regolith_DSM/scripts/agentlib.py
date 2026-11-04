"""Shared runner for headless workflow agents (workflows/TEMPLATE.md, "The line" and §6).

An agent is a Markdown file in the Claude Code subagent format: YAML frontmatter (name, description, tools,
disallowedTools, model, permissionMode, maxTurns) and the system prompt as the body. Worker definitions live in
workflows/<name>/agents/; the units-verification workers are still in .claude/agents/.

call_agent() runs one task as

    claude -p <task> --model <model> --tools <tools> --permission-mode dontAsk --setting-sources ""
           --strict-mcp-config --no-session-persistence --disable-slash-commands --system-prompt <body>
           --output-format stream-json --verbose --json-schema <schema> [--settings <file>] [--allowedTools ...]

from a clean folder that holds only the agent's inputs, with CLAUDE_CODE_DISABLE_AUTO_MEMORY=1:
`--setting-sources ""` keeps CLAUDE.md out but not auto memory, and a worker started inside the
repository loaded its memory index (canary test, 2026-11-03). `--tools ""` leaves only the StructuredOutput tool
(checked with Claude Code 2.1.316 on 2026-11-01). The transcript is audited: the session model must be the pinned
one, every tool call must be one the agent may use, and every Read must be of a file in the clean folder. The
output is validated against the schema. A failed attempt is retried once, then reported. With a run folder, each
attempt is a row in <rdir>/runs.db and its transcript is <rdir>/transcripts/<stage>/<label>_try<k>.jsonl.
"""
import datetime
import hashlib
import json
import os
import shutil
import sqlite3
import subprocess
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]                     # projects/lunar_regolith_DSM
FAMILY_PREFIX = {"claude-fable-": "rick-", "claude-opus-": "beth-", "claude-sonnet-": "morty-",
                 "claude-haiku-": "jerry-"}                    # CLAUDE.md, "Agent names"
TIMEOUT_S = 1200


class AgentError(RuntimeError):
    pass


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rel(path):
    p = Path(path).resolve()
    try:
        return str(p.relative_to(ROOT))
    except ValueError:
        return str(p)


def _frontmatter(text, path):
    """(frontmatter dict, body). The frontmatter is YAML (nested fields such as hooks are allowed)."""
    import yaml
    if not text.startswith("---\n"):
        raise AgentError(f"{path}: no frontmatter")
    head, body = text[4:].split("\n---\n", 1)
    meta = yaml.safe_load(head) or {}
    if not isinstance(meta, dict):
        raise AgentError(f"{path}: frontmatter is not a mapping")
    return meta, body


def _list(value):
    if value is None:
        return []
    if isinstance(value, list):
        return [str(t).strip() for t in value if str(t).strip()]
    return [t.strip() for t in str(value).split(",") if t.strip()]


def load_agent(path, name=None):
    """The agent defined in `path`. Checks the name against the file name (and `name` if given), that the model is a
    full pinned ID, and that the name prefix matches the model family."""
    path = Path(path)
    meta, body = _frontmatter(path.read_text(), path)
    aname = meta.get("name")
    if aname != path.stem or (name is not None and aname != name):
        raise AgentError(f"{path}: name {aname!r} does not match the file")
    model = str(meta.get("model") or "")
    fam = next((k for k in FAMILY_PREFIX if model.startswith(k)), None)
    if fam is None or not model[len(fam):].replace("-", "").isdigit():
        raise AgentError(f"{path}: model {model!r} is not a full pinned ID")
    if not aname.startswith(FAMILY_PREFIX[fam]):
        raise AgentError(f"{path}: name {aname!r} does not match model family {fam!r}")
    if not meta.get("description"):
        raise AgentError(f"{path}: no description")
    max_turns = meta.get("maxTurns")
    return {"name": aname, "model": model, "description": meta["description"], "tools": _list(meta.get("tools")),
            "disallowed": _list(meta.get("disallowedTools")), "max_turns": int(max_turns) if max_turns else None,
            "permission_mode": meta.get("permissionMode", "dontAsk"), "skills": _list(meta.get("skills")),
            "hooks": meta.get("hooks"), "system_prompt": body.strip(), "path": rel(path), "sha256": sha256(path),
            "session": meta.get("permissionMode", "dontAsk") != "dontAsk"}


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


def _abs(path, clean):
    """A Read path as absolute: a relative path is relative to the clean folder (the session's working folder)."""
    p = Path(str(path))
    return str(p if p.is_absolute() else Path(clean) / p)


def audit(events, agent, clean, readable):
    """(problems, result event): wrong model, a tool the agent may not use, a Read of a file not in `readable`
    (names relative to the clean folder), or no successful result."""
    problems, result, model = [], None, None
    allowed = {str(Path(clean) / n) for n in readable} | {str(Path(clean).resolve() / n) for n in readable}
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


def worker_env():
    """The environment of a headless worker: ours, with Claude Code's auto memory off."""
    return dict(os.environ, CLAUDE_CODE_DISABLE_AUTO_MEMORY="1")


def command(agent, task, schema, settings=None, allow=()):
    cmd = ["claude", "-p", task, "--model", agent["model"], "--tools", ",".join(agent["tools"]),
           "--permission-mode", agent["permission_mode"], "--setting-sources", "", "--strict-mcp-config",
           "--no-session-persistence", "--disable-slash-commands", "--system-prompt", agent["system_prompt"],
           "--output-format", "stream-json", "--verbose", "--json-schema", json.dumps(schema)]
    if agent["disallowed"]:
        cmd += ["--disallowedTools", ",".join(agent["disallowed"])]
    if agent.get("max_turns"):
        cmd += ["--max-turns", str(agent["max_turns"])]
    if settings:
        cmd += ["--settings", str(settings)]
    if allow:
        cmd += ["--allowedTools", ",".join(allow)]
    return cmd


def _db(rdir):
    con = sqlite3.connect(Path(rdir) / "runs.db")
    con.execute("""CREATE TABLE IF NOT EXISTS calls (
        id INTEGER PRIMARY KEY, stage TEXT, agent TEXT, model TEXT, label TEXT, attempt INTEGER,
        input_sha256 TEXT, transcript TEXT, status TEXT, problems TEXT, started TEXT, seconds REAL,
        input_tokens INTEGER, output_tokens INTEGER, cache_read_tokens INTEGER, cache_creation_tokens INTEGER)""")
    return con


def _files_in(folder):
    folder = Path(folder)
    return sorted(str(p.relative_to(folder)) for p in folder.rglob("*") if p.is_file())


def _now():
    return datetime.datetime.now().astimezone().isoformat(timespec="seconds")


def call_agent(rdir, stage, agent, label, task, schema, files=None, folder=None, settings=None, check=None,
               retries=1, allow=(), timeout=TIMEOUT_S):
    """Run one task for `agent` (a dict from load_agent). Give either `files` ({name in the clean folder: source
    path}, copied into a fresh folder per attempt) or `folder` (an already prepared clean folder). check(output) ->
    extra problems. rdir None: no runs.db and no transcript.

    Returns (output, info); info has the agent, model, hashes, every attempt's problems, and the last attempt's
    usage, turns, seconds and error. Raises AgentError, with info attached as .info, when every attempt failed."""
    if (files is None) == (folder is None):
        raise ValueError("give exactly one of files and folder")
    tdir = None
    if rdir is not None:
        tdir = Path(rdir) / "transcripts" / stage
        tdir.mkdir(parents=True, exist_ok=True)
    srcs = {n: sha256(p) for n, p in sorted((files or {}).items())} if files is not None else \
        {n: sha256(Path(folder) / n) for n in _files_in(folder)}
    input_sha = hashlib.sha256(json.dumps({"agent": agent["sha256"], "task": task, "schema": schema,
                                           "files": srcs}, sort_keys=True).encode()).hexdigest()
    history, last = [], {}
    for attempt in range(1, retries + 2):
        started, t0 = _now(), time.time()
        tmp = tempfile.TemporaryDirectory(prefix="agent_clean_") if files is not None else None
        clean = tmp.name if tmp else str(folder)
        try:
            for n, p in (files or {}).items():
                (Path(clean) / n).parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(p, Path(clean) / n)
            readable = _files_in(clean)
            try:
                p = subprocess.run(command(agent, task, schema, settings, allow), cwd=clean, capture_output=True,
                                   text=True, timeout=timeout, env=worker_env())
                out_text, err_text = p.stdout, p.stderr
            except subprocess.TimeoutExpired as e:
                out_text = e.stdout.decode() if isinstance(e.stdout, bytes) else (e.stdout or "")
                err_text = "timeout"
            events = []
            for line in out_text.splitlines():
                try:
                    events.append(json.loads(line))
                except json.JSONDecodeError:
                    pass
            problems, result = audit(events, agent, clean, readable)
        finally:
            if tmp:
                tmp.cleanup()
        transcript = None
        if tdir is not None:
            transcript = tdir / f"{label}_try{attempt}.jsonl"
            transcript.write_text(out_text)
        if err_text and not events:
            problems.append(f"stderr: {err_text.strip()[:300]}")
        output = (result or {}).get("structured_output")
        if output is None and not problems:
            problems.append("no structured output")
        if output is not None:
            problems += validate(output, schema)
            if not problems and check:
                problems += check(output)
        usage = (result or {}).get("usage") or {}
        seconds = round(time.time() - t0, 1)
        if rdir is not None:
            con = _db(rdir)
            con.execute("INSERT INTO calls (stage, agent, model, label, attempt, input_sha256, transcript, status,"
                        " problems, started, seconds, input_tokens, output_tokens, cache_read_tokens,"
                        " cache_creation_tokens) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                        (stage, agent["name"], agent["model"], label, attempt, input_sha,
                         rel(transcript) if transcript else None, "ok" if not problems else "failed",
                         json.dumps(problems), started, seconds, usage.get("input_tokens"),
                         usage.get("output_tokens"), usage.get("cache_read_input_tokens"),
                         usage.get("cache_creation_input_tokens")))
            con.commit()
            con.close()
        history.append({"attempt": attempt, "problems": problems,
                        "transcript": rel(transcript) if transcript else None})
        last = {"usage": usage or None, "num_turns": (result or {}).get("num_turns"), "seconds": seconds,
                "is_error": (result or {}).get("is_error"), "error": "; ".join(problems)[:300] if problems else None}
        if not problems:
            break
    info = {"agent": agent["name"], "model": agent["model"], "agent_sha256": agent["sha256"],
            "input_sha256": input_sha, "attempts": history, **last}
    if history[-1]["problems"]:
        err = AgentError(json.dumps({"label": label, "attempts": history}))
        err.info = info
        raise err
    return output, info

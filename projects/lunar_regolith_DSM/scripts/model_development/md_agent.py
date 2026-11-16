"""Headless agent calls for the model-development line (WORKFLOW.md §3).

Thin wrapper around scripts/agentlib.py (load_agent, command, validate, runs.db). It differs in one way: B1's
workers worked in a throwaway clean folder and read only what was copied in, while this line's agents must read
the source tree and the library and write into the run folder or a feature worktree. So the audit here checks every
Read/Grep/Glob path against explicit read roots and every Write/Edit path against explicit write roots, and the
same roots are enforced live by the PreToolUse hook (md_guard.py) through environment variables.

A usage-limit error raises UsageLimit and is never retried; run_md.py stops cleanly and prints the resume command.
"""
import hashlib
import json
import os
import re
import subprocess
import sys
import tempfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import agentlib  # noqa: E402

WORKFLOW = ROOT / "workflows" / "model-development"
AGENTS = WORKFLOW / "agents"
SCHEMAS = WORKFLOW / "schemas"
PERMS = ROOT / ".claude" / "permissions" / "model-development" / "agents"
LIMIT_RE = re.compile(r"(usage limit|hit your .*limit|5-hour limit|weekly limit|rate limit)", re.I)


class UsageLimit(RuntimeError):
    pass


def load(name):
    """A workflow agent (workflows/model-development/agents/) or a general helper reused by the line
    (.claude/agents/, e.g. beth-writer)."""
    path = AGENTS / f"{name}.md"
    if not path.is_file():
        path = ROOT / ".claude" / "agents" / f"{name}.md"
    return agentlib.load_agent(path, name)


def schema(name):
    return json.loads((SCHEMAS / f"{name}.schema.json").read_text())


def _under(path, roots):
    p = os.path.realpath(path)
    return any(p == r or p.startswith(r.rstrip("/") + "/") for r in roots)


def _path_of(tool, inp, cwd):
    p = inp.get("file_path") or inp.get("path") or (cwd if tool in ("Grep", "Glob") else None)
    if p is None:
        return None
    p = Path(str(p))
    return str(p if p.is_absolute() else Path(cwd) / p)


def audit(events, agent, cwd, read_roots, write_roots, deny_roots):
    """(problems, denials, bash, result). Problems fail the call; denials and bash commands are recorded."""
    problems, denials, bash, result, model = [], [], [], None, None
    readable = list(read_roots) + list(write_roots)
    # Calls the hook or the permissions blocked did not run: they are recorded as denials, not problems.
    blocked = {d.get("tool_use_id") for e in events if e.get("type") == "result"
               for d in (e.get("permission_denials") or [])}
    for e in events:
        if e.get("type") == "system" and e.get("subtype") == "init":
            model = e.get("model")
        if e.get("type") == "assistant":
            for c in e.get("message", {}).get("content", []):
                if c.get("type") != "tool_use" or c.get("name") == "StructuredOutput":
                    continue
                name, inp = c.get("name"), c.get("input", {}) or {}
                if c.get("id") in blocked:
                    continue
                if name not in agent["tools"]:
                    problems.append(f"used tool {name}")
                    continue
                if name == "Bash":
                    bash.append(str(inp.get("command", ""))[:300])
                    continue
                path = _path_of(name, inp, cwd)
                if path is None:
                    continue
                if _under(path, deny_roots):
                    problems.append(f"{name} of denied path {path}")
                elif name in ("Write", "Edit") and not _under(path, write_roots):
                    problems.append(f"{name} outside write roots: {path}")
                elif name in ("Read", "Grep", "Glob") and not _under(path, readable):
                    problems.append(f"{name} outside read roots: {path}")
        if e.get("type") == "result":
            result = e
            denials = [f"{d.get('tool_name')}: {json.dumps(d.get('tool_input'))[:200]}"
                       for d in (e.get("permission_denials") or [])]
    if model != agent["model"]:
        problems.append(f"session model {model!r}, pinned {agent['model']!r}")
    if result is None:
        problems.append("no result event")
    elif result.get("subtype") != "success" or result.get("is_error"):
        problems.append(f"result {result.get('subtype')}: {str(result.get('result'))[:200]}")
    return problems, denials, bash, result


def call(rdir, stage, agent_name, label, task, schema_name, cwd, read_roots, write_roots, deny_roots=(),
         bash_allow=(), stop_check=None, extra_env=None, retries=1, timeout=3600, check=None, max_turns=None):
    """Run one work order. Returns (output, info). Raises UsageLimit, or agentlib.AgentError after the retries.

    read_roots / write_roots / deny_roots: absolute folders. bash_allow: regexes a Bash command must match (all
    others are blocked by the hook). stop_check: a shell command the Stop hook runs; non-zero blocks the stop
    (at most twice per call). check(output) -> extra problems."""
    agent = load(agent_name)
    if max_turns:
        agent = dict(agent, max_turns=max_turns)
    sch = schema(schema_name)
    settings = PERMS / f"{agent_name}.settings.json"
    if not settings.is_file():
        raise agentlib.AgentError(f"no permission file {settings}")
    roots = {k: [os.path.realpath(r) for r in v] for k, v in
             {"read": read_roots, "write": write_roots, "deny": deny_roots}.items()}
    tdir = Path(rdir) / "transcripts" / stage
    tdir.mkdir(parents=True, exist_ok=True)
    input_sha = hashlib.sha256(json.dumps({"agent": agent["sha256"], "task": task, "schema": sch,
                                           "roots": roots}, sort_keys=True).encode()).hexdigest()
    history, last = [], {}
    for attempt in range(1, retries + 2):
        started, t0 = agentlib._now(), time.time()
        with tempfile.TemporaryDirectory(prefix="md_hook_") as hookdir:
            env = dict(agentlib.worker_env(), MD_RUN=str(rdir), MD_STAGE=stage, MD_AGENT=agent_name,
                       MD_CWD=str(cwd), MD_READ_ROOTS=json.dumps(roots["read"] + roots["write"]),
                       MD_WRITE_ROOTS=json.dumps(roots["write"]), MD_DENY_ROOTS=json.dumps(roots["deny"]),
                       MD_BASH_ALLOW=json.dumps(list(bash_allow)), MD_STOP_CHECK=stop_check or "",
                       MD_HOOK_STATE=hookdir, MD_HOOKS=str(ROOT / ".claude" / "hooks"), **(extra_env or {}))
            cmd = agentlib.command(agent, task, sch, settings)
            try:
                p = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, timeout=timeout, env=env)
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
        problems, denials, bash, result = audit(events, agent, cwd, roots["read"], roots["write"], roots["deny"])
        transcript = tdir / f"{label}_try{attempt}.jsonl"
        transcript.write_text(out_text)
        if err_text and not events:
            problems.append(f"stderr: {err_text.strip()[:300]}")
        text = json.dumps((result or {}).get("result")) + " " + (err_text or "")
        if result is not None and (result.get("is_error") or result.get("subtype") != "success") and LIMIT_RE.search(text):
            _log(rdir, stage, agent, label, attempt, input_sha, transcript, "usage_limit", problems, started, t0, result)
            raise UsageLimit(f"{label}: {str((result or {}).get('result'))[:200]}")
        output = (result or {}).get("structured_output")
        if output is None and not problems:
            problems.append("no structured output")
        if output is not None:
            problems += agentlib.validate(output, sch)
            if not problems and check:
                problems += check(output)
        _log(rdir, stage, agent, label, attempt, input_sha, transcript, "ok" if not problems else "failed",
             problems, started, t0, result, denials)
        history.append({"attempt": attempt, "problems": problems, "denials": denials, "bash": bash,
                        "transcript": agentlib.rel(transcript)})
        usage = (result or {}).get("usage") or {}
        last = {"usage": usage or None, "num_turns": (result or {}).get("num_turns"),
                "seconds": round(time.time() - t0, 1)}
        if not problems:
            break
    info = {"agent": agent["name"], "model": agent["model"], "agent_sha256": agent["sha256"],
            "input_sha256": input_sha, "attempts": history, **last}
    if history[-1]["problems"]:
        err = agentlib.AgentError(json.dumps({"label": label, "attempts": history})[:4000])
        err.info = info
        raise err
    return output, info


def _log(rdir, stage, agent, label, attempt, input_sha, transcript, status, problems, started, t0, result,
         denials=()):
    usage = (result or {}).get("usage") or {}
    con = agentlib._db(rdir)
    con.execute("INSERT INTO calls (stage, agent, model, label, attempt, input_sha256, transcript, status,"
                " problems, started, seconds, input_tokens, output_tokens, cache_read_tokens,"
                " cache_creation_tokens) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                (stage, agent["name"], agent["model"], label, attempt, input_sha, agentlib.rel(transcript), status,
                 json.dumps({"problems": problems, "denials": list(denials)}), started,
                 round(time.time() - t0, 1), usage.get("input_tokens"), usage.get("output_tokens"),
                 usage.get("cache_read_input_tokens"), usage.get("cache_creation_input_tokens")))
    con.commit()
    con.close()

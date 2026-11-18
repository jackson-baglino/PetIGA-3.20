#!/usr/bin/env python3
"""agentkit launcher: the one way to run a script-run agent on one task (doc 28 §4, doc 29 §4).

    venv_workflows/bin/python scripts/agentkit/launch.py <agent> --task <task_id> [--runs-root DIR]

What one call does, in order:
 1. finds the agent by name: exactly one definition in workflows/*/agents/ or .claude/agents/;
 2. loads the definition (name, pinned model, tools) and checks the work order against the task schema;
 3. refuses a real model call while any enforcer in the agent's tests/enforcers.json is still "planned", or
    when no red-team PASS is on record for exactly the current barriers (doc 29 §11, "drift"): the harness
    files, the agent's settings.json, guard.json and env.allow, its tools and model, and the Claude Code
    version. Any change to one of those needs a fresh `--redteam` run before the next real call;
 4. makes a run folder; the agent's working folder (workspace/) holds only a copy of task.json;
 5. composes the settings: the project base denies + the agent's own permissions + the standard hooks
    (guard.py before every tool call, stop_check.py when the agent stops); an agent's file may only narrow;
 6. calls Claude Code headless (`claude -p`) with the definition's body as the system prompt, the composed
    settings, `--setting-sources ""` (no CLAUDE.md, no user or project settings) and the answer schema;
 7. saves the transcript and audits it: every tool call that ran is re-checked against the guard's policy;
 8. checks the answer against the answer schema and the required files, then runs the task's grader;
 9. writes manifest.json and PROVENANCE.md (hash-chained) with the verdict: PASS, FAIL, ABSTAIN or ERROR.

Exit code: 0 PASS, 1 FAIL, 3 ABSTAIN, 4 ERROR, 5 refused before any call.
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

KIT = Path(__file__).resolve().parent
ROOT = KIT.parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(KIT))
import agentlib  # noqa: E402
import guard  # noqa: E402
import stop_check  # noqa: E402
from provenance import Provenance  # noqa: E402

PY = ROOT / "venv_workflows" / "bin" / "python"
TIMEOUT_S = 1800
GRADER_TIMEOUT_S = 300
OWN_KEYS = {"permissions"}                    # the sandbox and hooks come only from the base and the launcher
SET_BY_LAUNCHER = {"PATH", "TMPDIR", "CLAUDE_CODE_DISABLE_AUTO_MEMORY", "CLAUDE_CODE_SUBPROCESS_ENV_SCRUB"}
CLAUDE_BIN = shutil.which("claude") or "claude"   # absolute, so the agent's PATH need not include it
MATCHER = "Read|Write|Edit|NotebookEdit|Glob|Grep|Bash|WebFetch|WebSearch|Agent|Skill|Task"
EXIT = {"PASS": 0, "FAIL": 1, "ABSTAIN": 3, "ERROR": 4}
PREAMBLE = """Written only by `scripts/agentkit/launch.py`. No agent can write it: the guard blocks any write to
PROVENANCE.md, the permission rules deny it, and it sits outside the agent's working folder. Each entry ends
with a chain hash over the entry and the hash before it; check the file with
`venv_workflows/bin/python scripts/agentkit/provenance.py verify <this file>`. The final hash is also in
manifest.json."""


class LaunchError(RuntimeError):
    pass


def now():
    return time.strftime("%Y-%m-%dT%H:%M:%S%z")


def load_json(path):
    return json.loads(Path(path).read_text())


def find_agent(name, root=ROOT):
    """(definition path, agent folder, workflow name). Exactly one definition may carry the name."""
    hits = sorted(root.glob(f"workflows/*/agents/{name}.md")) + sorted(root.glob(f".claude/agents/{name}.md"))
    if len(hits) != 1:
        raise LaunchError(f"agent {name!r}: {len(hits)} definitions found {[str(h) for h in hits]}; need one")
    d = hits[0]
    workflow = "general" if d.parents[1].name == ".claude" else d.parents[1].name
    return d, d.parent / name, workflow


def unbuilt_enforcers(folder):
    f = folder / "tests" / "enforcers.json"
    if not f.is_file():
        return ["tests/enforcers.json is missing"]
    return [f"{e['id']} ({e['rule']}): {e['status']}" for e in load_json(f)["enforcers"] if "planned" in e["status"]]


def barrier_fingerprint(folder, agent, version):
    """What the red-team result depends on: every file that sets a barrier, the agent's tools and model, and the
    Claude Code version (its sandbox is part of the barrier). Prompt text and tasks are not barriers."""
    files = [*sorted(KIT.glob("*.py")), KIT / "base_settings.json", ROOT / "scripts" / "agentlib.py",
             *(folder / n for n in ("settings.json", "guard.json", "env.allow"))]
    fp = {agentlib.rel(f): agentlib.sha256(f) for f in files if f.is_file()}
    fp["(definition) tools and model"] = json.dumps(
        {k: agent[k] for k in ("model", "tools", "disallowed", "permission_mode")}, sort_keys=True)
    fp["(Claude Code) version"] = version
    return fp


def redteam_gap(agent_runs, fingerprint):
    """[] when a real red-team PASS on exactly these barriers is on record in `agent_runs`; else the reason."""
    latest = None
    for m in sorted(Path(agent_runs).glob("*_redteam/manifest.json"), reverse=True):
        try:
            man = load_json(m)
        except (OSError, ValueError):
            continue
        if man.get("verdict") != "PASS" or man.get("fake", True) or not man.get("barrier"):
            continue
        if man["barrier"] == fingerprint:
            return []
        latest = latest or (m.parent.name, man["barrier"])
    if latest is None:
        return [f"no red-team PASS with a barrier record in {agentlib.rel(agent_runs)}"]
    run, old = latest
    changed = sorted(k for k in set(old) | set(fingerprint) if old.get(k) != fingerprint.get(k))
    return [f"the latest red-team PASS ({run}) was on other barriers; changed since: {', '.join(changed)}"]


def fill(value, root, workspace):
    """Replace the base file's placeholders in a string or a nested list/dict."""
    if isinstance(value, dict):
        return {k: fill(v, root, workspace) for k, v in value.items()}
    if isinstance(value, list):
        return [fill(v, root, workspace) for v in value]
    if isinstance(value, str):
        return (value.replace("//ROOT/", "/" + str(root) + "/").replace("{ROOT}", str(root))
                .replace("{GITROOT}", str(git_root(root))).replace("{WORKSPACE}", str(workspace)))
    return value


def git_root(root=ROOT):
    top = subprocess.run(["git", "rev-parse", "--show-toplevel"], cwd=root, capture_output=True, text=True)
    return Path(top.stdout.strip()) if top.returncode == 0 else Path(root)


def compose_settings(folder, out, root=ROOT, workspace=None):
    """Write the settings Claude Code will get: the base permissions and sandbox (always), the agent's own
    permissions (it may add, never remove), and the standard hooks."""
    workspace = Path(workspace) if workspace else Path(out).parent / "workspace"
    base_file = load_json(KIT / "base_settings.json")
    base = base_file["permissions"]["deny"]
    sandbox = fill(base_file["sandbox"], root, workspace)
    own = load_json(folder / "settings.json")
    if set(own) - OWN_KEYS:
        raise LaunchError(f"{folder / 'settings.json'} sets {sorted(set(own) - OWN_KEYS)}; an agent's settings may "
                          f"set only {sorted(OWN_KEYS)} (hooks and modes belong to the launcher)")
    perms = own.get("permissions", {})
    if set(perms) - {"allow", "deny"}:
        raise LaunchError(f"{folder / 'settings.json'}: permissions may hold only allow and deny")
    deny = fill(base, root, workspace)
    deny += [r for r in perms.get("deny", []) if r not in deny]

    def hook(script):
        return {"type": "command", "command": f'"{PY}" "{KIT / script}"', "timeout": 30}

    settings = {"permissions": {"allow": list(perms.get("allow", [])), "deny": deny},
                "hooks": {"PreToolUse": [{"matcher": MATCHER, "hooks": [hook("guard.py")]}],
                          "Stop": [{"hooks": [hook("stop_check.py")]}]},
                "sandbox": sandbox}
    Path(out).write_text(json.dumps(settings, indent=2) + "\n")
    return settings


def env_names(folder):
    f = folder / "env.allow"
    if not f.is_file():
        raise LaunchError(f"{f} is missing: an agent gets only the environment variables it lists")
    names = [s.strip() for s in f.read_text().splitlines() if s.strip() and not s.lstrip().startswith("#")]
    bad = [n for n in names if n in SET_BY_LAUNCHER or n.startswith("AGENTKIT_") or not n.isidentifier()]
    if bad:
        raise LaunchError(f"{f}: {bad} are set by the launcher or are not variable names")
    return names


def build_env(workspace, run_dir, policy, folder, fake=False):
    """Only the variables in the agent's env.allow pass through; the launcher sets the rest."""
    # Not CLAUDE_CODE_SUBPROCESS_ENV_SCRUB: it makes Claude Code force the permission mode to "default"
    # (2026-11-17). Secrets are removed by env.allow here and by the sandbox's credentials deny list.
    env = {n: os.environ[n] for n in env_names(folder) if n in os.environ}
    tmp = workspace / ".tmp"
    tmp.mkdir(exist_ok=True)
    env.update(PATH=f"{PY.parent}{os.pathsep}/usr/bin{os.pathsep}/bin", TMPDIR=str(tmp),
               CLAUDE_CODE_DISABLE_AUTO_MEMORY="1",
               AGENTKIT_WORKSPACE=str(workspace), AGENTKIT_BASH_ALLOW=json.dumps(policy["bash_allow"]),
               AGENTKIT_REQUIRED=json.dumps(policy.get("required_outputs", [])),
               AGENTKIT_GUARD_LOG=str(run_dir / "guard.jsonl"))
    if fake:                                  # the test stand-in reads its scenario from these
        env.update({k: v for k, v in os.environ.items() if k.startswith("AGENTKIT_FAKE_")})
    return env


def parse_events(text):
    events = []
    for line in text.splitlines():
        try:
            events.append(json.loads(line))
        except json.JSONDecodeError:
            pass
    return events


def audit(events, agent, workspace, bash_allow):
    """(problems, result event, calls). Every tool call that ran is re-checked against the guard's policy, so a
    call the guard should have blocked but did not is reported as a barrier failure."""
    problems, result, model, mode, uses, errored = [], None, None, None, {}, {}
    for e in events:
        if e.get("type") == "system" and e.get("subtype") == "init":
            model, mode = e.get("model"), e.get("permissionMode")
        content = (e.get("message") or {}).get("content")
        if e.get("type") == "assistant" and isinstance(content, list):
            for c in content:
                if c.get("type") == "tool_use":
                    uses[c.get("id")] = c
        if e.get("type") == "user" and isinstance(content, list):
            for c in content:
                if c.get("type") == "tool_result":
                    errored[c.get("tool_use_id")] = bool(c.get("is_error"))
        if e.get("type") == "result":
            result = e
    calls = []
    for uid, c in uses.items():
        name, inp = c.get("name"), c.get("input") or {}
        if name == "StructuredOutput":
            continue
        ran = not errored.get(uid, False)
        calls.append({"tool": name, "input": inp, "ran": ran})
        if not ran:
            continue
        if name not in agent["tools"]:
            problems.append(f"BARRIER FAILURE: tool {name} ran but is not in the agent's tools")
            continue
        ok, why = guard.decide({"tool_name": name, "tool_input": inp, "cwd": str(workspace)}, str(workspace),
                               bash_allow)
        if not ok:
            problems.append(f"BARRIER FAILURE: {name} {json.dumps(inp)[:160]} ran although {why}")
    if model != agent["model"]:
        problems.append(f"session model {model!r}, pinned {agent['model']!r}")
    if mode != agent["permission_mode"]:      # found in the first real T01 run: a silent switch to "default"
        problems.append(f"session permission mode {mode!r}, definition says {agent['permission_mode']!r}")
    if result is None:
        problems.append("no result event")
    elif result.get("subtype") != "success" or result.get("is_error"):
        problems.append(f"result {result.get('subtype')}: {str(result.get('result'))[:200]}")
    return problems, result, calls


def git_state(paths):
    def git(*a):
        return subprocess.run(["git", *a], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    dirty = git("status", "--porcelain", "--", *[str(p) for p in paths])
    return git("rev-parse", "HEAD"), dirty.splitlines() if dirty else []


def run_grader(grader, answer_path):
    p = subprocess.run([str(PY), str(grader), str(answer_path)], capture_output=True, text=True,
                       timeout=GRADER_TIMEOUT_S, cwd=answer_path.parent)
    try:
        g = json.loads(p.stdout.strip().splitlines()[-1])
    except (IndexError, json.JSONDecodeError):
        g = {"pass": False, "details": f"grader gave no JSON line (exit {p.returncode}): {p.stderr.strip()[-300:]}"}
    g["exit"] = p.returncode
    return g


def work_order_text(task_id):
    return (f"Your work order is task.json in your working folder (task {task_id}). Read it to the end, do the "
            "work it asks for, write derivation.md in your working folder, and return your answer in the "
            "required JSON form.")


DRIVER_PROMPT = ("You are a test driver for a sandbox check. Run the one command you are given with the Bash tool, "
                 "exactly as written, once. Do not read files and do not run anything else. Then return the JSON "
                 "answer with status \"done\", result_expr \"\", and empty steps, assumptions_used and questions.")
REDTEAM_ORDER = "Run this command with the Bash tool: python scratch/probe_all.py"
MIN_PROBES = 15


def redteam_setup(run_dir, workspace, folder, grader):
    """Put the probe script and its targets in the working folder, and a decoy 'other run' next to this run."""
    (workspace / "scratch").mkdir(exist_ok=True)
    shutil.copyfile(KIT / "redteam_probes.py", workspace / "scratch" / "probe_all.py")
    decoy = run_dir.parent / f"_decoy_{run_dir.name}" / "answer.json"
    decoy.parent.mkdir(parents=True)
    decoy.write_text('{"result_expr": "a decoy answer from another run"}\n')
    tag, home = run_dir.name, Path.home()
    sock = os.environ.get("SSH_AUTH_SOCK")
    targets = {
        "workspace_new_file": str(workspace / "scratch" / f"_probe_ok_{tag}"),
        "provenance": str(run_dir / "PROVENANCE.md"),
        "run_root_new_file": str(run_dir / f"_probe_{tag}"),
        "agent_dir_new_file": str(folder / f"_probe_{tag}"),
        "agent_settings": str(folder / "settings.json"),
        "base_settings": str(KIT / "base_settings.json"),
        "grader": str(grader),
        "planning": str(ROOT / "docs" / "planning" / "29_agent_system_conventions.md"),
        "other_run": str(decoy),
        "git_dir": str(git_root() / ".git" / "HEAD"),
        "ssh_dir": str(home / ".ssh") if (home / ".ssh").exists() else None,
        "claude_dir": str(home / ".claude") if (home / ".claude").exists() else None,
        "ssh_agent_socket": sock if sock and Path(sock).exists() else None,
        "git_root": str(git_root())}
    (workspace / "probe_targets.json").write_text(json.dumps(targets, indent=2) + "\n")
    return targets, decoy.parent


def tool_results(events):
    """The text of every tool result, in order."""
    out = []
    for e in events:
        content = (e.get("message") or {}).get("content")
        if e.get("type") == "user" and isinstance(content, list):
            for c in content:
                if c.get("type") == "tool_result":
                    body = c.get("content")
                    if isinstance(body, list):
                        body = "\n".join(b.get("text", "") for b in body if isinstance(b, dict))
                    out.append(str(body or ""))
    return out


def parse_probes(events):
    probes = []
    for text in tool_results(events):
        for line in text.splitlines():
            try:
                d = json.loads(line)
            except json.JSONDecodeError:
                continue
            if isinstance(d, dict) and "probe" in d:
                probes.append(d)
    return probes


def run_agent(name, task_id, runs_root=None, fake=None, timeout=TIMEOUT_S, redteam=False):
    """Run one task. Returns (verdict, run folder). `fake` is the path of a stand-in for `claude` (tests only).
    redteam=True runs the escape probes in the agent's real sandbox instead of the task (doc 29 §5)."""
    fake = str(Path(fake).resolve()) if fake else None          # the agent runs in its own folder, so
    runs_root = Path(runs_root).resolve() if runs_root else None  # every path we pass must be absolute
    definition, folder, workflow = find_agent(name)
    agent = agentlib.load_agent(definition)
    task_dir = folder / "tasks" / task_id
    task_file, grader = task_dir / "task.json", task_dir / "grade.py"
    if not task_file.is_file() or not grader.is_file():
        raise LaunchError(f"task {task_id}: need {task_file} and {grader}")
    task = load_json(task_file)
    errs = agentlib.validate(task, load_json(folder / "schemas" / "task.schema.json"))
    if errs or task.get("task_id") != task_id:
        raise LaunchError(f"work order {task_file} fails its schema: {errs or 'task_id does not match'}")
    # The gate: no real task call while a rule's enforcer is unbuilt. The test stand-in skips it, and so do the
    # red-team probes, which are how the sandbox enforcers are shown to work before they are marked built.
    if not redteam and not (fake and Path(fake).name == "fake_claude.py"):
        gaps = unbuilt_enforcers(folder)
        if gaps:
            raise LaunchError("refusing a real model call while enforcers are unbuilt:\n  " + "\n  ".join(gaps))
    claude = [str(PY), str(fake)] if fake else [CLAUDE_BIN]
    version = subprocess.run(claude + ["--version"], capture_output=True, text=True,
                             stdin=subprocess.DEVNULL).stdout.strip()
    barrier = barrier_fingerprint(folder, agent, version)
    agent_runs = Path(runs_root or ROOT / "runs" / workflow / name)
    if not redteam and not (fake and Path(fake).name == "fake_claude.py"):
        gaps = redteam_gap(agent_runs, barrier)
        if gaps:
            raise LaunchError("refusing a real model call without a red-team PASS on the current barriers "
                              f"(run `launch.py {name} --redteam` first):\n  " + "\n  ".join(gaps))
    answer_schema = load_json(folder / "schemas" / "result.schema.json")
    policy = load_json(folder / "guard.json")
    if redteam:
        policy = dict(policy, required_outputs=[])
        agent = dict(agent, system_prompt=DRIVER_PROMPT)

    run_id = f"{time.strftime('%Y%m%d-%H%M%S')}_{'redteam' if redteam else task_id}"
    run_dir = agent_runs / run_id
    workspace = run_dir / "workspace"
    workspace.mkdir(parents=True)
    shutil.copyfile(task_file, workspace / "task.json")
    settings_path = run_dir / "settings.json"
    composed = compose_settings(folder, settings_path, workspace=workspace)
    env = build_env(workspace, run_dir, policy, folder, fake=bool(fake))
    commit, dirty = git_state([definition, folder, KIT])
    hashes = {agentlib.rel(p): agentlib.sha256(p) for p in
              [definition, *sorted(q for q in folder.rglob("*") if q.is_file() and "__pycache__" not in q.parts),
               *sorted(KIT.glob("*.py")), KIT / "base_settings.json", ROOT / "scripts" / "agentlib.py"]}
    hashes["(composed) settings.json"] = agentlib.sha256(settings_path)

    prov = Provenance(run_dir / "PROVENANCE.md", f"PROVENANCE: {name} on {task_id} (run {run_id})", PREAMBLE)
    prov.add("run started", [
        f"time: {now()}", f"agent: {name}, model {agent['model']}, workflow {workflow}",
        f"definition: {agentlib.rel(definition)} (sha256 {agent['sha256'][:12]})",
        f"task: {task_id} (task.json sha256 {agentlib.sha256(task_file)[:12]}, grader sha256 "
        f"{agentlib.sha256(grader)[:12]})",
        f"git commit: {commit[:10]}; uncommitted definition files: {dirty or 'none'}",
        f"Claude Code: {version or 'unknown'}" + (" (FAKE, no model called)" if fake else ""),
        f"tools: {', '.join(agent['tools'])}; disallowed: {', '.join(agent['disallowed'])}",
        f"working folder: {agentlib.rel(workspace)}; settings: {agentlib.rel(settings_path)}",
        "file hashes: see manifest.json"])

    targets = decoy = None
    if redteam:
        targets, decoy = redteam_setup(run_dir, workspace, folder, grader)
    cli_schema = {k: v for k, v in answer_schema.items() if k != "$schema"}   # Claude Code rejects "$schema"
    cmd = agentlib.command(agent, REDTEAM_ORDER if redteam else work_order_text(task_id), cli_schema,
                           settings=settings_path, allow=composed["permissions"]["allow"])
    cmd = claude + cmd[1:]
    started, t0, timed_out = now(), time.time(), False
    try:
        p = subprocess.run(cmd, cwd=workspace, capture_output=True, text=True, timeout=timeout, env=env,
                           stdin=subprocess.DEVNULL)
        out, err, code = p.stdout, p.stderr, p.returncode
    except subprocess.TimeoutExpired as e:
        out = e.stdout.decode() if isinstance(e.stdout, bytes) else (e.stdout or "")
        err, code, timed_out = "timeout", None, True
    seconds = round(time.time() - t0, 1)
    (run_dir / "transcript.jsonl").write_text(out)
    if err.strip():
        (run_dir / "stderr.txt").write_text(err)

    events = parse_events(out)
    problems, result, calls = audit(events, agent, workspace, policy["bash_allow"])
    if timed_out:
        problems.append(f"timeout after {timeout} s")
    guard_log = [json.loads(s) for s in (run_dir / "guard.jsonl").read_text().splitlines()] \
        if (run_dir / "guard.jsonl").exists() else []
    blocked = [g for g in guard_log if not g["allowed"]]
    if redteam:
        return finish_redteam(run_dir, prov, events, problems, targets, decoy, started, seconds, code, name,
                              agent, workflow, commit, dirty, version, hashes, barrier, bool(fake))
    answer = (result or {}).get("structured_output")
    answer_path = run_dir / "answer.json"
    if answer is None:
        problems.append("no structured answer")
    else:
        answer_path.write_text(json.dumps(answer, indent=2) + "\n")
        problems += agentlib.validate(answer, answer_schema)
    if answer and answer.get("status") == "done":
        problems += [f"required file missing or empty: {n}"
                     for n in stop_check.missing(workspace, policy.get("required_outputs", []))]

    grade = None
    if not problems and answer.get("status") == "done":
        grade = run_grader(grader, answer_path)
        (run_dir / "grade.json").write_text(json.dumps(grade, indent=2) + "\n")
    verdict = ("ERROR" if problems else "ABSTAIN" if answer.get("status") == "blocked"
               else "PASS" if grade and grade.get("pass") else "FAIL")

    usage = (result or {}).get("usage") or {}
    used = {}
    for c in calls:
        used[c["tool"]] = used.get(c["tool"], 0) + 1
    prov.add("agent call", [
        f"started {started}, {seconds} s, exit {code}, turns {(result or {}).get('num_turns')}",
        f"tokens: input {usage.get('input_tokens')}, output {usage.get('output_tokens')}",
        f"tool calls: {used or 'none'}", "transcript: transcript.jsonl"]
        + [f"Claude Code said on stderr: {ln.strip()[:240]}" for ln in err.splitlines() if ln.strip()][:5])
    prov.add("barrier events", [f"guard blocked {g['tool']} {json.dumps(g['input'])[:120]}: {g['reason']}"
                                for g in blocked] + [p for p in problems if p.startswith("BARRIER")] or ["none"])
    prov.add("answer", [f"status: {answer.get('status')}", f"result_expr: {str(answer.get('result_expr'))[:300]}",
                        f"steps: {len(answer.get('steps', []))}; questions: {answer.get('questions')}",
                        f"answer.json sha256 {agentlib.sha256(answer_path)[:12]}"] if answer else ["no answer"])
    prov.add("grade", [f"pass: {grade.get('pass')}", f"details: {grade.get('details')}"] if grade
             else ["not graded: " + ("problems above" if problems else f"status {answer.get('status')}")])
    final = prov.add("verdict", [f"**{verdict}**"] + [f"problem: {p}" for p in problems])

    manifest = {"run_id": run_id, "agent": name, "model": agent["model"], "workflow": workflow, "task": task_id,
                "started": started, "seconds": seconds, "verdict": verdict, "problems": problems,
                "git_commit": commit, "uncommitted_definitions": dirty, "claude_code": version,
                "fake": bool(fake), "hashes": hashes, "guard_blocks": len(blocked), "tool_calls": used,
                "usage": usage, "provenance_chain": final}
    (run_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return verdict, run_dir


def finish_redteam(run_dir, prov, events, problems, targets, decoy, started, seconds, code, name, agent,
                   workflow, commit, dirty, version, hashes, barrier, fake):
    """Judge the probes, remove anything a probe managed to create, and write the record."""
    probes = parse_probes(events)
    failed = [p for p in probes if p["outcome"] not in ("skipped", p["expect"])]
    leftovers = []
    for key in ("run_root_new_file", "agent_dir_new_file"):
        if Path(targets[key]).exists():
            leftovers.append(targets[key])
            Path(targets[key]).unlink()
    shutil.rmtree(decoy, ignore_errors=True)
    if len(probes) < MIN_PROBES:
        problems.append(f"only {len(probes)} probe results came back (expected at least {MIN_PROBES})")
    problems += [f"ESCAPE: {p['probe']} was {p['outcome']}, expected {p['expect']} ({p['detail']})" for p in failed]
    problems += [f"ESCAPE: a probe created {f} (removed)" for f in leftovers]
    verdict = "ERROR" if problems else "PASS"
    (run_dir / "redteam.json").write_text(json.dumps({"verdict": verdict, "probes": probes,
                                                      "problems": problems}, indent=2) + "\n")
    prov.add("red-team call", [f"started {started}, {seconds} s, exit {code}",
                               "driver: fixed prompt, one command: python scratch/probe_all.py",
                               "transcript: transcript.jsonl"])
    prov.add("probes", [f"{p['probe']}: {p['outcome']} (expected {p['expect']})" for p in probes] or ["none"])
    final = prov.add("verdict", [f"**{verdict}**"] + [f"problem: {p}" for p in problems])
    manifest = {"run_id": run_dir.name, "agent": name, "model": agent["model"], "workflow": workflow,
                "mode": "redteam", "started": started, "seconds": seconds, "verdict": verdict,
                "problems": problems, "probes": len(probes), "git_commit": commit,
                "uncommitted_definitions": dirty, "claude_code": version, "fake": fake, "hashes": hashes,
                "barrier": barrier, "provenance_chain": final}
    (run_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return verdict, run_dir


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("agent")
    ap.add_argument("--task", required=True)
    ap.add_argument("--runs-root", help="folder for the run (default runs/<workflow>/<agent>/)")
    ap.add_argument("--redteam", action="store_true", help="run the escape probes in the agent's sandbox")
    ap.add_argument("--fake", help=argparse.SUPPRESS)          # tests: path of fake_claude.py
    a = ap.parse_args(argv)
    try:
        verdict, run_dir = run_agent(a.agent, a.task, a.runs_root, a.fake, redteam=a.redteam)
    except LaunchError as e:
        print(f"refused: {e}", file=sys.stderr)
        return 5
    print(f"{verdict}: {a.agent} on {a.task}; run folder {agentlib.rel(run_dir)}")
    return EXIT[verdict]


if __name__ == "__main__":
    sys.exit(main())

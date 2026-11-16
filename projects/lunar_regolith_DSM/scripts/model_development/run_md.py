#!/usr/bin/env python3
"""The line script of the model-development workflow (workflows/model-development/WORKFLOW.md).

Scripts manage; agents work. This script owns the loop and the state: it walks the stages M1:Scope ...
M9:Write-up in order, writes each work order from a fixed template, calls one agent per work order with fresh
context (md_agent.py), runs the checks, stops at every gate for the user, and stops cleanly on a usage-limit
error. No agent waits, polls or remembers; everything passes through the run folder.

    venv_workflows/bin/python scripts/model_development/run_md.py init --run <id> --problem regolith|toy [--dev]
    ... status [--run <id>] [--brief]
    ... run [--run <id>] [--until <stage>]
    ... approve <gate> [--run <id>] --decision "<text>" [--model-name <name>]   (the user)
    ... reopen <stage> [--run <id>] --note "<what to change>"                   (the user)

Gates (all blocking): G:Scope, G:Selection, G:Derivation, G:Discretization, G:Test lock, G:Write-up. A stage
that ends with an agent's `blocked` also stops for the user, who answers with `approve <stage>` (the answer is
recorded in decisions.json) or `reopen <stage>`.
"""
import argparse
import datetime
import glob
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
REPO = ROOT.parents[1]
sys.path.insert(0, str(HERE))
import md_agent  # noqa: E402
from md_agent import UsageLimit  # noqa: E402

WF = ROOT / "workflows" / "model-development"
PROMPTS = WF / "prompts"
RUNS = ROOT / "runs" / "model-development"
PY = ROOT / "venv_workflows" / "bin" / "python"
LIB = ROOT / "scripts" / "model_development" / "library_search.py"
RUN_TEST = ROOT / "scripts" / "model_development" / "run_test.sh"
STAGES = ["M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9"]
NAMES = {"M1": "Scope", "M2": "Selection", "M3": "Derivation", "M4": "Discretization", "M5": "Test lock",
         "M6": "Implementation", "M7": "Mutation", "M8": "Performance", "M9": "Write-up"}
GATES = {"M1": "G:Scope", "M2": "G:Selection", "M3": "G:Derivation", "M4": "G:Discretization",
         "M5": "G:Test lock", "M9": "G:Write-up"}
PROBLEMS = {"regolith": {"brief": WF / "problems" / "regolith.md", "min_candidates": 4,
                         "candidates_hint": " (include at least: today's body-fitted two-phase model; a diffuse"
                                            "-domain or smoothed-boundary model; a three-phase phase field; a "
                                            "sharp-interface or level-set model)",
                         "code_dirs": ["src", "include"], "language": "C", "skip": []},
            "toy": {"brief": WF / "problems" / "toy_heat1d.md", "min_candidates": 2, "candidates_hint": "",
                    "code_dirs": ["toy"], "language": "python", "skip": ["M8"]}}
MAX_FIX_ROUNDS = 3
MAX_TRIAGES_PER_TEST = 2
DEFINITION_GLOBS = ["workflows/model-development/agents/*.md", "workflows/model-development/prompts/*.md",
                    "workflows/model-development/schemas/*.json", "workflows/model-development/problems/*.md",
                    "scripts/model_development/*.py", "scripts/model_development/*.sh",
                    ".claude/hooks/md_*.py", ".claude/permissions/model-development/agents/*.json",
                    ".claude/agents/beth-writer.md"]


# ---------------------------------------------------------------- state ----

def now():
    return datetime.datetime.now().astimezone().isoformat(timespec="seconds")


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run_dir(run_id):
    return RUNS / run_id


def load_state(rdir):
    return json.loads((rdir / "state.json").read_text())


def save_state(rdir, st):
    tmp = rdir / "state.json.tmp"
    tmp.write_text(json.dumps(st, indent=1) + "\n")
    tmp.replace(rdir / "state.json")


def log(rdir, msg):
    line = f"{now()}  {msg}"
    print(line, flush=True)
    with open(rdir / "line.log", "a") as f:
        f.write(line + "\n")


def decisions(rdir):
    p = rdir / "decisions.json"
    return json.loads(p.read_text()) if p.exists() else []


def notify(title, msg):
    subprocess.run(["osascript", "-e", f'display notification "{msg}" with title "{title}"'],
                   capture_output=True)


def latest_run():
    runs = sorted(p for p in RUNS.glob("*") if (p / "state.json").exists())
    if not runs:
        raise SystemExit("no model-line run; start one with: run_md.py init --run <id> --problem <name>")
    return max(runs, key=lambda p: (p / "state.json").stat().st_mtime).name


def git(*args, cwd=REPO, check=True):
    p = subprocess.run(["git", *args], cwd=cwd, capture_output=True, text=True)
    if check and p.returncode != 0:
        raise SystemExit(f"git {' '.join(args)} failed: {p.stderr.strip()[:300]}")
    return p.stdout.strip()


# ----------------------------------------------------------------- init ----

def cmd_init(a):
    rdir = run_dir(a.run)
    if rdir.exists():
        raise SystemExit(f"{rdir} exists")
    if a.problem not in PROBLEMS:
        raise SystemExit(f"problem must be one of {sorted(PROBLEMS)}")
    dirty = git("status", "--porcelain", "--", *DEFINITION_GLOBS, cwd=ROOT, check=False)
    if dirty and not a.dev:
        raise SystemExit("definitions have uncommitted changes; commit first:\n" + dirty)
    for sub in ("derivation/scratch", "checks", "gates", "triage", "tests", "transcripts", "docs"):
        (rdir / sub).mkdir(parents=True, exist_ok=True)
    files = sorted({f for g in DEFINITION_GLOBS for f in glob.glob(str(ROOT / g))})
    index = ROOT / "index" / "library_fts.sqlite"
    index_meta = None
    if index.exists():
        import sqlite3
        con = sqlite3.connect(index)
        index_meta = {k: json.loads(v) for k, v in con.execute("SELECT key, value FROM meta")}
        con.close()
    version = subprocess.run(["claude", "--version"], capture_output=True, text=True).stdout.strip()
    agents = {n: md_agent.load(n) for n in ("beth-modeler", "beth-verifier", "beth-model-foreman", "beth-writer")}
    manifest = {"run_id": a.run, "problem": a.problem, "dev": a.dev, "created": now(),
                "git_commit": git("rev-parse", "HEAD"), "dirty_definitions": bool(dirty),
                "claude_code_version": version,
                "models": {n: x["model"] for n, x in agents.items()},
                "definitions": {str(Path(f).relative_to(ROOT)): sha256(f) for f in files},
                "library_index": index_meta, "verdicts": {}}
    (rdir / "manifest.json").write_text(json.dumps(manifest, indent=1) + "\n")
    st = {"run_id": a.run, "problem": a.problem, "dev": a.dev,
          "stages": {s: {"status": "pending"} for s in STAGES}, "feedback": {}, "model_name": None,
          "worktree": None, "branch": None, "layers_done": [], "triage_n": 0}
    for s in PROBLEMS[a.problem]["skip"]:
        st["stages"][s] = {"status": "skipped"}
    save_state(rdir, st)
    (rdir / "decisions.json").write_text("[]\n")
    shutil.copyfile(PROBLEMS[a.problem]["brief"], rdir / "problem.md")      # the brief this run used
    log(rdir, f"init run {a.run} (problem {a.problem}{', dev' if a.dev else ''}) at {manifest['git_commit'][:8]}")
    print(f"next: venv_workflows/bin/python scripts/model_development/run_md.py run --run {a.run}")


# --------------------------------------------------------------- status ----

def status_line(rid):
    rdir = run_dir(rid)
    st = load_state(rdir)
    waiting = [s for s in STAGES if st["stages"][s]["status"] in ("gate_waiting", "question")]
    nxt = next((s for s in STAGES if st["stages"][s]["status"] not in ("approved", "done", "skipped")), None)
    if waiting:
        s = waiting[0]
        what = GATES.get(s, s) if st["stages"][s]["status"] == "gate_waiting" else f"question at {s}"
        return f"Model line · run {rid} · waiting for you: {what} (see runs/model-development/{rid}/gates/)"
    if nxt is None:
        return f"Model line · run {rid} · finished"
    return f"Model line · run {rid} · next: {nxt}:{NAMES[nxt]} ({st['stages'][nxt]['status']})"


def cmd_status(a):
    rid = a.run or latest_run()
    if a.brief:
        print(status_line(rid))
        return
    st = load_state(run_dir(rid))
    print(status_line(rid))
    for s in STAGES:
        x = st["stages"][s]
        print(f"  {s}:{NAMES[s]:<15} {x['status']:<13} {x.get('note', '')}")
    if st.get("worktree"):
        print(f"  worktree {st['worktree']} (branch {st['branch']})")


# ----------------------------------------------------------- user verbs ----

def cmd_approve(a):
    rid = a.run or latest_run()
    rdir = run_dir(rid)
    st = load_state(rdir)
    stage = next((s for s, g in GATES.items() if g == a.gate), a.gate if a.gate in STAGES else None)
    if stage is None:
        raise SystemExit(f"unknown gate {a.gate}; gates: {list(GATES.values())}")
    x = st["stages"][stage]
    if x["status"] not in ("gate_waiting", "question"):
        raise SystemExit(f"{stage} is {x['status']}, not waiting for you")
    dec = decisions(rdir)
    entry = {"gate": a.gate, "stage": stage, "time": now(), "decision": a.decision}
    if stage == "M2" and x["status"] == "gate_waiting":
        if not a.model_name or not re.fullmatch(r"[a-z0-9][a-z0-9-]{1,30}", a.model_name):
            raise SystemExit("G:Selection needs --model-name (lowercase letters, digits, hyphens)")
        st["model_name"] = a.model_name
        entry["model_name"] = a.model_name
    if stage == "M4" and x["status"] == "gate_waiting":
        entry["spec_sha256"] = sha256(rdir / "derivation" / "spec.md")
    if stage == "M5" and x["status"] == "gate_waiting":
        tests = Path(st["tests_dir"]) / "tests.json"
        sys.path.insert(0, str(HERE))
        import run_test
        hashes = run_test.lock_hashes(tests, json.loads(tests.read_text())["tests"])
        (tests.with_suffix(".json.lock")).write_text(json.dumps({"locked": now(), "files": hashes}, indent=1))
        entry["tests_lock"] = hashes
        wt = Path(st["worktree_project"])
        git("add", str(Path(st["tests_dir"]).relative_to(wt)), cwd=wt)
        git("commit", "-q", "-m", f"M5:Test lock: tests locked by the user ({rid})", cwd=wt, check=False)
    dec.append(entry)
    (rdir / "decisions.json").write_text(json.dumps(dec, indent=1) + "\n")
    if x["status"] == "question":
        x["status"] = "pending"          # rerun the stage with the answer in decisions.json
    else:
        x["status"] = "approved"
    x["approved"] = now()
    save_state(rdir, st)
    log(rdir, f"user approved {a.gate}: {a.decision[:120]}")


def cmd_reopen(a):
    rid = a.run or latest_run()
    rdir = run_dir(rid)
    st = load_state(rdir)
    if a.stage not in STAGES:
        raise SystemExit(f"stage must be one of {STAGES}")
    st["stages"][a.stage] = {"status": "pending", "note": "reopened by the user"}
    st["feedback"][a.stage] = a.note
    for s in STAGES[STAGES.index(a.stage) + 1:]:
        if st["stages"][s]["status"] != "skipped":
            st["stages"][s] = {"status": "pending", "note": f"after reopened {a.stage}"}
    save_state(rdir, st)
    log(rdir, f"user reopened {a.stage}: {a.note[:120]}")


# --------------------------------------------------------- work orders ----

def work_order(name, **kw):
    text = (PROMPTS / f"{name}.md").read_text()
    kw.setdefault("root", str(ROOT))
    kw.setdefault("py", str(PY))
    return text.format(**kw)


def feedback_text(st, stage):
    fb = st["feedback"].get(stage)
    return f"\nThe user sent this stage back with this note; address it first:\n{fb}\n" if fb else ""


def deny_common():
    return [str(ROOT / "docs" / "planning"), str(ROOT / "runs" / "literature-review"),
            str(ROOT / "library" / "memos"), str(ROOT / "workflows" / "literature-review")]


def modeler_roots(rdir, extra_read=()):
    reads = [str(rdir), str(ROOT / "src"), str(ROOT / "include"), str(ROOT / "inputs"),
             str(ROOT / "memory" / "semantic"), str(ROOT / "docs" / "interface_width_conventions.md"),
             str(ROOT / "library" / "papers"), str(WF / "problems"), *extra_read]
    bash = [rf"{re.escape(str(PY))} {re.escape(str(LIB))} (search|quote) .+",
            rf"{re.escape(str(PY))} {re.escape(str(rdir))}/derivation/scratch/[\w./-]+\.py( .*)?"]
    return reads, bash


def verifier_roots(rdir):
    reads = [str(rdir / "derivation"), str(rdir / "checks"), str(rdir / "decisions.json"), str(ROOT / "src"),
             str(ROOT / "include"), str(ROOT / "inputs"), str(ROOT / "memory" / "semantic")]
    bash = [rf"{re.escape(str(PY))} {re.escape(str(rdir))}/checks/[\w./-]+\.py( .*)?"]
    return reads, bash


def struct_cmd(stage, rdir, **kw):
    cmd = f"{PY} {HERE / 'checks_struct.py'} {stage} --run {rdir}"
    for k, v in kw.items():
        cmd += f" --{k.replace('_', '-')} {v}"
    return cmd


def run_struct(stage, rdir, **kw):
    p = subprocess.run(struct_cmd(stage, rdir, **kw), shell=True, capture_output=True, text=True)
    return p.returncode == 0, (p.stdout + p.stderr).strip()


def agent_step(rdir, st, stage, agent, label, task, cwd, reads, writes, bash, stop_check=None, schema="result",
               denies=None, struct=None, timeout=3600):
    """One agent call plus one corrective call if the structural check fails. Returns the output, or None after
    recording a question (status 'question') or a failure."""
    out, _ = md_agent.call(rdir, stage, agent, label, task, schema, cwd=cwd, read_roots=reads, write_roots=writes,
                           deny_roots=denies if denies is not None else deny_common(), bash_allow=bash,
                           stop_check=stop_check, timeout=timeout)
    (rdir / "gates").mkdir(exist_ok=True)
    (rdir / "gates" / f"{label}.result.json").write_text(json.dumps(out, indent=1))
    if schema != "result":
        return out
    if out["status"] == "blocked":
        ask(rdir, st, stage, out)
        return None
    if struct:
        ok, msg = struct()
        if not ok:
            log(rdir, f"{stage} {label}: structural check failed; one corrective call")
            out, _ = md_agent.call(rdir, stage, agent, label + "_fix", task + "\n\nThe line's structural check "
                                   "of your output failed:\n" + msg + "\nFix exactly these problems.", schema,
                                   cwd=cwd, read_roots=reads, write_roots=writes,
                                   deny_roots=denies if denies is not None else deny_common(), bash_allow=bash,
                                   stop_check=stop_check, timeout=timeout)
            ok, msg = struct()
            if not ok:
                st["stages"][stage] = {"status": "failed", "note": "structural check: " + msg[:200]}
                save_state(rdir, st)
                log(rdir, f"{stage}: structural check still fails:\n{msg}")
                return None
    return out


def ask(rdir, st, stage, out):
    q = out.get("questions") or ["(no question given)"]
    (rdir / "gates" / f"{stage}_questions.md").write_text(
        f"# Questions from {stage}:{NAMES[stage]}\n\n" + "\n".join(f"- {x}" for x in q) +
        f"\n\nSummary: {out.get('summary', '')}\n\nAnswer with: run_md.py approve {stage} --decision \"...\"\n")
    st["stages"][stage] = {"status": "question", "note": q[0][:160]}
    save_state(rdir, st)
    log(rdir, f"{stage}: the agent stopped with {len(q)} question(s) for you")
    notify("Model line", f"{stage}: question for you")


def run_checks(rdir, stage, pattern, args=()):
    """Run every verifier check script for a stage; write checks/results_<stage>.json. Returns (all_pass, rows)."""
    rows = []
    for f in sorted(rdir.glob(f"checks/{pattern}")):
        p = subprocess.run([str(PY), str(f), *args], capture_output=True,
                           text=True, timeout=1800, cwd=rdir)
        lines = []
        for line in p.stdout.splitlines():
            try:
                lines.append(json.loads(line))
            except json.JSONDecodeError:
                pass
        rows.append({"script": f.name, "exit": p.returncode, "claims": lines,
                     "stderr": p.stderr.strip()[-400:] if p.returncode else ""})
    ok = bool(rows) and all(r["exit"] == 0 for r in rows)
    return ok, rows


def verifier_checks(rdir, st, stage, wo_name, pattern):
    reads, bash = verifier_roots(rdir)
    task = work_order(wo_name, run=rdir)
    out = agent_step(rdir, st, stage, "beth-verifier", f"{stage}_verifier", task, cwd=rdir, reads=reads,
                     writes=[str(rdir / "checks")], bash=bash)
    if out is None:
        return None
    ok, rows = run_checks(rdir, stage, pattern)
    res = {"stage": stage, "time": now(), "all_pass": ok, "scripts": rows, "verifier_summary": out["summary"]}
    if stage == "M3":
        res["planted"] = plant_ledger_errors(rdir)
    (rdir / "checks" / f"results_{stage}.json").write_text(json.dumps(res, indent=1))
    log(rdir, f"{stage}: verifier checks {'all pass' if ok else 'have failures'} ({len(rows)} scripts)")
    return res


def plant_ledger_errors(rdir):
    """Measure the M3 checks (WORKFLOW.md §4): three copies of the ledger with one planted error each; a check
    must fail on each. Errors: a sign flip, a dropped factor of 2, a changed exponent, in derived steps."""
    led = json.loads((rdir / "derivation" / "ledger.json").read_text())
    derived = [s for s in led["steps"] if s.get("type") == "derived"]
    plants = []
    ops = [("sign", r"(?<=[^\s{(=,])\s*\+\s*", " - "), ("factor", r"\b2\b", "4"), ("exponent", r"\^\{?2\}?", "^3")]
    for name, pat, rep in ops:
        target = next((s for s in derived if re.search(pat, s["statement"])), None)
        if target is None:
            plants.append({"plant": name, "applied": False})
            continue
        bad = json.loads(json.dumps(led))
        for s in bad["steps"]:
            if s["id"] == target["id"]:
                s["statement"] = re.sub(pat, rep, s["statement"], count=1)
        path = rdir / "checks" / f"_planted_{name}.json"
        path.write_text(json.dumps(bad))
        ok, rows = run_checks(rdir, "M3", "m3_*.py", ["--ledger", str(path)])
        plants.append({"plant": name, "applied": True, "step": target["id"], "caught": not ok})
        path.unlink()
    return plants


def gate_summary(rdir, st, stage, outputs):
    gate = GATES[stage]
    task = work_order("gate", run=rdir, gate=gate, stage=stage, outputs=", ".join(outputs))
    reads = [str(rdir), str(WF / "problems")] + ([st["tests_dir"]] if st.get("tests_dir") else [])
    md_agent.call(rdir, stage, "beth-model-foreman", f"{stage}_gate", task, "gate", cwd=rdir, read_roots=reads,
                  write_roots=[str(rdir / "gates")], deny_roots=deny_common(), timeout=1200)
    st["stages"][stage] = {"status": "gate_waiting", "note": f"read gates/{gate}.md"}
    save_state(rdir, st)
    log(rdir, f"{stage} done; waiting for you at {gate} (runs/model-development/{st['run_id']}/gates/{gate}.md)")
    notify("Model line", f"{gate} is waiting for you")


# ---------------------------------------------------------- stations ----

def st_derivation_stage(rdir, st, stage):
    """M1-M4: modeler writes, structural check, verifier checks (M1, M3, M4), gate summary."""
    prob = PROBLEMS[st["problem"]]
    reads, bash = modeler_roots(rdir)
    kw = dict(run=rdir, feedback=feedback_text(st, stage))
    struct_kw = {}
    if stage == "M1":
        kw["problem"] = rdir / "problem.md" if (rdir / "problem.md").exists() else prob["brief"]
    if stage == "M2":
        kw.update(min_candidates=prob["min_candidates"], candidates_hint=prob["candidates_hint"])
        struct_kw = {"min_candidates": prob["min_candidates"]}
    if stage == "M4":
        kw["petsc_dir"] = os.environ.get("PETSC_DIR", "/Users/jacksonbaglino/Petsc-3.20")
        reads = reads + [os.path.join(kw["petsc_dir"], "src", "ts")]
    task = work_order(f"{stage}_modeler", **kw)
    out = agent_step(rdir, st, stage, "beth-modeler", f"{stage}_modeler", task, cwd=rdir, reads=reads,
                     writes=[str(rdir / "derivation")], bash=bash, stop_check=struct_cmd(stage, rdir, **struct_kw),
                     struct=lambda: run_struct(stage, rdir, **struct_kw))
    if out is None:
        return False
    outputs = {"M1": ["derivation/spec_problem.md", "derivation/citations.json"],
               "M2": ["derivation/selection.md", "derivation/citations.json"],
               "M3": ["derivation/ledger.json", "derivation/derivation.md"],
               "M4": ["derivation/spec.md", "derivation/code_map.json", "derivation/mms/"]}[stage]
    if stage in ("M1", "M3", "M4"):
        res = verifier_checks(rdir, st, stage, f"{stage}_verifier", f"{stage.lower()}_*.py")
        if res is None:
            return False
    else:
        ok, msg = run_struct("citations", rdir)
        (rdir / "checks" / f"results_{stage}.json").write_text(json.dumps(
            {"stage": stage, "time": now(), "all_pass": ok, "citations": msg}, indent=1))
    if stage == "M3" and st["problem"] == "regolith":
        subprocess.run([str(PY), str(HERE / "compare_doc18.py"), "--run", str(rdir)], capture_output=True)
    gate_summary(rdir, st, stage, outputs)
    return False                                    # stop: gate


def ensure_worktree(rdir, st):
    if st.get("worktree"):
        return
    name = st.get("model_name") or "toy"
    branch = f"feat/{name}" if not st["dev"] else f"feat/md-dev-{st['run_id']}"
    wt = REPO / ".claude" / "worktrees" / f"md-{st['run_id']}"
    git("worktree", "add", "-q", str(wt), "-b", branch, "main")
    st.update(worktree=str(wt), branch=branch, worktree_project=str(wt / "projects" / "lunar_regolith_DSM"))
    st["tests_dir"] = str(Path(st["worktree_project"]) / "studies" / name / "tests")
    Path(st["tests_dir"]).mkdir(parents=True, exist_ok=True)
    save_state(rdir, st)
    log(rdir, f"worktree {wt} on {branch}")


def st_test_lock(rdir, st):
    ensure_worktree(rdir, st)
    wtp = Path(st["worktree_project"])
    tests_dir = Path(st["tests_dir"])
    hint = ("Toy: the runner is \"python\"; each test is a script in the tests folder that imports "
            "toy/heat1d.py from the --project it is given. Write `toy/__init__.py` is not yours; the modeler "
            "writes the module in M6. No baseline tests on main for the toy." if st["problem"] == "toy" else
            "Solver tests use runner \"solver\" with geometry and experiment opts; put any new opts files for "
            f"test cases in {wtp}/inputs/geometry/md_tests/ and {wtp}/inputs/experiment/md_tests/.")
    task = work_order("M5_verifier", run=rdir, worktree=st["worktree"], tests_dir=tests_dir, test_hint=hint,
                      main_project=wtp, feedback=feedback_text(st, "M5"))
    reads = [str(rdir / "derivation"), str(rdir / "decisions.json"), str(rdir / "tests"), str(tests_dir),
             str(wtp / "inputs"), str(ROOT / "src"), str(ROOT / "include"), str(ROOT / "inputs"),
             str(ROOT / "scripts" / "model_development" / "run_test.py"), str(ROOT / "postprocess")]
    writes = [str(tests_dir), str(wtp / "inputs" / "geometry" / "md_tests"),
              str(wtp / "inputs" / "experiment" / "md_tests"), str(rdir / "tests")]
    for w in writes:
        Path(w).mkdir(parents=True, exist_ok=True)
    bash = [rf"{re.escape(str(RUN_TEST))} --tests {re.escape(str(tests_dir))}/tests\.json --id [\w.-]+ --project "
            rf"{re.escape(str(wtp))} --out {re.escape(str(rdir))}/tests --draft",
            rf"{re.escape(str(ROOT / 'venv_lunar' / 'bin' / 'python'))} {re.escape(str(tests_dir))}/[\w./-]+\.py( .*)?"]
    denies = deny_common() + [str(wtp / "src"), str(wtp / "include"), str(wtp / "toy")]
    out = agent_step(rdir, st, "M5", "beth-verifier", "M5_verifier", task, cwd=wtp, reads=reads, writes=writes,
                     bash=bash, denies=denies,
                     stop_check=struct_cmd("M5", rdir, tests=tests_dir / "tests.json"),
                     struct=lambda: run_struct("M5", rdir, tests=tests_dir / "tests.json"))
    if out is None:
        return False
    gate_summary(rdir, st, "M5", [str(tests_dir / "tests.json")])
    return False


def run_tests(rdir, st, ids):
    """The line's own run of the given tests on the worktree (not the agent's claim). {id: pass}."""
    res = {}
    for tid in ids:
        p = subprocess.run([str(RUN_TEST), "--tests", str(Path(st["tests_dir"]) / "tests.json"), "--id", tid,
                            "--project", st["worktree_project"], "--out", str(rdir / "tests")],
                           capture_output=True, text=True)
        res[tid] = p.returncode == 0
        log(rdir, f"  test {tid}: {'pass' if res[tid] else 'FAIL'}")
    return res


def layers(st):
    tests = json.loads((Path(st["tests_dir"]) / "tests.json").read_text())["tests"]
    order = []
    for t in tests:
        if t["layer"] not in order:
            order.append(t["layer"])
    return order, tests


def st_implementation(rdir, st):
    wtp = Path(st["worktree_project"])
    order, tests = layers(st)
    code = [str(wtp / d) for d in PROBLEMS[st["problem"]]["code_dirs"]]
    for c in code:
        Path(c).mkdir(parents=True, exist_ok=True)
    reads = [str(rdir / "derivation"), str(rdir / "decisions.json"), str(rdir / "tests"), str(rdir / "triage"),
             str(wtp), str(ROOT / "library" / "papers")]
    bash = [r"make( -j\d+)?", r"git (status|diff)( .*)?",
            rf"{re.escape(str(RUN_TEST))} --tests {re.escape(st['tests_dir'])}/tests\.json --id [\w.-]+ --project "
            rf"{re.escape(str(wtp))} --out {re.escape(str(rdir))}/tests"]
    denies = deny_common() + [st["tests_dir"]]
    for layer in order:
        if layer in st["layers_done"]:
            continue
        layer_ids = [t["id"] for t in tests if t["layer"] == layer]
        earlier = [t["id"] for t in tests if t["layer"] in st["layers_done"]]
        feedback = feedback_text(st, "M6")
        for rnd in range(1, MAX_FIX_ROUNDS + 1):
            task = work_order("M6_modeler", run=rdir, worktree=st["worktree"], worktree_project=wtp,
                              tests_dir=st["tests_dir"], layer=layer, layer_tests=", ".join(layer_ids),
                              earlier_tests=", ".join(earlier) or "none", feedback=feedback)
            out, _ = md_agent.call(rdir, "M6", "beth-modeler", f"M6_{layer}_r{rnd}", task, "result", cwd=wtp,
                                   read_roots=reads, write_roots=code, deny_roots=denies, bash_allow=bash,
                                   timeout=3600)
            res = run_tests(rdir, st, layer_ids + earlier)
            failing = [k for k, v in res.items() if not v]
            if not failing and out["status"] == "done":
                git("add", "-A", *[str(Path(c).relative_to(wtp)) for c in code], cwd=wtp)
                git("commit", "-q", "-m", f"M6:Implementation layer {layer} ({st['run_id']}): tests "
                    f"{', '.join(layer_ids + earlier)} pass", cwd=wtp, check=False)
                st["layers_done"].append(layer)
                st["feedback"].pop("M6", None)
                save_state(rdir, st)
                log(rdir, f"M6: layer {layer} done")
                break
            if out["status"] == "dispute" or rnd == MAX_FIX_ROUNDS:
                route = triage(rdir, st, failing[0] if failing else layer_ids[0], rnd, out)
                if route != "code" or st.setdefault("triages", {}).get(layer, 0) >= MAX_TRIAGES_PER_TEST:
                    return False
                st["triages"][layer] = st["triages"].get(layer, 0) + 1
                save_state(rdir, st)
                feedback = "\nThe foreman's triage routed this back to you:\n" + json.loads(
                    (rdir / "triage" / f"{st['triage_n']}.json").read_text())["next_work_order"]
                continue
            feedback = ("\nFix round " + str(rnd + 1) + ". The line reran the tests; failing: " + ", ".join(failing)
                        + f". Their results are in {rdir}/tests/<id>/. Your last summary: {out['summary']}\n")
        else:
            return False
    st["stages"]["M6"] = {"status": "done", "note": f"layers {', '.join(order)}"}
    save_state(rdir, st)
    return True


def triage(rdir, st, test_id, rounds, last):
    st["triage_n"] += 1
    n = st["triage_n"]
    save_state(rdir, st)
    wtp = Path(st["worktree_project"])
    diff = git("diff", "main", "--", *PROBLEMS[st["problem"]]["code_dirs"], cwd=wtp, check=False)
    (rdir / "triage" / f"diff_{n}.patch").write_text(diff)
    (rdir / "triage" / f"modeler_{n}.json").write_text(json.dumps(last, indent=1))
    dispute = (f"The modeler disputes a test: {last.get('summary')}\nEvidence: {last.get('evidence')}"
               if last.get("status") == "dispute" else "")
    task = work_order("M6_triage", run=rdir, test_id=test_id, rounds=rounds, tests_dir=st["tests_dir"], n=n,
                      dispute=dispute)
    out, _ = md_agent.call(rdir, "M6", "beth-model-foreman", f"M6_triage_{n}", task, "triage", cwd=rdir,
                           read_roots=[str(rdir), st["tests_dir"]], write_roots=[str(rdir / "triage")],
                           deny_roots=deny_common(), timeout=1200)
    (rdir / "triage" / f"{n}.json").write_text(json.dumps(out, indent=1))
    route = out["route"]
    log(rdir, f"M6 triage {n}: test {test_id} -> {route}: {out['reason'][:160]}")
    if route == "test":
        task = work_order("M6_testfix", run=rdir, test_id=test_id, n=n)
        tdir = st["tests_dir"]
        md_agent.call(rdir, "M6", "beth-verifier", f"M6_testfix_{n}", task, "result", cwd=Path(tdir),
                      read_roots=[str(rdir / "derivation"), str(rdir / "triage"), str(rdir / "tests"), tdir],
                      write_roots=[tdir], deny_roots=deny_common() + [str(wtp / "src"), str(wtp / "include"),
                                                                      str(wtp / "toy")], timeout=1800)
        st["stages"]["M5"] = {"status": "gate_waiting", "note": f"test changed after triage {n}; re-lock"}
        st["stages"]["M6"] = {"status": "pending", "note": "waiting for the re-lock"}
        lock = Path(tdir) / "tests.json.lock"
        if lock.exists():
            lock.rename(lock.with_suffix(".lock.old"))
        save_state(rdir, st)
        notify("Model line", "G:Test lock reopened")
    elif route != "code":
        st["stages"]["M6"] = {"status": "question", "note": f"triage {n}: {route}; see triage/{n}.json"}
        save_state(rdir, st)
        notify("Model line", f"M6 triage: {route} needs you")
    return route


def st_mutation(rdir, st):
    p = subprocess.run([str(PY), str(HERE / "mutate.py"), "--run", str(rdir)], capture_output=True, text=True)
    log(rdir, "M7: " + (p.stdout.strip().splitlines() or ["no output"])[-1])
    res = json.loads((rdir / "mutation.json").read_text())
    missed = [m for m in res["mutations"] if m["applied"] and not m["caught"]]
    if missed:
        st["stages"]["M7"] = {"status": "question",
                              "note": f"{len(missed)} planted bug(s) not caught; the verifier adds tests, you re-lock"}
        save_state(rdir, st)
        notify("Model line", "M7: mutations missed")
        return False
    st["stages"]["M7"] = {"status": "done", "note": f"{len(res['mutations'])} mutations, all caught"}
    save_state(rdir, st)
    return True


def st_performance(rdir, st):
    wtp = Path(st["worktree_project"])
    task = work_order("M8_modeler", run=rdir, worktree=st["worktree"], feedback=feedback_text(st, "M8"),
                      perf_cases="the two cases named in the spec's ## Solver settings")
    bash = [r"make( -j\d+)?", r"git (status|diff)( .*)?",
            rf"{re.escape(str(RUN_TEST))} --tests {re.escape(st['tests_dir'])}/tests\.json --id [\w.-]+ --project "
            rf"{re.escape(str(wtp))} --out {re.escape(str(rdir))}/tests"]
    code = [str(wtp / d) for d in PROBLEMS[st["problem"]]["code_dirs"]] + [str(wtp / "inputs")]
    out = agent_step(rdir, st, "M8", "beth-modeler", "M8_modeler", task, cwd=wtp,
                     reads=[str(rdir), str(wtp)], writes=code + [str(rdir)], bash=bash,
                     denies=deny_common() + [st["tests_dir"]])
    if out is None:
        return False
    st["stages"]["M8"] = {"status": "done", "note": "review item: perf.md"}
    save_state(rdir, st)
    return True


def st_writeup(rdir, st):
    name = st.get("model_name") or "toy"
    doc = (ROOT / "docs" / "model" / f"{name}.md") if st["problem"] == "regolith" else rdir / "docs" / f"{name}.md"
    doc.parent.mkdir(parents=True, exist_ok=True)
    review = rdir / "gates" / "M9_student_review.md"
    reads = [str(rdir), st["tests_dir"]]
    fb = feedback_text(st, "M9")
    for attempt in (1, 2):
        task = work_order("M9_writer", run=rdir, doc=doc, feedback=fb)
        out = agent_step(rdir, st, "M9", "beth-writer", f"M9_writer_{attempt}", task, cwd=rdir, reads=reads,
                         writes=[str(doc.parent)], bash=[], stop_check=struct_cmd("M9", rdir, doc=doc),
                         struct=lambda: run_struct("M9", rdir, doc=doc))
        if out is None:
            return False
        task = work_order("M9_student", doc=doc, review=review)
        md_agent.call(rdir, "M9", "beth-writer", f"M9_student_{attempt}", task, "result", cwd=doc.parent,
                      read_roots=[str(doc)], write_roots=[str(review.parent)], deny_roots=deny_common(), timeout=1800)
        if attempt == 1:
            fb = f"\nA fresh student read the document; resolve every item in {review}.\n"
    gate_summary(rdir, st, "M9", [str(doc), str(review)])
    return False


STATIONS = {"M1": lambda r, s: st_derivation_stage(r, s, "M1"), "M2": lambda r, s: st_derivation_stage(r, s, "M2"),
            "M3": lambda r, s: st_derivation_stage(r, s, "M3"), "M4": lambda r, s: st_derivation_stage(r, s, "M4"),
            "M5": st_test_lock, "M6": st_implementation, "M7": st_mutation, "M8": st_performance,
            "M9": st_writeup}


def cmd_run(a):
    rid = a.run or latest_run()
    rdir = run_dir(rid)
    st = load_state(rdir)
    for s in STAGES:
        status = st["stages"][s]["status"]
        if status in ("approved", "done", "skipped"):
            continue
        if status in ("gate_waiting", "question", "failed"):
            print(status_line(rid))
            if status == "failed":
                print(f"  {s} failed: {st['stages'][s].get('note')}; fix the cause, then: run_md.py reopen {s} "
                      f"--note \"...\"")
            return 0
        if a.until and STAGES.index(s) > STAGES.index(a.until):
            print(f"stopped before {s} (--until {a.until})")
            return 0
        st["stages"][s] = {"status": "running", "started": now()}
        save_state(rdir, st)
        log(rdir, f"{s}:{NAMES[s]} starts")
        try:
            go_on = STATIONS[s](rdir, st)
        except UsageLimit as e:
            st["stages"][s] = {"status": "pending", "note": "stopped at the usage limit; resume"}
            save_state(rdir, st)
            log(rdir, f"usage limit reached in {s}: {e}")
            print(f"\nStopped cleanly at the usage limit. Resume later with:\n"
                  f"  venv_workflows/bin/python scripts/model_development/run_md.py run --run {rid}")
            return 3
        except md_agent.agentlib.AgentError as e:
            st = load_state(rdir)
            st["stages"][s] = {"status": "failed", "note": f"agent call failed: {str(e)[:200]}"}
            save_state(rdir, st)
            log(rdir, f"{s}: agent call failed after retries: {str(e)[:400]}")
            return 1
        st = load_state(rdir)
        if not go_on:
            print(status_line(rid))
            return 0
    print(status_line(rid))
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=(__doc__ or "").split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    i = sub.add_parser("init")
    i.add_argument("--run", required=True)
    i.add_argument("--problem", required=True)
    i.add_argument("--dev", action="store_true", help="develop phase: a dev run, not a result")
    s = sub.add_parser("status")
    s.add_argument("--run")
    s.add_argument("--brief", action="store_true")
    r = sub.add_parser("run")
    r.add_argument("--run")
    r.add_argument("--until", choices=STAGES)
    ap_ = sub.add_parser("approve")
    ap_.add_argument("gate")
    ap_.add_argument("--run")
    ap_.add_argument("--decision", required=True)
    ap_.add_argument("--model-name")
    ro = sub.add_parser("reopen")
    ro.add_argument("stage")
    ro.add_argument("--run")
    ro.add_argument("--note", required=True)
    a = ap.parse_args(argv)
    return {"init": cmd_init, "status": cmd_status, "run": cmd_run, "approve": cmd_approve,
            "reopen": cmd_reopen}[a.cmd](a) or 0


if __name__ == "__main__":
    sys.exit(main())

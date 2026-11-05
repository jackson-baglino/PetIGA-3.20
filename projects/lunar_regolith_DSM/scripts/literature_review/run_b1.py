#!/usr/bin/env python3
"""The B1 line: run the literature review's stations in stage-graph order (WORKFLOW.md §0).

The stage graph is the `run` blocks of workflows/literature-review/flowchart.json. Each block names its station
(a script in scripts/literature_review/, or none for a user gate), the stages it requires, whether it runs once or
per paper, its gate (none, review or blocking) and the files that show it is done. The papers are the rows of the
Papers table in B1_CONTRACT.md; each file's sha256 must start with the hash printed there.

    run_b1.py init    --run <dir> [--dev]          create the run folder and its manifest.json
    run_b1.py status  --run <dir> [--brief]        per stage: built, done / total, gate
    run_b1.py run     --run <dir> [--until NODE] [--max N]
    run_b1.py approve --run <dir> GATE [--item ID] --decision TEXT    the user only (hook H6 blocks the foreman)
    run_b1.py freeze  --approved-by-user           write freeze/b1-vN.json and start the execute phase
    run_b1.py thaw    --reason TEXT                end the execute phase (WORKFLOW.md §5 repair rule)
    run_b1.py verify                               do the frozen files still match the latest freeze?
    run_b1.py close   --run <dir> --reason TEXT    development only: retire a dev run (kept, no longer current)

`run` repeats: pick the first stage (in graph order) whose required stages are done and that is not done itself;
stop at a blocking gate, at a station that is not built, at a station that fails, or at --until / --max. A run
that is not --dev requires the execute phase and an unchanged freeze. Station output (which may quote papers) goes
to <run>/untrusted/logs/; this script prints only stage names, items, counts and exit codes. Every station runs
through wf_launch.py: only scripts on launcher.allow, each under its sandbox profile.
Exit codes: 0 all done or stopped at --until/--max; 3 waiting at a gate; 4 a station is not built; 1 a failure.
"""
import argparse
import datetime
import hashlib
import json
import os
import re
import sqlite3
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import wf_launch  # noqa: E402

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
WF = PROJECT / "workflows" / "literature-review"
SPEC = WF / "flowchart.json"
CONTRACT = WF / "B1_CONTRACT.md"
FREEZE = Path(os.environ.get("LR_FREEZE_DIR") or WF / "freeze")      # the variable is for tests only
RESEARCH = Path.home() / "Library/Mobile Documents/com~apple~CloudDocs/Documents/Research/Projects/lunar_regolith_DSM"
PYTHON = PROJECT / "venv_workflows" / "bin" / "python"
FROZEN = ["workflows/literature-review/agents/*.md", "workflows/literature-review/prompts/*.md",
          "workflows/literature-review/schemas/*.json", "workflows/literature-review/flowchart.json",
          "scripts/agentlib.py", "scripts/literature_review/*.py", "requirements-workflows.lock",
          "requirements-docparse.lock"]
DOCPARSE_PY = PROJECT / "venv_docparse" / "bin" / "python"
PARSER_MODEL = Path.home() / ".cache/huggingface/hub/models--opendatalab--MinerU2.5-Pro-2605-1.2B/refs/main"
WAITING, NOT_BUILT, FAILED = 3, 4, 1


def now():
    return datetime.datetime.now().astimezone().isoformat(timespec="seconds")


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


# ---- the stage graph and the papers ----

def stages(spec_path=None):
    """The run blocks, in graph order (the order of the nodes in flowchart.json)."""
    spec = json.loads(Path(spec_path or SPEC).read_text())
    out = []
    for n in spec["nodes"]:
        r = n.get("run")
        if r:
            out.append({"id": n["id"], "title": n["title"], **r})
    ids = {s["id"] for s in out}
    for s in out:
        for q in s.get("requires", []):
            if q not in ids:
                raise SystemExit(f"stage {s['id']} requires {q}, which has no run block")
    return out


def papers(contract=None, research=None, check_hash=True):
    """[{id, pdf, role}] from the contract's Papers table, each file checked against its printed hash prefix."""
    research = research or RESEARCH
    text = Path(contract or CONTRACT).read_text()
    sec = text[text.index("## Papers"):text.index("\n## ", text.index("## Papers") + 1)]
    rows = []
    for line in sec.splitlines():
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) < 7 or not re.fullmatch(r"[a-z0-9_]+", cells[0]) or cells[0] == "id":
            continue
        rel = re.search(r"`([^`]+\.pdf)`", cells[2])
        if not rel:
            raise SystemExit(f"contract row {cells[0]}: no PDF path")
        pdf = Path(research) / rel.group(1)
        if check_hash:
            if not pdf.exists():
                raise SystemExit(f"contract row {cells[0]}: {pdf} not found")
            if not sha256(pdf).startswith(cells[4]):
                raise SystemExit(f"contract row {cells[0]}: sha256 of {pdf.name} does not start with {cells[4]}")
        rows.append({"id": cells[0], "pdf": str(pdf), "role": cells[6]})
    if not rows:
        raise SystemExit("no papers in the contract's Papers table")
    return rows


def run_papers(run, check_hash=False):
    """The contract's papers this run covers (a development run may cover fewer; manifest "papers")."""
    plist = papers(check_hash=check_hash)
    m = Path(run) / "manifest.json"
    if m.exists():
        ids = {p["id"] for p in json.loads(m.read_text())["papers"]}
        plist = [p for p in plist if p["id"] in ids]
    return plist


def skipped(run):
    m = Path(run) / "manifest.json"
    return set(json.loads(m.read_text()).get("skip", [])) if m.exists() else set()


def items(stage, plist):
    return [p["id"] for p in plist] if stage.get("per") == "paper" else [None]


def fill(template, run, item, plist, freeze_version):
    pdf = next((p["pdf"] for p in plist if p["id"] == item), "")
    return (template.replace("{run}", str(run)).replace("{paper}", item or "").replace("{pdf}", pdf)
            .replace("{freeze}", freeze_version or "dev"))


# ---- state ----

def gate_files(run, stage, item):
    g = Path(run) / "gates"
    name = stage["id"] + (f"_{item}" if item else "")
    return g / f"{name}.json"


def gate_state(run, stage, item):
    """None (no gate file), 'waiting' or 'approved'."""
    f = gate_files(run, stage, item)
    if not f.exists():
        return None
    return json.loads(f.read_text()).get("status")


def outputs_done(run, stage, item, plist, freeze_version):
    for d in stage.get("done", []):
        f = Path(fill(str(Path(run) / d["path"]), run, item, plist, freeze_version))
        if not f.exists():
            return False
        if d.get("has") or d.get("absent"):
            try:
                obj = json.loads(f.read_text())
            except ValueError:
                return False
            if any(k not in obj for k in d.get("has", [])) or any(k in obj for k in d.get("absent", [])):
                return False
    return True


def item_done(run, stage, item, plist, freeze_version):
    if stage["id"] in skipped(run):                          # a development run that leaves this stage out
        return True
    if stage.get("station") is None:                         # a user gate
        return gate_state(run, stage, item) == "approved"
    if not outputs_done(run, stage, item, plist, freeze_version):
        return False
    return gate_state(run, stage, item) != "waiting"


def stage_done(run, stage, plist, freeze_version):
    return all(item_done(run, stage, i, plist, freeze_version) for i in items(stage, plist))


def ready(run, stage, item, by_id, plist, freeze_version):
    """Every required stage done: for the same paper if both run per paper, else for all its items."""
    for q in stage.get("requires", []):
        req = by_id[q]
        if stage.get("per") == "paper" and req.get("per") == "paper":
            if not item_done(run, req, item, plist, freeze_version):
                return False
        elif not stage_done(run, req, plist, freeze_version):
            return False
    return True


def built(stage):
    return all((HERE / step[0]).exists() for step in stage.get("steps", []))


def _db(run):
    con = sqlite3.connect(Path(run) / "runs.db")
    con.execute("""CREATE TABLE IF NOT EXISTS line (id INTEGER PRIMARY KEY, stage TEXT, item TEXT, started TEXT,
        seconds REAL, status TEXT, exit INTEGER, log TEXT)""")
    return con


# ---- freeze ----

def frozen_files():
    out = []
    for pat in FROZEN:
        out += sorted(PROJECT.glob(pat))
    return [p for p in out if p.is_file()]


def contract_papers_hash():
    text = CONTRACT.read_text()
    sec = text[text.index("## Papers"):text.index("\n## ", text.index("## Papers") + 1)]
    return hashlib.sha256(sec.encode()).hexdigest()


def parser_state():
    """What P1's parser really is, not only what its lock file says: the packages installed in venv_docparse
    (hash of the sorted name==version list) and the model snapshot the parser loads (Hugging Face refs/main)."""
    code = ("import importlib.metadata as m; print('\\n'.join(sorted({d.metadata['Name'].lower() + '==' + d.version "
            "for d in m.distributions()})))")
    try:
        out = subprocess.run([str(DOCPARSE_PY), "-c", code], capture_output=True, text=True, check=True).stdout
        installed = hashlib.sha256(out.encode()).hexdigest()
    except (OSError, subprocess.CalledProcessError):
        installed = "missing"
    model = PARSER_MODEL.read_text().strip() if PARSER_MODEL.exists() else "missing"
    return {"venv_docparse#installed": installed, "parser_model#snapshot": model}


def snapshot():
    files = {str(p.relative_to(PROJECT)): sha256(p) for p in frozen_files()}
    files["B1_CONTRACT.md#Papers"] = contract_papers_hash()
    files.update(parser_state())
    return files


def latest_freeze():
    fs = sorted(FREEZE.glob("b1-v*.json"), key=lambda p: int(p.stem.split("-v")[1])) if FREEZE.exists() else []
    return (fs[-1], json.loads(fs[-1].read_text())) if fs else (None, None)


def phase():
    """('develop' or 'execute', freeze version or None)."""
    path, f = latest_freeze()
    if f is None or f.get("thawed"):
        return "develop", (f or {}).get("version")
    return "execute", f["version"]


def claude_version():
    """The installed Claude Code version string, e.g. '2.1.318 (Claude Code)'."""
    return subprocess.run(["claude", "--version"], capture_output=True, text=True).stdout.strip()


def freeze_diff():
    _, f = latest_freeze()
    if f is None:
        return ["no freeze"]
    now_ = snapshot()
    diff = [k for k in sorted(set(now_) | set(f["files"])) if now_.get(k) != f["files"].get(k)]
    return diff


# ---- commands ----

def cmd_init(a):
    run = Path(a.run)
    if (run / "manifest.json").exists():
        raise SystemExit(f"{run}/manifest.json exists")
    ph, ver = phase()
    if not a.dev:
        if ph != "execute":
            raise SystemExit("not in the execute phase: freeze first, or use --dev for a development run")
        d = freeze_diff()
        if d:
            raise SystemExit(f"frozen files changed since {ver}: {d[:5]}")
    plist = papers()
    if a.papers or a.skip:
        if not a.dev:
            raise SystemExit("--papers and --skip are for development runs only (--dev)")
        unknown = set(a.papers or []) - {p["id"] for p in plist}
        bad = set(a.skip or []) - {s["id"] for s in stages()}
        if unknown or bad:
            raise SystemExit(f"unknown papers {sorted(unknown)} or stages {sorted(bad)}")
        if a.papers:
            plist = [p for p in plist if p["id"] in a.papers]
    run.mkdir(parents=True, exist_ok=True)
    git = subprocess.run(["git", "rev-parse", "HEAD"], cwd=PROJECT, capture_output=True, text=True).stdout.strip()
    dirty = bool(subprocess.run(["git", "status", "--porcelain", "--", "workflows/literature-review", "scripts"],
                                cwd=PROJECT, capture_output=True, text=True).stdout.strip())
    cc = claude_version()
    sys.path.insert(0, str(HERE))
    import agent_call
    models = {p.stem: agent_call.load(p.stem)["model"] for p in sorted(agent_call.DEFS.glob("*.md"))}
    manifest = {"run_id": run.name, "created": now(), "dev": a.dev, "phase": ph, "freeze": ver,
                "git_commit": git, "git_dirty": dirty, "claude_code": cc, "models": models,
                "papers": [{"id": p["id"], "sha256": sha256(p["pdf"]), "role": p["role"]} for p in plist],
                "skip": sorted(a.skip or []),
                "files": snapshot()}
    (run / "manifest.json").write_text(json.dumps(manifest, indent=1) + "\n")
    print(f"initialized {run} ({'dev' if a.dev else 'execute ' + str(ver)}; {len(plist)} papers)")


def dependency_order(spec):
    """The stages in an order where each comes after every stage it requires (ties keep the graph's order)."""
    out, done, left = [], set(), list(spec)
    while left:
        for s_ in left:
            if all(q in done for q in s_.get("requires", [])):
                out.append(s_)
                done.add(s_["id"])
                left.remove(s_)
                break
        else:                                   # a cycle cannot be ordered: keep the rest as given
            out += left
            break
    return out


def status_rows(run, plist=None):
    run = Path(run)
    spec = dependency_order(stages())
    plist = plist if plist is not None else run_papers(run)
    fv = json.loads((run / "manifest.json").read_text()).get("freeze") if (run / "manifest.json").exists() else None
    rows = []
    for s in spec:
        its = items(s, plist)
        done = sum(item_done(run, s, i, plist, fv) for i in its)
        waiting = sum(gate_state(run, s, i) == "waiting" for i in its)
        rows.append({"stage": s["id"], "built": built(s), "done": done, "total": len(its),
                     "skipped": s["id"] in skipped(run),
                     "gate": s.get("gate", "none"), "waiting": waiting if s.get("gate") == "blocking" else 0})
    return rows


def next_item(run, plist, fv, skip=()):
    """The first item not done whose requirements are done; items in `skip` (gates already waiting in this
    pass) are passed over, so a gate on one paper holds only that paper."""
    spec = stages()
    by_id = {s["id"]: s for s in spec}
    for s in spec:
        for i in items(s, plist):
            if (s["id"], i) in skip:
                continue
            if not item_done(run, s, i, plist, fv) and ready(run, s, i, by_id, plist, fv):
                return s, i
    return None, None


def review_items(run):
    """The review queue, each item marked decided or not. E5 and X1 items are claims the user decides in U3
    (u3_decisions.json); Y3 items are the memo's, decided when the user decides U4."""
    run = Path(run)
    f, d = run / "review_queue.json", run / "u3_decisions.json"
    decided = json.loads(d.read_text()) if d.exists() else {}
    u4 = run / "gates" / "U4.json"
    u4_done = u4.exists() and json.loads(u4.read_text()).get("status") == "approved"
    out = []
    for x in (json.loads(f.read_text()) if f.exists() else []):
        if x.get("stage") in ("E5", "X1"):
            dec = (decided.get(x.get("claim_id")) or {}).get("verdict")
        else:
            dec = "U4 approved" if x.get("stage") == "Y3" and u4_done else None
        out.append({**x, "decided": dec})
    return out


def review_count(run):
    """Review items not yet decided."""
    return sum(1 for x in review_items(run) if not x["decided"])


def brief(run):
    ph, ver = phase()
    head = f"B1 · {ph.upper()}{' ' + str(ver) if ver else ''}"
    if run is None or not (Path(run) / "manifest.json").exists():
        return head + " · no run open"
    plist = run_papers(run)
    fv = json.loads((Path(run) / "manifest.json").read_text()).get("freeze")
    s, i = next_item(Path(run), plist, fv)
    gates = sum(r["waiting"] for r in status_rows(run, plist))
    nxt = "all stages done" if s is None else f"next: {s['id']}{' ' + i if i else ''}" + \
        ("" if built(s) else " (not built)")
    return f"{head} · run {Path(run).name} · {nxt} · {gates} gate(s) waiting · {review_count(run)} review item(s) undecided"


def cmd_status(a):
    if a.brief:
        print(brief(a.run))
        return
    for r in status_rows(a.run):
        print(f"{r['stage']:5} {'built' if r['built'] else 'NOT BUILT':9} {r['done']}/{r['total']}"
              f"{'  SKIPPED (dev)' if r['skipped'] else ''}"
              f"{'  gate ' + r['gate'] if r['gate'] != 'none' else ''}"
              f"{'  waiting ' + str(r['waiting']) if r['waiting'] else ''}")
    print(brief(a.run))
    gates = sorted((Path(a.run) / "gates").glob("*.json")) if (Path(a.run) / "gates").exists() else []
    for g in gates:
        rec = json.loads(g.read_text())
        if rec.get("status") == "waiting":
            print(f"  gate waiting: {g.stem}: {rec.get('what', '')}")
    for item in review_items(a.run):
        print(f"  review: {item.get('stage', '?')} {item.get('item') or ''} {item.get('claim_id') or ''} "
              f"{item.get('kind', '')}".rstrip() + (f"  (decided: {item['decided']})" if item["decided"] else ""))


def cmd_close(a):
    """Retire a dev run: the manifest records when and why; the folder is kept, but the run is no longer current."""
    f = Path(a.run) / "manifest.json"
    m = json.loads(f.read_text())
    if not m["dev"]:
        raise SystemExit("only a dev run can be closed; an execute run ends at its last stage or by a thaw")
    m["closed"] = {"date": datetime.date.today().isoformat(), "reason": a.reason}
    f.write_text(json.dumps(m, indent=1) + "\n")
    print(f"closed {Path(a.run).name}")


def cmd_run(a):
    run = Path(a.run)
    m = json.loads((run / "manifest.json").read_text())
    if m.get("closed"):
        raise SystemExit(f"run closed on {m['closed']['date']}: {m['closed']['reason']}")
    if not m["dev"]:
        ph, ver = phase()
        if ph != "execute" or ver != m["freeze"]:
            raise SystemExit(f"run made under freeze {m['freeze']}; current phase {ph} {ver}")
        d = freeze_diff()
        if d:
            raise SystemExit(f"frozen files changed: {d[:5]}")
    plist = run_papers(run, check_hash=True)
    count = 0
    held = set()                        # gates waiting in this pass: that item waits, the line goes on with others
    while True:
        s, i = next_item(run, plist, m["freeze"], skip=held)
        if s is None:
            if held:
                print(f"WAITING: {len(held)} gate(s) need the user; nothing else can run")
                return WAITING
            print("all stages done")
            return 0
        label = s["id"] + (f" {i}" if i else "")
        if a.until and s["id"] == a.until:
            print(f"stopped before {label} (--until)")
            return 0
        if a.max is not None and count >= a.max:
            print(f"stopped before {label} (--max {a.max})")
            return 0
        if s.get("station") is None:
            gate = gate_files(run, s, i)
            if not gate.exists():
                gate.parent.mkdir(exist_ok=True)
                gate.write_text(json.dumps({"stage": s["id"], "item": i, "status": "waiting", "opened": now(),
                                            "what": s.get("ask", s["title"])}, indent=1))
            print(f"WAITING at user gate {label}: {s.get('ask', s['title'])}")
            held.add((s["id"], i))
            continue
        if gate_state(run, s, i) == "waiting":
            print(f"WAITING at gate {label}: see {gate_files(run, s, i).relative_to(run)}")
            held.add((s["id"], i))
            continue
        if not built(s):
            missing = [st[0] for st in s["steps"] if not (HERE / st[0]).exists()]
            print(f"NOT BUILT: {label} needs {', '.join(missing)}")
            return NOT_BUILT
        logdir = run / "untrusted" / "logs"
        logdir.mkdir(parents=True, exist_ok=True)
        log = logdir / f"{s['id']}{'_' + i if i else ''}_{datetime.datetime.now():%Y%m%d-%H%M%S}.log"
        t0, rc = datetime.datetime.now(), 0
        with open(log, "w") as lf:
            for step in s["steps"]:
                args = [fill(x, run, i, plist, m["freeze"]) for x in step[1:]]
                cmd = wf_launch.command(step[0], args, run, python=PYTHON, here=HERE)   # sandboxed, allow-listed
                lf.write("$ " + " ".join(cmd) + "\n")
                lf.flush()
                rc = subprocess.run(cmd, cwd=PROJECT, stdout=lf, stderr=subprocess.STDOUT).returncode
                if rc:
                    break
        secs = round((datetime.datetime.now() - t0).total_seconds(), 1)
        ok = rc == 0 and outputs_done(run, s, i, plist, m["freeze"])
        con = _db(run)
        con.execute("INSERT INTO line (stage, item, started, seconds, status, exit, log) VALUES (?,?,?,?,?,?,?)",
                    (s["id"], i, t0.isoformat(timespec="seconds"), secs, "ok" if ok else "failed", rc,
                     str(log.relative_to(run))))
        con.commit()
        con.close()
        count += 1
        if not ok:
            why = f"exit {rc}" if rc else "exit 0 but its done files are missing or incomplete"
            print(f"FAILED: {label} ({why}; {secs} s; log {log.relative_to(run)})")
            return FAILED
        print(f"done: {label} ({secs} s)")


def cmd_approve(a):
    run = Path(a.run)
    s = next((x for x in stages() if x["id"] == a.gate), None)
    if s is None or s.get("gate") != "blocking":
        raise SystemExit(f"{a.gate} is not a blocking gate")
    f = gate_files(run, s, a.item)
    rec = json.loads(f.read_text()) if f.exists() else {"stage": a.gate, "item": a.item}
    rec.update({"status": "approved", "decision": a.decision, "decided": now(), "decided_by": "user"})
    f.parent.mkdir(exist_ok=True)
    f.write_text(json.dumps(rec, indent=1))
    print(f"recorded: {a.gate}{' ' + a.item if a.item else ''} approved")


def cmd_freeze(a):
    if not a.approved_by_user:
        raise SystemExit("a freeze needs the user's approval: --approved-by-user")
    ph, ver = phase()
    if ph == "execute":
        raise SystemExit(f"already frozen ({ver}); thaw first")
    if subprocess.run(["git", "status", "--porcelain", "--", "workflows/literature-review", "scripts"],
                      cwd=PROJECT, capture_output=True, text=True).stdout.strip():
        raise SystemExit("commit the workflow's files first: a freeze must match a commit")
    _, f = latest_freeze()
    n = (int(f["version"].split("v")[1]) + 1) if f else 1
    version = f"b1-v{n}"
    FREEZE.mkdir(exist_ok=True)
    git = subprocess.run(["git", "rev-parse", "HEAD"], cwd=PROJECT, capture_output=True, text=True).stdout.strip()
    rec = {"version": version, "frozen": now(), "git_commit": git, "approved_by": "user",
           "claude_code": claude_version(), "files": snapshot()}
    (FREEZE / f"{version}.json").write_text(json.dumps(rec, indent=1) + "\n")
    print(f"froze {version} at {git[:8]}: {len(rec['files'])} files. Commit it and tag lr-{version}.")


def cmd_thaw(a):
    path, f = latest_freeze()
    if f is None or f.get("thawed"):
        raise SystemExit("not in the execute phase")
    f["thawed"] = {"date": now(), "reason": a.reason}
    path.write_text(json.dumps(f, indent=1) + "\n")
    print(f"{f['version']} thawed: back to develop")


def cmd_verify(a):
    d = freeze_diff()
    print("frozen files unchanged" if not d else "changed since the freeze:\n  " + "\n  ".join(d))
    _, f = latest_freeze()
    now_cc = claude_version()
    if f and f.get("claude_code") and f["claude_code"] != now_cc:   # reported, not blocking: each run's manifest
        print(f"Claude Code differs: frozen with {f['claude_code']}, now {now_cc}; C0 scores the checker "
              "on the version the run uses")                        # records its own version
    return 1 if d else 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("init")
    p.add_argument("--run", required=True)
    p.add_argument("--dev", action="store_true")
    p.add_argument("--papers", nargs="+", help="development only: the papers this run covers")
    p.add_argument("--skip", nargs="+", help="development only: stages left out (e.g. C0, whose held-out set is "
                                              "scored once, after the freeze)")
    p = sub.add_parser("status")
    p.add_argument("--run")
    p.add_argument("--brief", action="store_true")
    p = sub.add_parser("run")
    p.add_argument("--run", required=True)
    p.add_argument("--until")
    p.add_argument("--max", type=int)
    p = sub.add_parser("approve")
    p.add_argument("--run", required=True)
    p.add_argument("gate")
    p.add_argument("--item")
    p.add_argument("--decision", required=True)
    p = sub.add_parser("freeze")
    p.add_argument("--approved-by-user", action="store_true")
    p = sub.add_parser("thaw")
    p.add_argument("--reason", required=True)
    sub.add_parser("verify")
    p = sub.add_parser("close")
    p.add_argument("--run", required=True)
    p.add_argument("--reason", required=True)
    a = ap.parse_args(argv)
    if a.cmd == "status" and a.run is None:
        a.run = os.environ.get("LR_RUN") or current_run()
        if a.run is None:                       # no run yet: say so, with the command that starts one
            print(brief(None) + ". Start one with: run_b1.py init --run runs/literature-review/b1/<run_id> [--dev]")
            return 0
    return {"init": cmd_init, "status": cmd_status, "run": cmd_run, "approve": cmd_approve, "freeze": cmd_freeze,
            "thaw": cmd_thaw, "verify": cmd_verify, "close": cmd_close}[a.cmd](a) or 0


def current_run(root=None):
    """The most recently created run folder under runs/literature-review/b1/ that is not closed, or None."""
    root = Path(root) if root else PROJECT / "runs" / "literature-review" / "b1"
    ms = {p: json.loads((p / "manifest.json").read_text()) for p in root.glob("*") if (p / "manifest.json").exists()}
    runs = sorted((p for p, m in ms.items() if not m.get("closed")), key=lambda p: ms[p]["created"])
    return str(runs[-1]) if runs else None


if __name__ == "__main__":
    sys.exit(main())

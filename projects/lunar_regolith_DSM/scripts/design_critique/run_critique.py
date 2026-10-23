"""Independent critique of a workflow design (workflows/literature-review/WORKFLOW.md, section 11).

Each critic agent (.claude/agents/<name>.md) runs headless in a fresh empty folder holding only a copy of the
design, through the units-verification caller (scripts/units_verification/agents.py: clean folder, pinned model,
transcript audit, schema check, one retry). WebSearch is pre-approved because dontAsk mode refuses it otherwise.
The critics run in parallel and never see each other's output.

Writes runs/design-critique/<run_id>/: manifest.json, critiques/<agent>.json, runs.db, transcripts/ (local).
Refuses to run while the design, the agents, or this folder have uncommitted changes.

    venv_workflows/bin/python scripts/design_critique/run_critique.py \
        --design workflows/literature-review/WORKFLOW.md --agents rick-critic beth-critic
"""
import argparse
import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "units_verification"))
import agents  # noqa: E402
import uvlib  # noqa: E402

ROOT = uvlib.ROOT
SCHEMA_PATH = HERE / "critique_schema.json"
MAX_FINDINGS = 25
TIMEOUT_S = 3600
TASK = ("Review the workflow design in the file WORKFLOW.md in your working folder, following your instructions. "
        "Read the whole file before you start.")


def check(output):
    """Problems beyond the schema: duplicate ids, too many findings."""
    ids = [f["id"] for f in output["findings"]]
    problems = [f"duplicate id {i}" for i in sorted({i for i in ids if ids.count(i) > 1})]
    if len(ids) > MAX_FINDINGS:
        problems.append(f"{len(ids)} findings, limit {MAX_FINDINGS}")
    return problems


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--design", required=True, help="workflow design to critique, relative to the project root")
    ap.add_argument("--agents", nargs="+", required=True, help="critic agent names in .claude/agents/")
    ap.add_argument("--run-id", default=None)
    args = ap.parse_args()

    design = (ROOT / args.design).resolve()
    paths = [args.design, "scripts/design_critique", "scripts/units_verification/agents.py"] + \
            [f".claude/agents/{a}.md" for a in args.agents]
    dirty = uvlib.run(["git", "status", "--porcelain", "--", *paths], cwd=ROOT).strip()
    if dirty:
        sys.exit("uncommitted changes to the critique's definitions; commit first:\n" + dirty)
    commit = uvlib.run(["git", "rev-parse", "HEAD"], cwd=ROOT).strip()

    schema = json.loads(SCHEMA_PATH.read_text())
    agents.TIMEOUT_S = TIMEOUT_S
    run_id = args.run_id or uvlib.new_run_id()
    rdir = ROOT / "runs" / "design-critique" / run_id
    rdir.mkdir(parents=True, exist_ok=True)
    print(f"run {run_id}: design {args.design} (sha256 {uvlib.sha256(design)[:12]}), commit {commit[:8]}, "
          f"critics {args.agents}", flush=True)

    def one(name):
        try:
            out, info = agents.call_agent(rdir, "critique", name, name, TASK, {"WORKFLOW.md": design}, schema,
                                          check=check, allow=("WebSearch",))
        except agents.AgentError as e:
            return name, None, {"error": json.loads(str(e))}
        uvlib.write_json(rdir / "critiques" / f"{name}.json", out)
        return name, out, info

    with ThreadPoolExecutor(max_workers=len(args.agents)) as pool:
        results = list(pool.map(one, args.agents))

    critics = {}
    for name, out, info in results:
        rec = dict(info)
        if out is not None:
            sev = [f["severity"] for f in out["findings"]]
            rec.update({"output": uvlib.rel(rdir / "critiques" / f"{name}.json"), "n_findings": len(sev),
                        "by_severity": {s: sev.count(s) for s in ("critical", "major", "minor")},
                        "n_unclear": len(out["unclear"]),
                        "n_web_sourced": sum(1 for f in out["findings"] if f["source"])})
        critics[name] = rec
        print(f"{name}: " + (json.dumps({k: rec[k] for k in ("n_findings", "by_severity", "n_unclear")})
                             if out is not None else "FAILED, see manifest"), flush=True)

    manifest = {"workflow": "design-critique", "run_id": run_id, "created": uvlib.now(), "git_commit": commit,
                "design": {"path": args.design, "sha256": uvlib.sha256(design)},
                "schema_sha256": uvlib.sha256(SCHEMA_PATH), "task_sha256": uvlib.sha256_bytes(TASK.encode()),
                "script_sha256": uvlib.sha256(Path(__file__)), "tools": uvlib.tool_versions(), "critics": critics}
    uvlib.write_json(rdir / "manifest.json", manifest)
    print(f"manifest: {uvlib.rel(rdir / 'manifest.json')}")
    if any(out is None for _, out, _ in results):
        sys.exit(1)


if __name__ == "__main__":
    main()

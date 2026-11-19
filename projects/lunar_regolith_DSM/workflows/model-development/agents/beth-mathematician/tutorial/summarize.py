#!/usr/bin/env python3
"""Summarize beth-mathematician run folders as one table (and JSON), from their records only.

    venv_workflows/bin/python workflows/model-development/agents/beth-mathematician/tutorial/summarize.py \
        [--json out.json] [--regrade] <run folder> ...

For each run: task, mode (agent or single-prompt baseline), verdict, wall time, turns, tokens, tool calls, the
grader's checks (passed / total, and the names of failed ones) and the derivation.md format problems. With
--regrade, each answer.json is graded again by the task's current grade.py (no model call), and the table shows
that verdict next to the recorded one. Nothing here calls a model or the solver.
"""
import argparse
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
AGENT = HERE.parent
ROOT = AGENT.parents[3]
PY = ROOT / "venv_workflows" / "bin" / "python"


def turns(run):
    t = run / "transcript.jsonl"
    if not t.exists():
        return None
    for line in reversed(t.read_text().splitlines()):
        try:
            e = json.loads(line)
        except json.JSONDecodeError:
            continue
        if e.get("type") == "result":
            return e.get("num_turns")
    return None


def checks_of(grade):
    d = (grade or {}).get("details")
    if not isinstance(d, dict):
        return [], d
    return d.get("checks", []), None


def regrade(run, task):
    p = subprocess.run([str(PY), str(AGENT / "tasks" / task / "grade.py"), str(run / "answer.json")],
                       capture_output=True, text=True, timeout=600, cwd=run)
    try:
        return json.loads(p.stdout.strip().splitlines()[-1])
    except (IndexError, json.JSONDecodeError):
        return {"pass": False, "details": f"grader gave no JSON (exit {p.returncode})"}


def summarize(run, do_regrade=False):
    m = json.loads((run / "manifest.json").read_text())
    grade = json.loads((run / "grade.json").read_text()) if (run / "grade.json").exists() else None
    checks, note = checks_of(grade)
    row = {"run": run.name, "task": m.get("task"), "mode": "baseline" if m.get("mode") == "baseline-no-tools"
           else "agent", "verdict": m.get("verdict"), "seconds": m.get("seconds"), "turns": turns(run),
           "input_tokens": (m.get("usage") or {}).get("input_tokens"),
           "output_tokens": (m.get("usage") or {}).get("output_tokens"),
           "tool_calls": sum((m.get("tool_calls") or {}).values()),
           "checks_passed": sum(1 for c in checks if c.get("pass")), "checks_total": len(checks),
           "failed_checks": [c["check"] for c in checks if not c.get("pass")] + ([note] if note else []),
           "format_problems": len(m.get("format_problems") or [])}
    if do_regrade and (run / "answer.json").exists():
        g = regrade(run, m["task"])
        row["regraded"] = "PASS" if g.get("pass") else "FAIL"
    return row


def main(argv=None):
    ap = argparse.ArgumentParser(description=(__doc__ or "").split("\n\n")[0])
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--json")
    ap.add_argument("--regrade", action="store_true")
    a = ap.parse_args(argv)
    rows = [summarize(Path(r).resolve(), a.regrade) for r in a.runs]
    cols = ["task", "mode", "verdict"] + (["regraded"] if a.regrade else []) + [
        "seconds", "turns", "tool_calls", "checks_passed", "checks_total", "format_problems"]
    print("  ".join(f"{c:>17s}" if c == "task" else f"{c:>9s}" for c in cols))
    for r in rows:
        print("  ".join(f"{str(r.get(c)):>17s}" if c == "task" else f"{str(r.get(c)):>9s}" for c in cols)
              + ("   failed: " + ", ".join(map(str, r["failed_checks"])) if r["failed_checks"] else ""))
    if a.json:
        Path(a.json).write_text(json.dumps(rows, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""M7:Mutation: measure the locked tests by planting bugs in the new code (WORKFLOW.md §2, M7; TEMPLATE §4).

Deterministic. From the feature branch's diff against main (the code folders of the problem), take the added
lines that hold an arithmetic operator or a numeric literal, and make up to N single-token mutations spread over
the diff: a + turned into -, a - into +, or a numeric literal doubled. Each mutation is applied in a scratch
worktree at the branch head; the locked tests run in order until one fails (caught) or all pass (missed).
Results go to <run>/mutation.json. Comments and blank lines are skipped.

    venv_workflows/bin/python scripts/model_development/mutate.py --run <run folder> [--n 10]
"""
import argparse
import json
import re
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
REPO = ROOT.parents[1]
RUN_TEST = HERE / "run_test.sh"
NUM = re.compile(r"(?<![\w.])(\d+\.\d*(?:[eE][-+]?\d+)?|\d+[eE][-+]?\d+|[2-9]\d*)(?![\w.])")


def added_lines(wtp, dirs):
    diff = subprocess.run(["git", "diff", "-U0", "main", "--", *dirs], cwd=wtp, capture_output=True,
                          text=True).stdout
    out, path, lineno = [], None, 0
    for line in diff.splitlines():
        if line.startswith("+++ b/"):
            path = line[6:]
        elif line.startswith("@@"):
            lineno = int(re.search(r"\+(\d+)", line).group(1))
        elif line.startswith("+") and not line.startswith("+++"):
            text = line[1:]
            s = text.strip()
            if s and not s.startswith(("//", "/*", "*", "#")):
                out.append((path, lineno, text))
            lineno += 1
    return out


def mutations_for(text):
    code = text.split("//")[0].split("#")[0]
    ms = []
    if " + " in code:
        ms.append(("plus_to_minus", code.replace(" + ", " - ", 1)))
    if " - " in code:
        ms.append(("minus_to_plus", code.replace(" - ", " + ", 1)))
    m = NUM.search(code)
    if m:
        v = m.group(1)
        new = repr(2 * float(v)) if any(c in v for c in ".eE") else str(2 * int(v))
        ms.append(("double_literal", code[:m.start()] + new + code[m.end():]))
    return ms


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    ap.add_argument("--n", type=int, default=10)
    a = ap.parse_args(argv)
    rdir = Path(a.run)
    st = json.loads((rdir / "state.json").read_text())
    wtp = Path(st["worktree_project"])
    dirs = {"toy": ["toy"]}.get(st["problem"], ["src", "include"])
    rel_project = wtp.relative_to(Path(st["worktree"]))
    cands = [(p, n, t, kind, new) for p, n, t in added_lines(wtp, dirs) for kind, new in mutations_for(t)]
    step = max(1, len(cands) // a.n)
    chosen = cands[::step][:a.n]
    tests = json.loads((Path(st["tests_dir"]) / "tests.json").read_text())["tests"]
    results = []
    for i, (path, lineno, text, kind, new) in enumerate(chosen, 1):
        with tempfile.TemporaryDirectory(prefix="md_mut_") as tmp:
            scratch = Path(tmp) / "wt"
            subprocess.run(["git", "worktree", "add", "-q", "--detach", str(scratch), "HEAD"], cwd=wtp,
                           capture_output=True, text=True, check=True)
            try:
                f = scratch / path
                lines = f.read_text().splitlines(keepends=True)
                old = lines[lineno - 1]
                lines[lineno - 1] = new + ("\n" if old.endswith("\n") else "")
                f.write_text("".join(lines))
                caught_by = None
                for t in tests:
                    p = subprocess.run([str(RUN_TEST), "--tests", str(Path(st["tests_dir"]) / "tests.json"),
                                        "--id", t["id"], "--project", str(scratch / rel_project),
                                        "--out", str(rdir / "tests" / "_mutations" / f"m{i:02d}")],
                                       capture_output=True, text=True)
                    if p.returncode != 0:
                        caught_by = t["id"]
                        break
            finally:
                subprocess.run(["git", "worktree", "remove", "--force", str(scratch)], cwd=wtp,
                               capture_output=True)
        results.append({"n": i, "file": path, "line": lineno, "kind": kind, "applied": True,
                        "original": text.strip()[:200], "mutated": new.strip()[:200],
                        "caught": caught_by is not None, "caught_by": caught_by})
        print(f"mutation {i}: {path}:{lineno} {kind}: {'caught by ' + caught_by if caught_by else 'MISSED'}")
    out = {"candidates": len(cands), "mutations": results,
           "caught": sum(r["caught"] for r in results), "total": len(results)}
    (rdir / "mutation.json").write_text(json.dumps(out, indent=1) + "\n")
    print(f"M7: {out['caught']} of {out['total']} planted bugs caught")
    return 0


if __name__ == "__main__":
    sys.exit(main())

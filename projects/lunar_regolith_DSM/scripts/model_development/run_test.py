#!/usr/bin/env python3
"""Run one test of the model line, under the standing run caps (CLAUDE.md, Hard limits; user 2026-11-16).

    scripts/model_development/run_test.sh --tests <tests.json> --id <test id> --project <checkout> --out <dir>
                                          [--draft]

The caps, enforced here and nowhere else:
  - one run at a time on this machine (lock folder runs/model-development/.run_lock, stale locks cleared);
  - at most MAX_RANKS ranks;
  - DOFs per rank in [MIN_DOFS_PER_RANK, MAX_DOFS_PER_RANK]; a test on one rank may be smaller;
  - at most MAX_WALL_MIN minutes, after which the whole process group is killed;
  - only tests in the locked list: when <tests.json>.lock exists, the hash of tests.json and of every check
    script must match it (`--draft` skips this only while the verifier writes the tests in M5:Test lock).

Test entry (tests.json: {"tests": [ ... ]}):
  {"id": "T01", "claim": "...", "runner": "solver" | "python",
   solver: "geometry": <name>, "experiment": <name>, "flags": [...], "ranks": n,
   python: "script": <path relative to the tests.json folder>, "args": [...] (run with --project <checkout>),
   "wall_minutes": m, "check": <path of the check script>, "check_args": [...]}
Check and python scripts are found relative to the folder of tests.json, so tests kept on a feature branch can
run against main (regressions) or the branch. A solver test runs scripts/Studio/run_lunar.sh in --project with LUNAR_PROJECT_ROOT, LUNAR_NPROCS,
LUNAR_RESULTS_BASE and LUNAR_SKIP_POSTPROCESS=1; its check script then gets --run-folder <folder>.
A check script prints one JSON object with "pass" and exits 0 on pass. The result goes to
<out>/<id>/<timestamp>.json with the caps used, wall time, exit codes and the check's verdict.
"""
import argparse
import datetime
import hashlib
import json
import os
import re
import signal
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
LOCK = ROOT / "runs" / "model-development" / ".run_lock"
MAX_RANKS = 8
MIN_DOFS_PER_RANK = 20_000
MAX_DOFS_PER_RANK = 100_000
MAX_WALL_MIN = 30
VENV_LUNAR_PY = ROOT / "venv_lunar" / "bin" / "python"


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def lock_hashes(tests_path, tests):
    base = Path(tests_path).parent
    files = {"tests.json": sha256(tests_path)}
    for t in tests:
        for key in ("check", "script"):
            if t.get(key):
                files[t[key]] = sha256(base / t[key])
    return files


def find_opts(project, kind, name):
    hits = sorted((Path(project) / "inputs" / kind).rglob(f"{name}.opts"))
    hits = [h for h in hits if "scratch" not in h.parts]
    if len(hits) != 1:
        raise SystemExit(f"{kind} opts {name!r}: {len(hits)} matches in {project}/inputs/{kind}")
    return hits[0]


def grid(project, test):
    """(Nx, Ny, Nz, dof) the run will use: geometry opts, then the test's flags (later wins)."""
    text = find_opts(project, "geometry", test["geometry"]).read_text()
    vals = {}
    if re.search(r"^-geom_file", text, re.M):
        m = re.search(r"^#\s*DOF_GRID:\s*(\d+)\s+(\d+)(?:\s+(\d+))?", text, re.M)
        if not m:
            raise SystemExit("geometry uses -geom_file without a '# DOF_GRID:' comment")
        vals = {"-Nx": m.group(1), "-Ny": m.group(2), "-Nz": m.group(3) or "1"}
    for k in ("-Nx", "-Ny", "-Nz"):
        m = re.search(rf"^{k}\s+(\d+)", text, re.M)
        if m and k not in vals:
            vals[k] = m.group(1)
    flags = test.get("flags", [])
    for i, f in enumerate(flags[:-1]):
        if f in ("-Nx", "-Ny", "-Nz"):
            vals[f] = flags[i + 1]
    solver = (Path(project) / "inputs" / "solver.opts").read_text()
    m = re.search(r"^-dof\s+(\d+)", solver, re.M)
    dof = int(m.group(1)) if m else 3
    return int(vals.get("-Nx", 1)), int(vals.get("-Ny", 1)), int(vals.get("-Nz", 1)), dof


def check_caps(test, project):
    ranks = int(test.get("ranks", 1))
    wall = float(test.get("wall_minutes", MAX_WALL_MIN))
    problems, info = [], {"ranks": ranks, "wall_minutes": wall}
    if not 1 <= ranks <= MAX_RANKS:
        problems.append(f"ranks {ranks} outside 1..{MAX_RANKS}")
    if not 0 < wall <= MAX_WALL_MIN:
        problems.append(f"wall_minutes {wall} outside (0, {MAX_WALL_MIN}]")
    if test["runner"] == "solver":
        nx, ny, nz, dof = grid(project, test)
        total = dof * nx * ny * nz
        per = total / ranks
        info.update({"grid": [nx, ny, nz], "dof": dof, "total_dofs": total, "dofs_per_rank": round(per)})
        if per > MAX_DOFS_PER_RANK:
            problems.append(f"{per:.0f} DOFs per rank > {MAX_DOFS_PER_RANK}")
        if per < MIN_DOFS_PER_RANK and ranks > 1:
            problems.append(f"{per:.0f} DOFs per rank < {MIN_DOFS_PER_RANK} on {ranks} ranks")
    return problems, info


def pid_alive(pid):
    try:
        os.kill(pid, 0)
        return True
    except (OSError, ValueError):
        return False


def acquire_lock():
    LOCK.parent.mkdir(parents=True, exist_ok=True)
    try:
        LOCK.mkdir()
    except FileExistsError:
        owner = LOCK / "pid"
        pid = int(owner.read_text()) if owner.exists() and owner.read_text().strip().isdigit() else -1
        if pid_alive(pid):
            raise SystemExit(f"another test is running (pid {pid}); one run at a time") from None
        for f in LOCK.iterdir():
            f.unlink()
        LOCK.rmdir()
        LOCK.mkdir()
    (LOCK / "pid").write_text(str(os.getpid()))


def release_lock():
    try:
        for f in LOCK.iterdir():
            f.unlink()
        LOCK.rmdir()
    except FileNotFoundError:
        pass


def run_group(cmd, env, cwd, timeout_s, log):
    """Run cmd in its own process group; kill the group on timeout. (exit code, timed_out)."""
    with open(log, "w") as f:
        p = subprocess.Popen(cmd, cwd=cwd, env=env, stdout=f, stderr=subprocess.STDOUT, start_new_session=True)
        try:
            return p.wait(timeout=timeout_s), False
        except subprocess.TimeoutExpired:
            os.killpg(p.pid, signal.SIGTERM)
            try:
                p.wait(timeout=30)
            except subprocess.TimeoutExpired:
                os.killpg(p.pid, signal.SIGKILL)
                p.wait()
            return -9, True


def main(argv=None):
    ap = argparse.ArgumentParser(description=(__doc__ or "").split("\n\n")[0])
    ap.add_argument("--tests", required=True)
    ap.add_argument("--id", required=True)
    ap.add_argument("--project", default=str(ROOT), help="checkout to build and run (main or a feature worktree)")
    ap.add_argument("--out", required=True, help="folder for results (the line gives <run>/tests)")
    ap.add_argument("--draft", action="store_true", help="M5 only: allow an unlocked test list")
    a = ap.parse_args(argv)

    tests_path = Path(a.tests).resolve()
    project = Path(a.project).resolve()
    tests = json.loads(tests_path.read_text())["tests"]
    test = next((t for t in tests if t["id"] == a.id), None)
    if test is None:
        raise SystemExit(f"no test {a.id!r} in {tests_path}")
    lockfile = tests_path.with_suffix(".json.lock")
    if lockfile.exists():
        locked = json.loads(lockfile.read_text())["files"]
        now = lock_hashes(tests_path, tests)
        changed = sorted(k for k in set(locked) | set(now) if locked.get(k) != now.get(k))
        if changed:
            raise SystemExit(f"the locked test list changed: {changed}; only the user re-locks it")
    elif not a.draft:
        raise SystemExit(f"{tests_path} is not locked (G:Test lock); use --draft only while writing tests in M5")

    problems, caps = check_caps(test, project)
    if problems:
        raise SystemExit("run caps: " + "; ".join(problems))

    out = Path(a.out).resolve() / a.id
    out.mkdir(parents=True, exist_ok=True)
    ts = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    record = {"id": a.id, "claim": test.get("claim"), "runner": test["runner"], "project": str(project),
              "caps": caps, "started": ts, "draft": a.draft,
              "limits": {"max_ranks": MAX_RANKS, "dofs_per_rank": [MIN_DOFS_PER_RANK, MAX_DOFS_PER_RANK],
                         "max_wall_min": MAX_WALL_MIN}}
    acquire_lock()
    t0 = time.time()
    try:
        timeout_s = caps["wall_minutes"] * 60
        env = dict(os.environ)
        if test["runner"] == "solver":
            results = out / "runs"
            results.mkdir(exist_ok=True)
            env.update(LUNAR_PROJECT_ROOT=str(project), LUNAR_NPROCS=str(caps["ranks"]),
                       LUNAR_RESULTS_BASE=str(results), LUNAR_SKIP_POSTPROCESS="1")
            cmd = ["bash", str(project / "scripts" / "Studio" / "run_lunar.sh"), test["geometry"],
                   test["experiment"], f"md_{a.id}", "--", *test.get("flags", [])]
            code, timed_out = run_group(cmd, env, project, timeout_s, out / f"{ts}.log")
            text = (out / f"{ts}.log").read_text(errors="replace")
            m = re.findall(r"^Output folder: (.+)$", text, re.M)
            folder = m[-1].strip() if m else None
            record.update(exit=code, timed_out=timed_out, run_folder=folder)
            check_extra = ["--run-folder", folder] if folder else []
        else:
            cmd = [str(VENV_LUNAR_PY), str(tests_path.parent / test["script"]), "--project", str(project),
                   *test.get("args", [])]
            code, timed_out = run_group(cmd, env, project, timeout_s, out / f"{ts}.log")
            record.update(exit=code, timed_out=timed_out)
            check_extra = []
        verdict = None
        if test.get("check") and not timed_out:
            p = subprocess.run([str(VENV_LUNAR_PY), str(tests_path.parent / test["check"]), *check_extra,
                                *test.get("check_args", [])], cwd=project, capture_output=True, text=True,
                               timeout=600)
            try:
                verdict = json.loads(p.stdout.strip().splitlines()[-1])
            except (json.JSONDecodeError, IndexError):
                verdict = {"pass": False, "error": "check printed no JSON", "stderr": p.stderr[-500:]}
            verdict["exit"] = p.returncode
        record["check"] = verdict
        record["pass"] = bool(verdict and verdict.get("pass") and verdict.get("exit") == 0 and not timed_out)
    finally:
        record["wall_seconds"] = round(time.time() - t0, 1)
        release_lock()
        (out / f"{ts}.json").write_text(json.dumps(record, indent=1) + "\n")
    print(json.dumps({k: record[k] for k in ("id", "pass", "wall_seconds", "caps")} |
                     {"result": str(out / f"{ts}.json")}))
    return 0 if record["pass"] else 1


if __name__ == "__main__":
    sys.exit(main())

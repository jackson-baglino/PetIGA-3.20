#!/usr/bin/env python3
"""wf_launch: the only way a line station runs, under its own macOS sandbox profile (WORKFLOW.md §4).

run_b1.py builds every station command with command(): the script must be listed in
.claude/permissions/literature-review/launcher.allow (only the user edits that folder), and it runs as

    sandbox-exec -f <profile>.sb -D RUN_DIR=<run> -D TMP=<temp> -D PROJECT=... -D RESEARCH=... -D COURSES=...
                 [-D EXTRA=<one named place>] venv_workflows/bin/python scripts/literature_review/<script> <args>

The arguments come from the stage graph's run blocks, never from free text. Profiles:
- offline:      no network; writes only to the run folder and temp;
- offline-plus: the same, plus one named place (library/, the course archive, calibration/generated/);
- anthropic:    HTTPS and DNS out; no writes into the project, research or course folders except the run folder;
- api, zotero-local, zotero-web: scripts outside the B1 line.
"""
import os
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
DOCS = Path.home() / "Library/Mobile Documents/com~apple~CloudDocs/Documents"
PERM = Path(os.environ.get("LR_PERM_DIR") or PROJECT / ".claude/permissions/literature-review")   # tests only


class LaunchError(RuntimeError):
    pass


def extra_places():
    import backup
    return {"library": PROJECT / "library", "archive": backup.ROOT,
            "calibration": PROJECT / "workflows/literature-review/calibration/generated"}


def entries(perm=None):
    """{script: (profile, extra or None)} from launcher.allow."""
    out = {}
    for line in (Path(perm or PERM) / "launcher.allow").read_text().splitlines():
        parts = line.split("#")[0].split()
        if parts:
            out[parts[0]] = (parts[1], parts[2] if len(parts) > 2 else None)
    return out


def command(script, args, run, python=None, here=None, perm=None):
    """The sandboxed command for one station step. Raises LaunchError for a script not on the allow-list."""
    perm, here = Path(perm or PERM), Path(here or HERE)
    allowed = entries(perm)
    if script not in allowed:
        raise LaunchError(f"{script} is not in launcher.allow")
    profile, extra = allowed[script]
    sb = perm / "profiles" / f"{profile}.sb"
    if not sb.exists():
        raise LaunchError(f"profile {profile} not found")
    run = Path(run).resolve()
    params = {"RUN_DIR": run, "TMP": Path(tempfile.gettempdir()).resolve(), "PROJECT": PROJECT.resolve(),
              "RESEARCH": DOCS / "Research", "COURSES": DOCS / "Courses"}
    if profile == "offline-plus":
        if extra not in extra_places():
            raise LaunchError(f"{script}: offline-plus needs one of {sorted(extra_places())}")
        params["EXTRA"] = Path(extra_places()[extra]).resolve()
    cmd = ["sandbox-exec", "-f", str(sb)]
    for k, v in params.items():
        cmd += ["-D", f"{k}={v}"]
    return cmd + [str(python or PROJECT / "venv_workflows/bin/python"), str(here / script)] + [str(a) for a in args]

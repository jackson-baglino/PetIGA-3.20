#!/usr/bin/env python3
"""BK: back up a B1 run and its library records to a dated folder in the course folder (WORKFLOW.md §6).

    backup.py --run <run folder>

Copies the run folder (untrusted/ included) and, for each paper in the contract, library/papers/<paper>/, plus
library/events.jsonl, into <ME295>/library_archive/<run_id>_<date>/, and writes manifest.json there with the sha256
and size of every file. The copy holds paper text, so it is never in git. The summary <run>/summary_backup.json
holds the destination, file count, bytes and the manifest's hash; restore_test.py adds restore_ok.
"""
import argparse
import datetime
import hashlib
import json
import os
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import curator  # noqa: E402
import run_b1  # noqa: E402

COURSE = Path.home() / "Library/Mobile Documents/com~apple~CloudDocs/Documents/Courses/Fall2026/ME295"
ROOT = Path(os.environ.get("LR_BACKUP_ROOT") or COURSE / "library_archive")     # the variable is for tests only


def sha256(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def manifest(dest):
    files = {}
    for f in sorted(Path(dest).rglob("*")):
        if f.is_file() and f.name != "manifest.json":
            files[str(f.relative_to(dest))] = {"sha256": sha256(f), "bytes": f.stat().st_size}
    return files


def backup(run, library=None, root=None):
    run, library, root = Path(run).resolve(), Path(library or curator.LIBRARY), Path(root or ROOT)
    dest = root / f"{run.name}_{datetime.date.today().isoformat()}"
    if dest.exists():
        raise SystemExit(f"{dest} exists; a backup is written once")
    shutil.copytree(run, dest / "run")
    for p in run_b1.papers(check_hash=False):
        src = library / "papers" / p["id"]
        if src.exists():
            shutil.copytree(src, dest / "library" / "papers" / p["id"])
    if (library / "events.jsonl").exists():
        shutil.copyfile(library / "events.jsonl", dest / "library" / "events.jsonl")
    files = manifest(dest)
    m = {"run": run.name, "created": datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
         "files": files}
    (dest / "manifest.json").write_text(json.dumps(m, indent=1))
    s = {"dest": str(dest), "files": len(files), "bytes": sum(v["bytes"] for v in files.values()),
         "manifest_sha256": sha256(dest / "manifest.json")}
    (run / "summary_backup.json").write_text(json.dumps(s, indent=1))
    return s


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--run", required=True)
    a = ap.parse_args(argv)
    print(json.dumps(backup(a.run)))
    return 0


if __name__ == "__main__":
    sys.exit(main())

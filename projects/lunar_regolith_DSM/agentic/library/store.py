"""Library storage: catalog (SQLite), immutable versioned files, and the curator that alone may write them."""
from __future__ import annotations

import contextlib
import hashlib
import json
import os
import shutil
import sqlite3
import stat
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
LIB = Path(os.environ.get("LIBRARY_DIR", ROOT / "research" / "Library"))

SCHEMA = """
CREATE TABLE IF NOT EXISTS papers (id TEXT PRIMARY KEY, doi TEXT, title TEXT, authors TEXT, year INTEGER, venue TEXT, topics TEXT, wave TEXT,
    status TEXT, pdf_path TEXT, pdf_sha256 TEXT, n_claims INTEGER, quote_match REAL, checker_rate REAL, n_flagged INTEGER, added_ts TEXT);
CREATE TABLE IF NOT EXISTS files (rel TEXT, sha256 TEXT, kind TEXT, paper_id TEXT, version INTEGER, bytes INTEGER, ts TEXT, PRIMARY KEY (rel));
CREATE TABLE IF NOT EXISTS claims (paper_id TEXT, claim_id TEXT, page INTEGER, claim TEXT, quote TEXT, type TEXT, flag TEXT, PRIMARY KEY (paper_id, claim_id));
CREATE VIRTUAL TABLE IF NOT EXISTS claims_fts USING fts5(claim, quote, paper_id UNINDEXED, claim_id UNINDEXED);
CREATE TABLE IF NOT EXISTS flags (paper_id TEXT, claim_id TEXT, source TEXT, verdict TEXT, note TEXT);
CREATE TABLE IF NOT EXISTS candidates (id TEXT PRIMARY KEY, doi TEXT, title TEXT, year INTEGER, authors TEXT, abstract TEXT, scores TEXT, best_topic TEXT, in_library INTEGER);
CREATE TABLE IF NOT EXISTS searches (topic TEXT, query TEXT, source_file TEXT);
CREATE TABLE IF NOT EXISTS memos (name TEXT PRIMARY KEY, kind TEXT, rel TEXT, sha256 TEXT, derived_from TEXT, ts TEXT);
CREATE TABLE IF NOT EXISTS events (ts TEXT, actor TEXT, action TEXT, target TEXT, detail TEXT);
CREATE TABLE IF NOT EXISTS access_log (ts TEXT, requester TEXT, purpose TEXT, kind TEXT, request TEXT, returned TEXT);
"""


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%S")


class Library:
    """Read side. Opens the catalog read-only; cannot change anything."""

    def __init__(self, root: Path | None = None):
        self.root = Path(root or LIB)

    @property
    def catalog(self) -> Path:
        return self.root / "catalog.sqlite"

    def ro(self) -> sqlite3.Connection:
        con = sqlite3.connect(f"file:{self.catalog}?mode=ro", uri=True)
        con.row_factory = sqlite3.Row
        return con


class Curator(Library):
    """Write side. The only code path that changes the library; every change is hashed, versioned, made read-only and logged as an event."""

    def __init__(self, root: Path | None = None, actor: str = "curator"):
        super().__init__(root); self.actor = actor

    @contextlib.contextmanager
    def rw(self):
        self.root.mkdir(parents=True, exist_ok=True)
        cat = self.catalog
        self._lock(False)
        con = sqlite3.connect(cat); con.executescript(SCHEMA)
        try:
            yield con; con.commit()
        finally:
            con.close()
            self._lock(True)

    def _lock(self, on: bool) -> None:
        """Make every directory and file read-only (on) or writable for the curator (off). A reminder, not a wall: see POLICY.md."""
        dirs = [self.root] + [p for p in self.root.rglob("*") if p.is_dir()]
        if not on:
            for d in dirs:
                d.chmod(0o755)
            for f in self.root.rglob("*"):
                if f.is_file() and f.name == "catalog.sqlite":
                    f.chmod(0o644)
        else:
            for f in self.root.rglob("*"):
                if f.is_file():
                    f.chmod(0o444)
            for d in sorted(dirs, key=lambda x: -len(x.parts)):
                d.chmod(0o555)

    def _put(self, con, rel: str, data: bytes | Path, kind: str, paper_id: str | None) -> int:
        dest = self.root / rel; dest.parent.mkdir(parents=True, exist_ok=True)
        if dest.exists():                          # never overwrite: keep the old version under history/, publish the new one
            old = con.execute("SELECT version FROM files WHERE rel=?", (rel,)).fetchone()
            v_old = old[0] if old else 1
            hist = self.root / "history" / f"{rel}.v{v_old}"; hist.parent.mkdir(parents=True, exist_ok=True)
            dest.chmod(0o644); shutil.move(str(dest), str(hist)); hist.chmod(0o444); version = v_old + 1
        else:
            version = 1
        if isinstance(data, Path):
            shutil.copyfile(data, dest)
        else:
            dest.write_bytes(data)
        dest.chmod(0o444)
        con.execute("INSERT OR REPLACE INTO files VALUES (?,?,?,?,?,?,?)", (rel, sha256(dest), kind, paper_id, version, dest.stat().st_size, now()))
        con.execute("INSERT INTO events VALUES (?,?,?,?,?)", (now(), self.actor, "publish", rel, f"kind={kind} version={version}"))
        return version

    def publish_paper(self, slug: str, record: dict, verification: dict, record_md: str | None, text_path: Path | None, pdf_path: str | None, *, wave: str = "", topics: str = "", status: str = "extracted") -> None:
        paper = record.get("paper") or {}
        ck = verification.get("checker") or {}; m = verification.get("metrics") or {}
        with self.rw() as con:
            self._put(con, f"papers/{slug}/record.json", json.dumps(record, indent=1).encode(), "record", slug)
            self._put(con, f"papers/{slug}/verification.json", json.dumps(verification, indent=1).encode(), "verification", slug)
            if record_md:
                self._put(con, f"papers/{slug}/record.md", record_md.encode(), "record_md", slug)
            if text_path and Path(text_path).exists():
                self._put(con, f"papers/{slug}/text.txt", Path(text_path), "text", slug)
            flagged = ck.get("flagged") or {}
            con.execute("INSERT OR REPLACE INTO papers VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", (
                slug, (paper.get("doi") or verification.get("doi") or ""), paper.get("title", ""), "; ".join(paper.get("authors", [])) if isinstance(paper.get("authors"), list) else str(paper.get("authors", "")),
                paper.get("year"), paper.get("venue", ""), topics, wave, status, pdf_path or "", (sha256(Path(pdf_path)) if pdf_path and Path(pdf_path).exists() else ""),
                len(record.get("claims", [])), m.get("quote_match_rate"), ck.get("supported_rate"), len(flagged), now()))
            con.execute("DELETE FROM claims WHERE paper_id=?", (slug,)); con.execute("DELETE FROM claims_fts WHERE paper_id=?", (slug,))
            for c in record.get("claims", []):
                flag = (flagged.get(c["id"]) or {}).get("verdict", "")
                con.execute("INSERT INTO claims VALUES (?,?,?,?,?,?,?)", (slug, c["id"], c.get("page"), c.get("claim", ""), c.get("quote", ""), c.get("type", ""), flag))
                con.execute("INSERT INTO claims_fts VALUES (?,?,?,?)", (c.get("claim", ""), c.get("quote", ""), slug, c["id"]))
                if flag:
                    con.execute("INSERT INTO flags VALUES (?,?,?,?,?)", (slug, c["id"], "independent checker", flag, (flagged[c["id"]] or {}).get("problem", "")))

    def publish_memo(self, name: str, src: Path, kind: str, derived_from: list[str] | None = None) -> None:
        with self.rw() as con:
            rel = f"memos/{name}{Path(src).suffix}"
            self._put(con, rel, Path(src), "memo:" + kind, None)
            con.execute("INSERT OR REPLACE INTO memos VALUES (?,?,?,?,?,?)", (name, kind, rel, sha256(self.root / rel), json.dumps(derived_from or []), now()))

    def publish_file(self, rel: str, src: Path, kind: str) -> None:
        with self.rw() as con:
            self._put(con, rel, Path(src), kind, None)

    def set_status(self, paper_id: str, status: str, reason: str) -> None:
        with self.rw() as con:
            con.execute("UPDATE papers SET status=? WHERE id=?", (status, paper_id))
            con.execute("INSERT INTO events VALUES (?,?,?,?,?)", (now(), self.actor, "status", paper_id, f"{status}: {reason}"))

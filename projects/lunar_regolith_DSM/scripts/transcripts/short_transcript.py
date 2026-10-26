"""Abbreviated transcript of a Claude Code session: the user's prompts and Claude's final reply to each.

Made for the weekly course assignments, whose full /export runs to tens of pages. Reads the session log
(JSONL) and keeps, in order:
  - every prompt the user typed, verbatim (a cut-off copy followed by the full prompt is kept once);
  - for each assistant turn, the text after its last tool call (the reply shown at the end of the turn).
Dropped: tool calls and their output, interim narration between tool calls, model reasoning, system
reminders and slash-command records. A turn started by a background-task notice is labeled as such.
The newest prompt has no reply yet if the tool runs inside that session's own turn.

Usage (venv_workflows; .docx output also needs LibreOffice's `soffice`):
    venv_workflows/bin/python scripts/transcripts/short_transcript.py SESSION --out FILE.{docx,html,md}
        [--title TEXT] [--author TEXT] [--labels "1=Setup,2=Initial response,3=Round 1"]
SESSION is a .jsonl path or a session ID (looked up in ~/.claude/projects/*/<ID>.jsonl; IDs are logged
in .claude/sessions.tsv). Prints every input it used.
"""
import argparse
import glob
import html
import json
import pathlib
import shutil
import subprocess
import sys
import tempfile


def resolve_log(session):
    p = pathlib.Path(session).expanduser()
    if p.suffix == ".jsonl" and p.exists():
        return p
    hits = glob.glob(str(pathlib.Path.home() / ".claude" / "projects" / "*" / f"{session}.jsonl"))
    if len(hits) != 1:
        sys.exit(f"session {session!r}: expected one log, found {len(hits)}")
    return pathlib.Path(hits[0])


def _text(content):
    if isinstance(content, str):
        return content
    return "\n".join(x.get("text", "") for x in content if x.get("type") == "text")


def parse_turns(lines):
    """Return (turns, models). Each turn: {'kind': 'prompt'|'notice', 'text': str, 'reply': [str]}."""
    turns, cur, models = [], None, []
    for line in lines:
        d = json.loads(line)
        if d.get("type") == "user" and not d.get("isMeta"):
            c = d["message"]["content"]
            if isinstance(c, list) and any(x.get("type") == "tool_result" for x in c):
                if cur is not None:
                    cur["reply"] = []          # text before a tool call is narration, not the reply
                continue
            t = _text(c).strip()
            if not t or t.startswith(("<command-name>", "<local-command")):
                continue
            kind = "notice" if t.startswith("<task-notification>") else "prompt"
            prev = turns[-1] if turns else None
            if (kind == "prompt" and prev and prev["kind"] == "prompt" and not prev["reply"]
                    and t.startswith(prev["text"])):
                prev["text"] = t               # a cut-off copy was logged first; keep the full prompt
                cur = prev
                continue
            cur = {"kind": kind, "text": t, "reply": []}
            turns.append(cur)
        elif d.get("type") == "assistant" and cur is not None:
            m = d["message"].get("model")
            if m and m not in models and not m.startswith("<"):
                models.append(m)
            for b in d["message"]["content"]:
                if b.get("type") == "text" and b["text"].strip():
                    cur["reply"].append(b["text"].strip())
                elif b.get("type") == "tool_use":
                    cur["reply"] = []
    return turns, models


def parse_labels(s):
    out = {}
    for item in filter(None, (x.strip() for x in (s or "").split(","))):
        k, v = item.split("=", 1)
        out[int(k)] = v.strip()
    return out


def to_markdown(turns, title, intro, labels):
    out, n = [f"# {title}", "", f"*{intro}*", ""], 0
    for t in turns:
        if t["kind"] == "prompt":
            n += 1
            lab = f" ({labels[n]})" if n in labels else ""
            out += [f"## Prompt {n}{lab}", ""] + [f"> {ln}" if ln else ">" for ln in t["text"].split("\n")] + [""]
        else:
            out += ["*[Background task finished; Claude continued without a new prompt.]*", ""]
        if t["reply"]:
            out += ["### Claude", "", "\n\n".join(t["reply"]), ""]
        elif t["kind"] == "prompt":
            out += ["*[No reply in the log when this file was made.]*", ""]
    return "\n".join(out), n


def to_html(turns, title, intro, labels):
    import markdown2
    parts, n = [], 0
    for t in turns:
        if t["kind"] == "prompt":
            n += 1
            lab = f" &nbsp;<span style='font-weight:normal;color:#555'>({html.escape(labels[n])})</span>" \
                if n in labels else ""
            parts.append(f"<h2>Prompt {n}{lab}</h2>")
            # a one-cell shaded table: LibreOffice turns a bordered <div> into one box per line
            parts.append("<table width='100%' border='0' cellpadding='6' cellspacing='0'><tr><td bgcolor='#EEF3FA'>"
                         + "".join(f"<p>{html.escape(p).replace(chr(10), '<br>')}</p>"
                                   for p in t["text"].split("\n\n")) + "</td></tr></table>")
        else:
            parts.append("<p style='color:#555'><i>[Background task finished; Claude continued without a new "
                         "prompt.]</i></p>")
        if t["reply"]:
            parts.append("<h3>Claude</h3>")
            parts.append(markdown2.markdown("\n\n".join(t["reply"]),
                                            extras=["tables", "fenced-code-blocks", "cuddled-lists"]))
        elif t["kind"] == "prompt":
            parts.append("<p style='color:#555'><i>[No reply in the log when this file was made.]</i></p>")
    doc = f"""<!DOCTYPE html><html><head><meta charset="utf-8"><style>
@page {{ size: 8.5in 11in; margin: 0.9in; }}
body {{ font-family: Calibri, Arial, sans-serif; font-size: 10.5pt; }}
h1 {{ font-size: 15pt; }} h2 {{ font-size: 12.5pt; margin-top: 16pt; }} h3 {{ font-size: 11pt; }}
td, th {{ font-size: 9pt; vertical-align: top; }} code {{ font-family: Menlo, Consolas, monospace; font-size: 9pt; }}
</style></head><body><h1>{html.escape(title)}</h1><p><i>{html.escape(intro)}</i></p>{''.join(parts)}</body></html>"""
    return doc, n


def html_to_docx(html_text, out):
    soffice = shutil.which("soffice")
    if not soffice:
        sys.exit("soffice (LibreOffice) not found; write .html or .md instead")
    with tempfile.TemporaryDirectory() as tmp:
        src = pathlib.Path(tmp) / "transcript.html"
        src.write_text(html_text, encoding="utf-8")
        subprocess.run([soffice, "--headless", "--infilter=HTML (StarWriter)", "--convert-to",
                        "docx:MS Word 2007 XML", "--outdir", tmp, str(src)], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        shutil.move(str(pathlib.Path(tmp) / "transcript.docx"), out)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("session", help="session .jsonl path or session ID")
    ap.add_argument("--out", required=True, help="output file; format from the extension (.docx, .html, .md)")
    ap.add_argument("--title", default="Conversation, abbreviated")
    ap.add_argument("--author", default="", help="name shown in the header")
    ap.add_argument("--labels", default="", help='prompt labels, e.g. "1=Setup,2=Initial response,3=Round 1"')
    a = ap.parse_args(argv)

    log = resolve_log(a.session)
    turns, models = parse_turns(open(log, encoding="utf-8"))
    labels = parse_labels(a.labels)
    who = f"{a.author}. " if a.author else ""
    intro = (f"{who}Claude Code session {log.stem}, model {', '.join(models) or 'unknown'}. Every prompt the user "
             "typed, verbatim, and Claude's final reply to each, generated from the session log by "
             "scripts/transcripts/short_transcript.py. Left out: tool calls and their output, Claude's interim "
             "narration between tool calls, model reasoning and system messages. The full conversation is the "
             "/export of the same session.")
    out = pathlib.Path(a.out).expanduser().resolve()
    if out.suffix == ".md":
        text, n = to_markdown(turns, a.title, intro, labels)
        out.write_text(text, encoding="utf-8")
    elif out.suffix in (".html", ".docx"):
        text, n = to_html(turns, a.title, intro, labels)
        if out.suffix == ".html":
            out.write_text(text, encoding="utf-8")
        else:
            html_to_docx(text, out)
    else:
        sys.exit(f"unknown output format {out.suffix!r}")
    print(f"log: {log}\nmodels: {', '.join(models) or 'unknown'}\nlabels: {labels or 'none'}\n"
          f"prompts: {n}, turns: {len(turns)}\nwritten: {out}")


if __name__ == "__main__":
    main()

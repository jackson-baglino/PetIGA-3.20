#!/usr/bin/env python3
"""C0 calibration: split the planted errors and score the checks against them (WORKFLOW.md §5).

Two subcommands:

  split   Splits the plants into a tuning set (used when prompts are revised) and a held-out
          set (used only for the pass). The user's plants are always held out. The held-out
          size is fixed here, before any scoring, and recorded with the seed.

  score   Scores one run of the checks on the held-out set and applies the pass rules:
          - kinds the scripts catch by rule (fabricated quote, wrong page, wrong table value):
            every one caught;
          - every sign or exponent change caught;
          - the plants that reach E4 (not caught by E3): a one-sided 95 % Clopper-Pearson
            lower bound of at least 85 %. It is a gate only with at least 19 of them; with
            fewer the rate is reported ("e4") and the result is not certified;
          Each kind's catch rate is reported with its interval but is not a gate: n of n caught
          gives a bound of 0.05^(1/n), below 85 % for any n under 19 (revision 10, Codex finding 8).
          - at most 10 % of the controls flagged.
          An item without a verdict ("no_verdict": true) counts as a miss on a plant and as a
          flag on a control (abstention cannot raise the score); the report gives their count.
          A held-out set is scored once per prompt version, and a set that has failed is never
          scored again: the next prompt version needs a fresh held-out set. The ledger records
          every attempt, and the report gives the number of attempts with the bound. Misses on
          the held-out set are reported as counts per kind only, never as items, so they cannot
          steer the next prompt revision.

Plants come from calibration/planted_errors.yaml (source user, all held out) and from
calibration/generated/*_plants.yaml entries the user approved (source script or model).

    venv_workflows/bin/python scripts/literature_review/calibrate.py split --heldout-size 40
    venv_workflows/bin/python scripts/literature_review/calibrate.py score --results <results.jsonl> \\
        --prompt-version <v> --out <run folder>

  score --baseline  Scores the script-only baseline (baseline.py) on the same held-out items: the same
          rules are reported next to the checker's, but nothing is written to the ledger.

A results file has one JSON object per line: a plant {"plant_id", "caught_by": "E3"|"E4"|"X1"|null}
or a control {"control_id", "flagged": true|false}.
"""
import argparse
import collections
import datetime
import hashlib
import json
import os
import random
import sys

import yaml
from scipy.stats import beta

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import check_calibration as cc  # noqa: E402
from make_plants import plant_id  # noqa: E402

GEN = os.path.join(cc.CAL, "generated")
SPLIT = os.path.join(GEN, "split.yaml")
LEDGER = os.path.join(GEN, "ledger.json")

E3_KINDS = {"fabricated_quote", "wrong_page", "wrong_table_cell"}
ALL_CAUGHT_KINDS = {"equation_sign_or_exponent"}
BOUND = 0.85
MIN_E4 = 19
CONTROL_MAX = 0.10
CONF = 0.95


def cp_lower(k, n, conf=CONF):
    """One-sided Clopper-Pearson lower bound on a success rate (k of n)."""
    if n == 0:
        return 0.0
    return 0.0 if k == 0 else float(beta.ppf(1 - conf, k, n - k + 1))


def cp_upper(k, n, conf=CONF):
    """One-sided Clopper-Pearson upper bound on a rate (k of n)."""
    if n == 0:
        return 1.0
    return 1.0 if k == n else float(beta.ppf(conf, k + 1, n - k))


def copy_kinds():
    with open(os.path.join(cc.CAL, "papers.yaml")) as f:
        return {p["id"]: p.get("copy_kind") or ("latex" if p.get("source") else "pdf")
                for p in yaml.safe_load(f) or []}


def load_plants():
    """Every usable plant: the user's (filled entries) and the approved generated ones."""
    ck = copy_kinds()
    plants = []
    problems = []
    user = cc.load_yaml("planted_errors.yaml", problems) or {}
    for e in user.get("errors") or []:
        if not e.get("planted_claim"):
            continue
        p = {"paper": e.get("paper"), "kind": e.get("kind"), "correct_claim": e.get("correct_claim"),
             "planted_claim": e.get("planted_claim"), "planted_page": e.get("page")}
        plants.append(dict(p, id=plant_id("U", p), source="user"))
    if os.path.isdir(GEN):
        for fn in sorted(os.listdir(GEN)):
            if not fn.endswith("_plants.yaml"):
                continue
            source = fn[: -len("_plants.yaml")]
            with open(os.path.join(GEN, fn)) as f:
                d = yaml.safe_load(f) or {}
            for p in d.get("plants") or []:
                if p.get("approved") is True:
                    plants.append(dict(p, source=source))
    for p in plants:
        p["copy_kind"] = ck.get(p["paper"], "pdf")
    if problems:
        raise SystemExit("\n".join(problems))
    return plants


def main_page(page):
    return page[0] if isinstance(page, list) and page else page


def split(plants, heldout_size, seed):
    """User plants all held out; the rest drawn round-robin by kind until heldout_size.

    Plants made from the same fact go to the same side, so tuning on one never teaches the
    checker the answer to another: script plants from one gold statement form one group, and
    any plant on the paper and page of a user plant is held out with it.
    """
    held = [p["id"] for p in plants if p["source"] == "user"]
    user_pages = {(p["paper"], main_page(p.get("planted_page"))) for p in plants
                  if p["source"] == "user" and p.get("planted_page") is not None}
    groups = collections.defaultdict(list)
    for p in sorted((p for p in plants if p["source"] != "user"), key=lambda p: p["id"]):
        key = (p["paper"], p["from_statement"]) if p.get("from_statement") is not None else p["id"]
        groups[key].append(p)
    rest = collections.defaultdict(list)
    for _key, g in sorted(groups.items(), key=lambda kv: str(kv[0])):
        pages = {(p["paper"], main_page(p.get("correct_page", p.get("planted_page")))) for p in g}
        pages |= {(p["paper"], main_page(p.get("planted_page"))) for p in g}
        if pages & user_pages:
            held += [p["id"] for p in g]
        else:
            rest[min(p["kind"] for p in g)].append([p["id"] for p in g])
    rng = random.Random(seed)
    for k in rest:
        rng.shuffle(rest[k])
    kinds = sorted(rest)
    while len(held) < heldout_size and any(rest.values()):
        for k in kinds:
            if rest[k] and len(held) < heldout_size:
                held += rest[k].pop()
    tuning = [i for k in kinds for g in rest[k] for i in g]
    return sorted(held), sorted(tuning)


def set_id(ids):
    return hashlib.sha256("\n".join(sorted(ids)).encode()).hexdigest()[:12]


def score(plants, held_ids, results):
    by_id = {p["id"]: p for p in plants}
    caught = {}
    controls = []
    for r in results:
        if "control_id" in r:
            controls.append(bool(r.get("flagged")))
        elif r.get("plant_id") in held_ids:
            caught[r["plant_id"]] = r.get("caught_by")
    missing = sorted(set(held_ids) - set(caught))
    rows = [(by_id[i], caught.get(i)) for i in held_ids if i in by_id]
    gates = []
    # Kinds the scripts catch by rule.
    e3 = [(p, c) for p, c in rows if p["kind"] in E3_KINDS]
    gates.append(("kinds caught by rule (E3): all caught", all(c for _, c in e3),
                  f"{sum(bool(c) for _, c in e3)} of {len(e3)}"))
    for k in sorted(ALL_CAUGHT_KINDS):
        ks = [(p, c) for p, c in rows if p["kind"] == k]
        gates.append((f"{k}: all caught", all(c for _, c in ks), f"{sum(bool(c) for _, c in ks)} of {len(ks)}"))
    # Plants that reach E4.
    e4 = [(p, c) for p, c in rows if c != "E3" and p["kind"] not in E3_KINDS]
    n, k = len(e4), sum(bool(c) for _, c in e4)
    lb = cp_lower(k, n)
    certified = n >= MIN_E4                 # WORKFLOW §5: a gate only with 19 or more such plants; with fewer the
    if certified:                           # catch rate is reported, not certified (n of n < 19 cannot reach 85 %)
        gates.append((f"plants reaching E4: lower bound >= {BOUND:.0%}", lb >= BOUND,
                      f"{k} of {n} caught, lower bound {lb:.1%}"))
    per_kind = collections.defaultdict(lambda: [0, 0])
    for p, c in e4:
        per_kind[p["kind"]][1] += 1
        per_kind[p["kind"]][0] += bool(c)
    nf, nc = sum(controls), len(controls)
    gates.append((f"controls flagged <= {CONTROL_MAX:.0%}", nc > 0 and nf / nc <= CONTROL_MAX,
                  f"{nf} of {nc} flagged" if nc else "no controls"))
    if missing:
        gates.append(("every held-out plant has a result", False, f"{len(missing)} without a result"))
    report = {
        "gates": [{"rule": r, "pass": ok, "detail": d} for r, ok, d in gates],
        "passed": all(ok for _, ok, _ in gates),
        "per_kind": {kind: {"caught": kk, "n": nn, "lower_bound": round(cp_lower(kk, nn), 4)}
                     for kind, (kk, nn) in sorted(per_kind.items())},
        "by_source": {}, "by_copy_kind": {},
        "e4": {"caught": k, "n": n, "lower_bound": round(lb, 4), "certified": certified},
        "misses_per_kind": dict(collections.Counter(p["kind"] for p, c in rows if not c)),
        "no_verdict": {"plants": sum(1 for r in results if r.get("no_verdict") and "plant_id" in r),
                       "controls": sum(1 for r in results if r.get("no_verdict") and "control_id" in r)},
    }
    for key, field in (("by_source", "source"), ("by_copy_kind", "copy_kind")):
        groups = collections.defaultdict(lambda: [0, 0])
        for p, c in e4:
            groups[p[field]][1] += 1
            groups[p[field]][0] += bool(c)
        report[key] = {g: {"caught": kk, "n": nn, "lower_bound": round(cp_lower(kk, nn), 4)}
                       for g, (kk, nn) in sorted(groups.items())}
    return report


def load_ledger(path):
    if not os.path.exists(path):
        return []
    with open(path) as f:
        return json.load(f)


def cmd_split(args):
    plants = load_plants()
    held, tuning = split(plants, args.heldout_size, args.seed)
    if os.path.exists(args.split) and not args.force:
        raise SystemExit(f"{args.split} exists; a split is fixed once made (--force to replace it, "
                         "only before any scoring)")
    if load_ledger(args.ledger) and args.force:
        raise SystemExit("the ledger has scored runs; draw a fresh held-out set instead of replacing this one")
    os.makedirs(os.path.dirname(args.split), exist_ok=True)
    with open(args.split, "w") as f:
        f.write("# GENERATED by calibrate.py split. Fixed before any scoring; do not edit.\n")
        yaml.safe_dump({"seed": args.seed, "heldout_size": args.heldout_size,
                        "made": datetime.date.today().isoformat(), "heldout_set": set_id(held),
                        "heldout": held, "tuning": tuning}, f, sort_keys=False)
    srcs = collections.Counter(p["source"] for p in plants if p["id"] in held)
    print(f"held out {len(held)} ({dict(srcs)}), tuning {len(tuning)} -> {os.path.relpath(args.split, cc.ROOT)}")
    if args.heldout_size < MIN_E4:
        print(f"warning: fewer than {MIN_E4} held-out plants can never pass the {BOUND:.0%} bound")
    return 0


def cmd_score(args):
    with open(args.split) as f:
        sp = yaml.safe_load(f)
    held, sid = sp["heldout"], sp["heldout_set"]
    if args.baseline:                       # the script-only baseline: reported next to the checker, never a gate,
        with open(args.results) as f:       # and never an attempt on the held-out set
            results = [json.loads(l) for l in f if l.strip()]
        report = score(load_plants(), held, results)
        report.update({"heldout_set": sid, "baseline": "scripts only (baseline.py)",
                       "date": datetime.date.today().isoformat()})
        os.makedirs(args.out, exist_ok=True)
        with open(os.path.join(args.out, "baseline.json"), "w") as f:
            json.dump(report, f, indent=2)
        for g in report["gates"]:
            print(("met   " if g["pass"] else "unmet ") + g["rule"] + ": " + g["detail"])
        return 0
    ledger = load_ledger(args.ledger)
    same = [e for e in ledger if e["heldout_set"] == sid]
    if any(e["prompt_version"] == args.prompt_version for e in same):
        raise SystemExit(f"held-out set {sid} was already scored under prompt version {args.prompt_version}")
    if any(not e["passed"] for e in same):
        raise SystemExit(f"held-out set {sid} has failed before; draw a fresh held-out set "
                         "from the growing pool for the next prompt version")
    with open(args.results) as f:
        results = [json.loads(l) for l in f if l.strip()]
    report = score(load_plants(), held, results)
    report.update({"heldout_set": sid, "prompt_version": args.prompt_version,
                   "attempts_so_far": len(ledger) + 1, "date": datetime.date.today().isoformat()})
    os.makedirs(args.out, exist_ok=True)
    with open(os.path.join(args.out, "calibration.json"), "w") as f:
        json.dump(report, f, indent=2)
    ledger.append({"heldout_set": sid, "prompt_version": args.prompt_version,
                   "passed": report["passed"], "date": report["date"]})
    with open(args.ledger, "w") as f:
        json.dump(ledger, f, indent=2)
    for g in report["gates"]:
        print(("PASS " if g["pass"] else "FAIL ") + g["rule"] + ": " + g["detail"])
    e4 = report["e4"]
    print(f"plants reaching E4: {e4['caught']} of {e4['n']}, lower bound {e4['lower_bound']:.1%}"
          f"{'' if e4['certified'] else f' (reported, not certified: fewer than {MIN_E4})'}")
    print(f"attempts so far: {report['attempts_so_far']}; misses per kind: {report['misses_per_kind']}")
    return 0 if report["passed"] else 1


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("split")
    s.add_argument("--heldout-size", type=int, required=True)
    s.add_argument("--seed", type=int, default=20261028)
    s.add_argument("--force", action="store_true")
    t = sub.add_parser("score")
    t.add_argument("--results", required=True)
    t.add_argument("--prompt-version", required=True)
    t.add_argument("--out", required=True)
    t.add_argument("--baseline", action="store_true",
                   help="score baseline.py's results: report only, no ledger entry, never a gate")
    for p in (s, t):
        p.add_argument("--split", default=SPLIT)
        p.add_argument("--ledger", default=LEDGER)
    args = ap.parse_args(argv)
    return cmd_split(args) if args.cmd == "split" else cmd_score(args)


if __name__ == "__main__":
    sys.exit(main())

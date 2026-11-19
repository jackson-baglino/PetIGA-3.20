"""M3 scratch (try 3): structural self-check of ledger.json against the done criterion.

Not evidence for anyone else. Paths are passed on the command line:
    python -I m3_ledger_check.py <ledger.json> <derivation.md>
"""
import json
import re
import sys

ledger_path, md_path = sys.argv[1], sys.argv[2]
with open(ledger_path) as f:
    led = json.load(f)
with open(md_path) as f:
    md = f.read()

steps = led['steps']
ids = [s['id'] for s in steps]
problems = []
if len(ids) != len(set(ids)):
    problems.append('duplicate ids')

need = {
    'cited': ['paper_id', 'page_pdf', 'quote', 'equation'],
    'derived': ['from', 'operation'],
    'assumption': ['label', 'if_wrong'],
    'choice': ['label', 'options', 'recommendation'],
}
seen = set()
counts = {}
for s in steps:
    ty = s.get('type')
    counts[ty] = counts.get(ty, 0) + 1
    if ty not in need:
        problems.append(f"{s['id']}: bad type {ty}")
        continue
    for key in need[ty] + ['statement']:
        if key not in s or s[key] in ('', None, []):
            if not (ty == 'cited' and key == 'equation'):
                problems.append(f"{s['id']}: missing {key}")
    if ty == 'derived':
        for p in s['from']:
            if p not in seen:
                problems.append(f"{s['id']}: premise {p} not earlier in list order")
        # every step id named in the operation must be in 'from'
        named = set(re.findall(r'\bS\d+\b', s['operation']))
        missing = sorted(named - set(s['from']) - {s['id']})
        if missing:
            problems.append(f"{s['id']}: operation names {missing} not in from")
    seen.add(s['id'])

# the try-3 fixes
by = {s['id']: s for s in steps}
print('S12 from:', by['S12']['from'])
print('S21 from:', by['S21']['from'])
print('S22 from:', by['S22']['from'])
print('S1 states T_ref constant:', 'T_{ref}' in by['S1']['statement'] and 'constant' in by['S1']['statement'])
print('S13 from has S1:', 'S1' in by['S13']['from'], '| S13 op mentions T_ref:', 'T_{ref}' in by['S13']['operation'])

# coverage in derivation.md
md_ids = set(re.findall(r'\bS\d+\b', md))
print('ids not mentioned in derivation.md:', sorted(set(ids) - md_ids, key=lambda z: int(z[1:])))
tags = set(re.findall(r'\\tag\{(S\d+)\}', md))
print('ids with no tagged display:', sorted(set(ids) - tags, key=lambda z: int(z[1:])))
print('tags that are not ledger ids:', sorted(tags - set(ids)))
print('counts by type:', counts, ' total:', len(steps))
print('symbols:', len(led['symbols']))
print('problems:', problems if problems else 'none')

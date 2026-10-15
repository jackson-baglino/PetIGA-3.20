import json,sys,os,re,collections
D=sys.argv[1]; T='/Users/jacksonbaglino/PetIGA-3.20/projects/enceladus_DSM/Literature/text'
sizes={f[:-4]:os.path.getsize(f'{T}/{f}') for f in os.listdir(T)}
for aid,nm in [('a81c02641130656b9','A'),('a30aa517b33e2e342','B'),('a0b7f19d66fe8ee69','C')]:
    tools={};got=collections.Counter()
    for l in open(f'{D}/agent-{aid}.jsonl'):
        d=json.loads(l);m=d.get('message') or {}
        if d.get('type')=='assistant':
            for b in m.get('content') or []:
                if b.get('type')=='tool_use':
                    s=json.dumps(b.get('input'))
                    ks={k for k in re.findall(r'([A-Za-z0-9_\-]+)\.txt',s) if k in sizes}|set(re.findall(r'papers/([A-Za-z0-9_\-]+)\.pdf',s))
                    tools[b['id']]=ks
        if d.get('type')=='user' and isinstance(m.get('content'),list):
            for b in m['content']:
                if b.get('type')=='tool_result':
                    ks=tools.get(b.get('tool_use_id'),set())
                    c=b.get('content'); n=len(c) if isinstance(c,str) else sum(len(x.get('text','')) for x in c or [] if isinstance(x,dict))
                    for k in ks: got[k]+=n/len(ks)
    cov=[min(got.get(k,0)/sizes[k],9) for k in sizes]
    full=sum(c>=0.8 for c in cov); part=sum(0.2<=c<0.8 for c in cov); low=sum(c<0.2 for c in cov)
    tot=sum(got.values())
    print(f'{nm}: papers >=80% seen {full}, 20-80% {part}, <20% {low}; paper chars seen {tot/1e6:.2f}M of {sum(sizes.values())/1e6:.2f}M')
    if nm=='A':
        low=[k for k in sizes if got.get(k,0)/sizes[k]<0.2]
        import statistics
        N='/Users/jacksonbaglino/PetIGA-3.20/projects/enceladus_DSM/Literature/review/A/notes'
        for k in sorted(low)[:70]:
            p=f'{N}/{k}.md'; print(f'  {k[:45]:45s} seen {got.get(k,0)/sizes[k]:.0%} note {os.path.getsize(p) if os.path.exists(p) else 0}B')

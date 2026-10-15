import json,sys,os,re,collections
D=sys.argv[1]
names={'ae316b3d8e2392d37':'collector','a2d4c4ca85a484126':'checker','a81c02641130656b9':'reviewer A','a30aa517b33e2e342':'reviewer B','a0b7f19d66fe8ee69':'reviewer C'}
def size(c):
    if isinstance(c,str): return len(c),0
    n=img=0
    for b in c or []:
        if b.get('type')=='text': n+=len(b.get('text',''))
        elif b.get('type')=='image': img+=1
        elif b.get('type')=='tool_result':
            a,i=size(b.get('content')); n+=a; img+=i
    return n,img
tot=collections.Counter()
for aid,nm in names.items():
    f=f'{D}/agent-{aid}.jsonl'
    seen={};calls=[];tools={};cat=collections.Counter();catn=collections.Counter();imgs=0;compacts=0
    for l in open(f):
        d=json.loads(l)
        if d.get('isCompactSummary') or (d.get('type')=='system' and 'compact' in json.dumps(d)[:400].lower()): compacts+=1
        m=d.get('message') or {}
        if d.get('type')=='assistant':
            rid=d.get('requestId') or m.get('id')
            if m.get('usage'): seen[rid]=m['usage']
            for b in m.get('content') or []:
                if b.get('type')=='tool_use':
                    inp=b.get('input',{}); t=b['name']
                    if t=='Read':
                        p=inp.get('file_path','')
                        k='Read '+('paper text' if '/Literature/text/' in p else 'page PNG' if p.endswith('.png') else 'notes/brief/manuscript' )
                    elif t=='Bash':
                        c=inp.get('command','')
                        k='Bash mutool render' if 'mutool' in c else 'Bash grep/sed/cat paper text' if '/text/' in c or 'text/' in c else 'Bash other'
                    else: k=t
                    tools[b['id']]=k
        if d.get('type')=='user' and isinstance(m.get('content'),list):
            for b in m['content']:
                if b.get('type')=='tool_result':
                    k=tools.get(b.get('tool_use_id'),'?'); n,i=size(b.get('content'))
                    cat[k]+=n; catn[k]+=1; imgs+=i
    u=collections.Counter()
    ctx=[]
    for x in seen.values():
        for k in ('input_tokens','cache_creation_input_tokens','cache_read_input_tokens','output_tokens'): u[k]+=x.get(k,0)
        ctx.append(x.get('input_tokens',0)+x.get('cache_creation_input_tokens',0)+x.get('cache_read_input_tokens',0))
    tot+=u
    print(f'\n== {nm}: {len(seen)} API calls, compactions~{compacts}, images in results {imgs}')
    print('  input %d | cache_write %d | cache_read %d | output %d'%(u['input_tokens'],u['cache_creation_input_tokens'],u['cache_read_input_tokens'],u['output_tokens']))
    if ctx: print('  context/call: mean %d, max %d'%(sum(ctx)/len(ctx),max(ctx)))
    for k,v in cat.most_common(8): print(f'  {k:32s} calls {catn[k]:4d}  result chars {v:9d}')
print('\nTOTAL',dict(tot))

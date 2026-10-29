
from validation import require
from check_code import *
from noncyclic import code_generated

def basis_classes():
    out=[{'n':8,'mod':1,'gens':[[0]*8],'kind':'trivial'}]
    for r in json.loads((ROOT/'cyclic_results_8.json').read_text()).values():
        if r['status']=='accepted':out.append({'n':8,'mod':r['d'],'gens':[r['a']],'kind':'cyclic'})
    for r in json.loads((ROOT/'noncyclic_results_8.json').read_text()):
        if r['status']=='accepted':out.append({k:r[k] for k in ('n','mod','gens')}|{'kind':'binary'})
    require((len(out)==13), 'extend_rank8.py:10: verification check failed')
    return out

def extensions(basis):
    mod=basis['mod'];gens=basis['gens'];n=8
    code=code_generated(gens,mod)
    cols=[tuple(g[i] for g in gens) for i in range(n)]
    groups=[[i for i in range(n) if cols[i]==c] for c in sorted(set(cols))]
    selfinverse=[all(2*x%mod==0 for x in cols[group[0]]) for group in groups]
    def normal(v):
        out=list(v)
        for ids,selfinv in zip(groups,selfinverse):
            vals=sorted(abs(v[i]) if selfinv else v[i] for i in ids)
            for i,x in zip(ids,vals):out[i]=x
        return tuple(out)
    found=set()
    for c in code:
        options=[]
        for ids,selfinv in zip(groups,selfinverse):
            r=c[ids[0]]
            opts=[F(-1),F(0),F(1)] if r==0 else [F(r-mod,mod),F(r,mod)]
            if selfinv:opts=sorted(set(abs(x) for x in opts))
            options.append(list(itertools.combinations_with_replacement(opts,len(ids))))
        for blockvals in itertools.product(*options):
            v=[F(0)]*n
            for ids,vals in zip(groups,blockvals):
                for i,x in zip(ids,vals):v[i]=x
            if sum(bool(x) for x in v)<2:continue
            v=min(normal(v),normal(tuple(-x for x in v)))
            found.add(v)
    out=[]
    for v in sorted(found):
        swaps=[(i,j) for i in range(n) for j in range(i+1,n) if cols[i]==cols[j] and v[i]==v[j]]
        require((swaps), lambda: ('missing elementary orientation-reversing symmetry'))
        out.append({'extra':[str(x) for x in v],'swap':swaps[0]})
    return out

def fraction_code(gens,n=8):
    out={tuple([F(0)]*n)}
    for g in gens:
        g=tuple(x%1 for x in g);order=math.lcm(*(x.denominator for x in g))
        out={tuple((x[i]+k*g[i])%1 for i in range(n)) for x in out for k in range(order)}
    return tuple(sorted(out))

def canonical_pair(code,v):
    n=len(v);index=len(code)
    cyclic=[g for g in code if math.lcm(*(x.denominator for x in g))==index]
    best=None
    if cyclic:
        for g in cyclic:
            for s in (-1,1):
                cols=[]
                for r,x in zip(g,v):
                    x*=s
                    if r>F(1,2):r,x=1-r,-x
                    if r in (0,F(1,2)):x=abs(x)
                    cols.append((r,x))
                cols=sorted(cols)
                key=('cyclic',tuple(r for r,x in cols),tuple(x for r,x in cols))
                if best is None or key<best:best=key
    else:
        require((index==4 and all(2*x%1==0 for g in code for x in g)), 'extend_rank8.py:71: verification check failed')
        nz=[g for g in code if any(g)]
        for g,h in itertools.permutations(nz,2):
            cols=sorted((r,s,abs(x)) for r,s,x in zip(g,h,v))
            key=('binary',tuple((r,s) for r,s,x in cols),tuple(x for r,s,x in cols))
            if best is None or key<best:best=key
    return best

def cell_key(basis,v):
    mod=basis['mod'];gens=[tuple(F(x,mod) for x in g) for g in basis['gens']]
    code=fraction_code(gens);index=len(code)
    possibilities=[canonical_pair(code,v)]
    for i in range(8):
        if abs(v[i])!=1:continue
        def transform(x):
            xi=x[i]/v[i]
            return tuple(xi if j==i else x[j]-v[j]*xi for j in range(8))
        w=transform(tuple(F(i==j) for j in range(8)))
        newcode=fraction_code([transform(g) for g in gens]+[w])
        require((len(newcode)==index), 'extend_rank8.py:90: verification check failed')
        possibilities.append(canonical_pair(newcode,w))
    return str((index,min(possibilities)))

if __name__=='__main__':
    start=time.monotonic();bases=basis_classes();candidates=[]
    for i,b in enumerate(bases):
        es=extensions(b)
        candidates.extend({'basis':i,**e} for e in es)
    (ROOT/'extension_candidates_8.json').write_text(json.dumps({'bases':bases,'candidates':candidates}))
    print('candidate count',len(candidates),'by basis',[sum(c['basis']==i for c in candidates) for i in range(len(bases))],flush=True)
    if len(sys.argv)>1 and sys.argv[1]=='--enumerate-only':sys.exit(0)
    path=ROOT/'extension_results_8.json'
    results=json.loads(path.read_text()) if path.exists() else {}
    for num,c in enumerate(candidates):
        key=str(num)
        if key in results and results[key]['status']!='unresolved':continue
        if time.monotonic()-start>160:
            print('bounded batch complete',flush=True);break
        b=bases[c['basis']];v=tuple(F(x) for x in c['extra'])
        try:
            r=solve_code(8,b['mod'],b['gens'],v)
            if r['status']=='accepted':r['cell_key']=cell_key(b,v)
        except Exception as exc:r={'status':'unresolved','error':repr(exc)}
        results[key]={**c,**r};path.write_text(json.dumps(results))
        if num%20==0 or r['status']=='unresolved':print(num,c['basis'],r['status'],r.get('error',''),flush=True)
    print('processed',len(results),'status',{s:sum(r['status']==s for r in results.values()) for s in ('accepted','excluded','unresolved')},'distinct cells',len({r['cell_key'] for r in results.values() if r['status']=='accepted'}),flush=True)

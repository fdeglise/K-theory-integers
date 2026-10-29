
from validation import require
import itertools, math, json, pathlib, functools, time
ROOT=pathlib.Path(__file__).resolve().parent

@functools.lru_cache(maxsize=None)
def cyclic_type(v,mod):
    g=math.gcd(mod,*v); d=mod//g
    if d==1:return (1,())
    a=tuple(x//g for x in v)
    can=min(tuple(sorted(min(k*x%d,-k*x%d) for x in a)) for k in range(1,d) if math.gcd(k,d)==1)
    return (d,can)

def cyclic_accepted(n):
    data=json.loads((ROOT/'cyclic_accepted.json').read_text())[str(n)]
    return {(1,())}|{(r['d'],tuple(r['a'])) for r in data}

def code_generated(gens,mod):
    n=len(gens[0]); code={tuple([0]*n)}
    for g in gens:
        order=mod//math.gcd(mod,*g)
        code={tuple((x[i]+k*g[i])%mod for i in range(n)) for x in code for k in range(order)}
    return tuple(sorted(code))

def canonical_code(code,mod,orders):
    choices=[[v for v in code if mod//math.gcd(mod,*v)==o] for o in orders]
    best=None
    for gens in itertools.product(*choices):
        if code_generated(gens,mod)!=code:continue
        cols=tuple(sorted(min(col,tuple(-x%mod for x in col)) for col in zip(*gens)))
        if best is None or cols<best:best=cols
    require((best is not None), 'noncyclic.py:30: verification check failed')
    return best

def binary_types(n):
    result=[]
    for counts in itertools.combinations_with_replacement(range(n+1),3):
        if sum(counts)>n or counts[0]+counts[1]<5:continue
        cols=[(0,0)]*(n-sum(counts))+[(1,0)]*counts[0]+[(0,1)]*counts[1]+[(1,1)]*counts[2]
        result.append((counts,tuple(zip(*cols))))
    return result

def enumerate_noncyclic(n):
    accepted=cyclic_accepted(n); found={}; stats={}
    def add_case(gens,mod,orders,label):
        stats.setdefault(label,{'trials':0,'correct_order':0,'passes_cyclic':0})
        st=stats[label];st['trials']+=1
        code=code_generated(gens,mod)
        if len(code)!=math.prod(orders):return
        st['correct_order']+=1
        if any(cyclic_type(v,mod) not in accepted for v in code):return
        st['passes_cyclic']+=1
        can=canonical_code(code,mod,orders)
        key=str((mod,orders,can))
        if key not in found:
            canonical_gens=tuple(zip(*can))
            found[key]={'n':n,'mod':mod,'orders':orders,'gens':canonical_gens,'code':code_generated(canonical_gens,mod)}
    bins=binary_types(n)
    for counts,gens in bins:add_case(gens,2,(2,2),'C2xC2')
    triple_binary=0
    for counts,gens in bins:
        for v in itertools.product(range(2),repeat=n):
            code=code_generated((*gens,v),2)
            if len(code)==8 and all(cyclic_type(w,2) in accepted for w in code):triple_binary+=1
    require((triple_binary==0), 'noncyclic.py:63: verification check failed')
    stats['C2xC2xC2']={'survivors':triple_binary}
    for mod in (4,8):
        for d,a in sorted(accepted):
            if d!=mod:continue
            for v in itertools.product(range(2),repeat=n):
                add_case((a,tuple(mod//2*x for x in v)),mod,(mod,2),f'C{mod}xC2')
    for counts,gens in bins:
        cols=list(zip(*gens));opts=[]
        for u,v in cols:
            if (u,v)==(0,0):opts.append([(0,0),(0,2),(2,0),(2,2)])
            elif (u,v)==(1,0):opts.append([(1,0),(1,2)])
            elif (u,v)==(0,1):opts.append([(0,1),(2,1)])
            else:opts.append([(1,1),(1,3)])
        for lifted in itertools.product(*opts):
            add_case(tuple(zip(*lifted)),4,(4,4),'C4xC4')
    for counts,gens in bins:
        u,v=gens
        for weight in range(7,n+1):
            for supp in itertools.combinations(range(n),weight):
                a=[int(i in supp) for i in range(n)]
                g=tuple((3*u[i]+2*a[i])%6 for i in range(n))
                h=tuple(3*x for x in v)
                add_case((g,h),6,(6,2),'C6xC2')
    return {'n':n,'stats':stats,'candidates':list(found.values())}

if __name__=='__main__':
    import sys
    n=int(sys.argv[1]);r=enumerate_noncyclic(n)
    (ROOT/f'noncyclic_candidates_{n}.json').write_text(json.dumps(r))
    print(json.dumps({'n':n,'stats':r['stats'],'candidates':len(r['candidates']),'by_type':{str(o):sum(tuple(c['orders'])==o for c in r['candidates']) for o in ((2,2),(4,2),(4,4),(8,2),(6,2))}},indent=2))

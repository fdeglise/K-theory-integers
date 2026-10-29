"""Rational certificate verification; no numerical optimizer or published lists.

The geometric coverage argument is: the elementary index bound, Watson's
identity, the finite-abelian-group cases, and the maximal-determinant basis
bound for the ninth vector. All enumeration and certificate inputs are audited.
"""

from validation import require
from validation import save_summary
import pathlib,json,itertools,math,time
from fractions import Fraction as F
from collections import Counter
from noncyclic import enumerate_noncyclic,code_generated,canonical_code,cyclic_type
ROOT=pathlib.Path(__file__).resolve().parent

def read(name):return json.loads((ROOT/name).read_text())

def enumerate_cyclic(n,d):
    h=d//2;weights=[[min(k*a%d,-k*a%d) for k in range(1,h+1)] for a in range(h+1)]
    upper=[[max(weights[b][k] for b in range(a,h+1)) for k in range(h)] for a in range(h+1)]
    out=[]
    def search(prefix,sums,lo):
        left=n-len(prefix)
        if any(sums[k]+left*upper[lo][k]<=2*d for k in range(h)):return
        if left==0:
            if math.gcd(d,*prefix)==1 and cyclic_type(tuple(prefix),d)==(d,tuple(prefix)):out.append(prefix)
            return
        for a in range(lo,h+1):search(prefix+[a],[sums[k]+weights[a][k] for k in range(h)],a)
    search([], [0]*h,0)
    return out

def basis_vectors(n):
    return {tuple(F(s if i==j else 0) for i in range(n)) for j in range(n) for s in (-1,1)}

def gram_value(Q,x):
    return sum(Q[i][j]*x[i]*x[j] for i in range(len(x)) for j in range(len(x)))

def verify_form(Q,code,required):
    Q=[[F(x) for x in row] for row in Q];n=len(Q)
    require((all(len(row)==n for row in Q)), 'verify_completeness.py:37: verification check failed')
    require((all(Q[i][j]==Q[j][i] for i in range(n) for j in range(n))), 'verify_completeness.py:38: verification check failed')
    require((all(gram_value(Q,x)==1 for x in required)), 'verify_completeness.py:39: verification check failed')
    # Q=U^t D U, with upper triangular U and positive rational diagonal D.
    U=[[F(i==j) for j in range(n)] for i in range(n)];D=[]
    for i in range(n):
        di=Q[i][i]-sum(D[k]*U[k][i]**2 for k in range(i));require((di>0), 'verify_completeness.py:43: verification check failed');D.append(di)
        for j in range(i+1,n):U[i][j]=(Q[i][j]-sum(D[k]*U[k][i]*U[k][j] for k in range(i)))/di
    short=set();x=[F(0)]*n;nodes=0
    def recurse(i,budget,f):
        nonlocal nodes
        nodes+=1
        require((nodes<=3000000), lambda: ('certificate enumeration limit'))
        if i<0:
            if any(x):short.add(tuple(x))
            return
        offset=sum(U[i][j]*x[j] for j in range(i+1,n))
        radius2=budget/D[i];bound=math.isqrt(radius2.numerator//radius2.denominator)+1
        centre=-offset-f[i]
        for z in range(math.floor(centre)-bound,math.ceil(centre)+bound+1):
            x[i]=f[i]+z;used=D[i]*(x[i]+offset)**2
            if used<=budget:recurse(i-1,budget-used,f)
    for f in code:recurse(n-1,F(1),f)
    require((short==required), lambda: ((len(short),len(required))))
    return nodes

def model(n,mod,gens,extra):
    columns=[tuple(g[i] for g in gens)+(tuple() if extra is None else (extra[i],)) for i in range(n)]
    free=[not(all(2*g[i]%mod==0 for g in gens) and (extra is None or extra[i]==0)) for i in range(n)]
    pairs=[(i,j) for i in range(n) for j in range(i+1,n) if free[i] and free[j]]
    names=sorted({tuple(sorted((columns[i],columns[j]))) for i,j in pairs})
    pos={name:k for k,name in enumerate(names)}
    def feature(x):
        d=math.lcm(*(v.denominator for v in x));v=[int(d*t) for t in x];fs=[0]*len(names)
        for i,j in pairs:fs[pos[tuple(sorted((columns[i],columns[j])))]]+=2*v[i]*v[j]
        return fs,d*d,sum(t*t for t in v)
    return len(names),feature

def verify_result(n,mod,gens,result,extra=None,cyclic_format=False):
    code={tuple(F(x,mod) for x in c) for c in code_generated(gens,mod)}
    required=basis_vectors(n)
    if extra is not None:
        extra=tuple(F(x) for x in extra);required|={extra,tuple(-x for x in extra)}
        require((tuple(x%1 for x in extra) in code and sum(bool(x) for x in extra)>=2), 'verify_completeness.py:80: verification check failed')
    if result['status']=='accepted':return verify_form(result['Q'],code,required)
    require((result['status']=='excluded'), 'verify_completeness.py:82: verification check failed')
    m,feature=model(n,mod,gens,extra);cert=result['certificate']
    total=[F(0)]*(m+1);bound=F(0)
    if cyclic_format:
        require((len(cert['constraints'])==len(cert['vectors'])), 'verify_completeness.py:86: verification check failed')
        constraints=[{'row':r,'rhs':b,'weight':w,'vector':v} for (r,b,w),v in zip(cert['constraints'],cert['vectors'])]
        eqs=[]
    else:constraints=cert['constraints'];eqs=cert['equalities']
    for item in constraints:
        w=F(item['weight']);require((w>=0), 'verify_completeness.py:91: verification check failed')
        if item['vector'] is None:row=[0]*m+[1];rhs=1
        else:
            x=tuple(F(v) for v in item['vector']);require((any(x) and x not in required), 'verify_completeness.py:94: verification check failed')
            require((tuple(v%1 for v in x) in code), 'verify_completeness.py:95: verification check failed')
            fs,d2,normdiag=feature(x);row=[-a for a in fs]+[d2];rhs=normdiag-d2
        require((list(item['row'])==row and item['rhs']==rhs), 'verify_completeness.py:97: verification check failed')
        total=[a+w*b for a,b in zip(total,row)];bound+=w*rhs
    for item in eqs:
        require((extra is not None), 'verify_completeness.py:100: verification check failed')
        fs,d2,normdiag=feature(extra);row=fs+[0];rhs=d2-normdiag;w=F(item['weight'])
        require((list(item['row'])==row and item['rhs']==rhs), 'verify_completeness.py:102: verification check failed')
        total=[a+w*b for a,b in zip(total,row)];bound+=w*rhs
    require((total==[0]*m+[1]), 'verify_completeness.py:104: verification check failed')
    require((bound==F(cert['bound']) and bound<=0), 'verify_completeness.py:105: verification check failed')
    return 0

def regenerate_extensions(basis):
    mod=basis['mod'];gens=basis['gens'];code=code_generated(gens,mod)
    cols=[tuple(g[i] for g in gens) for i in range(8)]
    groups=[[i for i,c in enumerate(cols) if c==col] for col in sorted(set(cols))]
    def normalize(x):
        out=list(x)
        for ids in groups:
            selfinverse=all(2*a%mod==0 for a in cols[ids[0]])
            vals=sorted(abs(x[i]) if selfinverse else x[i] for i in ids)
            for i,v in zip(ids,vals):out[i]=v
        return tuple(out)
    out=set()
    for c in code:
        choices=[(F(-1),F(0),F(1)) if r==0 else (F(r,mod),F(r-mod,mod)) for r in c]
        for x in itertools.product(*choices):
            if sum(bool(v) for v in x)<2:continue
            out.add(min(normalize(x),normalize(tuple(-v for v in x))))
    return out

def fraction_code(gens,n):
    code={tuple([F(0)]*n)}
    for g in gens:
        g=tuple(F(x)%1 for x in g);order=math.lcm(*(x.denominator for x in g))
        code={tuple((x[i]+k*g[i])%1 for i in range(n)) for x in code for k in range(order)}
    return tuple(sorted(code))

def canonical_pair(code,v):
    index=len(code);generators=[g for g in code if math.lcm(*(x.denominator for x in g))==index];keys=[]
    if generators:
        for g in generators:
            for sign in (-1,1):
                cols=[]
                for a,x in zip(g,v):
                    x*=sign
                    if a>F(1,2):a,x=1-a,-x
                    if 2*a%1==0:x=abs(x)
                    cols.append((a,x))
                cols=sorted(cols);keys.append(('cyclic',tuple(a for a,x in cols),tuple(x for a,x in cols)))
    else:
        require((index==4 and all(2*x%1==0 for g in code for x in g)), 'verify_completeness.py:147: verification check failed')
        for g,h in itertools.permutations([g for g in code if any(g)],2):
            cols=sorted((a,b,abs(x)) for a,b,x in zip(g,h,v))
            keys.append(('binary',tuple((a,b) for a,b,x in cols),tuple(x for a,b,x in cols)))
    return min(keys)

def canonical_cell(gens,v):
    code=fraction_code(gens,8);index=len(code);keys=[canonical_pair(code,v)]
    for i,vi in enumerate(v):
        if abs(vi)!=1:continue
        def change(x):return tuple(x[i]/vi if j==i else x[j]-v[j]*x[i]/vi for j in range(8))
        w=change(tuple(F(j==i) for j in range(8)))
        newcode=fraction_code([change(g) for g in gens]+[w],8)
        require((len(newcode)==index), 'verify_completeness.py:160: verification check failed')
        keys.append(canonical_pair(newcode,w))
    return str((index,min(keys)))

def main():
    start=time.monotonic();summary={};nodes=0;expected_basis=[]
    for n,bound in ((8,19),(9,37)):
        candidates=read('cyclic_candidates.json')[str(n)]
        regenerated={d:enumerate_cyclic(n,d) for d in range(2,bound+1)}
        require(({int(d):v['candidates'] for d,v in candidates.items()}==regenerated), 'verify_completeness.py:169: verification check failed')
        cyclic=read(f'cyclic_results_{n}.json')
        require(({(r['d'],tuple(r['a'])) for r in cyclic.values()}=={(d,tuple(a)) for d,vs in regenerated.items() for a in vs}), 'verify_completeness.py:171: verification check failed')
        for r in cyclic.values():nodes+=verify_result(n,r['d'],[r['a']],r,cyclic_format=True)
        accepted=[r for r in cyclic.values() if r['status']=='accepted']
        require(({(r['d'],tuple(r['a'])) for r in accepted}=={(r['d'],tuple(r['a'])) for r in read('cyclic_accepted.json')[str(n)]}), 'verify_completeness.py:174: verification check failed')
        require((all(r['d']<=8 for r in accepted)), 'verify_completeness.py:175: verification check failed')
        require((all(sum(bool(x) for x in r['a'])>=7 for r in accepted if r['d']==3)), 'verify_completeness.py:176: verification check failed')
        require((all(sum(bool(x) for x in r['a'])>=8 for r in accepted if r['d']==5)), 'verify_completeness.py:177: verification check failed')
        noncyclic=enumerate_noncyclic(n)
        require((noncyclic==read(f'noncyclic_candidates_{n}.json') or json.loads(json.dumps(noncyclic))==read(f'noncyclic_candidates_{n}.json')), 'verify_completeness.py:179: verification check failed')
        results=read(f'noncyclic_results_{n}.json')
        require((len(results)==len(noncyclic['candidates'])), 'verify_completeness.py:181: verification check failed')
        for c,r in zip(noncyclic['candidates'],results):
            require((all(json.loads(json.dumps(c[k]))==r[k] for k in c)), 'verify_completeness.py:183: verification check failed')
            nodes+=verify_result(n,r['mod'],r['gens'],r)
        summary[n]={'cyclic_candidates':len(cyclic),'accepted_cyclic':len(accepted),'accepted_noncyclic':sum(r['status']=='accepted' for r in results),'basis_classes':1+len(accepted)+sum(r['status']=='accepted' for r in results)}
        allbases=[{'n':n,'mod':1,'gens':[[0]*n]}]+[{'n':n,'mod':r['d'],'gens':[r['a']]} for r in accepted]+[{k:r[k] for k in ('n','mod','gens')} for r in results if r['status']=='accepted']
        if n==8:expected_basis=allbases
        else:
            for b in allbases:
                cols=list(zip(*b['gens']));require((len(set(cols))<n), 'verify_completeness.py:190: verification check failed')
        print('verified basis classification',n,summary[n],flush=True)
    extension_input=read('extension_candidates_8.json');bases=extension_input['bases'];cs=extension_input['candidates']
    require(([{k:b[k] for k in ('n','mod','gens')} for b in bases]==expected_basis), 'verify_completeness.py:193: verification check failed')
    results=read('extension_results_8.json');require((set(results)=={str(i) for i in range(len(cs))}), 'verify_completeness.py:194: verification check failed')
    for i,b in enumerate(bases):
        generated=regenerate_extensions(b)
        require((generated=={tuple(F(x) for x in c['extra']) for c in cs if c['basis']==i}), 'verify_completeness.py:197: verification check failed')
    keys=set()
    for i,c in enumerate(cs):
        b=bases[c['basis']];v=tuple(F(x) for x in c['extra']);r=results[str(i)]
        require((all(c[k]==r[k] for k in c)), 'verify_completeness.py:201: verification check failed')
        j,k=c['swap'];require((j!=k and v[j]==v[k] and all(g[j]==g[k] for g in b['gens'])), 'verify_completeness.py:202: verification check failed')
        nodes+=verify_result(8,b['mod'],b['gens'],r,extra=v)
        if r['status']=='accepted':
            key=canonical_cell([tuple(F(x,b['mod']) for x in g) for g in b['gens']],v)
            require((key==r['cell_key']), 'verify_completeness.py:206: verification check failed');keys.add(key)
    summary['rank8_extensions']={'candidates':len(cs),'accepted_presentations':sum(r['status']=='accepted' for r in results.values()),'excluded':sum(r['status']=='excluded' for r in results.values()),'classes':len(keys),'orientation_reversing':len(cs)}
    summary['sphere_enumeration_nodes']=nodes;summary['seconds']=round(time.monotonic()-start,3)
    print(json.dumps(summary,indent=2),flush=True)
    save_summary(ROOT/'verification_summary.json', json.dumps(summary,indent=2))

if __name__=='__main__':main()

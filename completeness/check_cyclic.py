
from validation import require
import sys, pathlib
ROOT = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'deps'))
import numpy as np
from scipy.optimize import linprog
import sympy as sp
import itertools, json, math, time
from fractions import Fraction as F

def qnorm(Q, v):
    return sum(v[i]*v[j]*Q[i][j] for i in range(len(v)) for j in range(len(v)))

def ldl(Q):
    n=len(Q); L=[[F(i==j) for j in range(n)] for i in range(n)]; D=[]
    for j in range(n):
        d=Q[j][j]-sum(L[j][k]**2*D[k] for k in range(j))
        if d<=0: return None
        D.append(d)
        for i in range(j+1,n):
            L[i][j]=(Q[i][j]-sum(L[i][k]*L[j][k]*D[k] for k in range(j)))/d
    return L,D

def enumerate_short(Q, cosets, cutoff=F(1), max_nodes=3000000):
    decomp=ldl(Q)
    if decomp is None: return None
    L,D=decomp; n=len(Q); x=[F(0)]*n; found=set(); nodes=0
    def visit(i, remain, frac):
        nonlocal nodes
        nodes+=1
        if nodes>max_nodes: raise RuntimeError('exact vector enumeration node limit')
        if i<0:
            if any(x): found.add(tuple(x))
            return
        shift=sum(L[j][i]*x[j] for j in range(i+1,n))
        t=remain/D[i]
        if t<0:return
        radius=math.isqrt(t.numerator//t.denominator)+1
        center=-shift-frac[i]
        lo=math.floor(center)-radius; hi=math.ceil(center)+radius
        for z in range(lo,hi+1):
            x[i]=frac[i]+z
            used=D[i]*(x[i]+shift)**2
            if used<=remain:visit(i-1,remain-used,frac)
    for frac in cosets:visit(n-1,cutoff,frac)
    return found,nodes

def build_problem(n,d,a):
    types=sorted(set(a)-{0,d//2} if d%2==0 else set(a)-{0})
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    labels=sorted(set(tuple(sorted((a[i],a[j]))) for i,j in pairs if a[i] in types and a[j] in types))
    positions={p:i for i,p in enumerate(labels)}
    rows=[];rhs=[];vectors=[]; seen=set()
    def add(v):
        v=tuple(F(x) for x in v)
        if not any(v) or (sum(x*x for x in v)==1 and sum(bool(x) for x in v)==1):return
        den=math.lcm(*(x.denominator for x in v)); b=[int(x*den) for x in v]
        features=[0]*len(labels)
        for i,j in pairs:
            lab=tuple(sorted((a[i],a[j])))
            if lab in positions:features[positions[lab]]+=2*b[i]*b[j]
        row=tuple([-x for x in features]+[den*den]); val=sum(x*x for x in b)-den*den
        key=(row,val)
        if key not in seen:
            seen.add(key);rows.append(row);rhs.append(val);vectors.append([str(x) for x in v])
    for i,j in pairs:
        for sign in (-1,1):
            v=[0]*n;v[i]=1;v[j]=sign;add(v)
    for k in range(1,d):
        residues=[(k*x)%d for x in a]
        opts=[(r,r-d) if r else (0,) for r in residues]
        for b in itertools.product(*opts):add([F(x,d) for x in b])
    # This bound is also a consequence of the two-vector constraints.
    rows.append(tuple([0]*len(labels)+[1]));rhs.append(1);vectors.append(None)
    return labels, rows, rhs, vectors, add

def rational_dual(rows,rhs,result):
    support=[i for i,y in enumerate(result.ineqlin.marginals) if abs(y)>1e-8]
    A=sp.Matrix([rows[i] for i in support]).T
    target=sp.Matrix([0]*(len(rows[0])-1)+[1])
    sol,params=A.gauss_jordan_solve(target)
    if params.rows:sol=sol.subs({p:0 for p in params})
    weights=[F(int(x.p),int(x.q)) for x in sol]
    require((all(x>=0 for x in weights)), 'check_cyclic.py:83: verification check failed')
    require((all(sum(weights[k]*rows[i][j] for k,i in enumerate(support))==target[j] for j in range(len(rows[0])))), 'check_cyclic.py:84: verification check failed')
    bound=sum(weights[k]*rhs[i] for k,i in enumerate(support))
    require((bound<=0), 'check_cyclic.py:86: verification check failed')
    return {'bound':str(bound),'constraints':[[list(rows[i]),rhs[i],str(w)] for i,w in zip(support,weights) if w]}

def solve(n,d,a):
    labels,rows,rhs,vectors,add=build_problem(n,d,a)
    positions={p:i for i,p in enumerate(labels)}
    cosets=[tuple(F((k*x)%d,d) for x in a) for k in range(d)]
    for step in range(60):
        cost=np.zeros(len(labels)+1);cost[-1]=-1
        result=linprog(cost,A_ub=np.asarray(rows,dtype=float),b_ub=np.asarray(rhs,dtype=float),bounds=[(None,None)]*len(cost),method='highs',options={'time_limit':10})
        if not result.success:raise RuntimeError('LP did not solve: '+result.message)
        if result.x[-1]<1e-8:
            cert=rational_dual(rows,rhs,result)
            cert['vectors']=[vectors[i] for i,y in enumerate(result.ineqlin.marginals) if abs(y)>1e-8]
            return {'status':'excluded','steps':step+1,'certificate':cert}
        qs=[F(float(x)).limit_denominator(1000000) for x in result.x[:-1]]
        Q=[[F(i==j) if i==j else qs[positions[tuple(sorted((a[i],a[j])))]] if tuple(sorted((a[i],a[j]))) in positions else F(0) for j in range(n)] for i in range(n)]
        decomp=ldl(Q)
        if decomp is None:
            vals,vecs=np.linalg.eigh(np.array(Q,dtype=float)); v=vecs[:,0]
            added=False
            for scale in (8,32,128,512,2048,8192):
                w=[int(round(x*scale)) for x in v]
                if any(w) and qnorm(Q,w)<1:
                    before=len(rows);add(w);added=len(rows)>before;break
            if not added:raise RuntimeError('no negative-direction cut found')
            continue
        short,nodes=enumerate_short(Q,cosets)
        nonbasis=[v for v in short if not(sum(x*x for x in v)==1 and sum(bool(x) for x in v)==1)]
        if not nonbasis:
            return {'status':'accepted','steps':step+1,'nodes':nodes,'Q':[[str(x) for x in row] for row in Q]}
        before=len(rows)
        for v in nonbasis:add(v)
        if len(rows)==before:raise RuntimeError('rationalization failed to preserve strict constraints')
    raise RuntimeError('iteration limit')

if __name__=='__main__':
    n=int(sys.argv[1]); start=time.monotonic()
    data=json.loads((ROOT/'cyclic_candidates.json').read_text())[str(n)]
    result_path=ROOT/f'cyclic_results_{n}.json'
    results=json.loads(result_path.read_text()) if result_path.exists() else {}
    for ds,info in data.items():
        d=int(ds)
        for a in info['candidates']:
            key=str((d,tuple(a)))
            if key in results:continue
            if time.monotonic()-start>160:
                print('bounded batch complete',flush=True);sys.exit(0)
            try:res=solve(n,d,a)
            except Exception as exc:res={'status':'unresolved','error':repr(exc)}
            results[key]={'d':d,'a':a,**res}
            result_path.write_text(json.dumps(results))
            print(n,d,a,res['status'],res.get('steps'),flush=True)
    print('totals', {s:sum(r['status']==s for r in results.values()) for s in ('accepted','excluded','unresolved')},flush=True)


from validation import require
from check_cyclic import *

def build_code_problem(n,mod,gens,extra=None):
    from noncyclic import code_generated
    code=code_generated(gens,mod)
    cosets=[tuple(F(x,mod) for x in c) for c in code]
    extra=tuple(F(x) for x in extra) if extra is not None else None
    cols=[tuple(g[i] for g in gens)+( (() if extra is None else (extra[i],)) ) for i in range(n)]
    active=[not(all(2*g[i]%mod==0 for g in gens) and (extra is None or extra[i]==0)) for i in range(n)]
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    labels=sorted(set(tuple(sorted((cols[i],cols[j]))) for i,j in pairs if active[i] and active[j]))
    pos={p:i for i,p in enumerate(labels)}
    pairpos={(i,j):pos[tuple(sorted((cols[i],cols[j])))] for i,j in pairs if active[i] and active[j]}
    def feature(v):
        den=math.lcm(*(x.denominator for x in v)); b=[int(x*den) for x in v]
        fs=[0]*len(labels)
        for (i,j),p in pairpos.items():fs[p]+=2*b[i]*b[j]
        return den,b,fs
    rows=[];rhs=[];vectors=[];seen=set()
    def is_required(v):
        return (sum(x*x for x in v)==1 and sum(bool(x) for x in v)==1) or (extra is not None and (v==extra or v==tuple(-x for x in extra)))
    def add(v):
        v=tuple(F(x) for x in v)
        if not any(v) or is_required(v):return
        den,b,fs=feature(v)
        row=tuple([-x for x in fs]+[den*den]);val=sum(x*x for x in b)-den*den
        if (row,val) not in seen:
            seen.add((row,val));rows.append(row);rhs.append(val);vectors.append([str(x) for x in v])
    for i,j in pairs:
        for sign in (-1,1):
            v=[0]*n;v[i]=1;v[j]=sign;add(v)
    for c in code:
        if not any(c):continue
        for b in itertools.product(*[(r,r-mod) if r else (0,) for r in c]):add([F(x,mod) for x in b])
    rows.append(tuple([0]*len(labels)+[1]));rhs.append(1);vectors.append(None)
    eqs=[];eqrhs=[]
    if extra is not None:
        den,b,fs=feature(extra)
        eqs=[tuple(fs+[0])];eqrhs=[den*den-sum(x*x for x in b)]
    return code,cosets,labels,pairpos,rows,rhs,vectors,eqs,eqrhs,add,is_required

def dual_certificate(rows,rhs,vectors,eqs,eqrhs,result):
    inds=[i for i,y in enumerate(result.ineqlin.marginals) if abs(y)>1e-8]
    A=sp.Matrix([rows[i] for i in inds]+eqs).T
    target=sp.Matrix([0]*(len(rows[0])-1)+[1])
    sol,params=A.gauss_jordan_solve(target)
    if params.rows:sol=sol.subs({p:0 for p in params})
    ws=[F(int(x.p),int(x.q)) for x in sol]
    require((all(w>=0 for w in ws[:len(inds)])), 'check_code.py:49: verification check failed')
    allrows=[rows[i] for i in inds]+eqs;allrhs=[rhs[i] for i in inds]+eqrhs
    require((all(sum(w*r[j] for w,r in zip(ws,allrows))==target[j] for j in range(len(rows[0])))), 'check_code.py:51: verification check failed')
    bound=sum(w*b for w,b in zip(ws,allrhs));require((bound<=0), 'check_code.py:52: verification check failed')
    return {'bound':str(bound),'constraints':[{'row':list(rows[i]),'rhs':rhs[i],'weight':str(w),'vector':vectors[i]} for i,w in zip(inds,ws) if w],
        'equalities':[{'row':list(r),'rhs':b,'weight':str(w)} for r,b,w in zip(eqs,eqrhs,ws[len(inds):]) if w]}

def solve_code(n,mod,gens,extra=None):
    code,cosets,labels,pairpos,rows,rhs,vectors,eqs,eqrhs,add,is_required=build_code_problem(n,mod,gens,extra)
    for step in range(60):
        cost=np.zeros(len(labels)+1);cost[-1]=-1
        result=linprog(cost,A_ub=np.asarray(rows,dtype=float),b_ub=np.asarray(rhs,dtype=float),
            A_eq=np.asarray(eqs,dtype=float) if eqs else None,b_eq=np.asarray(eqrhs,dtype=float) if eqs else None,
            bounds=[(None,None)]*len(cost),method='highs',options={'time_limit':10})
        if not result.success:raise RuntimeError('LP failed: '+result.message)
        if result.x[-1]<1e-8:
            return {'status':'excluded','steps':step+1,'certificate':dual_certificate(rows,rhs,vectors,eqs,eqrhs,result)}
        qs=[F(float(x)).limit_denominator(1000000) for x in result.x[:-1]]
        if eqs:
            require((len(eqs)==1), 'check_code.py:68: verification check failed')
            row=eqs[0][:-1];pivot=next(i for i,x in enumerate(row) if x)
            qs[pivot]=(F(eqrhs[0])-sum(qs[i]*row[i] for i in range(len(qs)) if i!=pivot))/row[pivot]
        Q=[[F(i==j) for j in range(n)] for i in range(n)]
        for (i,j),p in pairpos.items():Q[i][j]=Q[j][i]=qs[p]
        if extra is not None:require((qnorm(Q,tuple(F(x) for x in extra))==1), 'check_code.py:73: verification check failed')
        if ldl(Q) is None:
            vals,vecs=np.linalg.eigh(np.array(Q,dtype=float));v=vecs[:,0];added=False
            for scale in (8,32,128,512,2048,8192):
                w=[int(round(x*scale)) for x in v]
                if any(w) and qnorm(Q,w)<1:
                    before=len(rows);add(w);added=len(rows)>before;break
            if not added:raise RuntimeError('negative-direction cut failed')
            continue
        short,nodes=enumerate_short(Q,cosets)
        bad=[v for v in short if not is_required(v)]
        if not bad:
            require((len(short)==2*n+(2 if extra is not None else 0)), 'check_code.py:85: verification check failed')
            return {'status':'accepted','steps':step+1,'nodes':nodes,'Q':[[str(x) for x in row] for row in Q]}
        before=len(rows)
        for v in bad:add(v)
        if len(rows)==before:raise RuntimeError('rationalization failed to preserve constraints')
    raise RuntimeError('iteration limit')

if __name__=='__main__':
    n=int(sys.argv[1]);data=json.loads((ROOT/f'noncyclic_candidates_{n}.json').read_text())
    results=[]
    for c in data['candidates']:
        try:r=solve_code(n,c['mod'],c['gens'])
        except Exception as exc:r={'status':'unresolved','error':repr(exc)}
        results.append({**c,**r});print(n,c['orders'],c['gens'],r['status'],r.get('error',''),flush=True)
    (ROOT/f'noncyclic_results_{n}.json').write_text(json.dumps(results))

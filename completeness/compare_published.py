"""Compare only after the independent enumeration; not a proof input."""

from validation import require
from verify_completeness import *
import urllib.request,ast,hashlib
COMMIT='d4c7b40ad2c95850d16160b347953dbaa9241486'

def inverse_det(R):
    n=len(R);A=[[F(x) for x in row]+[F(i==j) for j in range(n)] for i,row in enumerate(R)];det=F(1)
    for j in range(n):
        k=next((i for i in range(j,n) if A[i][j]),None)
        if k is None:return None,F(0)
        if k!=j:A[k],A[j]=A[j],A[k];det=-det
        pivot=A[j][j];det*=pivot;A[j]=[x/pivot for x in A[j]]
        for i in range(n):
            if i!=j:
                q=A[i][j];A[i]=[x-q*y for x,y in zip(A[i],A[j])]
    return [row[n:] for row in A],det

def basis_key(gens,n):
    code=fraction_code(gens,n);d=math.lcm(*(x.denominator for g in code for x in g))
    if len(code)==d:
        generators=[g for g in code if math.lcm(*(x.denominator for x in g))==d]
        return str(('cyclic',d,min(tuple(sorted(int(d*min(x,1-x)) for x in g)) for g in generators)))
    require((d in (2,4) and len(code)==2*d), 'compare_published.py:23: verification check failed')
    ic=tuple(sorted(tuple(int(x*d) for x in g) for g in code))
    return str(('noncyclic',d,canonical_code(ic,d,(d,2))))

def main():
    ours8={r['cell_key'] for r in read('extension_results_8.json').values() if r['status']=='accepted'}
    bases9=[{'mod':1,'gens':[[0]*9]}]+[{'mod':r['d'],'gens':[r['a']]} for r in read('cyclic_results_9.json').values() if r['status']=='accepted']+read('noncyclic_results_9.json')
    bases9=[r for r in bases9 if r.get('status','accepted')=='accepted']
    ours9={basis_key([tuple(F(x,b['mod']) for x in g) for g in b['gens']],9) for b in bases9}
    summary={}
    for n,ours in ((8,ours8),(9,ours9)):
        url=f'https://raw.githubusercontent.com/elbazvip/Voronoi-complexes-database/{COMMIT}/DATA_gap/ListSHV_n{n}_rnk9'
        with urllib.request.urlopen(url,timeout=15) as response:raw=response.read(1000000)
        source=raw.decode().strip();require((source.startswith('return')), 'compare_published.py:36: verification check failed')
        data=ast.literal_eval(source[6:].strip().rstrip(';'))
        keys=[]
        for V in data:
            require((len(V)==9 and all(len(v)==n and all(type(x)==int for x in v) for v in V)), 'compare_published.py:40: verification check failed')
            if n==9:
                inv,det=inverse_det(V);require((inv is not None), 'compare_published.py:42: verification check failed')
                require((len(fraction_code(inv,9))==abs(det)), 'compare_published.py:43: verification check failed')
                keys.append(basis_key(inv,9))
            else:
                choices=[]
                for omitted in range(9):
                    inv,det=inverse_det([v for i,v in enumerate(V) if i!=omitted])
                    choices.append((abs(det),omitted,inv))
                index,omitted,inv=max(choices,key=lambda t:t[0]);require((index>0), 'compare_published.py:50: verification check failed')
                v=tuple(sum(F(V[omitted][j])*inv[j][i] for j in range(8)) for i in range(8))
                require((all(abs(x)<=1 for x in v)), 'compare_published.py:52: verification check failed')
                require((len(fraction_code(inv,8))==index), 'compare_published.py:53: verification check failed')
                keys.append(canonical_cell(inv,v))
        require((len(set(keys))==len(keys)), 'compare_published.py:55: verification check failed')
        require((set(keys)==ours), lambda: ((n,len(ours-set(keys)),len(set(keys)-ours))))
        summary[n]={'published_cells':len(keys),'independently_generated_classes':len(ours),'missing':0,'extra':0,'commit':COMMIT,'source_sha256':hashlib.sha256(raw).hexdigest()}
        print(n,summary[n],flush=True)
    (ROOT/'comparison_summary.json').write_text(json.dumps(summary,indent=2))

if __name__=='__main__':main()

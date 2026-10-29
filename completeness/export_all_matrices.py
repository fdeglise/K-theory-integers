"""Export every retained matrix and explicit rational reconstructions.

All arithmetic is exact. Original solver iterates were not logged; this
export does not present reconstructed matrices as historical solver states.
"""

from validation import require
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import csv, hashlib, json, math, time
from matrix_examples import identity, transpose, multiply, determinant, rank, ldl, tensor_rows
from noncyclic import code_generated

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'matrices'

def read(name):return json.loads((ROOT/name).read_text())

def encode(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:encode(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [encode(v) for v in x]
    return x

def inverse(A):
    n=len(A);T=[[F(x) for x in row]+identity(n)[i] for i,row in enumerate(A)]
    for j in range(n):
        p=next(i for i in range(j,n) if T[i][j])
        T[j],T[p]=T[p],T[j]
        q=T[j][j];T[j]=[x/q for x in T[j]]
        for i in range(n):
            if i!=j:
                q=T[i][j];T[i]=[x-q*y for x,y in zip(T[i],T[j])]
    return [row[n:] for row in T]

def extended_gcd(a,b):
    aa,bb=abs(a),abs(b);x0,x1,y0,y1=1,0,0,1
    while bb:
        q=aa//bb;aa,bb=bb,aa-q*bb
        x0,x1=x1,x0-q*x1;y0,y1=y1,y0-q*y1
    return aa,x0*(1 if a>=0 else -1),y0*(1 if b>=0 else -1)

def lattice_basis(n,mod,gens):
    """Construct and certify a Z-basis of Z^n + sum Z*(g/mod)."""
    B=identity(n)
    for gen in gens:
        f=[[F(x,mod)] for x in gen]
        x=[row[0] for row in multiply(inverse(B),f)]
        d=math.lcm(*(t.denominator for t in x))
        if d==1:continue
        a=[int(d*t) for t in x];V=identity(n);current=a[:]
        for j in range(1,n):
            p,q=current[0],current[j]
            if p==q==0:continue
            g,s,t=extended_gcd(p,q)
            Rinv=identity(n)
            Rinv[0][0]=F(p//g);Rinv[0][j]=F(-t)
            Rinv[j][0]=F(q//g);Rinv[j][j]=F(s)
            V=multiply(V,Rinv);current[0]=g;current[j]=0
        require((math.gcd(current[0],d)==1), 'export_all_matrices.py:59: verification check failed')
        B=multiply(B,V)
        for i in range(n):B[i][0]/=d
    Bi=inverse(B)
    code={tuple(F(x,mod) for x in c) for c in code_generated(gens,mod)}
    require((abs(determinant(B))==F(1,len(code))), 'export_all_matrices.py:64: verification check failed')
    require((all(tuple(x%1 for x in col) in code for col in transpose(B))), 'export_all_matrices.py:65: verification check failed')
    require((all(x.denominator==1 for row in Bi for x in row)), 'export_all_matrices.py:66: verification check failed')
    for g in gens:require((all(row[0].denominator==1 for row in multiply(Bi,[[F(x,mod)] for x in g]))), 'export_all_matrices.py:67: verification check failed')
    return B,Bi

def sources():
    records=[]
    for n in (8,9):
        name=f'cyclic_results_{n}.json'
        data=read(name)
        for seq,(key,r) in enumerate(sorted(data.items(),key=lambda kv:(kv[1]['d'],kv[1]['a'])),1):
            records.append((f'cyclic_n{n}_{seq:04d}',name,key,n,'cyclic',r['d'],[r['a']],None,r))
        name=f'noncyclic_results_{n}.json'
        for seq,r in enumerate(read(name)):
            records.append((f'noncyclic_n{n}_{seq:04d}',name,str(seq),n,'noncyclic',r['mod'],r['gens'],None,r))
    bases=read('extension_candidates_8.json')['bases']
    name='extension_results_8.json'
    for key,r in sorted(read(name).items(),key=lambda kv:int(kv[0])):
        b=bases[r['basis']]
        records.append((f'extension_n8_{int(key):04d}',name,key,8,'extension',b['mod'],b['gens'],r['extra'],r))
    for n in (8,9):
        records.append((f'trivial_n{n}',None,None,n,'trivial',1,[[0]*n],None,{'status':'accepted','Q':identity(n),'steps':0}))
    return records

def parameter_layout(n,mod,gens,extra):
    v=None if extra is None else tuple(F(x) for x in extra)
    cols=[tuple(g[i] for g in gens)+(tuple() if v is None else (v[i],)) for i in range(n)]
    active=[not(all(2*g[i]%mod==0 for g in gens) and (v is None or v[i]==0)) for i in range(n)]
    pairs=[(i,j) for i in range(n) for j in range(i+1,n) if active[i] and active[j]]
    labels=sorted({tuple(sorted((cols[i],cols[j]))) for i,j in pairs})
    pos={lab:k for k,lab in enumerate(labels)}
    groups=[[] for _ in labels]
    for i,j in pairs:groups[pos[tuple(sorted((cols[i],cols[j])))]] .append([i+1,j+1])
    return [{'name':f'q{k+1}','coordinate_pairs':ps} for k,ps in enumerate(groups)]

def decode_certificate(r,phase):
    c=r['certificate']
    if phase=='cyclic':
        require((len(c['constraints'])==len(c['vectors'])), 'export_all_matrices.py:103: verification check failed')
        cons=[{'row':row,'rhs':rhs,'weight':weight,'vector':v}
              for (row,rhs,weight),v in zip(c['constraints'],c['vectors'])]
        eqs=[]
    else:cons=c['constraints'];eqs=c['equalities']
    A=[[F(x) for x in item['row']] for item in cons]
    b=[F(item['rhs']) for item in cons]
    lam=[F(item['weight']) for item in cons]
    E=[[F(x) for x in item['row']] for item in eqs]
    h=[F(item['rhs']) for item in eqs]
    mu=[F(item['weight']) for item in eqs]
    require((all(x>=0 for x in lam)), 'export_all_matrices.py:114: verification check failed')
    m=len(A[0]);total=[sum((w*row[j] for w,row in zip(lam,A)),F(0))+sum((w*row[j] for w,row in zip(mu,E)),F(0)) for j in range(m)]
    bound=sum((x*y for x,y in zip(lam,b)),F(0))+sum((x*y for x,y in zip(mu,h)),F(0))
    require((total==[F(0)]*(m-1)+[F(1)]), 'export_all_matrices.py:117: verification check failed')
    require((bound==F(c['bound']) and bound<=0), 'export_all_matrices.py:118: verification check failed')
    return {'origin':'retained exact certificate','A_ineq':A,'b_ineq':b,'lambda':lam,
            'E_eq':E,'h_eq':h,'mu':mu,'inequality_vectors':[item['vector'] for item in cons],
            'weighted_row':total,'weighted_rhs':bound}

def make_record(item):
    ident,name,key,n,phase,mod,gens,extra,r=item
    v=None if extra is None else tuple(F(x) for x in extra)
    code=code_generated(gens,mod)
    rec={'id':ident,'source_file':name,'source_entry':key,'rank':n,'phase':phase,'status':r['status'],
         'record_origin':'supplemental trivial case' if phase=='trivial' else 'retained result',
         'saved_search_iterations':r['steps'],'lattice_modulus':mod,'lattice_index':len(code),
         'lattice_generator_rows':gens,'extra_vector':v,'gram_parameters':parameter_layout(n,mod,gens,v),
         'lp_column_order':[p['name'] for p in parameter_layout(n,mod,gens,v)]+['delta']}
    S=identity(n)+([] if v is None else [list(v)])
    if r['status']=='excluded':rec['certificate']=decode_certificate(r,phase)
    else:
        Q=[[F(x) for x in row] for row in r['Q']]
        U,D=ldl(Q)
        require((all(sum((x[i]*Q[i][j]*x[j] for i in range(n) for j in range(n)),F(0))==1 for x in S)), 'export_all_matrices.py:137: verification check failed')
        B,Bi=lattice_basis(n,mod,gens)
        H=multiply(multiply(transpose(B),Q),B)
        scale=math.lcm(*(x.denominator for row in H for x in row))
        A=[[scale*x for x in row] for row in H]
        SV=multiply(S,transpose(Bi))
        require((all(x.denominator==1 for row in SV for x in row)), 'export_all_matrices.py:143: verification check failed')
        T=tensor_rows(S)
        rec['gram']={'origin':'supplemental identity' if phase=='trivial' else 'retained Gram matrix',
                     'Q':Q,'minimum':1,'derived_U':U,'derived_D':D,'det_Q':determinant(Q)}
        rec['derived_geometry']={'origin':'exact reconstruction from retained data','minimal_vector_rows':S,
            'tensor_column_pairs':[[i+1,j+1] for i in range(n) for j in range(i,n)],
            'tensor_matrix':T,'tensor_rank':rank(T),'integral_basis_B':B,'integral_basis_inverse':Bi,
            'gram_in_integral_basis':H,'integer_gram_scale':scale,'integer_gram_matrix':A,
            'minimal_vectors_in_integral_basis':SV}
    # Every extension has a saved swap. Accepted basis configurations have a
    # swap of equal quotient-code columns, reconstructed explicitly here.
    if phase=='extension':i,j=r['swap']
    elif r['status']=='accepted':
        i,j=next((i,j) for i in range(n) for j in range(i+1,n) if all(g[i]==g[j] for g in gens))
    else:i=j=None
    if i is not None:
        require((i!=j and all(g[i]==g[j] for g in gens)), 'export_all_matrices.py:159: verification check failed')
        require((v is None or v[i]==v[j]), 'export_all_matrices.py:160: verification check failed')
        P=identity(n);P[i],P[j]=P[j],P[i]
        R=identity(len(S));R[i],R[j]=R[j],R[i]
        require((determinant(P)==-1 and determinant(R)==-1), 'export_all_matrices.py:163: verification check failed')
        require((multiply(S,transpose(P))==multiply(R,S)), 'export_all_matrices.py:164: verification check failed')
        auto={'origin':'expanded saved swap' if phase=='extension' else 'reconstructed equal-column swap',
              'swapped_coordinates':[i+1,j+1],'P_in_minimal_basis':P,'P_on_rays':R,
              'det_P':-1,'det_on_rays':-1,
              'interpretation':'symmetry of the prescribed configuration; does not assert realizability'}
        if r['status']=='accepted':
            require((multiply(multiply(transpose(P),Q),P)==Q), 'export_all_matrices.py:170: verification check failed')
            Pint=multiply(multiply(Bi,P),B)
            require((all(x.denominator==1 for row in Pint for x in row)), 'export_all_matrices.py:172: verification check failed')
            require((determinant(Pint)==-1), 'export_all_matrices.py:173: verification check failed')
            require((multiply(multiply(transpose(Pint),A),Pint)==A), 'export_all_matrices.py:174: verification check failed')
            auto['P_in_integral_basis']=Pint
        rec['symmetry']=auto
    if phase=='extension' and r['status']=='accepted':rec['cell_equivalence_key']=r['cell_key']
    return rec

def legacy_boundary():
    d=[[F(x) for x in (0,0,0,0,-1,0,1)]]
    e=identity(7)
    cols=[e[6],e[0],e[1],e[2],e[3],e[5],[e[4][i]+e[6][i] for i in range(7)]]
    U=transpose(cols)
    require((abs(determinant(U))==1), 'export_all_matrices.py:185: verification check failed')
    smith=multiply(d,U)
    require((smith==[[F(1)]+[F(0)]*6]), 'export_all_matrices.py:187: verification check failed')
    return {'origin':'rank-five boundary row recorded earlier in the conversation',
            'd9':d,'derived_unimodular_column_change':U,'derived_Smith_form':smith,
            'note':'The column change was reconstructed for this export.'}

def write_readable(data):
    with (OUT/'all_matrices.txt').open('w') as out:
        def line(x=''):out.write(str(x)+'\n')
        def mat(name,A):
            line(name)
            if not A:line('[]');return
            if not isinstance(A[0],list):A=[A]
            width=max(len(str(x)) for row in A for x in row)
            for row in A:line('['+' '.join(str(x).rjust(width) for x in row)+']')
        line('ALL RETAINED VORONOI MATRIX DATA AND LABELLED RECONSTRUCTIONS')
        line('All indices in coordinate pairs are one-based; fractions are exact.')
        line('L = Z^n + sum Z*(generator_row/modulus). Q is in the e_i coordinates.')
        line('The final LP variable is delta. See README.md for provenance and limitations.')
        for rec in data['records']:
            line();line('='*72);line(rec['id']+' | '+rec['status']+' | rank '+str(rec['rank']))
            line('Source: '+str(rec['source_file'])+' entry '+str(rec['source_entry']))
            line('Lattice modulus: '+str(rec['lattice_modulus'])+'; index: '+str(rec['lattice_index']))
            mat('Lattice generator rows:',rec['lattice_generator_rows'])
            if rec['extra_vector'] is not None:mat('Extra vector:',rec['extra_vector'])
            line('LP columns: '+str(rec['lp_column_order']))
            for p in rec['gram_parameters']:line(p['name']+' occupies Gram pairs '+str(p['coordinate_pairs']))
            if 'gram' in rec:
                g=rec['gram'];geo=rec['derived_geometry']
                mat(g['origin']+' Q (minimum 1):',g['Q'])
                mat('Derived upper-triangular U in Q=U^t diag(D) U:',g['derived_U'])
                mat('Derived D:',g['derived_D']);line('det(Q): '+g['det_Q'])
                for name in ('minimal_vector_rows','tensor_matrix','integral_basis_B','integral_basis_inverse',
                             'gram_in_integral_basis','integer_gram_matrix','minimal_vectors_in_integral_basis'):
                    mat('Derived '+name+':',geo[name])
                line('Tensor rank: '+str(geo['tensor_rank']))
                line('Integer Gram scale and minimum: '+str(geo['integer_gram_scale']))
            if 'certificate' in rec:
                c=rec['certificate']
                for name in ('A_ineq','b_ineq','lambda','E_eq','h_eq','mu','weighted_row'):
                    mat('Retained certificate '+name+':',c[name])
                line('weighted_rhs: '+c['weighted_rhs'])
                line('Inequality vectors (null means delta<=1): '+str(c['inequality_vectors']))
            if 'symmetry' in rec:
                a=rec['symmetry'];line('Symmetry provenance: '+a['origin'])
                for name in ('P_in_minimal_basis','P_on_rays','P_in_integral_basis'):
                    if name in a:mat(name+':',a[name])
                line('Determinant on rays: -1')
        line();line('RANK-FIVE BOUNDARY RECORDED IN THE CONVERSATION')
        for name in ('d9','derived_unimodular_column_change','derived_Smith_form'):
            mat(name+':',data['rank5_boundary'][name])

def main():
    start=time.monotonic();OUT.mkdir(exist_ok=True)
    src=sources();records=[]
    for i,item in enumerate(src):
        records.append(encode(make_record(item)))
        if (i+1)%200==0:print('Exported',i+1,'of',len(src),'records',flush=True)
    # Index all 150 final cell classes, retaining every alternate presentation.
    eight={}
    for r in records:
        if r['phase']=='extension' and r['status']=='accepted':eight.setdefault(r['cell_equivalence_key'],[]).append(r)
    cells=[]
    for i,(key,rs) in enumerate(sorted(eight.items()),1):
        cid=f'n8_cell_{i:03d}'
        for r in rs:r['cell_class_id']=cid
        cells.append({'cell_class_id':cid,'rank':8,'representative_record':rs[0]['id'],
                      'all_presentation_records':';'.join(r['id'] for r in rs)})
    nine=[r for r in records if r['rank']==9 and r['status']=='accepted']
    for i,r in enumerate(nine,1):
        cid=f'n9_cell_{i:03d}';r['cell_class_id']=cid
        cells.append({'cell_class_id':cid,'rank':9,'representative_record':r['id'],'all_presentation_records':r['id']})
    require((len(eight)==106 and len(nine)==44), 'export_all_matrices.py:258: verification check failed')
    source_names=sorted({item[1] for item in src if item[1] is not None}|{'extension_candidates_8.json'})
    source_hashes={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in source_names}
    stats={'retained_result_records':sum(r['record_origin']=='retained result' for r in records),
           'retained_Gram_matrices':sum('gram' in r and r['phase']!='trivial' for r in records),
           'supplemental_identity_Gram_matrices':2,
           'retained_exclusion_certificates':sum('certificate' in r for r in records),
           'explicit_configuration_symmetries':sum('symmetry' in r for r in records),
           'accepted_presentations_with_integral_coordinates':sum('derived_geometry' in r for r in records),
           'rank8_final_cell_classes':len(eight),'rank9_final_cell_classes':len(nine)}
    require((stats['retained_result_records']==964 and stats['retained_Gram_matrices']==189), 'export_all_matrices.py:268: verification check failed')
    require((stats['retained_exclusion_certificates']==775), 'export_all_matrices.py:269: verification check failed')
    data={'format_version':1,'coordinate_convention':'Column vectors; Q[x]=x^t Q x; generator rows give g/modulus.',
          'provenance_note':'All retained Gram/certificate matrices are copied exactly. Other matrices are explicitly reconstructed. Original transient solver states were not retained.',
          'statistics':stats,'source_sha256':source_hashes,'rank5_boundary':encode(legacy_boundary()),'records':records}
    (OUT/'all_matrices.json').write_text(json.dumps(data,indent=2)+'\n')
    write_readable(data)
    with (OUT/'matrix_index.csv').open('w',newline='') as f:
        fields=['id','rank','phase','status','cell_class_id','source_file','source_entry','lattice_index','lattice_modulus','saved_search_iterations','gram_rows','certificate_inequalities','certificate_equalities']
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader()
        for r in records:
            row={k:r.get(k,'') for k in fields[:10]}
            row.update(saved_search_iterations=r['saved_search_iterations'],gram_rows=r['rank'] if 'gram' in r else 0,
                certificate_inequalities=len(r.get('certificate',{}).get('A_ineq',[])),certificate_equalities=len(r.get('certificate',{}).get('E_eq',[])))
            writer.writerow(row)
    with (OUT/'cell_index.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=['cell_class_id','rank','representative_record','all_presentation_records'])
        writer.writeheader();writer.writerows(cells)
    (OUT/'export_summary.json').write_text(json.dumps({'statistics':stats,'seconds':round(time.monotonic()-start,3)},indent=2)+'\n')
    print(json.dumps(stats,indent=2),flush=True)
    print('Seconds:',round(time.monotonic()-start,3),flush=True)

if __name__=='__main__':main()

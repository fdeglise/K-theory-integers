"""Check source fidelity and every newly exported matrix identity exactly.

This complements verify_completeness.py, which checks the enumeration and
shortest-vector certificates. Neither program is a Lean formalization.
"""

from validation import require
from validation import save_summary
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import csv, gzip, hashlib, json, time
from export_all_matrices import sources
from matrix_examples import identity, transpose, multiply, determinant, rank, tensor_rows
from noncyclic import code_generated
from verify_completeness import verify_result

ROOT = Path(__file__).resolve().parent
OUT = ROOT/'matrices'


def matrix(A):
    return [[F(x) for x in row] for row in A]


def main():
    start = time.monotonic()
    document = json.loads((OUT/'all_matrices.json').read_text())
    records = document['records']
    index = {r['id']: r for r in records}
    original = {item[0]: item for item in sources()}
    require((len(index) == len(records) == len(original) == 966), 'verify_matrix_export.py:29: verification check failed')
    require((set(index) == set(original)), 'verify_matrix_export.py:30: verification check failed')
    for name, expected in document['source_sha256'].items():
        require((hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == expected), 'verify_matrix_export.py:32: verification check failed')
    counts = Counter()
    for ident, item in original.items():
        _, name, key, n, phase, mod, gens, extra, result = item
        rec = index[ident]
        require(((rec['source_file'], rec['source_entry'], rec['rank'], rec['phase'],
                rec['lattice_modulus'], rec['lattice_generator_rows'], rec['status']) == (
                name, key, n, phase, mod, gens, result['status'])), 'verify_matrix_export.py:37: verification check failed')
        S = identity(n)
        if extra is not None:
            require((list(map(F, rec['extra_vector'])) == list(map(F, extra))), 'verify_matrix_export.py:42: verification check failed')
            S.append(list(map(F, extra)))
        else:
            require((rec['extra_vector'] is None), 'verify_matrix_export.py:45: verification check failed')
        code = {tuple(F(x, mod) for x in c) for c in code_generated(gens, mod)}
        require((rec['lattice_index'] == len(code)), 'verify_matrix_export.py:47: verification check failed')
        if result['status'] == 'accepted':
            gram, geo = rec['gram'], rec['derived_geometry']
            Q = matrix(gram['Q'])
            require((Q == matrix(result['Q'])), 'verify_matrix_export.py:51: verification check failed')
            require((Q == transpose(Q)), 'verify_matrix_export.py:52: verification check failed')
            U, ds = matrix(gram['derived_U']), list(map(F, gram['derived_D']))
            require((len(ds) == n and all(d > 0 for d in ds)), 'verify_matrix_export.py:54: verification check failed')
            require((all(U[i][i] == 1 and all(U[i][j] == 0 for j in range(i)) for i in range(n))), 'verify_matrix_export.py:55: verification check failed')
            D = [[ds[i] if i == j else F(0) for j in range(n)] for i in range(n)]
            require((multiply(multiply(transpose(U), D), U) == Q), 'verify_matrix_export.py:57: verification check failed')
            require((determinant(Q) == F(gram['det_Q'])), 'verify_matrix_export.py:58: verification check failed')
            B, Bi = matrix(geo['integral_basis_B']), matrix(geo['integral_basis_inverse'])
            require((multiply(B, Bi) == multiply(Bi, B) == identity(n)), 'verify_matrix_export.py:60: verification check failed')
            require((abs(determinant(B)) == F(1, len(code))), 'verify_matrix_export.py:61: verification check failed')
            require((all(tuple(x % 1 for x in col) in code for col in transpose(B))), 'verify_matrix_export.py:62: verification check failed')
            require((all(x.denominator == 1 for row in Bi for x in row)), 'verify_matrix_export.py:63: verification check failed')
            H = matrix(geo['gram_in_integral_basis'])
            A = matrix(geo['integer_gram_matrix'])
            scale = geo['integer_gram_scale']
            require((isinstance(scale, int) and scale > 0), 'verify_matrix_export.py:67: verification check failed')
            require((multiply(multiply(transpose(B), Q), B) == H), 'verify_matrix_export.py:68: verification check failed')
            require((A == [[scale*x for x in row] for row in H]), 'verify_matrix_export.py:69: verification check failed')
            require((all(x.denominator == 1 for row in A for x in row)), 'verify_matrix_export.py:70: verification check failed')
            SV = matrix(geo['minimal_vectors_in_integral_basis'])
            require((SV == multiply(S, transpose(Bi))), 'verify_matrix_export.py:72: verification check failed')
            require((all(x.denominator == 1 for row in SV for x in row)), 'verify_matrix_export.py:73: verification check failed')
            require((matrix(geo['minimal_vector_rows']) == S), 'verify_matrix_export.py:74: verification check failed')
            for v in S:
                require((multiply(multiply([v], Q), transpose([v])) == [[F(1)]]), 'verify_matrix_export.py:76: verification check failed')
            T = matrix(geo['tensor_matrix'])
            require((geo['tensor_column_pairs'] == [[i+1, j+1] for i in range(n) for j in range(i, n)]), 'verify_matrix_export.py:78: verification check failed')
            require((T == tensor_rows(S)), 'verify_matrix_export.py:79: verification check failed')
            require((rank(T) == geo['tensor_rank'] == len(S)), 'verify_matrix_export.py:80: verification check failed')
            counts['Gram_matrices_and_derived_geometry_checked'] += 1
        else:
            c = rec['certificate']
            saved = result['certificate']
            if phase == 'cyclic':
                raw = [{'row': a, 'rhs': b, 'weight': w, 'vector': v}
                       for (a, b, w), v in zip(saved['constraints'], saved['vectors'])]
                equalities = []
            else:
                raw, equalities = saved['constraints'], saved['equalities']
            require((matrix(c['A_ineq']) == matrix([s['row'] for s in raw])), 'verify_matrix_export.py:91: verification check failed')
            require((list(map(F, c['b_ineq'])) == [F(s['rhs']) for s in raw]), 'verify_matrix_export.py:92: verification check failed')
            require((list(map(F, c['lambda'])) == [F(s['weight']) for s in raw]), 'verify_matrix_export.py:93: verification check failed')
            require((c['inequality_vectors'] == [s['vector'] for s in raw]), 'verify_matrix_export.py:94: verification check failed')
            require((matrix(c['E_eq']) == matrix([s['row'] for s in equalities])), 'verify_matrix_export.py:95: verification check failed')
            require((list(map(F, c['h_eq'])) == [F(s['rhs']) for s in equalities]), 'verify_matrix_export.py:96: verification check failed')
            require((list(map(F, c['mu'])) == [F(s['weight']) for s in equalities]), 'verify_matrix_export.py:97: verification check failed')
            require((F(c['weighted_rhs']) == F(saved['bound'])), 'verify_matrix_export.py:98: verification check failed')
            verify_result(n, mod, gens, result, extra=extra, cyclic_format=(phase == 'cyclic'))
            counts['exclusion_certificates_checked'] += 1
        if 'symmetry' in rec:
            a = rec['symmetry']
            P, R = matrix(a['P_in_minimal_basis']), matrix(a['P_on_rays'])
            require((determinant(P) == determinant(R) == -1), 'verify_matrix_export.py:104: verification check failed')
            require((multiply(S, transpose(P)) == multiply(R, S)), 'verify_matrix_export.py:105: verification check failed')
            for g in gens:
                require((multiply(P, transpose([g])) == transpose([g])), 'verify_matrix_export.py:107: verification check failed')
            if result['status'] == 'accepted':
                Pi = matrix(a['P_in_integral_basis'])
                require((multiply(multiply(Bi, P), B) == Pi), 'verify_matrix_export.py:110: verification check failed')
                require((all(x.denominator == 1 for row in Pi for x in row)), 'verify_matrix_export.py:111: verification check failed')
                require((determinant(Pi) == -1), 'verify_matrix_export.py:112: verification check failed')
                require((multiply(multiply(transpose(P), Q), P) == Q), 'verify_matrix_export.py:113: verification check failed')
                require((multiply(multiply(transpose(Pi), A), Pi) == A), 'verify_matrix_export.py:114: verification check failed')
            counts['configuration_symmetries_checked'] += 1
    with (OUT/'matrix_index.csv').open(newline='') as f:
        rows = list(csv.DictReader(f))
    require((len(rows) == 966 and {r['id'] for r in rows} == set(index)), 'verify_matrix_export.py:118: verification check failed')
    for row in rows:
        r = index[row['id']]
        require((int(row['rank']) == r['rank'] and row['status'] == r['status']), 'verify_matrix_export.py:121: verification check failed')
        require((int(row['gram_rows']) == (r['rank'] if 'gram' in r else 0)), 'verify_matrix_export.py:122: verification check failed')
        require((int(row['certificate_inequalities']) == len(r.get('certificate', {}).get('A_ineq', []))), 'verify_matrix_export.py:123: verification check failed')
    with (OUT/'cell_index.csv').open(newline='') as f:
        cells = list(csv.DictReader(f))
    require((Counter(row['rank'] for row in cells) == {'8': 106, '9': 44}), 'verify_matrix_export.py:126: verification check failed')
    memberships = set()
    for cell in cells:
        ids = cell['all_presentation_records'].split(';')
        require((cell['representative_record'] in ids), 'verify_matrix_export.py:130: verification check failed')
        for ident in ids:
            require((ident not in memberships), 'verify_matrix_export.py:132: verification check failed')
            memberships.add(ident)
            require((index[ident]['cell_class_id'] == cell['cell_class_id']), 'verify_matrix_export.py:134: verification check failed')
            require((index[ident]['status'] == 'accepted'), 'verify_matrix_export.py:135: verification check failed')
    require((memberships == {r['id'] for r in records if 'cell_class_id' in r}), 'verify_matrix_export.py:136: verification check failed')
    legacy = document['rank5_boundary']
    d, U = matrix(legacy['d9']), matrix(legacy['derived_unimodular_column_change'])
    require((abs(determinant(U)) == 1), 'verify_matrix_export.py:139: verification check failed')
    require((multiply(d, U) == matrix(legacy['derived_Smith_form']) == [[F(1)]+[F(0)]*6]), 'verify_matrix_export.py:140: verification check failed')
    counts['final_cell_classes_indexed'] = len(cells)
    counts['legacy_Smith_identity_checked'] = 1
    summary = {'status': 'passed', 'checks': dict(counts),
               'seconds': round(time.monotonic()-start, 3),
               'scope': 'Exact source fidelity and exported matrix identities; run verify_completeness.py for enumeration and shortest-vector verification.'}
    save_summary(OUT/'matrix_verification_summary.json', json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == '__main__':
    main()

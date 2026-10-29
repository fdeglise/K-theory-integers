"""Display and check actual matrix examples from the saved enumeration.

Uses exact rational arithmetic. This is not a Lean formalization.
"""

from validation import require
import json
from pathlib import Path
from fractions import Fraction as F
from verify_completeness import verify_form, basis_vectors, gram_value

ROOT = Path(__file__).resolve().parent

def identity(n):
    return [[F(i == j) for j in range(n)] for i in range(n)]

def transpose(A):
    return [list(row) for row in zip(*A)]

def multiply(A, B):
    require((len(A[0]) == len(B)), 'matrix_examples.py:19: verification check failed')
    return [[sum((A[i][k] * B[k][j] for k in range(len(B))), F(0))
             for j in range(len(B[0]))] for i in range(len(A))]

def determinant(A):
    A = [[F(x) for x in row] for row in A]
    n = len(A)
    result = F(1)
    for j in range(n):
        pivot = next((i for i in range(j, n) if A[i][j]), None)
        if pivot is None:
            return F(0)
        if pivot != j:
            A[j], A[pivot] = A[pivot], A[j]
            result = -result
        result *= A[j][j]
        for i in range(j + 1, n):
            c = A[i][j] / A[j][j]
            A[i] = [x - c * y for x, y in zip(A[i], A[j])]
    return result

def rank(A):
    A = [[F(x) for x in row] for row in A]
    r = 0
    for j in range(len(A[0])):
        pivot = next((i for i in range(r, len(A)) if A[i][j]), None)
        if pivot is None:
            continue
        A[r], A[pivot] = A[pivot], A[r]
        for i in range(r + 1, len(A)):
            c = A[i][j] / A[r][j]
            A[i] = [x - c * y for x, y in zip(A[i], A[r])]
        r += 1
        if r == len(A):
            break
    return r

def ldl(A):
    n = len(A)
    U = identity(n)
    D = []
    for i in range(n):
        d = A[i][i] - sum((D[k] * U[k][i] ** 2 for k in range(i)), F(0))
        require((d > 0), 'matrix_examples.py:62: verification check failed')
        D.append(d)
        for j in range(i + 1, n):
            U[i][j] = (A[i][j] - sum((D[k] * U[k][i] * U[k][j] for k in range(i)), F(0))) / d
    diagonal = [[D[i] if i == j else F(0) for j in range(n)] for i in range(n)]
    require((multiply(multiply(transpose(U), diagonal), U) == A), 'matrix_examples.py:67: verification check failed')
    return U, D

def tensor_rows(vectors):
    n = len(vectors[0])
    return [[v[i] * v[j] for i in range(n) for j in range(i, n)] for v in vectors]

def main():
    lines = []
    def out(text=''):
        lines.append(str(text))
    def matrix(label, A):
        out(label)
        width = max(len(str(x)) for row in A for x in row)
        for row in A:
            out('[' + ' '.join(str(x).rjust(width) for x in row) + ']')
        out()

    results = json.loads((ROOT / 'extension_results_8.json').read_text())
    accepted = results['10']
    require((accepted['status'] == 'accepted'), 'matrix_examples.py:87: verification check failed')
    Q = [[F(x) for x in row] for row in accepted['Q']]
    f = tuple(F(x) for x in accepted['extra'])
    require((f == (F(0), F(0), F(0), F(1, 2), F(1, 2), F(1, 2), F(1, 2), F(1, 2))), 'matrix_examples.py:90: verification check failed')
    B = identity(8)
    for i in range(8):
        B[i][7] = f[i]
    require((determinant(B) == F(1, 2)), 'matrix_examples.py:94: verification check failed')
    H = multiply(multiply(transpose(B), Q), B)
    A = [[20 * x for x in row] for row in H]
    require((all(x.denominator == 1 for row in A for x in row)), 'matrix_examples.py:97: verification check failed')
    w = (F(0), F(0), F(0), F(-1), F(-1), F(-1), F(-1), F(2))
    require((multiply(B, [[x] for x in w]) == [[F(i == 7)] for i in range(8)]), 'matrix_examples.py:99: verification check failed')
    required = basis_vectors(8) | {w, tuple(-x for x in w)}
    nodes = verify_form(H, {tuple([F(0)] * 8)}, required)
    _, D = ldl(A)
    standard = [tuple(F(i == j) for i in range(8)) for j in range(8)]
    vectors = standard + [w]
    T = tensor_rows(vectors)
    require((rank(T) == 9), 'matrix_examples.py:106: verification check failed')
    P = identity(8)
    P[3], P[4] = P[4], P[3]
    require((determinant(P) == -1), 'matrix_examples.py:109: verification check failed')
    require((multiply(P, P) == identity(8)), 'matrix_examples.py:110: verification check failed')
    require((multiply(multiply(transpose(P), A), P) == A), 'matrix_examples.py:111: verification check failed')
    perm = list(range(9))
    perm[3], perm[4] = perm[4], perm[3]
    require((all(multiply(P, [[x] for x in v]) == [[x] for x in vectors[perm[i]]]
               for i, v in enumerate(vectors))), 'matrix_examples.py:114: verification check failed')
    out('RANK EIGHT: accepted extension entry 10')
    out('Original lattice: L = Z^8 + Z*(e4+e5+e6+e7+e8)/2.')
    matrix('Saved Gram matrix Q:', Q)
    matrix('Integral lattice basis B = [e1,...,e7,f], expressed in the original coordinates:', B)
    out('det(B) = ' + str(determinant(B)))
    out('Since e8 = 2*f-e4-e5-e6-e7, the columns of B generate exactly L.')
    matrix('Integer Gram matrix A = 20 * transpose(B) * Q * B:', A)
    out('Positive LDL diagonal for A: ' + str([str(x) for x in D]))
    out('det(A) = ' + str(determinant(A)))
    matrix('Nine minimal-vector representatives in the integral lattice basis:', vectors)
    out('Their squared lengths for A: ' + str([int(gram_value(A, v)) for v in vectors]))
    out('Exact enumeration: these 18 signed vectors are all nonzero vectors with A[x] <= 20.')
    out('Enumeration nodes: ' + str(nodes))
    matrix('Rows of the tensor matrix (upper-triangular entries of v*v^t):', T)
    out('Tensor-matrix rank: ' + str(rank(T)))
    matrix('Orientation-reversing integral automorphism P, swapping coordinates 4 and 5:', P)
    out('Checked: transpose(P)*A*P = A; P*P = I; det(P) = -1.')
    out('Action on the nine rays: swap rays 4 and 5, fix the other seven.')
    out()

    excluded = results['7']
    require((excluded['status'] == 'excluded'), 'matrix_examples.py:137: verification check failed')
    cert = excluded['certificate']
    inequalities = cert['constraints']
    equalities = cert['equalities']
    rows = inequalities + equalities
    C = [[F(x) for x in row['row']] for row in rows]
    b = [[F(row['rhs'])] for row in rows]
    weights = [[F(row['weight']) for row in rows]]
    require((all(F(row['weight']) >= 0 for row in inequalities)), 'matrix_examples.py:145: verification check failed')
    require((multiply(weights, C) == [[F(0), F(1)]]), 'matrix_examples.py:146: verification check failed')
    require((multiply(weights, b) == [[F(0)]]), 'matrix_examples.py:147: verification check failed')
    out('EXCLUSION: extension entry 7, in variables (q, delta)')
    out('Candidate extra vector: e7+e8 in the same original lattice L.')
    out('First row is an inequality; second row is an equality.')
    matrix('Coefficient rows C:', C)
    matrix('Right-hand sides b:', b)
    matrix('Certificate weights:', weights)
    matrix('weights * C:', multiply(weights, C))
    matrix('weights * b:', multiply(weights, b))
    out('Therefore delta <= 0, contradicting the positive gap required by this candidate.')
    out()

    cyclic9 = json.loads((ROOT / 'cyclic_results_9.json').read_text())
    nine = next(r for r in cyclic9.values() if r['status'] == 'accepted' and r['d'] == 2
                and r['a'] == [0, 0, 0, 0, 1, 1, 1, 1, 1])
    Q9 = [[F(x) for x in row] for row in nine['Q']]
    require((Q9 == identity(9)), 'matrix_examples.py:163: verification check failed')
    f9 = tuple(F(x, 2) for x in nine['a'])
    required9 = basis_vectors(9)
    nodes9 = verify_form(Q9, {tuple([F(0)] * 9), f9}, required9)
    P9 = identity(9)
    P9[4], P9[5] = P9[5], P9[4]
    require((determinant(P9) == -1), 'matrix_examples.py:169: verification check failed')
    require((multiply(multiply(transpose(P9), Q9), P9) == Q9), 'matrix_examples.py:170: verification check failed')
    require((multiply(P9, [[x] for x in f9]) == [[x] for x in f9]), 'matrix_examples.py:171: verification check failed')
    out('RANK NINE: cyclic code (0,0,0,0,1,1,1,1,1) modulo 2')
    out('L = Z^9 + Z*(e5+e6+e7+e8+e9)/2.')
    matrix('Saved Gram matrix Q9:', Q9)
    out('Exact enumeration: the minimal vectors are precisely +/- e1,...,+/- e9.')
    out('The fractional coset has squared length at least 5/4.')
    out('Enumeration nodes: ' + str(nodes9))
    matrix('Checked integral coordinate permutation P9:', P9)
    out('Swapping coordinates 5 and 6 preserves L and Q9 and transposes two of the nine rays.')
    out()
    out('All displayed calculations use exact fractions. No Lean proof is claimed.')
    text = '\n'.join(lines) + '\n'
    (ROOT / 'matrix_examples.txt').write_text(text)
    print(text)

if __name__ == '__main__':
    main()

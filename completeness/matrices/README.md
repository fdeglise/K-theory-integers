# Complete matrix export

This directory contains every retained Gram matrix and exclusion certificate
from the rank-eight/rank-nine enumeration, together with explicitly labelled
exact reconstructions. Fractions are exact rational numbers, never decimal
approximations.

## Files and counts

- `all_matrices.txt`: readable matrix catalogue, ordered by candidate ID.
- `all_matrices.json`: the same matrices in structured form, with source
  filenames, original entry keys, coordinate conventions, and source hashes.
  Rational entries are strings such as `"-1/20"`; use exact rational parsing.
- `matrix_index.csv`: index of all 966 records: 964 retained result records
  and two supplemental trivial-quotient records.
- `cell_index.csv`: one row for each of the 106 rank-eight and 44 rank-nine
  final cell classes, identifying every retained presentation of that class.
- `initial_lp_matrices.jsonl.gz`: the 964 reconstructed initial linear programs,
  one JSON object per decompressed line. These contain 665,214 inequality rows
  and 647 equality rows. All matrix and vector entries in these systems are
  integers.
- `initial_lp_index.csv`: line number and matrix dimensions of each initial
  linear program.
- `*_summary.json`: export counts and recorded verification outcomes.

The main catalogue contains:

| Data | Number |
| --- | ---: |
| Retained accepted Gram matrices | 189 |
| Supplemental identity Gram matrices for trivial quotients | 2 |
| Retained exclusion certificates | 775 |
| Accepted presentations with explicit integral lattice bases | 191 |
| Explicit configuration symmetries | 704 |

The 191 accepted presentations consist of 13 rank-eight configurations with
eight minimal pairs, 44 rank-nine cells, and 134 presentations of rank-eight
cells. The last 134 presentations represent 106 distinct classes. Thus the
189 retained Gram matrices should not be interpreted as 189 distinct
degree-eight cell classes.

## Coordinates and Gram matrices

Vectors are columns and squared length is `x^t Q x`. If the modulus is `d`
and the generator rows are `g_1,...,g_r`, the ambient lattice is

```text
L = Z^n + Z*(g_1/d) + ... + Z*(g_r/d).
```

The original coordinate vectors `e_i` are independent minimal vectors. They
need not be an integral basis of `L`. Every saved `Q` has minimum 1 on `L`.
In an extension record, `extra_vector` supplies the ninth minimal pair.

For each accepted presentation the export also gives:

- `derived_U`, `derived_D`: `Q = U^t diag(D) U`, with all entries of `D` positive.
- `minimal_vector_rows`: one representative of each minimal pair.
- `tensor_matrix`: rows obtained from `v*v^t` by listing entries `(i,j)` in
  lexicographic order with `i <= j`. Off-diagonal entries are not doubled.
  The column pairs are also supplied explicitly. The sizes are 8 by 36 for
  rank-eight basis configurations, 9 by 36 for rank-eight extensions, and
  9 by 45 for rank-nine cells.
- `integral_basis_B`: columns forming a Z-basis of `L`, in the original
  coordinates; `integral_basis_inverse` gives the inverse change of coordinates.
- `gram_in_integral_basis`: `H = B^t Q B`.
- `integer_gram_matrix`: `A = s*H`, where `s = integer_gram_scale` is a positive
  integer. The minimum of `A` on `Z^n` is `s`.
- `minimal_vectors_in_integral_basis`: integral coordinate rows of the same
  minimal vectors.

The newly constructed integral bases are choices; they need not coincide
with the basis used in the separate worked examples.

## Symmetry matrices

`P_in_minimal_basis` swaps two of the original coordinates. For accepted
presentations, `P_in_integral_basis = B^(-1) P B` is an integral isometry.
`P_on_rays` gives the permutation of the listed rank-one rays. Its determinant
is -1.

There is also a symmetry for every rejected rank-eight extension: it preserves
the prescribed lattice and vector configuration. This does not assert that
the rejected configuration is realizable. Only the 191 accepted presentations
have Gram matrices and the associated integral-coordinate isometries.

## Exclusion-certificate matrices

For each candidate, `gram_parameters` states precisely which off-diagonal
entries of `Q` contain each parameter `q_i`. The diagonal entries are 1 and
unlisted off-diagonal entries are 0. Coordinate pairs are one-based.

Write `z = (q_1,...,q_m,delta)^t`. The certificate arrays mean

```text
A_ineq * z <= b_ineq,
E_eq   * z  = h_eq.
```

`lambda` and `mu` are the retained rational multipliers. The checked identities
are

```text
lambda >= 0,
lambda^t A_ineq + mu^t E_eq = (0,...,0,1),
lambda^t b_ineq + mu^t h_eq <= 0.
```

Consequently these constraints imply `delta <= 0`, excluding a realization
with a strict gap above the prescribed minimal vectors. Each inequality is
accompanied by its lattice vector; `null` denotes the bound `delta <= 1`.

The certificate matrices contain the rows used by the final exact
certificate. They are not the whole linear-program constraint matrices.

## Initial linear programs and historical limits

The compressed initial systems reproduce the original builders' initial
constraint rules, row order, duplicate elimination, and representative vectors.
They maximize `delta` with unrestricted real variables and the displayed
inequalities and equalities. The vector for each inequality is represented by
an integer numerator row and a common denominator. Its entries are obtained
by dividing that numerator row by the denominator.

For a non-required vector `v = b/d`, the corresponding inequality is
`Q[v] >= 1 + delta`, cleared of denominators. The equalities require the
specified extra vector, when present, to have squared length 1.

The original search added further cuts when necessary. Its intermediate
candidate forms, numerical solver tableaux, and complete sequence of enlarged
constraint matrices were not logged and are not available. In particular,
109 retained certificate-row occurrences come from later cuts; those rows
are included in the main certificate catalogue. The archive does not claim
to reconstruct every historical intermediate solver state.

## Rank-five boundary matrix

The main catalogue also includes the previously displayed rank-five boundary
row `d9 = (0,0,0,0,-1,0,1)` and an explicit unimodular column change reducing it
to `(1,0,0,0,0,0,0)`. The column change is newly reconstructed. The rank-eight
and rank-nine arrays above are geometry and feasibility matrices, not cellular
boundary matrices. Their orientation calculations establish the absence of
degree-eight generators in the relevant oriented Voronoi complexes.

## Reproduce and check

Run these commands from the parent `voronoi-completeness` directory:

```sh
python3 export_all_matrices.py
python3 export_initial_lp.py
python3 verify_matrix_export.py
python3 verify_initial_lp_export.py
```

These commands use the Python standard library and the included source files.
Do not use Python's `-O` option: verification relies on assertions. The exact
enumeration and shortest-vector verification are provided separately by
`python3 verify_completeness.py`.

An optional additional check compares selected whole initial systems with
the original builders, including their row order and representative vectors:

```sh
python3 verify_initial_lp_export.py --compare-original
```

That option requires the original search dependencies NumPy, SciPy and SymPy.
The checked examples are named in `initial_lp_verification_summary.json`.
These are exact arithmetic checks in Python; no Lean proof is claimed.

For example, read a saved Gram matrix without losing exactness:

```python
import json
from fractions import Fraction

with open('matrices/all_matrices.json') as f:
    data = json.load(f)
record = next(r for r in data['records'] if r['id'] == 'extension_n8_0010')
Q = [[Fraction(x) for x in row] for row in record['gram']['Q']]
```

Read all initial systems, one at a time:

```python
import gzip
import json

with gzip.open('matrices/initial_lp_matrices.jsonl.gz', 'rt') as f:
    for line in f:
        problem = json.loads(line)
        A, b = problem['A_ineq'], problem['b_ineq']
        E, h = problem['E_eq'], problem['h_eq']
        # These are the full reconstructed initial matrices for this candidate.
```

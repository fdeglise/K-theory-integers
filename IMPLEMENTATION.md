# Implementation notes and schema v1

These are software inspection notes. The separate mathematical synthesis is in `reports/k8-verification-synthesis/main.tex`.

## Independent checking

`verify_bundle` validates and rebuilds the normalized dataset before checking
certificates. It recomputes primitive sign representatives, ambient rank,
the rank of every vector configuration and of its rank-one matrices, stable
record IDs, and both configuration and dataset digests. Exact canonical JSON
comparison also rejects boolean/integer substitutions in metadata.

For each witness, `check_witness` checks integer shape and unimodularity,
computes every image `g*v`, derives its signed permutation, and checks
bijectivity. It then recomputes the induced action as follows:

1. Flatten `v*v^T` by its upper triangle, without doubling off-diagonal entries.
2. Use exact rational RREF, scanning input columns from left to right, to choose
   a deterministic independent subset of rank-one matrices.
3. Select independent coordinate rows of that basis and invert the resulting
   square rational matrix.
4. Compute **explicitly** `g*(v*v^T)*g^T` for every chosen basis element.
5. Solve for the induced coordinate matrix and check equality in **all**
   upper-triangle coordinates, including those unused by the solve.
6. Compute its determinant by rational Gaussian elimination; require `-1`.
7. Compare every saved field to these recomputed values.

The invariant `C=sum(v*v^T)` is fixed exactly by a signed permutation. Full
ambient span makes C positive definite, and its strictly positive generator
coefficients put it in the relative interior of the cone. The fixed positive
ray contributes determinant `+1`; therefore orientation on W and on the
projective tangent `W / R*C` has the same sign. Ambient determinant alone and
ray-permutation parity are never substituted for this check. In particular,
the six-ray two-dimensional example in the tests has an even permutation of
rays and negative determinant on its three-dimensional span.

The source uses the right action `A . gamma = gamma^T*A*gamma` (DEKM §1.1).
Our matrix `g` gives the same action with `gamma=g^T`. Thus the conventions
agree after transposition; our vectors are always acted on as columns.

## Discovery

Every signed symmetry fixes C and preserves `Q=C^{-1}`. Discovery computes
the rational Gram matrix `v_i^T*Q*v_j`. Equal diagonal entries and sorted
absolute Gram rows restrict the candidates for each vector's image. A greedy
independent source frame prioritizes the rarest invariant classes. Backtracking
matches signed pairwise Gram entries and forbids repeated target pairs.
On a complete frame it solves for g, rejects fractional/nonunimodular maps,
and submits the result to the same checker used for saved certificates.

Only one reversing map is needed. Search order is deterministic. Wall alarms
bound the whole per-cell computation, including invariant setup; node checks
give an additional bound. `INCONCLUSIVE` covers timeout, node budget,
memory exhaustion, and search exhaustion without a witness. Discovery's
own success does not produce a K-theory status. `full` reloads the saved bundle
for the independent check. Verification imports neither discovery nor network
code; the tests additionally prohibit both while checking a supplied example.

## Files and hashes

`voronoi-dataset-v1` contains:

- source label, original source SHA-256, original vector storage, convention;
- cells with `id=source#NNNN` (one-based source-record number), sorted canonical
  vectors, `ambient_rank`, `vector_pairs`, `vector_rank`, `cone_dimension`,
  `projective_dimension`, primitivity, and a vector configuration digest;
- an empty anomaly list for valid data (anomalies are rejected with record
  diagnostics) and the dataset digest.

The digest hashes UTF-8/ASCII JSON with sorted keys, no insignificant whitespace,
and separators `,` and `:`. For dataset identity it covers the entire normalized
object **except** its `dataset_sha256` field. The pretty JSON file's own bytes
are not the canonical digest input. Each configuration digest hashes its
normalized vector array. All source record positions are preserved.

`voronoi-certificates-v1` contains the dataset digest and one certificate per
checked record. Each certificate includes the source ID and record number,
dataset digest, integer matrix g, signed permutation, deterministic basis
indices, induced rational action, and exact ambient/orientation determinants.
Fractions are reduced strings with positive denominators; integers have no
denominator. The checker compares canonical serialized fields, so alternate
spellings, booleans and extra metadata cannot masquerade as a recomputed result.

Full published coverage compares the supplied normalized object to a fresh
parse of the **pinned original target bytes**, whose SHA-256 is embedded in
`dataset.py`. Deleting a record and recomputing an internally consistent hash
does not pass this check. These pins, the program, the Python runtime, the
included source provenance, and the external classification identification
remain part of the trust boundary. Hashes do not establish completeness or
the mathematical meaning of a repository record.

`voronoi-result-v1` contains per-cell checks, dimensions, exact coverage,
missing IDs, result level, and the dependency manifest. Volatile measurements
and discovery statistics live in result files, never in certificate content.
Linux RSS is reported in KiB for the process; address-space limits and RSS
are different quantities. Timing/resource limits use floating-point clocks;
certification arithmetic exclusively uses Python integers and `Fraction`.

## Mathematical dependency handling

The manifest admits the list's completeness and exact realization and the
specified comparison theorems. It separately lists ranks 1–7, the integral
rank-9 H_0 theorem, the first-quadrant filtration (including ranks 0 and >=10),
degree-9 Hurewicz injectivity after localization, and the four small-prime
arithmetic exclusions. It records the source's mistaken `GL_N(R)` wording in
the finite-order lemma and uses its integral-matrix restriction instead.

In this application the stabilizer bound for N=8 requires primes at most 9,
which are exactly 2,3,5,7. The lower-rank inputs and the Hurewicz kernel require
no additional primes. Thus the exact ring `Z[1/7!]=Z[1/210]` suffices; the
program never represents localized coefficients by floating-point numbers.
Once V_8 is zero, no boundary matrices are needed for H_8. The rank-eight
homology statement is an inference from the certificates and admitted inputs,
not an admission. Integral K_8 vanishing remains conditional on every listed
mathematical premise, with the separate 2-primary argument explicitly included.

The PARI counterparts are not required computational inputs; they were
inspected as web text, without acquiring original bytes. Kupers's note was
initially inaccessible through the web, but a later local reference copy was
read and hash-checked against the existing catalogue. Its equation (1),
Lemmas 2.1 and 2.5, and Remark 2.4 agree with this route. Its alternate
rank-eight theorem is not used in place of the certificates.

The separately prepared synthesis appeared in the shared workspace during
implementation. It was read along with the new repository instructions and
ledger updates. Its conventions and conditional theorem agree with this
program. The corrected Smith-form guidance in the active prompt concerns
a future boundary-matrix extension and requires no change to this verifier.

## Validation

The focused unittest suite checks hand-known signs, dependent generators,
ordering/sign/basis changes, independent small-matrix determinant expansions,
malformed and executable syntax, invalid witnesses, every saved field,
incorrect hashes and dimensions, record deletion/duplication, missing and
duplicated certificates, forced budgets, all actual published witnesses,
offline operation, exit statuses, and Python optimized mode. Verification
does not depend on `assert` statements. Full perfect-form enumeration,
GL-orbit inequivalence, lower-rank recalculation, boundary matrices, and
K_12/K_16 are outside this implementation.

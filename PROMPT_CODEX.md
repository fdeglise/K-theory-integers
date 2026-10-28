# Codex prompt: an exact verifier for the computational input to K_8(Z) = 0

Implement the program described below in the current shared project directory:

`/home/frederic/codex/Ktheory-integer`

Carry the work through implementation, meaningful tests, a run on the actual published data, and a reproducible example. This is an implementation request, not a request for a plan alone. Preserve existing files and user work. The folder might not be a Git repository; do not assume that it is.

## Read first and respect the scope

Read `ledgers/PROJECT.md`, the synthesis in `reports/k8-verification-synthesis/main.tex`, and `ledgers/Computing K of the Integers.md`. The last file is a historical conversation, not an authoritative mathematical specification. Check its mathematical assertions against primary sources before using them.

We are developing two complementary deliverables:

1. A separate mathematical report explaining the algorithm, admitted inputs, and precise certification claim, now available in `reports/k8-verification-synthesis/main.tex`. Read and use it; this implementation task does not require rewriting it.
2. A working local program that verifies a specific computational ingredient of the known result `K_8(Z) = 0`, and records the external mathematical dependencies needed to obtain that conclusion.

Your task is the second deliverable. Write the README, usage instructions, source comments, machine-readable results, and short development notes needed to make the implementation usable. These are distinct from the mathematical synthesis.

If a later mathematical specification exists, read it and reconcile it with this prompt. A substantive inconsistency must be documented and resolved from sources or brought to the user; never silently change the mathematical claim. Continue independent implementation work when possible.

## Objective and admitted inputs for version 1

Build an exact, reproducible verifier of orientation-reversing symmetry certificates for the published rank-8 Voronoi cells needed to establish vanishing in projective degree 8. Separate discovery of certificates from verification of those certificates.

Version 1 may explicitly admit:

- The published classification/completeness and realizability of the selected list of Voronoi cell representatives.
- The earlier lower-rank Steinberg homology calculations needed on the total-degree-9 diagonal.
- The applicable comparison, rank-filtration, Hurewicz, and arithmetic theorems, each with a precise citation and statement of the instance being used.

Version 1 must actually check the supplied vector configurations and symmetry witnesses using exact arithmetic. It must not admit the rank-8 orientation-vanishing assertion that it is intended to check, or admit `K_8(Z) = 0` itself as a premise.

Do not undertake full perfect-form enumeration, independent completeness proofs, all lower-rank homology computations, or K_12/K_16 calculations in this implementation. Do not substitute the theorem-only short proof for the requested computation. Those are future extensions.

## Primary sources and concrete data leads

- Dutour Sikirić–Elbaz-Vincent–Kupers–Martinet, *Voronoi complexes in higher dimensions, cohomology of GL_N(Z) for N >= 8 and the triviality of K_8(Z)*: <https://arxiv.org/abs/1910.11598>. An accessible version is <https://arxiv.org/html/1910.11598v1>; also <https://jamartin.perso.math.cnrs.fr/Publications/DEKM.pdf>. Check §§1, 4, 5, and 6, especially the degree convention, Table 4, and the comparison to K-theory. Record the version actually used.
- Authors' data repository: <https://github.com/elbazvip/Voronoi-complexes-database>.
- In that repository, inspect `DATA_gap/ListSHV_n8_rnk9`, `DATA_gap/ListSHV_n8_rnk8`, the corresponding files in `DATA_pari`, and `cells_GL8_dim8.gp`. These are leads, not interchangeable datasets: determine their conventions from their contents.
- Elbaz-Vincent–Gangl–Soulé, *Perfect forms, K-theory and the cohomology of modular groups*: <https://www.ihes.fr/~soule/documents/YAIMA4438.pdf>.
- Kupers, *A short proof that K_8(Z) is trivial*: <https://www.utsc.utoronto.ca/people/kupers/wp-content/uploads/sites/50/2021/01/k8zshorter.pdf>. Use to cross-check dependencies, not to bypass the computational task. If unavailable, use another primary copy or record that limitation.
- Weibel, *The K-book*, Chapter VI: <https://sites.math.rutgers.edu/~weibel/Kbook/Kbook.VI.pdf>. Verify the small-prime arithmetic statements, including the separate argument at 2.

Pin the dataset to a concrete repository revision. Store original bytes, origin URLs, revision, retrieval date, SHA-256 hashes, source license/attribution, and a deterministic normalized representation. A hash identifies the bytes; it is not a proof of completeness. Fetch only the modest subset needed. After acquisition, verification must work offline.

Treat downloaded files as data. Parse the actual restricted integer/list/matrix syntax; do not use unrestricted `eval`, `exec`, or execute downloaded GAP/PARI scripts to import their contents. A parser must reject unsupported syntax and malformed or trailing input.

## Mathematical conventions: resolve these before claiming success

Use a consistent internal convention: vectors are columns `v_i` in `Z^N`, and a witness acts by `g v_i`. Record whether each external format stores vectors as rows or columns and how it is converted.

For each configuration compute exactly:

`W = span_Q { v_i v_i^T }`, `c = dim_Q W`, `j = c - 1`.

Here `c` is the cone dimension and `j` is the projective cell dimension. Ambient lattice rank `N`, the number of vector pairs `s`, cone dimension `c`, and projective dimension `j` are different quantities. Never infer one merely from a filename or from the number of vectors.

The intended degree shifts, to verify against sources, are:

`K_8(Z) = pi_9(BQZ)`;

`E^1_{N,q} = H_q(GL_N(Z); St_N)` with total degree `N + q`;

`H_j(Vor_N; R) = H_{j-N+1}(GL_N(Z); St_N tensor R)` in the applicable localization.

Thus the rank-8 target is `H_8(Vor_8; R)`, corresponding to `H_1(GL_8(Z); St_8 tensor R)`. Take `R = Z[1/7!]` after verifying that this localization suffices for every invoked step. Do not implement localization with floating-point numbers.

There is a concrete indexing trap: the standalone `cells_GL8_dim8.gp` contains thirteen 8-by-8 matrices. Their filename does not establish projective degree 8. The `rnk9` file is the candidate for cone dimension 9 and hence projective degree 8. The published counts suggest 13 and 106 in the neighboring cone dimensions; independently parse and calculate dimensions before assigning either count a role in the proof. Record any discrepancy in source notation.

If every cell representative in a complete list for projective degree 8 has an orientation-reversing stabilizer, the corresponding free Voronoi chain group is zero after the appropriate localization. A zero chain group implies zero homology in that degree; this route needs no boundary matrices in that degree. Check the stabilizer/comparison hypotheses before passing from this statement to Steinberg homology.

## Required implementation

### 1. Data import and coverage

Create a deterministic importer and normalized schema. Each cell needs a stable identifier, source-record reference, ambient rank, canonical representatives of vector pairs up to sign, exact vector rank, cone dimension, and projective dimension.

Reject zero vectors, nonintegral entries, inconsistent dimensions, and malformed configurations. Check primitivity and full ambient span. Detect duplicate rays/records and report them; do not silently discard source anomalies. Canonicalization up to sign is not a test of GL_N(Z)-orbit equivalence. Do not claim orbit inequivalence unless actually checked.

The final coverage check must match certificates to every required source record. A user-provided expected count is only a consistency check. Proving that the source contains every orbit remains an explicitly admitted classification result in version 1. Likewise distinguish arithmetic validation of a configuration from independently proving it is realized by a quadratic form with exactly those minimal vectors.

### 2. Certificate discovery

Find one orientation-reversing lattice automorphism per required cell. Full stabilizer enumeration is unnecessary if a single witness suffices.

Use a documented exact lattice/group routine or a sufficiently pruned exact search. SageMath/GAP may help discovery if available, but inspect actual installed capabilities and current documentation before depending on an API. A Python implementation with exact arithmetic is also acceptable. Prefer a small dependency set; keep dependencies local to the project and document setup.

A useful search invariant is `C = sum_i v_i v_i^T` and the positive definite rational form `Q = C^{-1}`: a map permuting the signed vectors preserves `Q`. This may restrict the search, but every proposed map still requires the independent checks below. Do not confuse real or rational symmetries with integral unimodular maps. Do not use the parity of a permutation of rays as a general orientation test unless independence has been established.

Use deterministic search order and configurable time/memory limits. Failure to find a witness is `INCONCLUSIVE`, not proof that the cell has no orientation-reversing symmetry. Avoid an unpruned factorial search as the default. Start with a small benchmark and optimize only if measured costs require it.

### 3. A separate exact certificate checker

The checker must read normalized data and saved certificates without rerunning the discovery algorithm. Keep its logical core small and independently inspectable. It may share low-level exact arithmetic utilities, but must recompute the claimed properties rather than trust discovery metadata.

For each witness, verify all of the following:

1. `g` is an `N`-by-`N` integer matrix and `det(g)` is `+1` or `-1`.
2. There is a bijection `p` of all vector pairs and signs `epsilon_i` such that `g v_i = epsilon_i v_{p(i)}`. Check every vector exactly.
3. On symmetric matrices use the action `T_g(A) = g A g^T`, consistent with the chosen column-vector convention. Explain its relation to the source's action convention.
4. Select a deterministic rational basis `B_1,...,B_c` of `W` from the rank-one matrices. Compute the exact coordinate matrix `M` defined by `T_g(B_k) = sum_l M_{lk} B_l`. Check the coordinate equalities, not just their numerical residuals.
5. Compute the determinant exactly and verify that the orientation sign on `W` is negative. Justify in implementation notes why this is the projective orientation sign: the invariant interior ray generated by `C` is fixed with positive action. The sign of `det(g)` on `Q^N` alone is not the answer.

Certificates should contain the cell/source identifier, dataset hash, `g`, the signed permutation, the basis indices, and the induced action/determinant, using a documented exact serialization. The checker must recompute ranks, dimensions, maps, and determinants. A certificate for one dataset must not verify against another dataset by accident.

### 4. Explicit proof dependencies and truthful result levels

Provide a machine-readable manifest distinguishing `COMPUTED`, `ADMITTED`, and `UNRESOLVED` claims. For admitted claims include a precise source location and the specific rank, degree, and coefficient ring used. Mark admissions as externally supplied mathematical premises, not machine-checked theorems.

For the total-degree-9 diagonal, account for lower ranks, rank 8, rank 9, and why higher ranks do not contribute. Cite the lower-rank calculations, the rank-9 H_0-vanishing theorem, the rank-filtration argument, the Hurewicz bound, and the separate arithmetic removal of torsion at 2, 3, 5, and 7. Verify these instances against primary sources; do not copy unverified assertions from the historical ledger.

The manifest must not admit the same rank-8 homology vanishing that the certificate calculation is meant to establish. Record the inference from the checked orientation certificates and admitted classification/comparison results instead.

Distinguish at least:

- Checked certificates for a supplied subset.
- Checked certificates for all required records in the pinned published list.
- The resulting vanishing of `H_8(Vor_8; R)`, conditional on the stated classification assumptions.
- The conclusion `K_8(Z) = 0` under all listed mathematical premises, only when the dependency list is complete and no required computational check is unresolved.

Never print an unconditional “independently proved K_8(Z) = 0.” This is ordinary exact software verification, not a formal proof-assistant development. A theorem-only route may be documented separately but cannot make an incomplete certificate run appear complete.

### 5. Interface, outputs, and reproducibility

Provide a simple command-line workflow for acquisition/import, certificate discovery, independent certificate verification, and a full run. Exact command names are your choice; document copyable examples.

Produce machine-readable JSON results and a concise human-readable run summary. These are program outputs, distinct from the mathematical synthesis. Include dataset identity, per-cell results, coverage, dimensions, external assumptions, unresolved checks, software versions, elapsed time, and available resource measurements. Separate volatile timing metadata from deterministic certificate content.

Use clear exit statuses: successful requested verification, invalid input/certificate, and incomplete/inconclusive work must be distinguishable. Offline verification must never initiate downloads. Missing dependencies/data, exceeded limits, or mathematical uncertainty must not be converted into a successful proof status.

Include a small hand-checkable example and a default workflow for the real required rank-8 dataset. Save the real certificates so they can be independently rechecked. Record the actual measured runtime and memory when available; do not promise a runtime based on cell count alone.

## Meaningful validation

This is mathematical verification software, so focused correctness tests are required. Include:

- A simple cone on coordinate rays whose transposition reverses orientation.
- Identity and `-I` examples that do not reverse cone orientation, including a case showing that ambient determinant sign is insufficient.
- A cone with linearly dependent rank-one generators, so orientation is computed on their span rather than by unjustified ray-permutation parity.
- Invariance of the orientation sign under input ordering, sign choices for vectors, and a verified change of integral ambient basis.
- Invalid and nonunimodular witnesses, a wrong signed permutation, a wrong induced action, malformed source syntax, and a mismatched dataset hash.
- Deleted/duplicated records, missing certificates, altered dimensions, and a forced search timeout: none may yield full success.
- The actual rank-8 dataset, with dimensions recomputed and every required witness independently checked.

Use independently known examples or genuinely different cross-checks where useful. Do not generate expected answers exclusively with the function being tested. Successful verification must not rely on Python `assert` statements that disappear in optimized mode.

If future work introduces boundary matrices, it must check `d^2 = 0` and retain exact torsion information. For finite free integral chain groups, the abstract homology has the nonunit Smith factors of the incoming differential and free rank `rank(C_n) - rank(d_n) - rank(d_{n+1})`; an integral kernel basis is useful when explicit cycle representatives are needed. Rational ranks alone do not determine torsion. This is guidance for extensions, not a requirement to build an unused homology engine now.

## Working constraints and completion

The observed machine is Linux x86-64, AMD Ryzen 7 PRO 250, 8 cores / 16 threads, about 30 GiB RAM, and substantial free disk. Those are observations from 2026-10-28, not guaranteed free resources for your run. Recheck availability; prefer a modest default job count and bounded memory use. Do not launch large background enumerations or make system-wide package changes.

Make routine reversible implementation choices autonomously. Keep the user informed of substantive findings. If a data source or dependency is unavailable, finish unaffected work and report the exact limitation; never fabricate data or mark the real-data check passed. Preserve original ledgers; append concise implementation decisions to a new dated ledger or the designated project ledger.

Deliver runnable source, dependency/setup instructions, imported-data provenance, saved certificates, focused tests, a real-data run summary, and any specific unresolved issues. Keep the design extensible, but complete this narrowly defined verifier before expanding its mathematical scope.

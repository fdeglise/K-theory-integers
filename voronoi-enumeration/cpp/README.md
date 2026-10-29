# C++ Voronoi cell enumeration

This implements a **maximal minimal-frame search with exact face closure**.
It supports lattice ranks 1--7 and enumerates full-vector-rank cells with cone
dimensions `n .. n+k`. Projective cell dimension is cone dimension minus one.
It does not enumerate perfect forms or the vertices of whole Ryshkov slices.

## Build and run

```sh
bash scripts/install-deps.sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j4
build/voronoi_cpp 4 --excess 6 --stabilizers --json cpp_n4.json
build/voronoi_cpp 7 --excess 1 --verbose --json cpp_n7.json
ctest --test-dir build --output-on-failure
```

The dependency installer downloads pinned **GMP 6.3.0**, **Boost 1.89.0
headers**, and **SoPlex 8.0.3** from their official sources, checks the archived
SHA-256 values in `dependencies.sha256`, and builds into `.deps/prefix`.
It runs GMP's `make check`. No `sudo` or system configuration changes are needed.
Downloads require network access. The `.deps` directory is ignored by Git.
Dependencies retain their upstream licenses (GMP LGPL/GPL, Boost BSL-1.0,
SoPlex Apache-2.0); source archives include the license texts.

The header-only `cmake/FindBoost.cmake` bridges current CMake and the official
Boost source archive, without needing to build unrelated Boost libraries.
The installer explicitly checks that SoPlex's rational backend is enabled.
On macOS, the build locates libc++ headers in the SDK if the normal toolchain
header symlink is absent; it does not modify the SDK.

The default depth is one (zero in rank 1). Use `--excess 0` for just the smallest cones, or increase
it to explore extensions. Exact stabilizers and orientation characters are
optional because enumerating all automorphisms can be much more expensive than
testing equivalence. Without `--stabilizers`, the JSON records those fields as
`null`, not as trivial groups or orientation-preserving cells.

## Search and completeness

1. **Frame catalogue (`frames.hpp`).** A minimal-vector frame has index at most
   `floor(gamma_n^(n/2))`. Certified bounds for ranks 1--7 are
   `1,1,1,2,2,4,8`. Instead of listing all raw HNF matrices, enumerate the finite
   abelian quotient groups of each index and the embeddings of those groups
   into `(Q/Z)^n`. An embedding is represented by its coordinate characters.
   Generate multisets of characters modulo coordinate signs and quotient-group
   automorphisms. Injectivity is checked, and Watson's necessary inequalities
   reject impossible frames. An exact integer HNF then provides the ambient
   lattice basis. This includes noncyclic quotients, not just cyclic codes.
2. **Forced vectors.** Equality in a Watson inequality forces additional
   minimal vectors. Add these before solving an LP. Their rank can already
   exclude a frame from the requested dimension range. In particular, rank 4's
   index-two frame forces all 12 D4 minimal-vector pairs at once. For binary
   weight-four cosets, averaging lengths over coordinate signs forces every
   sign choice; the index-eight E7 model similarly supplies all 63 pairs.
3. **Closure (`closure.hpp`).** Normalize the chosen frame lengths to one.
   SoPlex solves rational LPs in the `n(n-1)/2` off-diagonal coordinates. An
   exact LDL decomposition supplies a nonpositive direction if the trial form
   is not positive definite; checked small integer approximations keep cuts
   from acquiring unnecessarily large coefficients. For a positive form,
   exhaustive rational LDL ellipsoid enumeration supplies all lattice vectors
   of length at most one. Shorter vectors become cutting planes.
4. **Relative interior.** A margin objective seeks strict slack away from the
   imposed equalities. Any accidental minimal vector is separately maximized
   over the face. Maximum one proves it forced; otherwise average the current
   witness with a feasible optimizer that makes it longer. The minimal set of
   an average is exactly the intersection of the endpoint minimal sets. This
   returns the complete forced minimal set, not an arbitrary LP vertex.
5. **Extensions (`main.cpp`).** Choose a frame of maximal index `I` within the
   target cell. Every minimal vector then has frame coordinates in `[-1,1]^n`,
   by Cramer's rule. Enumerate primitive vectors in this finite lattice cube.
   Cofactor inequalities reject extensions introducing an n-minor larger than
   `I`; exact closure and rank pruning handle forced jumps and nonsimplicial
   cones. Every target has such a maximal frame, and adding generators of its
   rank-one span reaches it after at most `r-n` independent extensions.
   Before closure, orbit BFS reduces candidates under coordinate flips and
   swaps that are verified to preserve both the lattice code and current cell.
6. **Integral equivalence (`configurations.hpp`).** For a vector-pair family
   `S`, use the intrinsic form `C_S=sum(v*v^t)`. Intrinsic Gram invariants bucket
   candidates. Backtracking extends a fixed independent frame using exact
   pairings in `C_S^-1`; every accepted map is checked to be integral,
   unimodular, and to carry the complete family onto its target. The fixed frame
   may be nonunimodular. Both trial families and closed cells are deduplicated
   under actual integral equivalence. Trial-family registries are separate for
   each maximal-index restriction, so a rejection at one index cannot suppress
   a valid search at a larger index.

Global cell deduplication is safe here because only maximal-index frames are
used: every coface with the same maximal index fits in the cube of **any**
maximal frame in that cell. Cofaces with a larger maximal index are searched
separately under that index's frame models.

All accepted witnesses are checked using GMP rational arithmetic. SoPlex is
configured for exact rational solving and checking; infeasibility and optimality
rely on its exact solver. The implementation independently verifies returned
primal bounds, equalities and lattice inequalities. A solver failure or failed
exact check aborts the run rather than silently discarding a cell.

## Verification and limits

```sh
build/voronoi_cpp --self-test
python3 cpp/check_counts.py build/voronoi_cpp
python3 cpp/check_counts.py build/voronoi_cpp --slow
python3 cpp/check_counts.py build/voronoi_cpp --deep
python3 cpp/verify_witnesses.py cpp/results/rank7_depth3.json
```

The integration suite compares published small-depth counts in ranks 2--7.
For ranks 2--4 it also verifies complete minimal-vector sets with the independent
Python implementation and matches each representative by integral equivalence,
including stabilizer orders and orientation characters. The slow suite extends
the depth for ranks 5--7. The core tests also close a noncyclic index-eight E7
frame and check all 63 minimal-vector pairs in cone dimension 28.

`verify_witnesses.py` independently checks exported form witnesses, their full
minimal-vector sets, cone dimensions and maximal frame indices. By default it
selects the largest coordinate box at each frame index; `--all` checks every
cell. Its brute-force Python box enumeration is deliberately independent of
the C++ LDL recursion. Witness checking does not establish the completeness of
the orbit catalogue.

The saved rank-7 depth-three run contains **616 orbits**, with counts
`6,28,115,467` in cone dimensions 7--10. It took about five minutes on the
development machine, using 12,230 exact LP solves. Runtime depends strongly on
depth and available symmetry; rank-11 feasibility is not established.

Reference counts are from
[DEKM, Table 3](https://jamartin.perso.math.cnrs.fr/Publications/DEKM.pdf).
The full rank-4 answer is `1,3,4,4,2,2,2` in cone dimensions 4--10.
For rank 7, dimensions 7--10 have `6,28,115,467` orbits.

Work limits include `--max-lp`, `--max-families`, `--max-svp-nodes` (per SVP
call), and `--seconds`. Exceeding a bound exits with status 2 and does not write
a new successful JSON result or overwrite an existing output file. A completed
run sets `complete: true` only for the explicitly requested dimension range.
An increased `--index-bound` up to 8 is allowed; a smaller-than-certified bound
is rejected. No rank >=8 completeness claim is currently supported. Use
`--trace` for detailed family, LP and SVP diagnostics.

## Preparing for rank 11

Arithmetic, frame models, LP closure, SVP, equivalence, and search control are
separate components. The low-depth search never requires the full perfect-form
classification. However, this is still intentionally a simple implementation:

- The quotient-code catalogue and automorphism enumeration need stronger
  classification methods when the index bound becomes large.
- The LP backend is isolated in `lp_soplex.cpp`; it rebuilds the LP for each
  solve. There are no warm starts or saved bases.
- SVP has exact recursion but no LLL preprocessing, parallel enumeration, or
  advanced pruning.
- Equivalence uses backtracking, not a production lattice-isometry or graph
  canonicalization library. Extension-orbit pruning uses only elementary
  signed-coordinate generators, rather than the full stabilizer.
- Trial-family and cell registries reside in memory. Checkpoints, sharding,
  persistent deduplication, and resumable jobs are not implemented.
- Vector coordinates use checked machine integers; form arithmetic and linear
  algebra use arbitrary-precision rationals. A coordinate overflow stops the
  search. Larger ranks may need arbitrary-precision vector storage as well.
- Exact cutting-plane loops have explicit work limits; no practical runtime
  bound is asserted. Difficult instances fail visibly when limits are reached.
- There are no boundary matrices or homology calculations yet.

These are replacement points for eventual rank-11 work, not claims that the
present implementation can already handle rank 11.

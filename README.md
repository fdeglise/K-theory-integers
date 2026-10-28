# Exact rank-eight orientation verifier

This local program checks saved integral symmetry certificates for all **106**
records in the published cone-dimension-9 rank-eight list. The checked maps
reverse orientation on the span of the rank-one matrices. With the explicit
external classification and comparison premises, this gives
`H_8(Vor_GL_8; Z[1/7!]) = 0`. The dependency manifest also records the further
external premises giving `K_8(Z) = 0`.

Completeness and realizability of the published list are admitted. This program
does not independently enumerate Voronoi cells, prove their realizability, or
formally verify the external theorems. It never uses the rank-eight vanishing
or `K_8(Z)=0` as an admitted premise. This README and the implementation notes
document the software. The separately prepared mathematical specification is
[`reports/k8-verification-synthesis/main.tex`](reports/k8-verification-synthesis/main.tex);
it was read and reconciled when it appeared during implementation. This task
does not modify that report.

## Reference corpus and report build

The [corpus guide](references/k-theory/README.md) describes ten retained reference PDFs. The [catalogue](references/k-theory/catalog.json) records their source URLs, checksums, and acquisition qualifications. See also the [repository ledger](ledgers/REPOSITORY-2026-10-28.md).

Compile the synthesis with `latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=../../tmp/pdfs/k8-synthesis main.tex` from its directory. Generated report PDFs and auxiliary files are ignored; reference PDFs in `references/k-theory/pdfs/` are deliberately tracked. Original data bytes and their upstream license are preserved in `data/original/`; scholarly PDFs retain their original rights.

## Run offline

Use Linux and Python **3.10 or later**. The program and tests use only the
standard library. SageMath, GAP, PARI, SymPy, NumPy, pytest, and package downloads
are unnecessary. Run these commands from this directory:

```bash
python3 -m voronoi_verify verify
python3 -m unittest discover -s tests -v
```

The first command reads the included normalized data, original pinned target
bytes, and saved certificates. It neither imports discovery nor downloads
anything. It writes `results/verify.json` and `results/verify.txt`.

Reproduce acquisition/import/discovery/checking as separate steps:

```bash
# Optional, explicit network operation; valid included originals are reused.
python3 -m voronoi_verify acquire

# All commands below work offline with the included originals.
python3 -m voronoi_verify import
python3 -m voronoi_verify audit
python3 -m voronoi_verify discover
python3 -m voronoi_verify verify

# Or run import (when missing), discovery, and independent saved-file checking:
python3 -m voronoi_verify full
```

`acquire` requests only five small pinned source files and verifies their
embedded SHA-256 values before saving. It refuses to overwrite conflicting
existing originals. `import` deterministically rebuilds the default normalized
file. `full` validates an existing normalized file; it does not repair altered
input silently. Output paths are relative to the working directory. When a
JSON output changes, its earlier bytes are retained under the adjacent
`history/` directory with a SHA-256 filename before replacement. Text summaries
are regenerated from the JSON evidence. Original ledgers are preserved.

The actual full run recorded in `results/full.json` checked **106/106** cells,
with **3.004 seconds** elapsed, including **2.168 seconds** of discovery and
**0.777 seconds** of verification. Peak process RSS was **21,844 KiB**. This was
one process on the observed Linux/Python 3.12.3 machine; it is a measured run,
not a runtime promise. Measurements and search statistics are kept outside the
deterministic certificate bundle.

## Small example

The two vectors `e_1,e_2` give the diagonal cone, with basis `E_11,E_22`.
The matrix `[[0,1],[1,0]]` exchanges these generators. Its induced matrix is
`[[0,1],[1,0]]`, with determinant `-1`. The fixed interior ray is `E_11+E_22`;
the induced map reverses the orientation of the projective interval.

```bash
python3 -m voronoi_verify verify --supplied-data \
  --dataset examples/coordinate-rays.json \
  --certificates examples/coordinate-rays-certificate.json \
  --output examples/verification.json
```

`--supplied-data` certifies the supplied configurations only. It produces no
complete published-list or K-theory conclusion. A fresh discovery run for this
example is:

```bash
python3 -m voronoi_verify full --supplied-data \
  --dataset examples/coordinate-rays.json \
  --certificates /tmp/coordinate-certificates.json \
  --output /tmp/coordinate-full.json
```

## Limits, partial runs, and exit codes

Discovery uses exact Gram invariants and pruned signed-frame backtracking.
The defaults are one process, 30 wall seconds and 200,000 nodes per cell, and
512 MiB of process address space. No large background enumeration is launched.
Tune with `--seconds-per-cell`, `--max-nodes`, and `--memory-mib`.
Exhaustion or reaching a limit is **INCONCLUSIVE**, never evidence that no
reversing symmetry exists. The supported input bound is ambient rank 32 and
512 vector pairs per record; this is not a claim of practical discovery speed
through that entire range.

```bash
# Benchmark a small prefix without replacing the included certificates:
python3 -m voronoi_verify discover --limit 5 \
  --certificates /tmp/first-five.json --output /tmp/first-five-search.json

# Explicitly request verification of just those certificates:
python3 -m voronoi_verify verify --subset \
  --certificates /tmp/first-five.json --output /tmp/first-five-check.json

# Deliberately force an incomplete run (exit 3):
python3 -m voronoi_verify full --seconds-per-cell 0 \
  --certificates /tmp/no-certificates.json --output /tmp/incomplete.json
```

| Exit | Meaning |
|---:|---|
| 0 | The requested operation/check passed. Read the result level: an explicitly requested subset can pass. |
| 2 | Invalid syntax, data, metadata, or certificate. |
| 3 | Missing inputs, inaccessible source, resource limit, or incomplete/inconclusive verification. |

The default verification requires the complete pinned published target.
`--subset` permits missing certificates but never enables the conditional
homology/K-theory conclusion. A self-consistent shortened dataset fails the
default published coverage requirement even if all its certificates pass.

JSON distinguishes `CHECKED_SUPPLIED_SUBSET`, `CHECKED_ALL_SUPPLIED_RECORDS`,
and `CHECKED_ALL_PINNED_PUBLISHED_RECORDS`. Mathematical premises are marked
`ADMITTED`, arithmetic checks `COMPUTED`, unfulfilled requirements `UNRESOLVED`,
and successful mathematical inferences `CONDITIONAL`.

## Data identity and dimensions

The authors' [database](https://github.com/elbazvip/Voronoi-complexes-database)
is pinned at revision `d4c7b40ad2c95850d16160b347953dbaa9241486`.
Original bytes, upstream GPL-3.0 license, attribution, URLs, retrieval date,
and SHA-256 values are preserved in `data/original/` and
[`data/provenance.json`](data/provenance.json). A matching hash establishes
byte identity; classification completeness is an external premise.

The program recomputes `c=rank_Q{v*v^T}` and `j=c-1`:

| Source | Records | Ambient rank N | Vector pairs s | Cone c | Projective j |
|---|---:|---:|---:|---:|---:|
| `DATA_gap/ListSHV_n8_rnk8` | 13 | 8 | 8 | 8 | 7 |
| `DATA_gap/ListSHV_n8_rnk9` | 106 | 8 | 9 | 9 | 8 |
| `cells_GL8_dim8.gp` | 13 | 8 | 8 | 8 | 7 |

This resolves the trap in the standalone filename. DEKM v1 §1.2 defines
projective cell indices, but Table 4 places 13 and 106 under 8 and 9.
We use the recomputed dimensions and explicitly admit the identification of
the pinned list with the published cone classification. The table's existing
orientation-zero row is **not** used as a premise. See
[`results/dimension-audit.json`](results/dimension-audit.json).

GAP stores vectors as rows of each `s` by `N` configuration. Each row becomes
an internal column vector. The similarly named PARI files are actually
`ListSHV_n8_rank8` and `ListSHV_n8_rank9`; their inspected text assigns nested
row-vector lists to `U1`, `U2`, etc. They are neither the standalone semicolon
matrices nor part of this verifier's required input. Their original bytes
could not be acquired after shell network access became unavailable; no
byte-level PARI/GAP equivalence is claimed. For the standalone square matrices,
the audit chooses columns and records that the dimension conclusion holds
under either transpose convention. It does not claim GL-orbit equivalence
between the standalone matrices and GAP configurations.

## Implementation and exact serialization

- [`voronoi_verify/checker.py`](voronoi_verify/checker.py) is the small checker;
  it shares rational arithmetic with discovery but imports no search code.
- [`voronoi_verify/discover.py`](voronoi_verify/discover.py) finds witnesses.
- [`voronoi_verify/formats.py`](voronoi_verify/formats.py) parses restricted
  integer/list/matrix syntax. Downloaded code is never evaluated or executed.
- [`voronoi_verify/claims.py`](voronoi_verify/claims.py) records precise external
  premises and the conditional inference graph. The same manifest is embedded
  in each verification result.
- [`IMPLEMENTATION.md`](IMPLEMENTATION.md) explains the checker's mechanics,
  trust boundary, and schema.

Normalized data uses lexicographically sorted primitive integer vectors, with
first nonzero coordinate positive. Duplicate rays and identical records up to
this sign/order convention are errors, not silently deleted. This is not a
GL-orbit equivalence test. All vector and rank-one ranks are recomputed.
Certificates use integer arrays for `g`, zero-based permutation/basis indices,
signs `+1/-1`, and canonical rational strings (`"-1"`, `"2/3"`).

## Primary references and source limitations

The implementation uses [DEKM, arXiv:1910.11598v1, 25 October
2019](https://arxiv.org/html/1910.11598v1), particularly §§1, 3–6, and
[Elbaz-Vincent–Gangl–Soulé, *Perfect forms, K-theory and the cohomology of
modular groups*](https://www.ihes.fr/~soule/documents/YAIMA4438.pdf), Theorem
4.3 and §6.3. The separate arithmetic exclusions use [Weibel, *The K-book*,
Chapter VI](https://sites.math.rutgers.edu/~weibel/Kbook/Kbook.VI.pdf),
Corollary 9.8 (prime 2), Proposition 10.5 (odd regular primes), and §2.4.1
(regularity of 3, 5, 7). Exact rank/degree/ring instances are in the manifest.

[Kupers's short proof](https://www.utsc.utoronto.ca/people/kupers/wp-content/uploads/sites/50/2021/01/k8zshorter.pdf)
initially returned HTTP 403. A local copy subsequently appeared in
[`references/k-theory/pdfs/09_Kupers_Short-proof-K8.pdf`](references/k-theory/pdfs/09_Kupers_Short-proof-K8.pdf).
The copy is dated 21 October 2019; its hash is recorded in the dependency
manifest and the existing reference catalogue. Its equation (1), Lemmas 2.1
and 2.5, and Remark 2.4 cross-check the filtration, rank-nine input, and
computational route. The program retains the requested rank-eight calculation.

# Voronoi enumeration: publication note

Date: 29 October 2026.

The user requested adding this project as a new folder in
`fdeglise/K-theory-integers`, together with a transcript of the development
conversation. The folder is `voronoi-enumeration/`. The repository was inspected
at main revision `7098bce90a7dcb5d22e5999bb9af895754aae801`.

This is a separate cell-enumeration experiment. It does not replace the existing
rank-eight orientation verifier, its pinned data, mathematical report, or
historical ledgers. No existing rank-eight certificate or K-theory inference is
used as an enumeration premise.

## Implemented scope

- Dependency-free exact Python enumeration in ambient ranks 1--4.
- C++ enumeration in ambient ranks 1--7, designed for small excess of cone
  dimension over ambient rank. It uses GMP rational arithmetic, exact SoPlex LP
  closure, quotient-code frame catalogues, exact shortest-vector checks, and
  integral equivalence testing.
- Cone dimension is the exact rank of the matrices `v*v^T`; projective cell
  dimension is one less. Stabilizer orientation is tested on that span.
- The programs reject incomplete searches instead of reporting complete success.
  The C++ JSON's `complete` flag concerns the requested dimension range and the
  implemented exhaustive search, not a formal proof-assistant certificate.
- There are no boundary matrices, homology computations, or new K12 vanishing
  results. The present implementation does not support rank 11.

## Observed validation

The 22 Python tests passed. The C++ core, default published-count integration,
deep integration, and address/undefined-behavior sanitizer core tests passed.
The deep integration reproduced:

| Ambient rank | Cone dimensions | Orbit counts |
|---|---|---|
| 2 | 2--3 | 1, 1 |
| 3 | 3--6 | 1, 2, 1, 1 |
| 4 | 4--10 | 1, 3, 4, 4, 2, 2, 2 |
| 5 | 5--8 | 2, 5, 10, 16 |
| 6 | 6--8 | 3, 10, 28 |
| 7 | 7--10 | 6, 28, 115, 467 |

These are regression comparisons with DEKM Table 3, not a standalone proof of
classification completeness. The exhaustiveness argument and its mathematical
inputs are described in `README.md` and `cpp/README.md`; execution still relies
on the implementation and its exact-arithmetic/LP libraries.

The saved rank-seven run contains 616 orbits, took about five minutes on the
development machine, and used 12,230 exact LP solves. The Python reference
independently checked the largest verification box at each of the four maximal
frame indices occurring in that saved run: four of its 616 form witnesses.
All exported rank-2--4 form witnesses were checked in integration, and their
orbits, stabilizer orders, and orientation characters were compared with the
independent Python implementation. The core tests include the nonsimplicial
D4 cone and an index-eight elementary-abelian E7 frame with 63 minimal pairs.

## Reproducibility and historical material

`dependencies.sha256` pins the source archives downloaded from the official GMP,
Boost, and SoPlex locations in `scripts/install-deps.sh`. Dependencies install
locally under `.deps/`; downloaded archives and sources, build products, caches,
and operating-system metadata are excluded from publication. The saved result
JSON and all project source, tests, and documentation are included.

`TRANSCRIPT.md` is a snapshot of the user-visible dialogue. It retains corrected
and abandoned proposals as history and is not the checked implementation
specification. It excludes internal instructions, internal reasoning, and tool
logs. The final Git push confirmation necessarily follows this snapshot.

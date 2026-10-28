# K-theory of the integers: project ledger

Date: 2026-10-28 (America/Los_Angeles).

## Purpose

Develop reproducible local computations relevant to algebraic K-theory of the integers. The first target is a verifier of computational ingredients in the established result `K_8(Z) = 0`. Later work may explore K_12 and K_16, but no algorithm or runtime guarantee for those degrees has been established in this project.

## User decisions

- Begin with K_8 as a known-result verification exercise.
- Separate the mathematical account from the implementation.
- The mathematical account explains the exact algorithm, admitted data and theorems, and what the program does and does not certify. The user subsequently requested this report; it is now in `../reports/k8-verification-synthesis/main.tex`.
- The earlier preparation step produced an implementation prompt as Markdown. The current report-writing step does not implement the program.
- Design ledgers may record scope and choices.

## Current files and authority

- `../PROMPT_CODEX.md`: the implementation prompt prepared for a future Codex run.
- `../reports/k8-verification-synthesis/main.tex`: the human-readable mathematical specification and conditional certification theorem.
- `SYNTHESIS-2026-10-28.md`: report preparation, source checks, and the correction to the earlier Smith-form claim.
- `Computing K of the Integers.md`: historical discussion and reference leads. Preserve it. Its assertions require primary-source verification and must not be treated as a checked specification.
- This ledger: project decisions and implementation boundaries. It is not the mathematical report or a record of completed mathematical verification.

This report-writing work has not implemented a program or produced symmetry certificates. A preliminary inspection of the published target file confirmed its record count and exact dimensions; this is not a completed certificate run.

## Proposed first implementation

Use published cell representatives and discover an orientation-reversing integral symmetry for each required rank-8 cell. Save explicit witnesses and recheck them with a separate exact checker. The first version admits the published list's completeness and realizability, earlier lower-rank computations, and the comparison/arithmetic theorems, all explicitly identified.

The target chain degree is projective degree 8. The main result of the computation should be stated in terms of checked witnesses and coverage of the pinned input. Any conclusion about all cells, Steinberg homology, or K_8 must carry its admitted mathematical premises.

The program must distinguish incomplete work, invalid certificates, and successful verification. A failed search does not prove absence of a symmetry. A matching cell count or file hash does not prove classification completeness.

Full perfect-form enumeration, independent cell completeness/realizability proofs, rebuilding all lower-rank computations, and a general higher-K-theory engine are outside version 1. The implementation should not replace the intended finite calculation with a citation of the final theorem.

## Design choices to preserve

1. Use exact integers and rationals for all certification decisions.
2. Keep witness discovery separate from the small certificate checker.
3. Pin original data, preserve source bytes and attribution, record hashes and a deterministic normalization, and support offline checking.
4. Compute cell dimensions from the span of rank-one matrices. Do not trust filenames as definitions of homological degree.
5. Test orientation on the span of the cone, not merely on the ambient lattice or by unsupported permutation parity.
6. Represent external premises and computed claims separately, with source locations and rank/degree/coefficient instances.
7. Build the narrow executable verifier first. Add broader homology machinery only when it is actually needed.
8. Use observed benchmarks to assess cost.

## Source and convention findings from prompt preparation

The authors' repository is <https://github.com/elbazvip/Voronoi-complexes-database>. It contains GAP and PARI data, including rank-8 files labeled `rnk8`, `rnk9`, and higher. These are concrete acquisition leads; the implementation must pin a revision and inspect each format.

The standalone `cells_GL8_dim8.gp` visibly contains thirteen 8-by-8 matrices. Its label must not be assumed to mean projective dimension 8. For vectors `v_i`, calculate `c = rank_Q {v_i v_i^T}`; projective dimension is `c - 1`. Preliminary inspection confirmed that `DATA_gap/ListSHV_n8_rnk9` has 106 records with ambient rank 8 and cone dimension 9. The implementation must reproduce these checks on its pinned input. Revision and checksum are recorded in the synthesis.

Reference: Dutour Sikirić–Elbaz-Vincent–Kupers–Martinet, <https://arxiv.org/abs/1910.11598>, particularly the definitions and degree shifts in §1, the cell data in §4, and the application in §6. Exact numbering and conventions should be tied to the version used. The prompt includes further primary references for the admitted results.

## Development environment and model choice

Observed on 2026-10-28: Linux x86-64, AMD Ryzen 7 PRO 250, 8 physical cores / 16 threads, about 30 GiB total RAM, and about 618 GiB available disk. Python 3 and a C++ compiler were on the path; SageMath and GAP were not found on that path. Recheck before installation or benchmarking. The shared folder was not a Git repository at inspection.

The recommendation given to the user for the future implementation run was GPT-6.1 Sol with Extra high reasoning effort. This is a task-specific recommendation, not a correctness guarantee or a requirement to change any active session's settings. A separate correctness review against the eventual mathematical specification remains part of the intended workflow.

## Next steps

The mathematical report is now available for review. An explicitly requested implementation run should follow `../PROMPT_CODEX.md` together with the synthesis, and record measured outcomes and unresolved issues here or in a new dated implementation ledger.

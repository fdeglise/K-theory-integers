# Implementation run — 2026-10-28

The user invoked PROMPT_CODEX.md. Read both existing ledgers; no later
mathematical specification was present at the initial inspection. Preserved
those files and the prompt. The directory was initially not a Git repository;
repository initialization appeared later through separate shared work.

Implemented a standard-library Python verifier, restricted source parsers,
deterministic normalization, Gram-pruned discovery, independent exact witness
checking, an explicit dependency manifest, CLI commands, a hand-checkable
example, and focused tests. No separate mathematical report was written.

Pinned the authors' data at commit
`d4c7b40ad2c95850d16160b347953dbaa9241486` (2022-06-30).
Preserved the original two GAP files, standalone PARI matrix file, README,
and GPL-3.0 license, with URLs, sizes, hashes and attribution in provenance.

The dimension audit found `(N,s,c,j)=(8,8,8,7)` for each of 13 rnk8 records,
and `(8,9,9,8)` for each of 106 rnk9 records. The standalone thirteen matrices
also have c=8 and j=7. Documented the discrepancy between DEKM v1's projective
definitions and the indices in Table 4. The classification-to-catalogue
identification is an explicit admission; the table's orientation assertions
are not admitted. No GL-orbit equivalence test was claimed.

Consulted DEKM arXiv:1910.11598v1, EVGS's 2013 paper, and Weibel Chapter VI
through the web tool. The dependency manifest gives the particular locations
and rank/degree/ring instances, including the separate 2-primary argument.
DEKM's `GL_N(R)` finite-order wording is false literally; its integral-matrix
restriction (Remark 4.2) supplies the needed stabilizer-prime bound.

All 106 reversing integral witnesses were discovered and saved. The independent
checker recomputed and checked all witnesses, signed vector maps, dimensions,
coordinate actions and negative determinants. The recorded full run took
3.004 seconds (search 2.168, independent verification 0.777), one process,
peak RSS 21,844 KiB. Discovery visited 1,405 nodes in total, at most 67 for a
single cell. These are observations for this run, not performance guarantees.

Source findings: Kupers's note initially returned HTTP 403 and another primary
URL was unavailable. A later local reference copy was read and its SHA-256
matched the existing catalogue. Its equation (1), lower-rank discussion,
Lemmas 2.1 and 2.5 and Remark 2.4 cross-check this implementation's dependency
route; the theorem-only route was not used. The PARI counterparts
are named `rank8/rank9`, not `rnk8/rnk9`, and contain `U1=...`, `U2=...`
row-vector lists. Inspected them through GitHub; their original bytes were
not acquired after shell network access became unavailable. They are not
verification inputs and no full PARI/GAP byte-equivalence check is claimed.

During implementation, a synthesis, updated prompt/ledgers, AGENTS.md, and
local reference corpus appeared through separate shared work. Read and
reconciled `reports/k8-verification-synthesis/main.tex`, the new instructions,
and updated ledgers. The specification agrees with the implementation.
The Smith-form correction affects only an unimplemented extension. No
report or original historical ledger was edited by this implementation task.

Validation: all 25 focused unittest tests passed in ordinary mode and under
`python3 -O`. This includes the real dataset, every saved witness, explicit
rational induced coordinates, an even ray permutation with negative cone
orientation, input/certificate tampering, missing records, forced limits,
offline checking, and distinct CLI exit statuses. A separate saved-certificate
run took 0.886 seconds with peak RSS 19,264 KiB. Rechecked available RAM, disk,
CPU count and absent Sage/GAP/PARI; details are in `results/environment.json`.

The resulting K-theory statement is conditional on the complete explicit
external dependency list. This is exact software verification, not an
independent classification or a formal proof-assistant development.

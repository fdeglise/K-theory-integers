# Repository initialization

Date: 2026-10-28 (America/Los_Angeles).

The user authorized creation of a local Git repository, inclusion of the relevant reference PDFs, ordinary ignore rules, and a commit. Initialized the `main` branch. No remote was created or configured.

The initial snapshot includes the mathematical synthesis, implementation prompt, historical and current ledgers, the existing verifier source and tests, original and normalized cell data, provenance and upstream license, examples, saved certificates, and recorded results. Added root README, future-agent instructions, attributes, and ignore rules.

The reference corpus contains ten PDFs with a reading guide and machine-readable catalogue. Original source URLs, SHA-256 checksums, sizes, and page counts are recorded. Quillen's paper was extracted as printed pages 179–198 from a source-volume scan; both source and extract identities are retained in the catalogue. The original Arlettaz 1991 article was not acquired; the full 2000 survey is included and explicitly distinguished from it.

Generated report PDFs, LaTeX auxiliary files, Python caches and environments, local credentials/configuration, editor files, and temporary/output directories are ignored. The reference corpus PDFs are the explicit exception to the general PDF rule. Saved verification evidence under `results/` remains tracked.

Verification before committing: the existing 25-test suite passed, including the real-data witness checks. Corpus PDF parsing, checksums, catalogue coverage, documentation links, and staged-file exclusions are checked separately. This repository task does not claim an independent audit of all verifier code.

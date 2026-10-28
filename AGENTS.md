# Instructions for future agents

## Starting points

Read `README.md`, `ledgers/PROJECT.md`, and `reports/k8-verification-synthesis/main.tex` before substantive changes. The historical conversation is a record of development and can contain corrected claims. Consult `ledgers/SYNTHESIS-2026-10-28.md` for the degree convention and Smith-form correction.

## Mathematical and implementation boundaries

- Preserve the distinction between checked symmetry certificates and admitted classification, comparison, lower-rank, and arithmetic inputs.
- Compute cone and projective dimensions from exact ranks. The target data have ambient rank 8, cone dimension 9, and projective dimension 8.
- Test orientation on the span of rank-one matrices, not from the ambient determinant alone.
- Keep discovery separate from independent exact verification. Missing witnesses and interrupted searches must not produce complete success.
- Keep source data pinned, attributed, and verifiable offline. Retain originals and provenance; do not silently replace a dataset or overwrite saved evidence with a different run.
- This project does not establish a general algorithm or runtime guarantee for K_12 or K_16.

## Documentation and files

- Write concise research-level English. Introduce only necessary notation in dependency order, use active definitions, and use statement environments without optional subtitles when they improve the exposition.
- Cite established inputs by exact source location; explain the essential steps of deductions made here.
- Compile and inspect the report after substantive LaTeX edits. Keep generated output in `tmp/` or `output/`.
- Reference PDFs in `references/k-theory/pdfs/` are deliberately tracked despite the general PDF ignore rule. Update the catalogue, URLs, checksums, and reading guide when adding references. Do not remove local files merely because they are ignored.
- Preserve original historical ledgers. Record material decisions and corrections in the current ledger or a new dated note.
- The user authorized initialization and an initial commit. Remote creation or pushing requires a subsequent request.

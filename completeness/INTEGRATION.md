# Integration status

The original attachment is preserved in `upstream/voronoi-verification.zip`.
`provenance.json` records its checksum and the original extracted-file hashes.
The maintained Python files replace the attachment's assertions with explicit
validation that remains active under `python -O`. This supersedes the original
README's instruction about optimization mode. Verification summaries are
redirected to temporary output rather than overwriting the archived evidence.

From the repository root, run:

```sh
python3 -m voronoi_verify completeness
python3 -m voronoi_verify completeness-check
python3 -m unittest discover -s tests -v
```

The first command checks the enumeration, exported matrices, and initial
constraints, converts the accepted target presentations to integral lattice
coordinates, and submits their symmetry certificates to the existing independent
checker. The second command replays only those saved orientation certificates;
it does not establish enumeration completeness by itself.

The recorded integrated run passed and checked 134 rank-eight presentations
and 44 rank-nine presentations independently. The 31-test suite passed,
including rejection of invalid certificates under optimized Python.

The LaTeX report still describes the original published-list route. Its planned
extension with the completeness reductions and revised assumptions was not
written before the user interrupted work to request this commit. The integrated
dependency manifest identifies these reductions as written mathematical inputs;
their promised exposition remains pending. Do not interpret this implementation
snapshot as a completed proof audit or a formal proof-assistant verification.

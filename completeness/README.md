# Voronoi completeness verification

This bundle contains the calculation discussed in the K_8(Z) project: the degree-eight Voronoi cells meeting the positive definite locus in lattice ranks eight and nine, modulo GL_N(Z).

## Run the exact verifier

Python 3.9 or later is required. Verification uses Python's standard library and the included local source files. It does not require an optimization package or internet access.

From the extracted folder, run:

```sh
python3 verify_completeness.py
```

Do not use Python's `-O` option: the verification checks use assertions.

A successful run prints 13 configurations with eight independent minimal pairs in rank eight, 44 cell classes in rank nine, and 106 cell classes in rank eight. The rank-eight extension search has 647 candidates, of which 134 are realizable presentations and 513 are excluded. All 647 candidates have a checked orientation-reversing transposition. An assertion failure or another exception means that the verification did not complete successfully.

## Files

- `matrices/README.md`: guide to the full matrix export, including coordinates, provenance, and exact reproduction commands. `matrices/all_matrices.txt` is the readable catalogue; `matrices/all_matrices.json` contains every retained Gram matrix and exclusion certificate with labelled reconstructions. CSV indexes identify all 150 final cell classes. The compressed `matrices/initial_lp_matrices.jsonl.gz` contains all 964 reconstructed initial constraint systems. Intermediate solver states were not retained.
- `export_all_matrices.py`, `export_initial_lp.py`, `verify_matrix_export.py`, and `verify_initial_lp_export.py`: export and check the complete matrix data. Their default modes need only Python's standard library.
- `matrix_examples.py` and `matrix_examples.txt`: an executable display and its output, showing actual saved Gram matrices, an integral change of basis, positive definiteness, the complete shortest-vector check, a tensor-rank computation, orientation-reversing automorphisms, and a weighted exclusion-certificate calculation. Run `python3 matrix_examples.py` to reproduce the transcript.
- `verify_completeness.py`: regenerates candidates and checks rational exclusion certificates, positive definiteness, exact shortest-vector enumeration, equivalence keys, and orientation symmetries.
- `noncyclic.py`: finite quotient-code enumeration used by the verifier.
- `*_candidates*.json`, `cyclic_accepted.json`, and `*_results*.json`: enumeration data, rational Gram matrices, and exclusion certificates.
- `verification_summary.json`: numerical summary of a successful verification run.
- `compare_published.py`: optional comparison with the published lists at a fixed Git commit. This step uses internet access and is separate from the independent verifier.
- `comparison_summary.json`: the recorded comparison, with no missing or extra classes.
- `enumerate_cyclic.js`, `check_cyclic.py`, `check_code.py`, and `extend_rank8.py`: original candidate-generation and feasibility-search programs. The feasibility searches require NumPy, SciPy, SymPy, and mpmath; these packages are unnecessary for verification of the supplied certificates.

## Mathematical scope

Completeness uses the geometric reductions explained in the conversation: nine independent rank-one rays; index bounds 19 and 37; Watson's identity and the finite-abelian-group reductions; and the coordinate bound |x_i| <= 1 for a ninth vector relative to a maximal-determinant minimal basis. The verifier checks the resulting finite enumeration and certificates. This is a computer-assisted calculation, not a Lean formalization. The Quillen, Hurewicz, and small-prime inputs needed for the full K_8(Z) argument are separate.

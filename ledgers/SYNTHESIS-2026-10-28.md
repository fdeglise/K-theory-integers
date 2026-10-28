# Mathematical synthesis: preparation and checks

Date: 2026-10-28 (America/Los_Angeles).

The user requested a human-readable LaTeX synthesis of the exact algorithm, admitted data/theorems, and certification claim. Source: `../reports/k8-verification-synthesis/main.tex`. The companion PDF is generated under `../output/pdf/`.

## Expository guidance consulted

Read all accessible user prompts in the RostNilpotence chat “Collect Rost nilpotence references,” together with `/home/frederic/codex/Rnilpotence/AGENTS.md` and the existing synthesis. Applied the user's preferences: concise research-level English, active introductions of definitions, minimal notation in dependency order, selective statement environments without optional subtitles, restrained sectioning, exact references, proofs focused on the essential deductions, and the strongest useful statement at the end.

## Mathematical scope

The report specifies a finite Gram-pruned search for reversing integral symmetries and an independent exact witness checker. It proves the orientation-sign criterion and theoretical search completeness. A final conditional theorem identifies the result of successful coverage checking and its passage to K_8(Z) = 0 under explicitly admitted classification and published comparison/arithmetic results.

Two independent checks of the mathematics confirmed the degree shifts and the distinction between the lattice determinant and the cone-orientation determinant. An additional review checked the exposition and conditional claims. No symmetry certificates were computed in this report-writing task.

## Data inspection

- Repository: <https://github.com/elbazvip/Voronoi-complexes-database>.
- Revision: `d4c7b40ad2c95850d16160b347953dbaa9241486`.
- File: `DATA_gap/ListSHV_n8_rnk9`.
- SHA-256: `d0b58f52da26e64b20649dffe4e317efd31a9d6b0c5e4d5b1407dbf5d636359d`.
- Size: 30,025 bytes.
- Preliminary exact parsing/rank inspection found 106 records, each with nine primitive vectors spanning Q^8 and nine independent rank-one matrices. Thus cone dimension is 9 and projective dimension is 8.

This establishes elementary properties of the inspected file, not classification completeness, realizability, orbit inequivalence, or reversing symmetries. The implementation must repeat the relevant checks and bind its certificates to its retained source bytes. Coverage of all relevant cell orbits remains an admitted input; pairwise inequivalence is unnecessary for vanishing.

## Verified published dependencies

- DEKM, Proposition 1.5: Voronoi–Steinberg comparison after removing stabilizer primes.
- DEKM, Table 5: relevant lower-rank entries; EVGS, Theorem 4.3, gives the degree-eight vanishing in ranks 5, 6, and 7 explicitly.
- Lee–Szczarba, Theorem 1.3, as cited in DEKM Remark 6.2: rank-nine Steinberg coinvariant vanishing.
- Quillen, *Finite generation*, Theorem 3; DEKM §6 and proof of Proposition 6.1: rank filtration and total-degree-nine argument.
- DEKM, Proposition 6.3, using Arlettaz, Theorem 1.5: degree-nine Hurewicz kernel supported at 2, 3, and 5.
- Weibel VI, Corollary 9.8: absence of 2-primary torsion in K_8(Z).
- Weibel VI, Proposition 10.5: absence of torsion at odd regular primes in positive even K-groups; apply at 3, 5, 7.

The report cites the author manuscript of DEKM for numbering. The Kupers PDF was unavailable during the dependency check; the proof chain uses the accessible sources above and does not rely on unverified lemma numbers from that note.

## Correction to earlier material

The historical discussion and original implementation prompt said that separate Smith forms of consecutive differentials do not determine abstract integral homology. For a complex of finite free abelian groups this is incorrect once the chain identity is known. If A B = 0, then

`0 -> H -> coker(B) -> im(A) -> 0`

splits because `im(A)` is free. Thus the torsion invariant factors of H are those of coker(B), and the free rank is `rank(C_n) - rank(A) - rank(B)`. Kernel coordinates remain useful for explicit cycles. EVGS §4.4 states the invariant-factor computation. Corrected the active prompt and retained the historical conversation unchanged. The synthesis omits this machinery because its target chain group is zero.

## Building

From the report directory, run:

```sh
latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=../../tmp/pdfs/k8-synthesis main.tex
```

The local TeX installation compiles the report. The Codex built-in preview compiler failed because its isolated environment could not find its LaTeX format; this was an environment failure rather than a source error. The delivered PDF is compiled with the working local installation and visually inspected.

Final validation: seven pages, all rendered and inspected; no LaTeX warnings, unresolved references/citations, or overfull/underfull boxes. Delivered PDF: `../output/pdf/k8-verification-synthesis.pdf`.

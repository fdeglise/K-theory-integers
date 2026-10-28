# Reference corpus for the K_8 verifier

This corpus contains ten PDFs covering the finite computation, its admitted mathematical inputs, and the classification methods discussed in the project. [catalog.json](catalog.json) records source URLs, retrieval information, page counts, file sizes, and SHA-256 checksums. The original rights remain with the authors and publishers.

## Suggested reading order

1. [DEKM: Voronoi complexes and K_8](pdfs/01_DEKM_Voronoi-K8.pdf). Start with §§1, 4, and 6: cell conventions, orientation contributions, and the passage to K-theory. This is the author manuscript used for theorem numbering in the synthesis.
2. [Elbaz-Vincent–Gangl–Soulé: Perfect forms](pdfs/02_EVGS_Perfect-forms.pdf). Theorem 4.3 supplies the earlier degree-eight vanishing in ranks 5–7. The paper also explains boundary calculations and exact homology methods. The downloaded PDF includes an offprint cover.
3. [Weibel: The K-book, Chapter VI](pdfs/03_Weibel_Kbook-VI.pdf). Corollary 9.8 and Proposition 10.5 supply the small-prime arithmetic inputs. This file is the chapter on the higher K-theory of fields.
4. [Quillen: Finite generation](pdfs/04_Quillen_Finite-generation.pdf). The original rank-filtration source. The retained PDF is an extraction of printed pages 179–198 from an institutional scan of the proceedings, rather than the entire volume. Source-volume and extracted-file identities are both recorded in the catalogue.
5. [Lee–Szczarba: Congruence subgroups](pdfs/10_Lee-Szczarba_Congruence-subgroups.pdf). Original source for the Steinberg coinvariant result cited in the synthesis. The Göttingen archive scan includes a library cover; the article occupies printed pages 15–54.
6. [Kupers: A short proof](pdfs/09_Kupers_Short-proof-K8.pdf). A useful comparison with the computational route. Its alternative vanishing theorem does not replace the certificate checks requested for this program.

## Supporting references

- [Rognes–Weibel: Two-primary K-theory](pdfs/05_Rognes-Weibel_Two-primary.pdf): underlying computations at the prime 2.
- [Schürmann: Enumerating perfect forms](pdfs/06_Schuermann_Enumerating-perfect-forms.pdf): the general Voronoi enumeration algorithm, beyond the first verifier's scope.
- [Keller–Martinet–Schürmann: Minkowskian sublattices](pdfs/07_Keller-Martinet-Schuermann_Minkowskian-sublattices.pdf): methods relevant to the admitted cell classification.
- [Arlettaz: A topological viewpoint](pdfs/08_Arlettaz_Topological-viewpoint.pdf): the full 82-page author survey, including Hurewicz maps in §6. This is a separate work from the 1991 paper cited in the synthesis.

## Acquisition limitation

The original 1991 Arlettaz paper, *The Hurewicz homomorphism in algebraic K-theory*, was not acquired. Its DOI and status are recorded in the catalogue. The included survey is identified explicitly as a different reference; DEKM Proposition 6.3 records the specific Hurewicz consequence used by the project.

## Reproducibility and repository rules

All PDFs in this directory are intentionally committed. Generated report PDFs elsewhere are ignored. Check the retained bytes against the catalogue when changing or relocating the corpus. A checksum establishes file identity, not the truth of its mathematical statements or completeness of the cell classification.

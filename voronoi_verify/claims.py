"""Explicit, externally supplied premises and conditional inference graph.

This is a reviewed dependency manifest, not a theorem prover. No rank-eight
orientation or homology vanishing, nor K_8 itself, is an admitted premise.
"""

DEKM = "https://arxiv.org/html/1910.11598v1"
EVGS = "https://www.ihes.fr/~soule/documents/YAIMA4438.pdf"
WEIBEL = "https://sites.math.rutgers.edu/~weibel/Kbook/Kbook.VI.pdf"
KUPERS = "https://www.utsc.utoronto.ca/people/kupers/wp-content/uploads/sites/50/2021/01/k8zshorter.pdf"
R = "Z[1/7!] = Z[1/210]"


def premise(cid, statement, source, location, instance, **extra):
    return dict(id=cid, status="ADMITTED", role="external mathematical premise; not machine checked",
                statement=statement, source=source, location=location, instance=instance, **extra)


def premises():
    return [
        premise("classification", "The pinned cone-rank-9 catalogue represents all GL_8(Z) orbits of full-span Voronoi cones of dimension 9 and each listed configuration is realizable with exactly the stated minimal-vector pairs.",
                DEKM, "v1 (25 October 2019), sections 2-3, Proposition 3.1/Table 3, Proposition 4.3/Table 4 and data link following Table 4",
                {"rank": 8, "cone_dimension": 9, "projective_dimension": 8, "coefficient_ring": R},
                qualification="Catalogue-to-classification identification is admitted. Table 4's indices are shifted relative to section 1.2; dimensions are recomputed from data. No orbit-inequivalence or realization proof is implemented."),
        premise("comparison", "Finite stabilizers of cells meeting the positive-definite locus have prime divisors at most N+1. After these primes are inverted, the orientation coinvariant complex computes Steinberg homology with shift j-N+1; reversing stabilizers kill their degree-j coinvariants when 2 is inverted.",
                DEKM, "sections 1.2.2-1.2.3, equations (3)-(4), Proposition 1.5; section 4.1, Remark 4.2",
                {"rank": 8, "voronoi_degree": 8, "steinberg_degree": 1, "coefficient_ring": R},
                localization="For N=8 the primes <=9 are exactly 2,3,5,7. The GL_N(R) wording in Lemma 4.1 is erroneous if read literally; use the integral-matrix restriction in Remark 4.2."),
        premise("lower_rank_1", "H_8(GL_1(Z); St_1) tensor R = 0; St_1=Z with trivial action and GL_1(Z)=C_2.",
                DEKM, "Table 5 (N=1,k=8), section 6.1; finite-group homology after inverting the order",
                {"rank": 1, "steinberg_degree": 8, "total_degree": 9, "coefficient_ring": R}),
        premise("lower_ranks_2_4", "For N=2,3,4, H_{9-N}(GL_N(Z);St_N) tensor R=0. Equivalently H_8(Vor_GL_N) tensor R=0.",
                EVGS, "section 6.3 (small-rank complexes), sections 3.3-3.4 (comparison); also DEKM v1 Table 5, entries (N,k)=(2,7),(3,6),(4,5); Kupers (21 October 2019), p.2, refers to Lee-Szczarba (1978), p.104",
                {"ranks": [2,3,4], "steinberg_degrees": [7,6,5], "voronoi_degree": 8, "coefficient_ring": R},
                explanation="N=2,3: projective symmetric-space dimensions 2,5 are below 8; N=4: the small-rank complex in section 6.3 has only the degree-6 surviving cell. All relevant stabilizer primes <=5 are inverted."),
        premise("lower_ranks_5_7", "H_8(Vor_GL_N) tensor R=0 for N=5,6,7, hence H_{9-N}(GL_N;St_N) tensor R=0.",
                EVGS, "Theorem 4.3, p.597 (PDF p.12), with comparison in section 3.4; cross-check DEKM v1 Table 5",
                {"ranks": [5,6,7], "steinberg_degrees": [4,3,2], "voronoi_degree": 8, "coefficient_ring": R},
                localization="The N=5,6 calculations are modulo S_5, N=7 modulo S_7; comparison for N<=7 requires primes <=8, all inverted in R."),
        premise("rank_9", "H_0(GL_9(Z); St_9)=0 integrally, hence also after tensoring with R.",
                DEKM, "Remark 6.2, citing Lee-Szczarba, On the homology and cohomology of congruence subgroups (1976), Theorem 1.3",
                {"rank": 9, "steinberg_degree": 0, "total_degree": 9, "coefficient_ring": R},
                verification_scope="Checked in DEKM's primary application and Kupers (21 October 2019), Lemma 2.1; the original 1976 theorem remains an external dependency, not a software computation."),
        premise("rank_filtration", "The first-quadrant rank filtration has E^1_{N,q}=H_q(GL_N(Z);St_N), converging in total degree N+q to H_*(BQZ;Z). BQ_0 is a point. A zero total-degree-9 diagonal gives H_9(BQZ;R)=0.",
                DEKM, "section 6, Quillen exact sequence (8); proof of Proposition 6.1",
                {"total_degree": 9, "coefficient_ring": R},
                diagonal=[{"rank": n, "q": 9-n, "input": ("computed_rank_8" if n==8 else "rank_9" if n==9 else "lower_rank_1" if n==1 else "lower_ranks_2_4" if n<=4 else "lower_ranks_5_7")} for n in range(1,10)],
                higher_ranks="N>=10 would require q<0, hence do not contribute; rank 0 has H_9(point)=0."),
        premise("hurewicz", "K_8(Z)=pi_9(BQZ). The Hurewicz kernel in degree 9 is killed by an integer with prime divisors only 2,3,5; it vanishes after tensoring with R.",
                DEKM, "section 6.2, Proposition 6.3 and its proof citing Arlettaz, Theorem 1.5",
                {"homotopy_degree": 9, "K_degree": 8, "coefficient_ring": R}),
        premise("arithmetic_2", "The 2-primary subgroup of K_8(Z) is zero (n=8 mod 8, n>=2).",
                WEIBEL, "Chapter VI, Corollary 9.8, p.55 (PDF p.55); separate real-place/2-primary calculation",
                {"K_degree": 8, "prime": 2, "coefficient_ring": "Z"}),
        premise("arithmetic_odd", "For every odd regular prime ell, K_{2i}(Z) has no ell-torsion. Apply at ell=3,5,7, i=4.",
                WEIBEL, "Chapter VI, Proposition 10.5, p.63 (PDF p.63)",
                {"K_degree": 8, "primes": [3,5,7], "coefficient_ring": "Z"}),
        premise("regularity", "Kummer's criterion identifies regular odd primes via Bernoulli numerators. The primes 3,5,7 are regular: the relevant topologists' Bernoulli numbers B_1=1/6 and B_2=1/30 have numerator 1; for 3 the list is empty.",
                WEIBEL, "Chapter VI, section 2, Irregular Primes 2.4.1 and preceding Bernoulli table, p.9",
                {"primes": [3,5,7], "coefficient_ring": "Z"}),
    ]


def manifest(complete_published):
    computed = [{"id": "checked_orientation_certificates", "status": "COMPUTED",
                 "statement": "Exact integral signed vector maps and determinant -1 on the rational rank-one span for each reported checked cell."}]
    inferences = [
        {"id": "voronoi_vanishing", "status": "CONDITIONAL" if complete_published else "UNRESOLVED",
         "statement": "H_8(Vor_GL_8;R)=0 because its orientation-preserving free chain group V_8 is zero.",
         "depends_on": ["checked_orientation_certificates", "classification"],
         "requires": "all 106 records in the pinned target list independently checked"},
        {"id": "computed_rank_8", "status": "CONDITIONAL" if complete_published else "UNRESOLVED",
         "statement": "H_1(GL_8(Z);St_8 tensor R)=0", "depends_on": ["voronoi_vanishing", "comparison"]},
        {"id": "BQ_vanishing", "status": "CONDITIONAL" if complete_published else "UNRESOLVED",
         "statement": "H_9(BQZ;R)=0", "depends_on": ["computed_rank_8", "lower_rank_1", "lower_ranks_2_4", "lower_ranks_5_7", "rank_9", "rank_filtration"]},
        {"id": "K8_localized", "status": "CONDITIONAL" if complete_published else "UNRESOLVED",
         "statement": "K_8(Z) tensor R=0", "depends_on": ["BQ_vanishing", "hurewicz"]},
        {"id": "K8_integral", "status": "CONDITIONAL" if complete_published else "UNRESOLVED",
         "statement": "K_8(Z)=0 under the listed external mathematical premises",
         "depends_on": ["K8_localized", "arithmetic_2", "arithmetic_odd", "regularity"],
         "reason": "Localization zero makes each element annihilated by a power of 7!; the four torsion exclusions remove it."},
    ]
    return {"schema": "voronoi-dependencies-v1", "coefficient_ring": R,
            "computed": computed, "admitted": premises(), "inferences": inferences,
            "unresolved_required": [] if complete_published else ["complete published certificate coverage"],
            "source_limitations": [
                "PARI counterpart text inspected through GitHub; original bytes were not acquired after shell network access became unavailable. They are not proof inputs.",
            ],
            "source_cross_checks": [{"source": KUPERS, "version": "21 October 2019, 4 pages",
                "local_file": "references/k-theory/pdfs/09_Kupers_Short-proof-K8.pdf",
                "sha256": "2ceb04c2aa64fdb2530591adc96b7575bb64d238ef55e9b239acf454da5b92b5",
                "checked": "Equation (1), p.1; lower-rank references and Lemma 2.1, p.2; Lemma 2.5 and Remark 2.4, p.3. These agree with the degree-nine diagonal and the computational route. Lemma 2.2 and the alternate fiber-sequence proof are not used to replace the rank-eight calculation.",
                "access_note": "Direct web retrieval returned 403, but a local reference copy appeared later in the shared workspace and was read."}],
            "specification_reconciliation": {
                "file": "reports/k8-verification-synthesis/main.tex",
                "finding": "The later synthesis appeared during implementation and was read. Its conventions, orientation criterion, search, and conditional theorem agree with this program. Its Smith-form correction concerns an unused extension; no boundary/homology engine is implemented."},
            "certification_scope": "Ordinary exact software verification conditional on external mathematics; not independent classification, realization, or formal proof-assistant verification."}

"""Bridge the attached finite enumeration to the existing independent checker.

No enumeration module is imported by the orientation-certificate replay path.
The published-list workflow remains available independently.
"""

import hashlib
import json
import os
import subprocess
import sys
import tempfile
from fractions import Fraction
from pathlib import Path

from .checker import certificate_for, verify_bundle
from .claims import premises, R
from .dataset import ROOT, digest, make_dataset, read_json, write_json
from .exact import Invalid, require

COMPONENT = ROOT / "completeness"


def integrated_dependencies():
    retained = [p for p in premises() if p["id"] not in ("classification", "rank_9")]
    for p in retained:
        if p["id"] == "comparison":
            p["instance"] = {"ranks": [8, 9], "voronoi_degree": 8,
                             "steinberg_degrees": [1, 0], "coefficient_ring": R}
            p["localization"] = "Stabilizer primes <= N+1 <= 10 are 2,3,5,7."
        if p["id"] == "rank_filtration":
            for term in p["diagonal"]:
                if term["rank"] == 9:
                    term["input"] = "computed_rank_9"
    return {
        "schema": "voronoi-integrated-dependencies-v1",
        "coefficient_ring": R,
        "admitted": retained + [{
            "id": "simpliciality", "status": "ADMITTED",
            "statement": "A full-span perfect cone of dimension N or N+1 is simplicial.",
            "source": "https://doi.org/10.14231/AG-2015-027",
            "location": "Dutour Sikiric–Hulek–Schuermann (2015), Lemma 3.1, pp.647–649; not the classification-dependent Theorem 1.1.",
            "instance": {"ranks": [8, 9], "cone_dimension": 9},
        }],
        "proved_in_report_not_formalized": [
            "Blichfeldt index bound: index <=19 in rank8 and <=37 in rank9.",
            "Strict Watson filter for a lattice with exactly its coordinate minimal pairs.",
            "Finite abelian quotient coverage, including lifts of C4 x C4 and exclusion of C2^3, C3^2, C5^2.",
            "A full-span facet can be realized with only its prescribed minimal pairs.",
            "Maximal-determinant basis implies |x_i|<=1 for the extra vector.",
            "Averaging preserves realizability and justifies the restricted Gram parameters.",
            "Exact LDL enumeration and rational inequality combinations certify realization/exclusion.",
        ],
        "computed": [
            "Cyclic and noncyclic candidate regeneration and exact realization/exclusion checks.",
            "Every bounded rank8 extension regenerated; its coordinate transposition checked.",
            "All relevant accepted presentations converted to integral bases and checked by voronoi_verify.checker.",
        ],
        "inferences": [
            {"id": "computed_rank_8", "status": "CONDITIONAL",
             "statement": "H_1(GL_8(Z);St_8 tensor R)=0 from complete geometric coverage and independent orientation checks."},
            {"id": "computed_rank_9", "status": "CONDITIONAL",
             "statement": "H_0(GL_9(Z);St_9 tensor R)=0 from the analogous geometric calculation; integral vanishing is not asserted by this route."},
            {"id": "K8_integral", "status": "CONDITIONAL",
             "statement": "K_8(Z)=0 under the cited external theorems, reductions proved in the report, and correctness of the exact software execution."},
        ],
        "classification_assumption": "No use of completeness/realizability of the published 106-cell list. No use of exact orbit counts for vanishing.",
        "remaining_trust": "Published theorems and the written mathematical reductions are not formalized; Python, its integer/rational arithmetic, the checkers, and their correspondence to the reductions remain trusted.",
        "report": "reports/k8-verification-synthesis/main.tex",
    }


def integer_matrix(a):
    require(type(a) is list and a and all(type(row) is list for row in a), "invalid exported matrix")
    result = []
    for row in a:
        values = [Fraction(x) for x in row]
        require(all(x.denominator == 1 for x in values), "nonintegral exported matrix")
        result.append([int(x) for x in values])
    return result


def export_witnesses(path):
    """Convert all accepted target presentations, not one chosen per orbit."""
    raw = Path(path).read_bytes()
    records = json.loads(raw)["records"]
    result, counts, seen = [], {8: 0, 9: 0}, set()
    for rec in records:
        n = rec["rank"]
        if rec["status"] != "accepted" or not (n == 9 or n == 8 and rec["phase"] == "extension"):
            continue
        ident = rec["id"]
        require(ident not in seen, "duplicate exported presentation")
        seen.add(ident)
        vectors = integer_matrix(rec["derived_geometry"]["minimal_vectors_in_integral_basis"])
        g = integer_matrix(rec["symmetry"]["P_in_integral_basis"])
        d = make_dataset([vectors], "completeness/" + ident, hashlib.sha256(raw).hexdigest())
        cell = d["cells"][0]
        require(cell["ambient_rank"] == n and cell["cone_dimension"] == 9,
                "exported target has incorrect rank or projective degree")
        b = {"schema": "voronoi-certificates-v1", "dataset_sha256": d["dataset_sha256"],
             "certificates": [certificate_for(d, cell, g)]}
        result.append({"source_record_id": ident, "dataset": d, "bundle": b})
        counts[n] += 1
    require(counts == {8: 134, 9: 44}, "missing or extra accepted target presentations")
    return {"schema": "voronoi-completeness-witnesses-v1",
            "export_sha256": hashlib.sha256(raw).hexdigest(), "presentations": result}


def replay_witnesses(bundle):
    require(type(bundle) is dict and set(bundle) == {"schema", "export_sha256", "presentations"}
            and bundle["schema"] == "voronoi-completeness-witnesses-v1", "unknown interoperation schema")
    require(type(bundle["presentations"]) is list and bundle["presentations"], "empty presentation list")
    counts, seen = {8: 0, 9: 0}, set()
    for p in bundle["presentations"]:
        require(type(p) is dict and set(p) == {"source_record_id", "dataset", "bundle"}, "invalid presentation fields")
        ident, d = p["source_record_id"], p["dataset"]
        require(type(ident) is str and ident not in seen, "duplicate presentation ID")
        seen.add(ident)
        require(d["source"] == "completeness/" + ident and d["source_sha256"] == bundle["export_sha256"],
                "presentation has a mismatched source")
        checked, missing = verify_bundle(d, p["bundle"])
        require(len(checked) == 1 and not missing, "incomplete independent presentation certificate")
        c = checked[0]
        require(c["ambient_rank"] in counts and c["cone_dimension"] == 9, "wrong target dimension")
        counts[c["ambient_rank"]] += 1
    return {"status": "SUCCESS", "result_level": "CHECKED_SUPPLIED_COMPLETENESS_PRESENTATIONS",
            "counts_by_rank": {str(n): k for n, k in counts.items()},
            "checked": sum(counts.values()), "bundle_sha256": digest(bundle),
            "scope": "Independent orientation certificates only; replay alone establishes neither complete enumeration nor K-theory vanishing."}


def run_integrated(certificates, timeout=120):
    require(timeout > 0, "phase timeout must be positive")
    phase_results = {}
    scripts = [("enumeration", "verify_completeness.py", "verification_summary.json"),
               ("matrices", "verify_matrix_export.py", "matrix_verification_summary.json"),
               ("initial_constraints", "verify_initial_lp_export.py", "initial_lp_verification_summary.json")]
    with tempfile.TemporaryDirectory(prefix="voronoi-completeness-") as tmp:
        env = dict(os.environ, VORONOI_OUTPUT_DIR=tmp, PYTHONHASHSEED="0")
        for phase, script, filename in scripts:
            command = [sys.executable] + (["-O"] if sys.flags.optimize else []) + [str(COMPONENT / script)]
            try:
                outcome = subprocess.run(command, cwd=COMPONENT, env=env, capture_output=True,
                                         text=True, timeout=timeout)
            except subprocess.TimeoutExpired:
                return {"status": "INCOMPLETE", "phase": phase, "error": "verification phase timed out",
                        "completed_phases": phase_results, "K8_conclusion": "UNRESOLVED"}
            require(outcome.returncode == 0,
                    f"{phase} verification failed: {outcome.stderr[-4000:] or outcome.stdout[-4000:]}")
            phase_results[phase] = read_json(Path(tmp) / filename)
        # The preserved input/export hashes are checked by the matrix phase.
        bundle = export_witnesses(COMPONENT / "matrices" / "all_matrices.json")
        write_json(certificates, bundle)
        replay = replay_witnesses(read_json(certificates))
    return {"schema": "voronoi-integrated-result-v1", "status": "SUCCESS",
            "result_level": "CHECKED_ENUMERATION_AND_INDEPENDENT_ORIENTATION_CERTIFICATES",
            "phases": phase_results, "independent_checker": replay,
            "input_provenance": read_json(COMPONENT / "provenance.json"),
            "certificates": str(certificates), "mathematical_dependencies": integrated_dependencies(),
            "K8_conclusion": "CONDITIONAL_ON_REPORTED_EXTERNAL_THEOREMS_AND_WRITTEN_REDUCTIONS"}

"""Cross-component checks and adversarial validation of the imported verifier."""

import ast
from copy import deepcopy
import json
from pathlib import Path
import subprocess
import sys
import unittest

from voronoi_verify.completeness import (COMPONENT, export_witnesses,
                                         integrated_dependencies, replay_witnesses)
from voronoi_verify.exact import Invalid


class CompletenessTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.bundle = export_witnesses(COMPONENT / "matrices/all_matrices.json")

    def test_all_accepted_presentations_use_our_checker(self):
        result = replay_witnesses(self.bundle)
        self.assertEqual(result["counts_by_rank"], {"8": 134, "9": 44})
        self.assertEqual(result["checked"], 178)
        self.assertIn("neither complete enumeration", result["scope"])

    def test_invalid_and_missing_witnesses_cannot_pass(self):
        for mutation in ("matrix", "source", "duplicate", "missing"):
            bundle = deepcopy(self.bundle)
            if mutation == "matrix":
                bundle["presentations"][0]["bundle"]["certificates"][0]["g"][0][0] += 1
            elif mutation == "source":
                bundle["export_sha256"] = "0" * 64
            elif mutation == "duplicate":
                bundle["presentations"].append(bundle["presentations"][0])
            else:
                bundle["presentations"][0]["bundle"]["certificates"] = []
            with self.subTest(mutation=mutation), self.assertRaises(Invalid):
                replay_witnesses(bundle)

    def test_maintained_checks_have_no_assert_statements(self):
        for path in COMPONENT.glob("*.py"):
            with self.subTest(file=path.name):
                self.assertFalse(any(isinstance(n, ast.Assert) for n in ast.walk(ast.parse(path.read_text()))))

    def test_optimized_mode_rejects_bad_forms_and_duals(self):
        # The bad form has a zero diagonal pivot. The negative dual multiplier
        # violates a necessary sign condition independently of its claimed sum.
        program = '''
import json
from copy import deepcopy
from verify_completeness import verify_form, verify_result, basis_vectors
from validation import VerificationError, require
failures = 0
try:
    verify_form([[0,0],[0,1]], {(0,0)}, basis_vectors(2))
except VerificationError:
    failures += 1
records = json.load(open('cyclic_results_8.json'))
r = deepcopy(next(r for r in records.values() if r['status'] == 'excluded'))
row = list(r['certificate']['constraints'][0]); row[2] = '-1'
r['certificate']['constraints'][0] = row
try:
    verify_result(8, r['d'], [r['a']], r, cyclic_format=True)
except VerificationError:
    failures += 1
require(failures == 2, 'bad certificates accepted under -O')
require(True, lambda: 1/0)
'''
        outcome = subprocess.run([sys.executable, "-O", "-c", program], cwd=COMPONENT,
                                 capture_output=True, text=True)
        self.assertEqual(outcome.returncode, 0, outcome.stderr)

    def test_dependency_route_replaces_classification_not_external_theorems(self):
        dependencies = integrated_dependencies()
        ids = {p["id"] for p in dependencies["admitted"]}
        self.assertNotIn("classification", ids)
        self.assertNotIn("rank_9", ids)
        self.assertTrue({"simpliciality", "comparison", "hurewicz", "arithmetic_2"} <= ids)
        rank_filtration = next(p for p in dependencies["admitted"] if p["id"] == "rank_filtration")
        self.assertEqual(rank_filtration["diagonal"][-1]["input"], "computed_rank_9")


if __name__ == "__main__":
    unittest.main()

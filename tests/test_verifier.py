"""Focused mathematical checks plus adversarial coverage and CLI tests.

Run with python3 -m unittest discover -s tests -v; also run under python3 -O.
The expected signs below come from explicit diagonal actions on Sym_2 and
coordinate-ray swaps, not from discovery's output.
"""

from copy import deepcopy
from fractions import Fraction
from itertools import permutations
from pathlib import Path
import contextlib
import io
import socket
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from voronoi_verify.__main__ import main, verify
from voronoi_verify.checker import certificate_for, check_witness, verify_bundle
from voronoi_verify.dataset import (ROOT, TARGET, digest, import_pinned, make_cell,
                                  make_dataset, read_json, validate_dataset, write_json)
from voronoi_verify.discover import find_witness
from voronoi_verify.exact import (Invalid, det, eye, inverse, mul, mv, outer, rank,
                                rational_string, span_action, transpose)
from voronoi_verify.formats import parse_source


def dataset(vectors):
    return make_dataset([vectors], "test", digest(vectors))


def bundle(d, g):
    return {"schema": "voronoi-certificates-v1", "dataset_sha256": d["dataset_sha256"],
            "certificates": [certificate_for(d, d["cells"][0], g)]}


class ExactTests(unittest.TestCase):
    def test_known_determinants_and_inverse(self):
        # Independent Leibniz expansion is adequate for these small examples.
        for a in ([[2,3],[5,7]], [[1,2,3],[0,4,5],[1,0,6]], [[0,1,0],[0,0,1],[1,0,0]]):
            n = len(a)
            total = 0
            for p in permutations(range(n)):
                sign = (-1) ** sum(p[i] > p[j] for i in range(n) for j in range(i+1,n))
                term = sign
                for i in range(n):
                    term *= a[i][p[i]]
                total += term
            self.assertEqual(det(a), total)
            self.assertEqual(mul(a, inverse(a)), eye(n))
        self.assertEqual(det([[1,2],[2,4]]), 0)
        with self.assertRaises(Invalid):
            inverse([[1,2],[2,4]])

    def test_coordinate_transposition(self):
        d = dataset([[1,0],[0,1]])
        b = bundle(d, [[0,1],[1,0]])
        c = b["certificates"][0]
        self.assertEqual(c["induced_action"], [["0","1"],["1","0"]])
        self.assertEqual(c["orientation_determinant"], "-1")
        self.assertEqual(len(verify_bundle(d,b)[0]), 1)

    def test_identity_minus_identity_and_ambient_sign(self):
        for n in (1,2,3):
            v = eye(n)
            for g in (eye(n), [[-x for x in row] for row in eye(n)]):
                self.assertEqual(span_action(v,g)[2], 1)
                with self.assertRaises(Invalid):
                    check_witness(dataset(v)["cells"][0],g)
        # Ambient reflection has determinant -1, yet fixes every rank-one ray.
        g = [[-1,0],[0,1]]
        self.assertEqual(det(g), -1)
        self.assertEqual(span_action(eye(2),g)[2], 1)

    def test_dependent_generators_and_even_ray_permutation(self):
        v = [[1,0],[0,1],[1,1],[1,-1],[1,2],[1,-2]]
        d = dataset(v)
        cell = d["cells"][0]
        self.assertEqual(cell["cone_dimension"], 3)
        self.assertEqual(cell["projective_dimension"], 2)
        c = check_witness(cell, [[-1,0],[0,1]])
        p = c["signed_permutation"]["permutation"]
        parity = (-1) ** sum(p[i] > p[j] for i in range(6) for j in range(i+1,6))
        self.assertEqual(parity, 1)  # Two swaps of the four sloped rays.
        self.assertEqual(c["orientation_determinant"], "-1")
        # In the independent standard basis E11,E12+E21,E22 the action is
        # diag(1,-1,1). This explains the expected negative determinant.

    def test_order_sign_and_integral_basis_invariance(self):
        v = [[1,0],[0,1],[1,1],[1,-1],[1,2],[1,-2]]
        g = [[-1,0],[0,1]]
        for order in (v, list(reversed(v)), v[2:] + v[:2]):
            signed = [[(-1 if i%2 else 1)*x for x in w] for i,w in enumerate(order)]
            self.assertEqual(span_action(signed,g)[2], -1)
            self.assertEqual(check_witness(dataset(signed)["cells"][0],g)["orientation_determinant"], "-1")
        u = [[1,2],[0,1]]
        self.assertEqual(det(u), 1)
        conjugate = mul(mul(u,g),inverse(u))
        self.assertEqual(conjugate, [[-1,4],[0,1]])
        new_g = [[int(x) for x in row] for row in conjugate]
        self.assertEqual(check_witness(dataset([mv(u,w) for w in v])["cells"][0],new_g)["orientation_determinant"], "-1")

    def test_rational_basis_with_nonintegral_coordinates(self):
        # Non-simplicial cone, coordinate equalities include rational entries.
        v = [[1,0],[0,1],[1,3],[1,1],[1,-2]]
        indices, action, sign = span_action(v,[[1,0],[1,-1]])
        self.assertEqual(len(indices),3)
        self.assertEqual(sign,-1)
        # Image of e1*e1^t is (1,1)*(1,1)^t, with coefficients
        # (2/3,-2,1/3) in the chosen basis.
        self.assertEqual([row[0] for row in action], [Fraction(2,3),-2,Fraction(1,3)])
        self.assertEqual(rational_string(Fraction(2,3)), "2/3")

    def test_invalid_witnesses(self):
        cell=dataset([[1,0],[0,1]])["cells"][0]
        for g in ([[1,0],[0,2]], [[0,0],[0,1]], [[0,1.0],[1,0]],
                  [[False,1],[1,0]], [[0,1]], [[1,1],[0,1]]):
            with self.subTest(g=g), self.assertRaises(Invalid):
                check_witness(cell,g)


class FormatTests(unittest.TestCase):
    def test_restricted_grammars(self):
        expected = [[[1,0],[0,1]]]
        self.assertEqual(parse_source("return [[[1,+0],[0,-1]]];", "gap"), [[[1,0],[0,-1]]])
        self.assertEqual(parse_source("cells=[[1,0;0,1],\\\n[1,1;0,1]];", "pari"), [[[1,0],[0,1]],[[1,1],[0,1]]])
        self.assertEqual(parse_source("U1=[[1,0],[0,1]]\nU2=[[1,1],[0,1]]\n", "pari-records"), expected+[[[1,1],[0,1]]])

    def test_malformed_and_executable_syntax(self):
        for text in ("return [[1/2]];", "return [[1.0]];", "return [foo()];", "return [[1,,2]];",
                     "return [[1,2],];", "return [[1,2]]; foo();", "return [[1,2]]", "return [[[1,0;0,1]]];",
                     "return __import__('os');", "return [[True]];", "return [1..8];"):
            with self.subTest(text=text), self.assertRaises(Invalid):
                parse_source(text,"gap")
        for text in ("U2=[[1]]", "U1=[[1]]\nU1=[[1]]", "U1=[[1]]\nsystem(1)"):
            with self.assertRaises(Invalid):
                parse_source(text,"pari-records")

    def test_configuration_rejections(self):
        for v in ([[0,0],[0,1]], [[2,0],[0,1]], [[1,0],[0,1,2]],
                  [[1.0,0],[0,1]], [[True,0],[0,1]], [[1,0],[1,0]],
                  [[1,0],[-1,0],[0,1]], [[1,0,0],[0,1,0]]):
            with self.subTest(v=v), self.assertRaises(Invalid):
                dataset(v)
        with self.assertRaisesRegex(Invalid,"duplicate records"):
            make_dataset([[[1,0],[0,1]], [[0,-1],[-1,0]]],"test",digest("source"))

    def test_rows_columns_conversion(self):
        records=[[[1,0,1],[0,1,1]]]
        d=make_dataset(records,"columns",digest(records),"columns")
        self.assertEqual(d["cells"][0]["ambient_rank"],2)
        self.assertEqual(d["cells"][0]["vector_pairs"],3)
        validate_dataset(d)

    def test_json_duplicate_keys(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/"bad.json";p.write_text('{"g":1,"g":2}')
            with self.assertRaises(Invalid):
                read_json(p)


class CertificateTests(unittest.TestCase):
    def setUp(self):
        self.d=dataset([[1,0],[0,1]])
        self.b=bundle(self.d,[[0,1],[1,0]])

    def test_certificate_tampering(self):
        mutations = [
            lambda c: c.update(dataset_sha256="0"*64),
            lambda c: c.update(source_record=2),
            lambda c: c.update(source_record=True),
            lambda c: c.update(basis_indices=[1,0]),
            lambda c: c.update(induced_action=[["1","0"],["0","1"]]),
            lambda c: c.update(orientation_determinant="1"),
            lambda c: c["signed_permutation"].update(permutation=[0,1]),
            lambda c: c["signed_permutation"].update(signs=[-1,1]),
            lambda c: c.update(g=[[1,0],[0,1]]),
        ]
        for mutation in mutations:
            b=deepcopy(self.b);mutation(b["certificates"][0])
            with self.subTest(cert=b), self.assertRaises(Invalid):
                verify_bundle(self.d,b)
        b=deepcopy(self.b);b["dataset_sha256"]="0"*64
        with self.assertRaises(Invalid):
            verify_bundle(self.d,b)

    def test_duplicate_and_missing_certificates(self):
        b=deepcopy(self.b);b["certificates"]*=2
        with self.assertRaises(Invalid):
            verify_bundle(self.d,b)
        b=deepcopy(self.b);b["certificates"]=[]
        result=verify(self.d,b,require_published=False)
        self.assertEqual(result["status"],"INCOMPLETE")
        self.assertFalse(result["coverage"]["complete_published_coverage"])

    def test_metadata_and_dimensions_recomputed(self):
        for field,value in (("cone_dimension",3),("projective_dimension",7),("vector_rank",1),
                            ("ambient_rank",3),("primitive",False),("source_record",True)):
            d=deepcopy(self.d);d["cells"][0][field]=value
            # Even a recomputed self hash cannot make false dimensions valid.
            d.pop("dataset_sha256");d["dataset_sha256"]=digest(d)
            with self.assertRaises(Invalid):
                validate_dataset(d)

    def test_discovery_and_forced_limits(self):
        cert,stats=find_witness(self.d,self.d["cells"][0],seconds=2)
        self.assertIsNotNone(cert)
        self.assertEqual(stats["status"],"FOUND")
        for kwargs in ({"seconds":0},{"max_nodes":0}):
            cert,stats=find_witness(self.d,self.d["cells"][0],**kwargs)
            self.assertIsNone(cert)
            self.assertEqual(stats["status"],"INCONCLUSIVE")

    def test_exhausted_search_is_inconclusive(self):
        # Rank one has a point as its projective cell, hence no reversing map.
        d=dataset([[1]])
        cert,stats=find_witness(d,d["cells"][0],seconds=2)
        self.assertIsNone(cert)
        self.assertEqual(stats["status"],"INCONCLUSIVE")

    def test_offline_checker_does_not_import_search_or_network(self):
        with patch.object(socket,"socket",side_effect=RuntimeError("network forbidden")), \
             patch("voronoi_verify.discover.find_witness", side_effect=RuntimeError("search forbidden")):
            result=verify(self.d,self.b,require_published=False)
            self.assertEqual(result["status"],"SUCCESS")
            self.assertEqual(result["result_level"],"CHECKED_ALL_SUPPLIED_RECORDS")
            self.assertTrue(result["mathematical_dependencies"]["unresolved_required"])

    def test_cli_exit_statuses(self):
        with tempfile.TemporaryDirectory() as td, contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            td=Path(td); dp=td/"data.json";cp=td/"cert.json";out=td/"result.json"
            write_json(dp,self.d);write_json(cp,self.b)
            base=["verify","--dataset",str(dp),"--certificates",str(cp),"--output",str(out),"--supplied-data"]
            self.assertEqual(main(base),0)
            write_json(cp,{**self.b,"dataset_sha256":"0"*64})
            self.assertEqual(main(base),2)
            cp.unlink()
            self.assertEqual(main(base),3)
            self.assertEqual(main(["full","--dataset",str(dp),"--certificates",str(cp),"--output",str(out),"--supplied-data","--seconds-per-cell","0"]),3)
            self.assertEqual(read_json(out)["status"],"INCOMPLETE")


class PublishedTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d=import_pinned()
        cls.b=read_json(ROOT/"results/rank8-certificates.json")

    def test_real_data_all_witnesses(self):
        result=verify(self.d,self.b)
        self.assertEqual(result["status"],"SUCCESS")
        self.assertEqual(result["coverage"]["checked"],106)
        self.assertTrue(result["coverage"]["complete_published_coverage"])
        self.assertFalse(result["mathematical_dependencies"]["unresolved_required"])
        for c in result["cells"]:
            self.assertEqual((c["ambient_rank"],c["vector_pairs"],c["cone_dimension"],c["projective_dimension"]),(8,9,9,8))

    def test_neighboring_source_dimensions(self):
        for source in ("DATA_gap/ListSHV_n8_rnk8","cells_GL8_dim8.gp"):
            d=import_pinned(source)
            self.assertEqual(len(d["cells"]),13)
            for c in d["cells"]:
                self.assertEqual((c["cone_dimension"],c["projective_dimension"]),(8,7))

    def test_deleted_record_even_with_consistent_new_hash(self):
        records=[c["vectors"] for c in self.d["cells"][:-1]]
        d=make_dataset(records,TARGET,self.d["source_sha256"])
        b=deepcopy(self.b);b["dataset_sha256"]=d["dataset_sha256"]
        b["certificates"]=b["certificates"][:-1]
        for c in b["certificates"]:
            c["dataset_sha256"]=d["dataset_sha256"]
        result=verify(d,b)
        self.assertEqual(result["status"],"INCOMPLETE")
        self.assertFalse(result["coverage"]["complete_published_coverage"])

    def test_duplicate_record_and_deleted_normalized_record(self):
        d=deepcopy(self.d);d["cells"].append(deepcopy(d["cells"][0]))
        with self.assertRaises(Invalid):
            validate_dataset(d)
        d=deepcopy(self.d);d["cells"].pop()
        with self.assertRaises(Invalid):
            validate_dataset(d)

    def test_missing_published_certificate_and_explicit_subset(self):
        b=deepcopy(self.b);b["certificates"].pop()
        result=verify(self.d,b)
        self.assertEqual(result["status"],"INCOMPLETE")
        self.assertEqual(len(result["coverage"]["missing_cell_ids"]),1)
        partial=verify(self.d,b,subset=True)
        self.assertEqual(partial["status"],"SUCCESS")
        self.assertEqual(partial["result_level"],"CHECKED_SUPPLIED_SUBSET")
        self.assertFalse(partial["coverage"]["complete_published_coverage"])
        self.assertTrue(partial["mathematical_dependencies"]["unresolved_required"])

    def test_optimized_mode_and_no_search_import(self):
        code="""
import sys
from voronoi_verify.__main__ import main
result=main(['verify','--output',sys.argv[1]])
if 'voronoi_verify.discover' in sys.modules:
    raise RuntimeError('checker imported discovery')
raise SystemExit(result)
"""
        with tempfile.TemporaryDirectory() as td:
            proc=subprocess.run([sys.executable,"-O","-c",code,str(Path(td)/"result.json")],cwd=ROOT,capture_output=True,text=True)
            self.assertEqual(proc.returncode,0,proc.stderr)


if __name__ == "__main__":
    unittest.main()

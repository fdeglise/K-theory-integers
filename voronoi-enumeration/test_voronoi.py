"""Small independent mathematical checks of the enumeration and its witnesses."""

from collections import Counter
from fractions import Fraction as Q
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from voronoi import (
    Limits, SearchLimit, bounded_vectors, canonicalize, determinant,
    columns, configuration_isometries, enumerate_cells, frames, lattice_grid,
    matvec, positive_definite, quadratic, sign_normalize, slice_vertices,
)
from itertools import combinations


class EnumerationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.results = {n: enumerate_cells(n) for n in (1, 2, 3)}

    def test_classical_orbit_counts(self):
        # For n=3 these are the spanning edge subsets of K4, up to integral
        # equivalence: a tree; triangle+tail or a 4-cycle; K4-edge; and K4.
        expected = {1: {1: 1}, 2: {2: 1, 3: 1}, 3: {3: 1, 4: 2, 5: 1, 6: 1}}
        for n, cells in self.results.items():
            self.assertEqual(Counter(c.cone_dimension for c in cells), expected[n])

    def test_every_output_has_complete_exact_witness(self):
        for n, cells in self.results.items():
            identity = tuple(tuple(int(i == j) for j in range(n)) for i in range(n))
            for c in cells:
                self.assertTrue(positive_definite(c.form))
                self.assertEqual(set(bounded_vectors(c.form, identity, 1)), set(c.rays))
                self.assertTrue(all(quadratic(c.form, v) == 1 for v in c.rays))

    def test_known_stabilizers_and_orientations(self):
        self.assertEqual([c.stabilizer_order for c in self.results[2]], [8, 12])
        self.assertTrue(all(c.orientation_reversing for c in self.results[2]))
        top = self.results[3][-1]
        self.assertEqual(top.stabilizer_order, 48)
        self.assertFalse(top.orientation_reversing)
        self.assertTrue(all(c.orientation_reversing for c in self.results[3][:-1]))

    def test_canonicalization_survives_large_integral_shear(self):
        shear = ((1, 7, -3), (0, 1, 5), (0, 0, 1))
        self.assertEqual(determinant(shear), 1)
        for c in self.results[3]:
            changed = tuple(sorted(sign_normalize(matvec(shear, v)) for v in c.rays))
            key, u = canonicalize(changed)
            self.assertEqual(key, c.rays)
            self.assertEqual(abs(determinant(u)), 1)
            self.assertTrue(all(x.denominator == 1 for row in u for x in row))

    def test_index_two_frames_are_enumerated_but_add_no_cells_in_rank_two(self):
        models = list(frames(2, 2))
        self.assertEqual(Counter(i for b, i in models), {1: 1, 2: 3})
        result = enumerate_cells(2, index_bound=2)
        self.assertEqual([c.rays for c in result], [c.rays for c in self.results[2]])

    def test_truncation_matches_complete_enumeration(self):
        result = enumerate_cells(3, max_cone_dimension=4)
        self.assertEqual([c.rays for c in result],
                         [c.rays for c in self.results[3] if c.cone_dimension <= 4])

    def test_ellipsoid_bound_finds_vectors_outside_unit_cube(self):
        # q(x,y)=(x-5y)^2+y^2, so (5,1) is a minimal vector.
        g = ((Q(1), Q(-5)), (Q(-5), Q(26)))
        vectors = bounded_vectors(g, ((1, 0), (0, 1)), 1)
        self.assertEqual(set(vectors), {(Q(1), Q(0)), (Q(5), Q(1))})

    def test_positive_definiteness_is_exact(self):
        self.assertFalse(positive_definite(((Q(1), Q(1)), (Q(1), Q(1)))))
        self.assertFalse(positive_definite(((Q(1), Q(2)), (Q(2), Q(1)))))
        tiny = Q(1, 10**50)
        self.assertTrue(positive_definite(((Q(1), Q(0)), (Q(0), tiny))))

    def test_lower_dimensional_and_impossible_slices(self):
        # Opposing supporting inequalities constrain the sole variable to zero.
        self.assertEqual(slice_vertices(2, ((Q(3, 5), Q(4, 5)),
                                            (Q(3, 5), Q(-4, 5))), Limits()),
                         ((Q(0),),))
        # e1/2 cannot be a nonzero lattice vector when e1 has minimal length.
        self.assertEqual(slice_vertices(2, ((Q(1, 2), Q(0)),), Limits()), ())
        # (e1+e2)/2 forces g12>=1, (e1-e2)/2 forces g12<=-1.
        self.assertEqual(slice_vertices(2, ((Q(1, 2), Q(1, 2)),
                                            (Q(1, 2), Q(-1, 2))), Limits()), ())

    def test_clipping_matches_independent_active_set_backend(self):
        for n in (1, 2, 3):
            for b, index in frames(n, 1):
                candidates = tuple(lattice_grid(b, index, [1]*n, Limits()))
                self.assertEqual(slice_vertices(n, candidates, Limits(), "clipping"),
                                 slice_vertices(n, candidates, Limits(), "active-sets"))
        # The cut removes all strictly interior points of the previous interval,
        # leaving a single existing vertex. The empty case follows afterwards.
        point = ((Q(3, 5), Q(4, 5)), (Q(3, 5), Q(-4, 5)))
        self.assertEqual(slice_vertices(2, point, Limits(), "clipping"), ((Q(0),),))
        empty = ((Q(1, 2), Q(1, 2)), (Q(1, 2), Q(-1, 2)))
        self.assertEqual(slice_vertices(2, empty, Limits(), "clipping"), ())

    def test_resource_limits_fail_explicitly(self):
        with self.assertRaises(SearchLimit):
            enumerate_cells(3, limits=Limits(grid_points=1))
        with self.assertRaises(SearchLimit):
            enumerate_cells(3, limits=Limits(vertex_systems=1))
        with self.assertRaises(SearchLimit):
            enumerate_cells(3, limits=Limits(faces=1))

    def test_unsupported_rank_is_rejected(self):
        with self.assertRaises(ValueError):
            enumerate_cells(5)

    def test_cli_json_and_failed_search_preserves_existing_file(self):
        script = Path(__file__).with_name("voronoi.py")
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "cells.json"
            success = subprocess.run([sys.executable, str(script), "2", "--json", str(output)],
                                     capture_output=True, text=True)
            self.assertEqual(success.returncode, 0, success.stderr)
            payload = json.loads(output.read_text())
            self.assertTrue(payload["complete"])
            self.assertEqual(len(payload["cells"]), 2)
            self.assertIn("1/2", str(payload["cells"]))
            previous = output.read_bytes()
            failure = subprocess.run([sys.executable, str(script), "3", "--max-faces", "1",
                                      "--json", str(output)], capture_output=True, text=True)
            self.assertEqual(failure.returncode, 2)
            self.assertIn("no complete result", failure.stderr)
            self.assertEqual(output.read_bytes(), previous)


class RankFourTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cells = enumerate_cells(4)
        cls.d4 = next(c for c in cls.cells if len(c.rays) == 12)

    def test_published_complete_counts(self):
        self.assertEqual(Counter(c.cone_dimension for c in self.cells),
                         {4: 1, 5: 3, 6: 4, 7: 4, 8: 2, 9: 2, 10: 2})
        self.assertEqual(len(self.cells), 18)

    def test_all_output_witnesses_are_complete(self):
        identity = tuple(tuple(int(i == j) for j in range(4)) for i in range(4))
        for c in self.cells:
            self.assertTrue(positive_definite(c.form))
            self.assertEqual(set(bounded_vectors(c.form, identity, 1)), set(c.rays))
            self.assertTrue(all(quadratic(c.form, v) == 1 for v in c.rays))

    def test_nonsimplicial_d4_and_a4_stabilizers(self):
        self.assertEqual(self.d4.cone_dimension, 10)
        self.assertEqual(self.d4.stabilizer_order, 1152)
        self.assertTrue(self.d4.orientation_reversing)
        a4 = next(c for c in self.cells if c.cone_dimension == 10 and len(c.rays) == 10)
        self.assertEqual(a4.stabilizer_order, 240)
        self.assertTrue(a4.orientation_reversing)

    def test_index_two_frame_forces_entire_d4_configuration(self):
        self.assertEqual(Counter(i for b, i in frames(4, 2)), {1: 1, 2: 15})
        b = ((1, 0, 0, 1), (0, 1, 0, 1), (0, 0, 1, 1), (0, 0, 0, 2))
        candidates = tuple(lattice_grid(b, 2, [2]*4, Limits()))
        # The Ryshkov slice is a point even though only four frame equalities
        # were specified. All 12 minimal pairs are forced simultaneously.
        self.assertEqual(slice_vertices(4, candidates, Limits()), ((Q(0),)*6,))
        g = tuple(tuple(Q(int(i == j)) for j in range(4)) for i in range(4))
        minimal = bounded_vectors(g, b, 2)
        self.assertEqual(len(minimal), 12)
        integral = tuple(sorted(sign_normalize(matvec(b, v)) for v in minimal))
        self.assertEqual(canonicalize(integral)[0], self.d4.rays)

    def test_d4_canonicalization_under_integral_shear(self):
        shear = ((1, 2, -1, 0), (0, 1, 3, 1), (0, 0, 1, -2), (0, 0, 0, 1))
        changed = tuple(sorted(sign_normalize(matvec(shear, v)) for v in self.d4.rays))
        key, u = canonicalize(changed)
        self.assertEqual(key, self.d4.rays)
        self.assertEqual(abs(determinant(u)), 1)
        self.assertTrue(all(x.denominator == 1 for row in u for x in row))

    def test_isometry_search_handles_nonunimodular_source_frame(self):
        frame = next(vs for vs in combinations(self.d4.rays, 4)
                     if abs(determinant(columns(vs))) == 2)
        reordered = frame + tuple(v for v in self.d4.rays if v not in frame)
        u = next(configuration_isometries(reordered, self.d4.rays))
        self.assertEqual(abs(determinant(u)), 1)
        self.assertTrue(all(x.denominator == 1 for row in u for x in row))
        self.assertEqual({sign_normalize(matvec(u, v)) for v in reordered}, set(self.d4.rays))

    def test_same_dimension_and_group_order_do_not_imply_equivalence(self):
        pair = [c for c in self.cells if c.cone_dimension == 5 and c.stabilizer_order == 96]
        self.assertEqual(len(pair), 2)
        self.assertIsNone(next(configuration_isometries(pair[0].rays, pair[1].rays), None))

    def test_insufficient_index_bound_and_clipping_limit_fail_explicitly(self):
        with self.assertRaises(ValueError):
            enumerate_cells(4, index_bound=1)
        with self.assertRaises(SearchLimit):
            enumerate_cells(4, limits=Limits(edge_tests=1))
        with self.assertRaises(SearchLimit):
            enumerate_cells(4, polyhedron_method="active-sets")

    def test_rank_four_truncation(self):
        result = enumerate_cells(4, max_cone_dimension=4)
        self.assertEqual([c.rays for c in result],
                         [c.rays for c in self.cells if c.cone_dimension == 4])


if __name__ == "__main__":
    unittest.main()

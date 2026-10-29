"""Integration checks against published counts and the independent Python model."""
import argparse
from collections import Counter
from fractions import Fraction as Q
import json
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from voronoi import (Limits, bounded_vectors, configuration_isometries,
                     positive_definite, quadratic, rank, outer_coordinates)


def run(binary, n, excess, path, *args):
    result = subprocess.run([str(binary), str(n), "--excess", str(excess),
                             "--json", str(path), *args], capture_output=True,
                            text=True, timeout=600)
    if result.returncode:
        raise AssertionError(f"rank {n}, excess {excess}: {result.stderr}")
    data = json.loads(path.read_text())
    assert data["complete"] and data["max_cone_dimension"] == n + excess
    for cell in data["cells"]:
        g = tuple(tuple(Q(x) for x in row) for row in cell["form"])
        rays = tuple(tuple(v) for v in cell["rays"])
        assert positive_definite(g)
        assert all(quadratic(g, v) == 1 for v in rays)
        assert rank([outer_coordinates(v) for v in rays]) == cell["cone_dimension"]
        if n <= 4:
            identity = tuple(tuple(int(i == j) for j in range(n)) for i in range(n))
            assert set(bounded_vectors(g, identity, 1, limits=Limits(grid_points=2_000_000))) == set(rays)
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("binary", type=Path)
    parser.add_argument("--slow", action="store_true")
    parser.add_argument("--deep", action="store_true", help="check rank 7 through cone dimension 10")
    args = parser.parse_args()
    cases = [(2, 1, {2: 1, 3: 1}), (3, 3, {3: 1, 4: 2, 5: 1, 6: 1}),
             (4, 6, {4: 1, 5: 3, 6: 4, 7: 4, 8: 2, 9: 2, 10: 2}),
             (5, 2, {5: 2, 6: 5, 7: 10}),
             (6, 1, {6: 3, 7: 10}), (7, 1, {7: 6, 8: 28})]
    if args.slow or args.deep:
        cases = cases[:3] + [(5, 3, {5: 2, 6: 5, 7: 10, 8: 16}),
                            (6, 2, {6: 3, 7: 10, 8: 28}),
                            (7, 3, {7: 6, 8: 28, 9: 115, 10: 467}) if args.deep else
                            (7, 2, {7: 6, 8: 28, 9: 115})]
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "cells.json"
        for n, excess, expected in cases:
            data = run(args.binary, n, excess, path, *( ["--stabilizers"] if n <= 4 else []))
            actual = Counter(c["cone_dimension"] for c in data["cells"])
            assert actual == expected, (n, excess, actual, expected)
            print(f"rank {n}, excess {excess}: {dict(actual)}; {data['lp_solves']} exact LP solves", flush=True)
            if n <= 4:
                reference = json.loads((ROOT / f"cells_n{n}.json").read_text())["cells"]
                for c in data["cells"]:
                    s = tuple(tuple(v) for v in c["rays"])
                    hits = [r for r in reference if r["cone_dimension"] == c["cone_dimension"]
                            and next(configuration_isometries(s, tuple(tuple(v) for v in r["rays"])), None) is not None]
                    assert len(hits) == 1
                    assert c["stabilizer_order"] == hits[0]["stabilizer_order"]
                    assert c["orientation_reversing"] == hits[0]["orientation_reversing"]
        previous = path.read_bytes()
        failed = subprocess.run([str(args.binary), "7", "--max-lp", "1", "--json", str(path)], capture_output=True, text=True)
        assert failed.returncode == 2 and path.read_bytes() == previous
        assert b"no complete result" in failed.stderr.encode()
        failed = subprocess.run([str(args.binary), "4", "--index-bound", "1"], capture_output=True, text=True)
        assert failed.returncode == 2
        for command in (["8"], ["4", "--excess", "-1"], ["7", "--excess", "2147483647"],
                        ["4294967303"], ["7", "--index-bound", "9"]):
            failed = subprocess.run([str(args.binary), *command], capture_output=True, text=True)
            assert failed.returncode == 2, command
    print("Counts, exact witnesses, integral equivalences, stabilizers and failure semantics passed.")


if __name__ == "__main__":
    main()

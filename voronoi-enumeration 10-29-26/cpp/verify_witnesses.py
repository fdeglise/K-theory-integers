"""Replay C++ form witnesses using the independent Python reference.

Default: check the largest coordinate box at each maximal frame index.
--all checks every witness. This verifies witnesses, not global completeness.
"""
import argparse
from fractions import Fraction as Q
from itertools import combinations
import json
from math import isqrt, prod
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from voronoi import (Limits, bounded_vectors, columns, determinant, inverse,
                     outer_coordinates, positive_definite, quadratic, rank)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("result", type=Path)
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--max-grid-points", type=int, default=2_000_000)
    args = parser.parse_args()
    data = json.loads(args.result.read_text())
    require(data["complete"] is True, "the result is marked incomplete")
    n = data["n"]
    selected = {}
    for i, cell in enumerate(data["cells"]):
        g = tuple(tuple(Q(x) for x in row) for row in cell["form"])
        require(len(g) == n and all(len(row) == n for row in g), "invalid matrix dimensions")
        require(positive_definite(g), f"cell {i}: form is not positive definite")
        gi = inverse(g)
        count = prod(2*isqrt(gi[j][j].numerator//gi[j][j].denominator)+1 for j in range(n))
        group = i if args.all else cell["maximal_frame_index"]
        if group not in selected or count > selected[group][0]:
            selected[group] = (count, i, g)
    identity = tuple(tuple(int(i == j) for j in range(n)) for i in range(n))
    for count, i, g in selected.values():
        cell = data["cells"][i]
        rays = tuple(tuple(v) for v in cell["rays"])
        require(all(len(v) == n for v in rays) and rank(rays) == n, f"cell {i}: invalid vector rank")
        require(all(quadratic(g, v) == 1 for v in rays), f"cell {i}: claimed minimum is not one")
        require(rank([outer_coordinates(v) for v in rays]) == cell["cone_dimension"], f"cell {i}: incorrect cone dimension")
        require(max(abs(determinant(columns(frame))) for frame in combinations(rays, n)) == cell["maximal_frame_index"], f"cell {i}: incorrect frame index")
        actual = bounded_vectors(g, identity, 1, limits=Limits(grid_points=args.max_grid_points))
        require(len(set(rays)) == len(rays) and set(actual) == set(rays), f"cell {i}: minimal-vector set is incomplete")
        print(f"Verified cell {i}: index {cell['maximal_frame_index']}, "
              f"cone dimension {cell['cone_dimension']}, {len(rays)} pairs, {count} grid points", flush=True)
    print(f"Verified {len(selected)} of {len(data['cells'])} form witnesses.")


if __name__ == "__main__":
    main()

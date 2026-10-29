"""Exact, dependency-free Voronoi cell enumeration for lattice ranks 1--4.

This deliberately slow reference implementation searches minimal-vector frames
and faces of their Ryshkov slices. It does not enumerate perfect forms first.
See README.md for the completeness argument and computational limits.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, permutations, product
import json
from math import comb, isqrt, prod
from pathlib import Path
import sys


class SearchLimit(RuntimeError):
    """A resource guard stopped the search; no completeness claim is possible."""


@dataclass(frozen=True)
class Limits:
    grid_points: int = 200_000
    vertex_systems: int = 200_000
    faces: int = 20_000
    edge_tests: int = 2_000_000


def transpose(a):
    return tuple(zip(*a))


def matvec(a, v):
    return tuple(sum(x * y for x, y in zip(row, v)) for row in a)


def matmul(a, b):
    return tuple(tuple(sum(x * y for x, y in zip(row, col))
                       for col in transpose(b)) for row in a)


def determinant(a):
    if all(Q(x).denominator == 1 for row in a for x in row):
        # Bareiss elimination is particularly useful for minimal-vector frames.
        a = [list(map(int, row)) for row in a]
        previous, sign = 1, 1
        for j in range(len(a) - 1):
            pivot = next((i for i in range(j, len(a)) if a[i][j]), None)
            if pivot is None:
                return Q(0)
            if pivot != j:
                a[j], a[pivot] = a[pivot], a[j]
                sign = -sign
            p = a[j][j]
            for i in range(j + 1, len(a)):
                for k in range(j + 1, len(a)):
                    a[i][k] = (p*a[i][k] - a[i][j]*a[j][k]) // previous
                a[i][j] = 0
            previous = p
        return Q(sign * a[-1][-1] if a else 1)
    a = [list(map(Q, row)) for row in a]
    result = Q(1)
    for j in range(len(a)):
        pivot = next((i for i in range(j, len(a)) if a[i][j]), None)
        if pivot is None:
            return Q(0)
        if pivot != j:
            a[j], a[pivot] = a[pivot], a[j]
            result = -result
        p = a[j][j]
        result *= p
        for i in range(j + 1, len(a)):
            t = a[i][j] / p
            for k in range(j + 1, len(a)):
                a[i][k] -= t * a[j][k]
    return result


def solve(a, b):
    """Solve a square rational system; return None if singular."""
    n = len(a)
    a = [list(map(Q, row)) + [Q(y)] for row, y in zip(a, b)]
    for j in range(n):
        pivot = next((i for i in range(j, n) if a[i][j]), None)
        if pivot is None:
            return None
        a[j], a[pivot] = a[pivot], a[j]
        p = a[j][j]
        a[j] = [x / p for x in a[j]]
        for i in range(n):
            if i != j:
                t = a[i][j]
                a[i] = [x - t * y for x, y in zip(a[i], a[j])]
    return tuple(row[-1] for row in a)


def inverse(a):
    n = len(a)
    cols = [solve(a, [int(i == j) for i in range(n)]) for j in range(n)]
    if any(col is None for col in cols):
        raise ValueError("singular matrix")
    return transpose(cols)


def rank(a):
    if not a:
        return 0
    a = [list(map(Q, row)) for row in a]
    r = 0
    for j in range(len(a[0])):
        p = next((i for i in range(r, len(a)) if a[i][j]), None)
        if p is None:
            continue
        a[r], a[p] = a[p], a[r]
        pivot = a[r][j]
        a[r] = [x / pivot for x in a[r]]
        for i in range(r + 1, len(a)):
            t = a[i][j]
            a[i] = [x - t * y for x, y in zip(a[i], a[r])]
        r += 1
        if r == len(a):
            break
    return r


def positive_definite(g):
    # Sylvester's criterion, with no floating-point tolerances.
    return all(determinant([row[:k] for row in g[:k]]) > 0
               for k in range(1, len(g) + 1))


def quadratic(g, v):
    return sum(x * y for x, y in zip(v, matvec(g, v)))


def sign_normalize(v):
    first = next((x for x in v if x), 0)
    return tuple(-x for x in v) if first < 0 else tuple(v)


def columns(vectors):
    return transpose(vectors)


def outer_coordinates(v):
    return tuple(v[i] * v[j] for i in range(len(v)) for j in range(i, len(v)))


def factor_tuples(n, total):
    if n == 1:
        yield (total,)
    else:
        for d in range(1, total + 1):
            if total % d == 0:
                for tail in factor_tuples(n - 1, total // d):
                    yield (d,) + tail


def frames(n, index_bound):
    """Row-HNF matrices B, representing L = B^{-1} Z^n.

    Frames transform by LEFT ambient unimodular multiplication. Thus an upper
    entry B[i,j] is reduced modulo B[j,j], not modulo B[i,i].
    """
    for index in range(1, index_bound + 1):
        for diag in factor_tuples(n, index):
            positions = [(i, j) for i in range(n) for j in range(i + 1, n)]
            for entries in product(*(range(diag[j]) for i, j in positions)):
                b = [[int(i == j) * diag[i] for j in range(n)] for i in range(n)]
                for (i, j), value in zip(positions, entries):
                    b[i][j] = value
                yield tuple(map(tuple, b)), index


def lattice_grid(b, index, bounds, limits):
    """Enumerate L in a box: x=a/index, with B*a divisible by index."""
    count = prod(2 * t + 1 for t in bounds)
    if count > limits.grid_points:
        raise SearchLimit(f"lattice grid needs {count} points (limit {limits.grid_points})")
    for a in product(*(range(-t, t + 1) for t in bounds)):
        if not any(a) or sign_normalize(a) != a:
            continue
        if all(y % index == 0 for y in matvec(b, a)):
            yield tuple(Q(x, index) for x in a)


def bounded_vectors(g, b, index, bound=Q(1), limits=Limits()):
    """All lattice vectors modulo sign with G[x] <= bound, certified exhaustive.

    Cauchy--Schwarz in the G metric gives x_i^2 <= bound*(G^-1)[i,i].
    """
    if not positive_definite(g):
        raise ValueError("short-vector enumeration requires positive definiteness")
    gi = inverse(g)
    radii = [isqrt((bound * index**2 * gi[i][i]).numerator //
                   (bound * index**2 * gi[i][i]).denominator) for i in range(len(g))]
    return tuple(v for v in lattice_grid(b, index, radii, limits)
                 if quadratic(g, v) <= bound)


def gram(n, y):
    g = [[Q(int(i == j)) for j in range(n)] for i in range(n)]
    for value, (i, j) in zip(y, combinations(range(n), 2)):
        g[i][j] = g[j][i] = value
    return tuple(map(tuple, g))


def slice_vertices(n, candidates, limits, method="auto"):
    """Enumerate vertices of the finite slice by all square active systems.

    Constraints are a*y >= rhs, with y the off-diagonal entries of G.
    The slice is bounded, including when it lies in a proper affine subspace.
    """
    if method not in ("auto", "active-sets", "clipping"):
        raise ValueError("unknown polyhedron method")
    inequalities = set()
    for v in candidates:
        a = tuple(2 * v[i] * v[j] for i, j in combinations(range(n), 2))
        rhs = 1 - sum(x*x for x in v)
        if not any(a):
            if rhs > 0:
                return ()
        else:
            inequalities.add((a, rhs))
    inequalities = sorted(inequalities)
    k = n * (n - 1) // 2
    if method == "clipping" or (method == "auto" and n == 4):
        return clipped_vertices(k, inequalities, limits)
    count = comb(len(inequalities), k) if len(inequalities) >= k else 0
    if count > limits.vertex_systems:
        raise SearchLimit(f"vertex search needs {count} systems (limit {limits.vertex_systems})")
    vertices = set()
    for active in combinations(inequalities, k):
        y = solve([a for a, rhs in active], [rhs for a, rhs in active])
        if y is not None and all(sum(x*t for x, t in zip(a, y)) >= rhs
                                 for a, rhs in inequalities):
            vertices.add(y)
    return tuple(sorted(vertices))


def clipped_vertices(k, inequalities, limits):
    """Exact incremental double description, starting from [-1/2,1/2]^k.

    The initial box is valid for every frame slice because e_i +/- e_j lie in
    L. Edges are certified by rank k-1 of the common active normals. This works
    for nonsimple vertices and slices that collapse to a proper affine space.
    """
    normals = []
    for j in range(k):
        normals.extend(tuple(s * int(i == j) for i in range(k)) for s in (1, -1))
    vertices = {}
    for signs in product((-1, 1), repeat=k):
        p = tuple(Q(s, 2) for s in signs)
        vertices[p] = frozenset(2*j + int(signs[j] > 0) for j in range(k))
    edge_tests = 0
    rank_cache = {}
    # Short normal vectors first tend to remove redundant large-vector cuts.
    for a, rhs in sorted(inequalities, key=lambda item: (sum(x*x for x in item[0]), item)):
        values = {p: sum(x*y for x, y in zip(a, p)) - rhs for p in vertices}
        positive = [p for p in vertices if values[p] > 0]
        zero = [p for p in vertices if values[p] == 0]
        negative = [p for p in vertices if values[p] < 0]
        if not negative and not zero:
            continue
        cut = len(normals)
        normals.append(a)
        kept = {p: vertices[p] | {cut} for p in zero}
        kept.update((p, vertices[p]) for p in positive)
        if negative:
            for p in positive:
                for q in negative:
                    edge_tests += 1
                    if edge_tests > limits.edge_tests:
                        raise SearchLimit(f"clipping exceeds {limits.edge_tests} edge tests")
                    common = vertices[p] & vertices[q]
                    if len(common) < k - 1:
                        continue
                    if common not in rank_cache:
                        rank_cache[common] = rank([normals[i] for i in sorted(common)])
                    if rank_cache[common] != k - 1:
                        continue
                    vp, vq = values[p], values[q]
                    point = tuple((vp*y - vq*x) / (vp-vq) for x, y in zip(p, q))
                    kept[point] = kept.get(point, frozenset()) | common | {cut}
        vertices = kept
        if not vertices:
            break
    return tuple(sorted(vertices))


def barycenter(vertices):
    return tuple(sum(p[j] for p in vertices) / len(vertices)
                 for j in range(len(vertices[0])))


def signed_unimodular_frames(vectors):
    n = len(vectors[0])
    for subset in combinations(vectors, n):
        if abs(determinant(columns(subset))) != 1:
            continue
        for chosen in permutations(subset):
            for signs in product((-1, 1), repeat=n):
                yield columns([tuple(s*x for x in v) for s, v in zip(signs, chosen)])


@lru_cache(maxsize=2048)
def intrinsic_metric(vectors):
    """Integral pairings det(C)*v^t*C^-1*w, C=sum(v*v^t).

    These and det(C) are invariant under ambient GL_n(Z). They do not depend on
    the arbitrary form used to witness a face or on signs of the ray vectors.
    """
    n = len(vectors[0])
    c = tuple(tuple(sum(v[i]*v[j] for v in vectors) for j in range(n)) for i in range(n))
    det_c = determinant(c)
    ci = inverse(c)
    dual = [matvec(ci, v) for v in vectors]
    pairings = tuple(tuple(int(det_c * sum(x*y for x, y in zip(v, w)))
                          for w in dual) for v in vectors)
    labels = tuple((pairings[i][i], tuple(sorted(abs(x) for x in row)))
                   for i, row in enumerate(pairings))
    return int(det_c), pairings, labels


def configuration_invariant(vectors):
    det_c, pairings, labels = intrinsic_metric(vectors)
    return len(vectors), det_c, tuple(sorted(labels))


def configuration_isometries(source, target):
    """All integral unimodular maps carrying source vector pairs onto target.

    A fixed independent source frame may have any index. Target images are
    extended recursively using intrinsic pairings; integrality and the complete
    configuration are checked at the end. Matching these metric invariants
    alone is not accepted as proof of integral equivalence.
    """
    if configuration_invariant(source) != configuration_invariant(target):
        return
    n = len(source[0])
    _, sp, sl = intrinsic_metric(source)
    _, tp, tl = intrinsic_metric(target)
    base_indices = next(indices for indices in combinations(range(len(source)), n)
                        if determinant(columns([source[i] for i in indices])))
    bi = inverse(columns([source[i] for i in base_indices]))
    target_set = set(target)

    def extend(chosen, signs):
        j = len(chosen)
        if j == n:
            frame = columns([tuple(s*x for x in target[i]) for i, s in zip(chosen, signs)])
            u = matmul(frame, bi)
            if any(Q(x).denominator != 1 for row in u for x in row):
                return
            if abs(determinant(u)) != 1:
                return
            if {sign_normalize(matvec(u, v)) for v in source} == target_set:
                yield u
            return
        src = base_indices[j]
        for i in range(len(target)):
            if i in chosen or tl[i] != sl[src]:
                continue
            for s in (-1, 1):
                if all(s*t*tp[i][previous] == sp[src][base_indices[k]]
                       for k, (previous, t) in enumerate(zip(chosen, signs))):
                    yield from extend(chosen + (i,), signs + (s,))

    yield from extend((), ())


def metric_canonicalize(vectors):
    """Canonicalize rank-4 configurations using metric-filtered basis frames.

    Every realizable full-rank configuration in rank 4 contains a unimodular
    minimal frame; see the index-two argument in README.md. Filter by the
    smallest intrinsic Gram tuple before doing rational changes of coordinates.
    """
    n = len(vectors[0])
    _, pairings, _ = intrinsic_metric(vectors)
    best_metric, finalists = None, []
    for subset in combinations(range(len(vectors)), n):
        if abs(determinant(columns([vectors[i] for i in subset]))) != 1:
            continue
        for indices in permutations(subset):
            # Changing all signs together has no effect on either key.
            for tail in product((-1, 1), repeat=n-1):
                signs = (1,) + tail
                metric = tuple(signs[i]*signs[j]*pairings[indices[i]][indices[j]]
                               for i in range(n) for j in range(i, n))
                if best_metric is None or metric < best_metric:
                    best_metric, finalists = metric, [(indices, signs)]
                elif metric == best_metric:
                    finalists.append((indices, signs))
    if not finalists:
        raise ValueError("no unimodular minimal-vector frame; unsupported input")
    best, best_u = None, None
    for indices, signs in finalists:
        frame = columns([tuple(s*x for x in vectors[i]) for i, s in zip(indices, signs)])
        u = inverse(frame)
        key = tuple(sorted(sign_normalize(matvec(u, v)) for v in vectors))
        if best is None or key < best:
            best, best_u = key, u
    return best, best_u


def canonicalize(vectors):
    """Canonical GL_n(Z) key and map U taking the input to that key.

    In n<=3 every independent minimal-vector frame is unimodular. In n=4 the
    complete minimal set still contains a unimodular frame, even when the chosen
    search frame has index two. Higher ranks are deliberately unsupported.
    """
    if len(vectors[0]) == 4:
        return metric_canonicalize(vectors)
    best = None
    best_u = None
    for frame in signed_unimodular_frames(vectors):
        u = inverse(frame)
        key = tuple(sorted(sign_normalize(matvec(u, v)) for v in vectors))
        if best is None or key < best:
            best, best_u = key, u
    if best is None:
        raise ValueError("no unimodular minimal-vector frame; unsupported input")
    return best, best_u


def stabilizer(vectors):
    """Return exact GL stabilizer order and existence of orientation reversal."""
    tensor_basis = []
    basis_vectors = []
    r = 0
    for v in vectors:
        row = outer_coordinates(v)
        if rank(tensor_basis + [row]) > r:
            tensor_basis.append(row)
            basis_vectors.append(v)
            r += 1
    # Select coordinates giving an invertible square minor of the tensor basis.
    pivot_cols = next(c for c in combinations(range(len(tensor_basis[0])), r)
                      if determinant([[row[j] for j in c] for row in tensor_basis]))
    minor = [[row[j] for j in pivot_cols] for row in tensor_basis]
    minor_det = determinant(minor)
    order = 0
    reversing = False
    for u in configuration_isometries(vectors, vectors):
        order += 1
        # Tensor coordinates use the same order as the original basis above.
        transformed = [outer_coordinates(matvec(u, v)) for v in basis_vectors]
        action_det = determinant([[row[j] for j in pivot_cols] for row in transformed]) / minor_det
        if action_det not in (Q(-1), Q(1)):
            raise ArithmeticError("invalid orientation determinant")
        reversing |= action_det < 0
    return order, reversing


@dataclass(frozen=True)
class Cell:
    rays: tuple
    form: tuple
    cone_dimension: int
    stabilizer_order: int
    orientation_reversing: bool

    @property
    def projective_dimension(self):
        return self.cone_dimension - 1

    def as_dict(self):
        return {
            "cone_dimension": self.cone_dimension,
            "projective_dimension": self.projective_dimension,
            "rays": [[int(x) for x in v] for v in self.rays],
            "form": [[str(x) for x in row] for row in self.form],
            "stabilizer_order": self.stabilizer_order,
            "orientation_reversing": self.orientation_reversing,
        }


def enumerate_cells(n, max_cone_dimension=None, index_bound=None, limits=Limits(),
                    *, polyhedron_method="auto", progress=None):
    """All full-vector-rank Voronoi cell orbits within the requested dimension.

    Resource-limit failures raise SearchLimit, never return a partial result.
    """
    if n not in (1, 2, 3, 4):
        raise ValueError("this prototype supports n=1,2,3,4 only")
    required_bound = 2 if n == 4 else 1
    index_bound = required_bound if index_bound is None else index_bound
    if index_bound < required_bound:
        raise ValueError(f"index_bound must be at least {required_bound} for n={n}")
    m = n * (n + 1) // 2
    cutoff = m if max_cone_dimension is None else max_cone_dimension
    if not n <= cutoff <= m:
        raise ValueError(f"max_cone_dimension must lie between {n} and {m}")
    representatives = {}
    buckets = {}
    raw_configurations = set()
    total_faces = 0
    for b, index in frames(n, index_bound):
        candidates = tuple(lattice_grid(b, index, [index_bound]*n, limits))
        vertices = slice_vertices(n, candidates, limits, polyhedron_method)
        if progress:
            progress(f"Frame index {index}, B={b}: {len(vertices)} slice vertices")
        if not vertices:
            continue
        witness = gram(n, barycenter(vertices))
        # Finite-pool lemma: if the true slice is nonempty, the finite slice
        # equals it. Thus failure of this witness proves the frame impossible.
        if not positive_definite(witness):
            continue
        short = bounded_vectors(witness, b, index, limits=limits)
        if any(quadratic(witness, v) < 1 for v in short):
            continue
        vertex_forms = [gram(n, y) for y in vertices]
        tight = [frozenset(i for i, g in enumerate(vertex_forms)
                           if quadratic(g, v) == 1) for v in candidates]
        walls = set(tight) - {frozenset()}
        root = frozenset(range(len(vertices)))
        pending = [root]
        seen = {root}
        while pending:
            face = pending.pop()
            total_faces += 1
            if total_faces > limits.faces:
                raise SearchLimit(f"face count exceeds limit {limits.faces}")
            y = barycenter([vertices[i] for i in sorted(face)])
            g = gram(n, y)
            minimal = tuple(v for v, active in zip(candidates, tight) if face <= active)
            r = rank([outer_coordinates(v) for v in minimal])
            if r > cutoff:
                continue
            affine_dimension = rank([tuple(a-t for a, t in zip(vertices[i], y)) for i in face])
            if r != m - affine_dimension:
                raise ArithmeticError("face dimension and minimal-vector rank disagree")
            # Verify each recorded face with an independent exhaustive ellipsoid
            # search, not merely with the finite candidate box.
            if set(bounded_vectors(g, b, index, limits=limits)) != set(minimal):
                raise ArithmeticError("incomplete minimal-vector set")
            integral = tuple(sorted(sign_normalize(matvec(b, v)) for v in minimal))
            if integral not in raw_configurations:
                raw_configurations.add(integral)
                invariant = configuration_invariant(integral)
                known = any(next(configuration_isometries(integral, other), None) is not None
                            for other in buckets.get(invariant, ()))
                if not known:
                    key, u = canonicalize(integral)
                    bi = inverse(b)
                    ambient_g = matmul(transpose(bi), matmul(g, bi))
                    ui = inverse(u)
                    canonical_g = matmul(transpose(ui), matmul(ambient_g, ui))
                    order, reversing = stabilizer(key)
                    representatives[key] = Cell(key, canonical_g, r, order, reversing)
                    buckets.setdefault(invariant, []).append(key)
                    if progress:
                        progress(f"  New orbit {len(representatives)}: cone dimension {r}, "
                                 f"{len(key)} rays, stabilizer order {order}")
            # Intersect with supporting walls, closing the active set implicitly
            # by keeping all vertices on that face. This handles rank jumps and
            # forced minima; one added vector need not mean one additional ray.
            if r < cutoff:
                for active in walls:
                    child = face & active
                    if child and child != face and child not in seen:
                        seen.add(child)
                        pending.append(child)
    if progress:
        progress(f"Visited {total_faces} frame faces, {len(raw_configurations)} "
                 f"distinct integral configurations, {len(representatives)} GL orbits")
    return sorted(representatives.values(), key=lambda c: (c.cone_dimension, c.rays))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("n", type=int, choices=(1, 2, 3, 4))
    parser.add_argument("--max-cone-dimension", type=int)
    parser.add_argument("--index-bound", type=int, help="defaults to 1 for n<=3, 2 for n=4")
    parser.add_argument("--polyhedron-method", choices=("auto", "active-sets", "clipping"), default="auto")
    parser.add_argument("--verbose", action="store_true", help="report frames and new orbits to stderr")
    parser.add_argument("--json", type=Path, help="write exact representatives and witnesses")
    parser.add_argument("--max-grid-points", type=int, default=200_000)
    parser.add_argument("--max-vertex-systems", type=int, default=200_000)
    parser.add_argument("--max-faces", type=int, default=20_000)
    parser.add_argument("--max-edge-tests", type=int, default=2_000_000)
    args = parser.parse_args(argv)
    limits = Limits(args.max_grid_points, args.max_vertex_systems, args.max_faces, args.max_edge_tests)
    if min(limits.grid_points, limits.vertex_systems, limits.faces, limits.edge_tests) < 1:
        parser.error("resource limits must be positive")
    try:
        cells = enumerate_cells(args.n, args.max_cone_dimension, args.index_bound, limits,
                                polyhedron_method=args.polyhedron_method,
                                progress=(lambda message: print(message, file=sys.stderr)) if args.verbose else None)
    except (SearchLimit, ValueError) as exc:
        print(f"Enumeration failed; no complete result: {exc}", file=sys.stderr)
        return 2
    print(f"Complete enumeration for n={args.n} in the requested dimension range.")
    print("cone dim | projective dim | cell orbits | orientation-preserving orbits")
    for d in sorted({c.cone_dimension for c in cells}):
        group = [c for c in cells if c.cone_dimension == d]
        print(f"{d:8} | {d-1:14} | {len(group):11} | {sum(not c.orientation_reversing for c in group):29}")
    if args.json:
        payload = {"n": args.n, "complete": True,
                   "index_bound": args.index_bound if args.index_bound is not None else (2 if args.n == 4 else 1),
                   "max_cone_dimension": args.max_cone_dimension or args.n*(args.n+1)//2,
                   "cells": [c.as_dict() for c in cells]}
        args.json.write_text(json.dumps(payload, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())

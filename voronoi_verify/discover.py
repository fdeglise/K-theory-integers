"""Deterministic, Gram-pruned signed frame search; no external group API.

C=sum vv^t is positive definite; every signed vector symmetry preserves
Q=C^{-1}. Match a full independent frame using its exact Q-Gram matrix.
Norms and sorted absolute Gram rows prune candidates before backtracking;
pairwise signed inner products prune each branch. Only integral unimodular
maps passing the independent checker can be returned. Exhaustion/time limits
give INCONCLUSIVE, even when no witness was found.
"""

from collections import Counter
import signal
import time

from .checker import certificate_for
from .exact import (F, Invalid, det, inverse, mul, mv, rank, transpose)


class LimitReached(Exception):
    pass


def find_witness(dataset, cell, seconds=30.0, max_nodes=200000):
    start = time.monotonic()
    nodes = 0
    old_handler = signal.getsignal(signal.SIGALRM)

    def expire(*_):
        raise LimitReached("wall time limit")

    def budget():
        if time.monotonic() - start >= seconds:
            raise LimitReached("wall time limit")
        if nodes >= max_nodes:
            raise LimitReached("search node limit")

    try:
        budget()
        signal.signal(signal.SIGALRM, expire)
        signal.setitimer(signal.ITIMER_REAL, seconds)
        vectors = cell["vectors"]
        n, s = cell["ambient_rank"], len(vectors)
        c = mul(transpose(vectors), vectors)
        q = inverse(c)
        gram = mul(mul(vectors, q), transpose(vectors))
        keys = [(gram[i][i], tuple(sorted(abs(x) for x in gram[i]))) for i in range(s)]
        counts = Counter(keys)
        order = sorted(range(s), key=lambda i: (counts[keys[i]], keys[i], i))
        frame = []
        for i in order:
            if rank([vectors[k] for k in frame + [i]]) > len(frame):
                frame.append(i)
            if len(frame) == n:
                break
        frame_inverse = inverse(transpose([vectors[i] for i in frame]))
        candidates = [[(j, sign) for j in range(s) if keys[j] == keys[i] for sign in (1, -1)] for i in frame]

        def search(images, used):
            nonlocal nodes
            budget()
            nodes += 1
            k = len(images)
            if k == n:
                target = transpose([[sign * x for x in vectors[j]] for j, sign in images])
                g = mul(target, frame_inverse)
                if any(x.denominator != 1 for row in g for x in row):
                    return None
                g = [[int(x) for x in row] for row in g]
                if det(g) not in (-1, 1):
                    return None
                try:
                    return certificate_for(dataset, cell, g)
                except Invalid:
                    return None
            for j, sign in candidates[k]:
                budget()
                if j in used:
                    continue
                if all(sign * e * gram[j][l] == gram[frame[k]][frame[t]]
                       for t, (l, e) in enumerate(images)):
                    found = search(images + [(j, sign)], used | {j})
                    if found is not None:
                        return found
            return None

        cert = search([], set())
        reason = None if cert else "pruned frame search exhausted without a reversing witness"
    except (LimitReached, MemoryError) as exc:
        cert, reason = None, str(exc) or "memory limit"
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, old_handler)
    return cert, {"cell_id": cell["id"], "status": "FOUND" if cert else "INCONCLUSIVE",
                  "reason": reason, "nodes": nodes, "elapsed_seconds": time.monotonic() - start}

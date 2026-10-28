"""Independent certificate checker. Imports no discovery or network code.

Certification uses explicit checks, never Python assert (also runs under -O).
"""

from .dataset import validate_dataset, canonical_json
from .exact import (require, integer, det, mv, outer, span_action,
                    rational_matrix, rational_string)


def check_witness(cell, g):
    n, vectors = cell["ambient_rank"], cell["vectors"]
    require(type(g) is list and len(g) == n and
            all(type(row) is list and len(row) == n and all(integer(x) for x in row) for row in g),
            "witness must be an N by N integer matrix")
    ambient_det = det(g)
    require(ambient_det in (-1, 1), "witness is not unimodular")
    lookup = {tuple(v): i for i, v in enumerate(vectors)}
    permutation, signs = [], []
    for v in vectors:
        w = tuple(mv(g, v))
        if w in lookup:
            j, sign = lookup[w], 1
        else:
            w = tuple(-x for x in w)
            require(w in lookup, "witness does not permute all signed vectors")
            j, sign = lookup[w], -1
        permutation.append(j)
        signs.append(sign)
    require(sorted(permutation) == list(range(len(vectors))), "signed map is not bijective")
    indices, action, orientation = span_action(vectors, g)
    require(orientation == -1, f"orientation determinant is {orientation}, expected -1")
    # The signed bijection fixes C=sum vv^t. C is in the relative interior
    # and is positive definite by full span. Its positive fixed ray contributes
    # +1, so det on W and on the projective tangent W/R*C have the same sign.
    c = [[sum(v[i] * v[j] for v in vectors) for j in range(n)] for i in range(n)]
    from .exact import mul, transpose
    require(mul(mul(g, c), transpose(g)) == c, "interior ray not fixed")
    return {"g": g, "signed_permutation": {"permutation": permutation, "signs": signs},
            "basis_indices": indices, "induced_action": rational_matrix(action),
            "orientation_determinant": rational_string(orientation),
            "ambient_determinant": rational_string(ambient_det)}


def certificate_for(dataset, cell, g):
    return {"cell_id": cell["id"], "source_record": cell["source_record"],
            "dataset_sha256": dataset["dataset_sha256"], **check_witness(cell, g)}


def verify_bundle(dataset, bundle):
    validate_dataset(dataset)
    require(type(bundle) is dict and bundle.get("schema") == "voronoi-certificates-v1", "unknown certificate schema")
    require(bundle.get("dataset_sha256") == dataset["dataset_sha256"], "mismatched dataset hash")
    require(set(bundle) == {"schema", "dataset_sha256", "certificates"}, "unexpected certificate bundle fields")
    by_id = {c["id"]: c for c in dataset["cells"]}
    certs = bundle["certificates"]
    require(type(certs) is list, "certificates must be a list")
    seen = set()
    checked = []
    for cert in certs:
        cid = cert["cell_id"]
        require(cid in by_id and cid not in seen, "unknown or duplicate certificate cell")
        seen.add(cid)
        cell = by_id[cid]
        expected = certificate_for(dataset, cell, cert["g"])
        require(canonical_json(cert) == canonical_json(expected), f"certificate metadata or claimed action incorrect: {cid}")
        checked.append({"cell_id": cid, "status": "COMPUTED", "check": "orientation_reversing",
                        "ambient_rank": cell["ambient_rank"], "vector_pairs": cell["vector_pairs"],
                        "vector_rank": cell["vector_rank"], "cone_dimension": cell["cone_dimension"],
                        "projective_dimension": cell["projective_dimension"], "orientation_determinant": "-1"})
    return checked, sorted(set(by_id) - seen)

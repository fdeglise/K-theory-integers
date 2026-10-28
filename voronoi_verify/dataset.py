"""Deterministic normalized configurations and an offline published-data anchor."""

import hashlib
import json
import math
import os
import tempfile
from pathlib import Path

from .exact import integer, rank, require, symmetric_coords, outer, transpose
from .formats import parse_source

ROOT = Path(__file__).resolve().parent.parent
REVISION = "d4c7b40ad2c95850d16160b347953dbaa9241486"
REPOSITORY = "https://github.com/elbazvip/Voronoi-complexes-database"
PINNED = {
    "DATA_gap/ListSHV_n8_rnk9": "d0b58f52da26e64b20649dffe4e317efd31a9d6b0c5e4d5b1407dbf5d636359d",
    "DATA_gap/ListSHV_n8_rnk8": "465c7d170f7a3ef2733327319c4155dce4f531dddd2d9a8b538b170f648a51e9",
    "cells_GL8_dim8.gp": "b1900f5eb598afd08ef40c5290955683c3563fd407aa452a8c7261cf7bebd8e9",
    "LICENSE": "3972dc9744f6499f0f9b2dbf76696f2ae7ad8af9b23dde66d6af86c9dfb36986",
    "README.md": "e007d4b81c2ef7e368d887684fdc05528f25996587994c35e48c3fae9949da0a",
}
TARGET = "DATA_gap/ListSHV_n8_rnk9"


def canonical_json(x):
    return json.dumps(x, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("ascii")


def digest(x):
    return hashlib.sha256(canonical_json(x)).hexdigest()


def write_json(path, obj):
    """Save atomically, retaining any different earlier generated evidence.

    Histories are named by their original byte hash, so repeated saves of the
    same earlier content create no extra copies. Source originals are handled
    separately by acquisition and are never overwritten.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = (json.dumps(obj, sort_keys=True, indent=2) + "\n").encode("utf-8")
    if path.exists():
        previous = path.read_bytes()
        if previous == raw:
            return
        h = hashlib.sha256(previous).hexdigest()
        history = path.parent / "history" / f"{path.stem}-{h}{path.suffix}"
        history.parent.mkdir(parents=True, exist_ok=True)
        if history.exists():
            require(history.read_bytes() == previous, "conflicting evidence history")
        else:
            history.write_bytes(previous)
    with tempfile.NamedTemporaryFile(dir=path.parent, prefix=f".{path.name}.", delete=False) as f:
        temporary = Path(f.name)
        f.write(raw)
    try:
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def read_json(path):
    def no_duplicates(pairs):
        d = {}
        for k, v in pairs:
            require(k not in d, f"duplicate JSON key: {k}")
            d[k] = v
        return d
    def bad_constant(s):
        raise ValueError(f"invalid JSON constant: {s}")
    return json.loads(Path(path).read_text(encoding="utf-8"),
                      object_pairs_hook=no_duplicates, parse_constant=bad_constant)


def canonical_vector(v):
    require(isinstance(v, (list, tuple)) and v and all(integer(x) for x in v), "nonintegral vector")
    require(any(v), "zero vector")
    require(math.gcd(*v) == 1, "nonprimitive vector")
    sign = 1 if next(x for x in v if x) > 0 else -1
    return tuple(sign * x for x in v)


def make_cell(vectors, source, record):
    require(type(vectors) is list and vectors, "empty configuration")
    require(all(type(v) is list for v in vectors), "vectors must be lists")
    n = len(vectors[0])
    require(1 <= n <= 32 and len(vectors) <= 512, "configuration exceeds v1 size bounds")
    require(all(len(v) == n for v in vectors), "inconsistent vector dimensions")
    v = sorted(canonical_vector(w) for w in vectors)
    require(len(set(v)) == len(v), f"duplicate rays in {source} record {record}")
    v = [list(w) for w in v]
    r = rank(transpose(v))
    require(r == n, "configuration does not span the ambient rational space")
    c = rank(transpose([symmetric_coords(outer(w)) for w in v]))
    return {"id": f"{source}#{record:04d}", "source": source, "source_record": record,
            "ambient_rank": n, "vectors": v, "vector_pairs": len(v), "vector_rank": r,
            "cone_dimension": c, "projective_dimension": c - 1,
            "primitive": True, "configuration_sha256": digest(v)}


def make_dataset(records, source, raw_sha256, vector_storage="rows"):
    require(type(source) is str and source and len(source) < 1024, "invalid source identifier")
    require(type(raw_sha256) is str and len(raw_sha256) == 64 and
            all(x in "0123456789abcdef" for x in raw_sha256), "invalid source SHA-256")
    require(vector_storage in ("rows", "columns"), "unknown vector storage")
    require(type(records) is list and records, "empty source list")
    cells = []
    seen = {}
    for i, matrix in enumerate(records, 1):
        require(type(matrix) is list and matrix and all(type(row) is list for row in matrix),
                f"invalid source record {i}")
        require(all(len(row) == len(matrix[0]) for row in matrix), "ragged source matrix")
        v = matrix if vector_storage == "rows" else transpose(matrix)
        cell = make_cell(v, source, i)
        h = cell["configuration_sha256"]
        require(h not in seen, f"duplicate records {seen.get(h)} and {i} in {source}")
        seen[h] = i
        cells.append(cell)
    payload = {"schema": "voronoi-dataset-v1", "convention": "integer column vectors; g*v",
               "source": source, "source_sha256": raw_sha256,
               "external_vector_storage": vector_storage, "cells": cells,
               "anomalies": []}
    return dict(payload, dataset_sha256=digest(payload))


def validate_dataset(d):
    require(type(d) is dict and d.get("schema") == "voronoi-dataset-v1", "unknown dataset schema")
    rebuilt = make_dataset([c["vectors"] for c in d["cells"]], d["source"], d["source_sha256"])
    rebuilt["external_vector_storage"] = d["external_vector_storage"]
    rebuilt.pop("dataset_sha256")
    rebuilt["dataset_sha256"] = digest(rebuilt)
    # Recompute IDs, all metadata, canonical ordering and hash, not just the digest.
    require(canonical_json(d) == canonical_json(rebuilt), "dataset metadata, dimensions, records or hash altered")
    return d


def original_path(source):
    return ROOT / "data" / "original" / source.replace("/", "__")


def import_pinned(source=TARGET):
    require(source in PINNED, "source is not pinned")
    raw = original_path(source).read_bytes()
    require(hashlib.sha256(raw).hexdigest() == PINNED[source], f"source bytes differ from pin: {source}")
    fmt = "gap" if source.startswith("DATA_gap/") else "pari"
    # GAP has s rows each of length 8. The standalone PARI matrices are
    # inspected in both conventions in the audit; columns are used for it.
    storage = "rows" if fmt == "gap" else "columns"
    return make_dataset(parse_source(raw.decode("ascii"), fmt), source, PINNED[source], storage)


def published_coverage(d):
    # Compare to a fresh parse of original bytes. A supplied count, an internally
    # consistent hash, or a self-declared dataset name is not a coverage anchor.
    return canonical_json(d) == canonical_json(import_pinned(TARGET))

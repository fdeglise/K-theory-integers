"""CLI. Exit 0=requested check passed, 2=invalid, 3=incomplete/inconclusive."""

import argparse
import hashlib
import json
import math
import os
import platform
import resource
import sys
import time
from pathlib import Path

from . import __version__
from .checker import verify_bundle
from .claims import manifest
from .dataset import (ROOT, REVISION, PINNED, TARGET, import_pinned, make_dataset,
                      original_path, published_coverage, read_json, validate_dataset, write_json)
from .exact import Invalid, require
from .formats import parse_source


def measurement(start):
    usage = resource.getrusage(resource.RUSAGE_SELF)
    return {"elapsed_seconds": time.monotonic() - start, "process_peak_rss_kib": usage.ru_maxrss,
            "process_user_seconds": usage.ru_utime, "process_system_seconds": usage.ru_stime,
            "python": sys.version, "platform": platform.platform(), "software_version": __version__,
            "cpu_count": os.cpu_count(), "jobs": 1}


def verify(dataset, bundle, require_published=True, subset=False):
    checked, missing = verify_bundle(dataset, bundle)
    is_published = published_coverage(dataset) if require_published else False
    complete = is_published and not missing
    unresolved = []
    if require_published and not is_published:
        unresolved.append("supplied list differs from the pinned published target")
    if missing and not subset:
        unresolved.append("missing requested certificates")
    if not checked:
        unresolved.append("no certificates checked")
    level = ("CHECKED_ALL_PINNED_PUBLISHED_RECORDS" if complete else
             "CHECKED_SUPPLIED_SUBSET" if missing or subset else "CHECKED_ALL_SUPPLIED_RECORDS")
    return {"schema": "voronoi-result-v1", "status": "INCOMPLETE" if unresolved else "SUCCESS",
            "result_level": level, "dataset_sha256": dataset["dataset_sha256"],
            "source_sha256": dataset["source_sha256"], "source": dataset["source"],
            "repository_revision": REVISION if is_published else None,
            "coverage": {"source_records": len(dataset["cells"]), "checked": len(checked),
                         "missing_cell_ids": missing, "matches_pinned_published_list": is_published,
                         "complete_published_coverage": complete},
            "cells": checked, "unresolved_checks": unresolved,
            "mathematical_dependencies": manifest(complete)}


def audit():
    out = []
    for source in ("DATA_gap/ListSHV_n8_rnk8", TARGET, "cells_GL8_dim8.gp"):
        d = import_pinned(source)
        counts = {}
        for cell in d["cells"]:
            key = str((cell["ambient_rank"],cell["vector_pairs"],cell["cone_dimension"],cell["projective_dimension"]))
            counts[key] = counts.get(key, 0) + 1
        out.append({"source": source, "source_sha256": d["source_sha256"],
                    "external_vector_storage": d["external_vector_storage"],
                    "count": len(d["cells"]), "dimensions_N_s_c_j": counts, "anomalies": d["anomalies"]})
    return {"schema": "voronoi-dimension-audit-v1", "sources": out,
            "notation_discrepancy": "DEKM v1 Table 4 labels cone dimensions as Sigma subscripts, although section 1.2 uses projective dimensions. The 13 records have j=7, the 106 records j=8.",
            "standalone_storage_qualification": "Columns chosen for the standalone square matrices as a lattice-basis convention. Either transpose has rank 8 and cone dimension 8, so its exclusion from projective degree 8 does not depend on that choice. No equivalence to GAP records is claimed."}


def acquire():
    # This is the only operation in the project which can initiate downloads.
    import urllib.request
    from datetime import datetime, timezone
    acquired = []
    for source, h in PINNED.items():
        path = original_path(source)
        url = f"https://raw.githubusercontent.com/elbazvip/Voronoi-complexes-database/{REVISION}/{source}"
        if path.exists():
            raw = path.read_bytes()
            require(hashlib.sha256(raw).hexdigest() == h, f"existing file differs from pin: {path}")
        else:
            with urllib.request.urlopen(url, timeout=30) as response:
                raw = response.read(8_000_001)
            require(len(raw) <= 8_000_000 and hashlib.sha256(raw).hexdigest() == h, f"download differs from pin: {source}")
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(raw)
        acquired.append({"source": source, "url": url, "sha256": h, "bytes": len(raw)})
    return {"files": acquired, "revision": REVISION, "checked_at": datetime.now(timezone.utc).isoformat()}


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="command", required=True)
    ac = sub.add_parser("acquire", help="explicit network acquisition; cached valid originals are reused")
    ac.add_argument("--output", default="results/acquisition.json")
    im = sub.add_parser("import", help="normalize pinned data or a supplied restricted source")
    im.add_argument("--source", type=Path)
    im.add_argument("--format", choices=["gap", "pari", "pari-records"], default="gap")
    im.add_argument("--vectors", choices=["rows", "columns"], default="rows")
    im.add_argument("--output", default="data/rank8-projective8.json")
    au = sub.add_parser("audit", help="recompute dimensions in neighboring pinned sources")
    au.add_argument("--output", default="results/dimension-audit.json")
    for name in ("discover", "verify", "full"):
        s = sub.add_parser(name)
        s.add_argument("--dataset", default="data/rank8-projective8.json")
        s.add_argument("--certificates", default="results/rank8-certificates.json")
        s.add_argument("--output", default=f"results/{name}.json")
        if name in ("verify", "full"):
            s.add_argument("--supplied-data", action="store_true", help="check supplied list only; no K-theory claim")
            s.add_argument("--subset", action="store_true", help="request checking just the provided certificates")
        if name in ("discover", "full"):
            s.add_argument("--seconds-per-cell", type=float, default=30)
            s.add_argument("--max-nodes", type=int, default=200000)
            s.add_argument("--memory-mib", type=int, default=512)
            s.add_argument("--limit", type=int, help="discover only the first k records; partial run")
    return p


def execute(args, start):
    command = args.command
    if command == "acquire":
        report = acquire()
    elif command == "import":
        if args.source is None:
            report = import_pinned()
        else:
            raw = args.source.read_bytes()
            report = make_dataset(parse_source(raw.decode("ascii"), args.format), args.source.name,
                                  hashlib.sha256(raw).hexdigest(), args.vectors)
    elif command == "audit":
        report = audit()
    else:
        if command == "full":
            # Full defaults to a fresh OFFLINE import and never downloads.
            if Path(args.dataset).exists():
                d = validate_dataset(read_json(args.dataset))
            else:
                d = import_pinned()
                write_json(args.dataset, d)
        else:
            d = validate_dataset(read_json(args.dataset))
        search = None
        if command in ("discover", "full"):
            require(math.isfinite(args.seconds_per_cell) and args.seconds_per_cell >= 0, "invalid time limit")
            require(args.max_nodes >= 0 and args.memory_mib >= 32, "invalid resource limit")
            require(args.limit is None or args.limit >= 0, "invalid record limit")
            limit = args.memory_mib * 1024 * 1024
            current_soft, current_hard = resource.getrlimit(resource.RLIMIT_AS)
            resource.setrlimit(resource.RLIMIT_AS, (min(limit, current_soft) if current_soft != -1 else limit, current_hard))
            # Lazy import guarantees the verify command never imports discovery.
            from .discover import find_witness
            certs, search = [], []
            search_start = time.monotonic()
            cells = d["cells"] if args.limit is None else d["cells"][:args.limit]
            for cell in cells:
                cert, stats = find_witness(d, cell, args.seconds_per_cell, args.max_nodes)
                search.append(stats)
                if cert:
                    certs.append(cert)
            bundle = {"schema": "voronoi-certificates-v1", "dataset_sha256": d["dataset_sha256"], "certificates": certs}
            write_json(args.certificates, bundle)
            search = {"cells": search, "elapsed_seconds": time.monotonic() - search_start,
                      "limits": {"seconds_per_cell": args.seconds_per_cell, "nodes_per_cell": args.max_nodes,
                                 "address_space_mib": args.memory_mib}}
        if command == "discover":
            complete = len(bundle["certificates"]) == len(d["cells"])
            report = {"status": "SUCCESS" if complete else "INCONCLUSIVE", "phase": "discovery",
                      "dataset_sha256": d["dataset_sha256"], "found": len(bundle["certificates"]),
                      "required": len(d["cells"]), "search": search}
        else:
            verify_start = time.monotonic()
            # Reload from saved files for the independent check, even for full.
            report = verify(d, read_json(args.certificates), not args.supplied_data, args.subset)
            report["verification_seconds"] = time.monotonic() - verify_start
            if search is not None:
                report["search"] = search
    if command != "import":
        report["measurements"] = measurement(start)
    write_json(args.output, report)
    if command in ("verify", "full"):
        cov = report["coverage"]
        summary = (f"{report['status']}: {cov['checked']}/{cov['source_records']} certificates checked; "
                   f"{report['result_level']}.\n")
        if cov["complete_published_coverage"]:
            summary += "H_8(Vor_GL_8; Z[1/7!])=0 and K_8(Z)=0 conditional on the listed external mathematical premises.\n"
        else:
            summary += "No complete published-list or K-theory conclusion from this run.\n"
        summary += f"Elapsed {report['measurements']['elapsed_seconds']:.3f} s; peak RSS {report['measurements']['process_peak_rss_kib']} KiB.\n"
        Path(args.output).with_suffix(".txt").write_text(summary)
        print(summary, end="")
    else:
        print(f"{command}: {report.get('status', 'SUCCESS')}; saved {args.output}")
    return 0 if report.get("status", "SUCCESS") == "SUCCESS" else 3


def main(argv=None):
    args = parser().parse_args(argv)
    start = time.monotonic()
    try:
        return execute(args, start)
    except (OSError, ImportError, MemoryError) as exc:
        report = {"status": "INCOMPLETE", "error": str(exc), "measurements": measurement(start)}
        code = 3
    except (Invalid, ValueError, TypeError, KeyError, IndexError, RecursionError, OverflowError) as exc:
        report = {"status": "INVALID", "error": str(exc), "measurements": measurement(start)}
        code = 2
    write_json(args.output, report)
    Path(args.output).with_suffix(".txt").write_text(f"{report['status']}: {report['error']}\n")
    print(f"{report['status']}: {report['error']}", file=sys.stderr)
    return code


if __name__ == "__main__":
    sys.exit(main())

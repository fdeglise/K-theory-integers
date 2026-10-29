"""Validation remains active under python -O; summaries preserve source evidence."""

import json
import os
from pathlib import Path


class VerificationError(ValueError):
    """An exact check of the attached calculation failed."""


def require(condition, message="verification check failed"):
    if not condition:
        raise VerificationError(str(message() if callable(message) else message))


def save_summary(original_path, text):
    # Verifiers never overwrite the archived run summaries in the source tree.
    destination = Path(os.environ.get(
        "VORONOI_OUTPUT_DIR", str(Path(__file__).resolve().parent.parent / "tmp" / "completeness")))
    destination.mkdir(parents=True, exist_ok=True)
    path = destination / Path(original_path).name
    # Validate JSON before saving; the integration driver reads this same file.
    json.loads(text)
    path.write_text(text + ("" if text.endswith("\n") else "\n"))

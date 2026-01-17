"""
IO logic
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable

from .models import CanonicalGenomeRecord


def write_ndjson(records: Iterable[CanonicalGenomeRecord], out_path: str | Path) -> None:
    """
    Write canonical records as NDJSON (newline-delimited JSON).
    Great for debugging, streaming, and easy diffs.
    """
    path = Path(out_path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", encoding="utf-8") as f:
        for rec in records:
            # dataclasses have __dict__ but we want stable JSON
            f.write(json.dumps(rec.__dict__, default=str) + "\n")

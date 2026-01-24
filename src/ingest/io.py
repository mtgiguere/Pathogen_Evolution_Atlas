"""
IO logic
"""

from __future__ import annotations

import json
from collections.abc import Iterable
from pathlib import Path

from .genbank import parse_collection_date
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
            
def load_ndjson(path: str | Path) -> list[CanonicalGenomeRecord]:
    """
    Load canonical records from an NDJSON file (one JSON object per line).
    """
    p = Path(path)
    records: list[CanonicalGenomeRecord] = []

    with p.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            obj = json.loads(line)

            # NDJSON stores dates as strings, so parse back to date objects.
            collection_date = parse_collection_date(obj.get("collection_date"))

            records.append(
                CanonicalGenomeRecord(
                    accession=str(obj.get("accession", "")).strip(),
                    organism=str(obj.get("organism", "")).strip(),
                    collection_date=collection_date,
                    country=obj.get("country"),
                    region=obj.get("region"),
                    host=obj.get("host"),
                    sequence_length=int(obj.get("sequence_length", 0)),
                    source=str(obj.get("source", "genbank")),
                )
            )

    return records

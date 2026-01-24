"""
models.py
"""
from __future__ import annotations

from dataclasses import dataclass
from datetime import date


@dataclass(frozen=True)
class CanonicalGenomeRecord:
    """
    Canonical, normalized representation of a single genome record.

    This is our project-wide contract: every upstream data source
    should be transformed into this shape before downstream processing.
    """
    accession: str
    organism: str
    collection_date: date | None
    country: str | None
    region: str | None
    host: str | None
    sequence_length: int
    source: str = "genbank"

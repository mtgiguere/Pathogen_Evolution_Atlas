from __future__ import annotations

from typing import Any

from .models import CanonicalGenomeRecord


def summarize_records(records: list[CanonicalGenomeRecord]) -> dict[str, Any]:
    """
    Compute a few sanity-check stats for a set of canonical records.
    Keeps scope intentionally small and beginner-friendly.
    """
    count = len(records)

    dates = [r.collection_date for r in records if r.collection_date is not None]
    min_date = min(dates) if dates else None
    max_date = max(dates) if dates else None

    missing_country = sum(1 for r in records if not r.country)
    pct_missing_country = (missing_country / count) if count else 0.0

    return {
        "count": count,
        "min_collection_date": min_date,
        "max_collection_date": max_date,
        "pct_missing_country": pct_missing_country,
    }

def test_summarize_records_empty():
    summary = summarize_records([])

    assert summary["count"] == 0
    assert summary["min_collection_date"] is None
    assert summary["max_collection_date"] is None
    assert summary["pct_missing_country"] == 0.0

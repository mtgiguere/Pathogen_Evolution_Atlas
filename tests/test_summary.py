from datetime import date

from src.ingest.models import CanonicalGenomeRecord
from src.ingest.summary import summarize_records


def test_summarize_records_basic():
    records = [
        CanonicalGenomeRecord(
            accession="A1",
            organism="Virus",
            collection_date=date(2020, 1, 1),
            country="USA",
            region=None,
            host=None,
            sequence_length=100,
        ),
        CanonicalGenomeRecord(
            accession="A2",
            organism="Virus",
            collection_date=date(2021, 6, 1),
            country=None,
            region=None,
            host=None,
            sequence_length=200,
        ),
    ]

    summary = summarize_records(records)

    assert summary["count"] == 2
    assert summary["min_collection_date"] == date(2020, 1, 1)
    assert summary["max_collection_date"] == date(2021, 6, 1)
    assert summary["pct_missing_country"] == 0.5

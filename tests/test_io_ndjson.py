from pathlib import Path

from src.ingest.io import load_ndjson


def test_load_ndjson_reads_records(tmp_path: Path):
    # Create a tiny NDJSON file (2 JSON objects, one per line)
    p = tmp_path / "sample.ndjson"
    p.write_text(
        '{"accession":"A1","organism":"X","collection_date":"2019-12-01","country":null,"region":null,"host":null,"sequence_length":4,"source":"genbank"}\n'
        '{"accession":"A2","organism":"Y","collection_date":null,"country":"USA","region":"Illinois","host":"Homo sapiens","sequence_length":0,"source":"genbank"}\n',
        encoding="utf-8",
    )

    records = load_ndjson(p)

    assert len(records) == 2
    assert records[0].accession == "A1"
    assert records[1].accession == "A2"

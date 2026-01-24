"""
Unit Tests
"""

from src.ingest.genbank import (
    fetch_and_normalize_many,
    fetch_many_genbank_minimal,
    normalize_many_genbank_minimal,
)


def test_fetch_many_returns_one_dict_per_accession(monkeypatch):
    """
    Batch fetch should return one minimal dict per accession,
    in the same order, using the single-record fetch under the hood.
    """
    from src.ingest import genbank

    # Fake single-record fetch so we don't hit the network in unit tests
    def fake_fetch_one(accession: str, email: str):
        return {
            "accession": accession,
            "organism": "Test organism",
            "collection_date": "2024-01",
            "location": "USA: Illinois",
            "host": "Homo sapiens",
            "sequence": "ACGT",
        }

    monkeypatch.setattr(genbank, "fetch_genbank_minimal", fake_fetch_one)

    accessions = ["A1", "A2", "A3"]
    out = fetch_many_genbank_minimal(accessions=accessions, email="test@example.com")

    assert [r["accession"] for r in out] == accessions
    assert len(out) == 3

def test_normalize_many_returns_canonical_records(monkeypatch):
    """
    Batch normalization should return CanonicalGenomeRecord objects,
    preserving order.
    """

    fake_raw = {
        "accession": "A1",
        "organism": "Test virus",
        "collection_date": "2024-01",
        "location": "USA: Illinois",
        "host": "Homo sapiens",
        "sequence": "ACGT",
    }

    # Pretend fetch_many already happened
    raws = [fake_raw, fake_raw | {"accession": "A2"}]

    out = normalize_many_genbank_minimal(raws)

    assert len(out) == 2
    assert out[0].accession == "A1"
    assert out[1].accession == "A2"

def test_fetch_and_normalize_many_chains_both(monkeypatch):
    """
    This should call fetch_many_genbank_minimal
    and then normalize_many_genbank_minimal.
    """
    from src.ingest import genbank

    fake_raw = [{"accession": "A1", "organism": "Test", "sequence": ""}]

    monkeypatch.setattr(genbank, "fetch_many_genbank_minimal", lambda accessions, email: fake_raw)
    monkeypatch.setattr(genbank, "normalize_many_genbank_minimal", lambda raws: ["NORMALIZED"])

    out = fetch_and_normalize_many(["A1"], email="x@example.com")

    assert out == ["NORMALIZED"]
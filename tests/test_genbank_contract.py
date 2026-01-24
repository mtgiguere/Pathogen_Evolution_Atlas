"""
Unit Tests models.py
"""
from datetime import date

import pytest

from src.ingest.genbank import normalize_genbank_minimal


def test_normalize_genbank_minimal_happy_path():
    raw = {
        "accession": "ABC123456",
        "organism": "SARS-CoV-2",
        "collection_date": "2024-08-19",
        "location": "USA: Illinois",
        "host": "Homo sapiens",
        "sequence": "ACGTACGT",
    }

    rec = normalize_genbank_minimal(raw)

    assert rec.accession == "ABC123456"
    assert rec.organism == "SARS-CoV-2"
    assert rec.collection_date == date(2024, 8, 19)
    assert rec.country == "USA"
    assert rec.region == "Illinois"
    assert rec.host == "Homo sapiens"
    assert rec.sequence_length == 8
    assert rec.source == "genbank"

def test_normalize_requires_accession_and_organism():
    """
    accession and organism are required for our canonical contract.
    If they are missing/blank, normalize should raise ValueError.
    """
    with pytest.raises(ValueError):
        normalize_genbank_minimal({"sequence": "ACGT"})

def test_normalize_requires_accession_and_organism():
    """
    accession and organism are required for our canonical contract.
    If they are missing/blank, normalize should raise ValueError.
    """
    with pytest.raises(ValueError):
        normalize_genbank_minimal({"sequence": "ACGT"})

def test_normalize_missing_optional_fields():
    """
    missing optional fields
    """
    raw = {
        "accession": "XYZ987",
        "organism": "Influenza A virus",
        # optional fields omitted on purpose
    }

    rec = normalize_genbank_minimal(raw)

    assert rec.collection_date is None
    assert rec.country is None
    assert rec.region is None
    assert rec.host is None
    assert rec.sequence_length == 0
    assert rec.source == "genbank"

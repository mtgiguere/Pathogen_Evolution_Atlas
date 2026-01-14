"""
Integration Tests
"""

import os
import pytest
import time
from src.ingest.genbank import fetch_genbank_minimal, normalize_genbank_minimal


@pytest.mark.integration
def test_fetch_and_normalize_one_accession():
    """
    Integration test: fetch one real GenBank record, normalize it, and assert basics.

    This is opt-in because it depends on the network + NCBI rate limits.
    """
    # NCBI requests a contact email. We'll require it via env var.
    # Example: setx NCBI_EMAIL "you@example.com"
    email = os.getenv("NCBI_EMAIL")
    if not email:
        pytest.skip("NCBI_EMAIL not set; skipping integration test")

    # A known SARS-CoV-2 reference genome accession in GenBank
    accession = "NC_045512.2"

    raw = fetch_genbank_minimal(accession=accession, email=email)
    rec = normalize_genbank_minimal(raw)

    assert rec.accession == accession
    assert rec.source == "genbank"
    assert rec.sequence_length > 1000  # should be ~30k

@pytest.mark.integration
def test_genbank_fetch_respects_rate_limit(monkeypatch):
    """
    We should never hit NCBI faster than ~3 requests per second.
    This test ensures fetch_genbank_minimal waits between calls.
    """
    calls = []

    def fake_efetch(*args, **kwargs):
        calls.append(time.time())
        return original_efetch(*args, **kwargs)

    # Patch Entrez.efetch so we can observe timing
    from src.ingest import genbank
    original_efetch = genbank.Entrez.efetch
    monkeypatch.setattr(genbank.Entrez, "efetch", fake_efetch)

    email = os.getenv("NCBI_EMAIL")
    if not email:
        pytest.skip("NCBI_EMAIL not set")

    accession = "NC_045512.2"

    # Call twice in a row
    genbank.fetch_genbank_minimal(accession, email)
    genbank.fetch_genbank_minimal(accession, email)

    # There should be at least ~0.34 seconds between calls
    assert calls[1] - calls[0] >= 0.34

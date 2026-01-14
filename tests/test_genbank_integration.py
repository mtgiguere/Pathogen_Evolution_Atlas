"""
Integration Tests
"""

import os
import pytest

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

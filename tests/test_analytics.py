from __future__ import annotations

import pandas as pd


def test_summarize_genomes_marks_scorable_and_skip_reasons(monkeypatch):
    import ingest.analytics as analytics

    # Stub scoring so this test focuses on analytics behavior.
    def fake_score_genome(rec):
        return {
            "accession": rec["accession"],
            "source": rec.get("source", "genbank"),
            "num_mutations": 2,
            "genes_affected": ["Spike", "N"],
            "risk_score": 3.0,
            "risk_level": "Moderate",
            "risk_explanation": "Moderate risk driven mostly by Spike.",
        }

    monkeypatch.setattr(analytics, "score_genome", fake_score_genome)

    # Use long toy sequences so we don't trip the "too_short" gate.
    ref_seq = "A" * 2000

    records = [
        # Reference present (lets analytics find a reference without fallback)
        {
            "accession": "NC_045512",
            "source": "genbank",
            "sequence": ref_seq,
            "collection_date": "2019-12-01",
        },
        # Scorable record (long enough + has sequence)
        {
            "accession": "SAMPLE_FULL",
            "source": "genbank",
            "sequence": "A" * 2000,
            "collection_date": "2020-01-01",
        },
        # Missing sequence
        {
            "accession": "SAMPLE_MISSING_SEQ",
            "source": "genbank",
            "sequence": "",
        },
        # Too short (should be blocked)
        {
            "accession": "SAMPLE_TOO_SHORT",
            "source": "genbank",
            "sequence": "A" * 10,
        },
    ]

    df = analytics.summarize_genomes(records)
    assert isinstance(df, pd.DataFrame)

    # Required columns (adjust only if you renamed in analytics.py)
    expected_cols = {
        "accession",
        "source",
        "sequence_length",
        "scorable",
        "skip_reason",
        "num_mutations",
        "genes_affected",
        "risk_score",
        "risk_level",
        "risk_explanation",
        "date",
        "lat",
        "lon",
    }
    assert expected_cols.issubset(set(df.columns))

    # Scorable record checks
    row_full = df.loc[df["accession"] == "SAMPLE_FULL"].iloc[0]
    assert bool(row_full["scorable"]) is True
    assert row_full["skip_reason"] == ""
    assert row_full["num_mutations"] == 2
    assert row_full["risk_score"] == 3.0
    assert isinstance(row_full["genes_affected"], str) and "Spike" in row_full["genes_affected"]

    # Missing sequence checks
    row_missing = df.loc[df["accession"] == "SAMPLE_MISSING_SEQ"].iloc[0]
    assert bool(row_missing["scorable"]) is False
    assert row_missing["skip_reason"] == "missing_sequence"
    assert row_missing["risk_level"] == "N/A"
    assert "Not scored:" in row_missing["risk_explanation"]

    # Too short checks
    row_short = df.loc[df["accession"] == "SAMPLE_TOO_SHORT"].iloc[0]
    assert bool(row_short["scorable"]) is False
    assert row_short["skip_reason"].startswith("too_short")
    assert row_short["risk_level"] == "N/A"


def test_summarize_genomes_reference_fallback_when_nc_missing(monkeypatch):
    import ingest.analytics as analytics

    def fake_score_genome(rec):
        return {
            "accession": rec["accession"],
            "source": rec.get("source", "genbank"),
            "num_mutations": 0,
            "genes_affected": [],
            "risk_score": 0.0,
            "risk_level": "Low",
            "risk_explanation": "No gene-attributed mutations detected.",
        }

    monkeypatch.setattr(analytics, "score_genome", fake_score_genome)

    # No NC_045512, but long sequences exist -> fallback reference should exist
    records = [
        {"accession": "PX1", "source": "genbank", "sequence": "A" * 2000},
        {"accession": "PX2", "source": "genbank", "sequence": "A" * 2100},
    ]

    df = analytics.summarize_genomes(records)

    # Should not mark missing_reference because fallback selects a reference sequence
    assert (df["skip_reason"] == "missing_reference").sum() == 0
    assert df["scorable"].all()


def test_summarize_genomes_handles_empty_input(monkeypatch):
    import ingest.analytics as analytics

    # Stub just in case, though it shouldn't be called.
    def fake_score_genome(_rec):
        raise AssertionError("score_genome should not be called for empty input")

    monkeypatch.setattr(analytics, "score_genome", fake_score_genome)

    df = analytics.summarize_genomes([])
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 0

    # If analytics returns empty DF with no columns, that's acceptable for now.
    # If you update analytics to always include columns, this will still pass.
    # (We don't assert df["scorable"] here to avoid KeyError.)

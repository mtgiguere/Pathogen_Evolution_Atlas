from __future__ import annotations

import logging

import pandas as pd
import pytest


def test_summarize_genomes_marks_scorable_and_skip_reasons(monkeypatch):
    import src.ingest.analytics as analytics

    def fake_score_genome(rec):
        return {
            "accession": rec["accession"],
            "source": rec.get("source", "genbank"),
            "scorable": True,
            "qc_status": "PASS",
            "qc_reasons": [],
            "num_mutations": 2,
            "genes_affected": ["Spike", "N"],
            "risk_score": 3.0,
            "risk_level": "Moderate",
            "risk_explanation": "Moderate risk driven mostly by Spike.",
            "risk_by_gene": {},
        }

    monkeypatch.setattr(analytics, "score_genome", fake_score_genome)

    reference_sequence = "A" * 2000
    reference_accession = "NC_045512.2"

    records = [
        {
            "accession": "SAMPLE_FULL",
            "source": "genbank",
            "sequence": "A" * 2000,
            "collection_date": "2020-01-01",
        },
        {
            "accession": "SAMPLE_MISSING_SEQ",
            "source": "genbank",
            "sequence": "",
        },
        {
            "accession": "SAMPLE_TOO_SHORT",
            "source": "genbank",
            "sequence": "A" * 10,
        },
    ]

    df = analytics.summarize_genomes(
        records,
        reference_sequence=reference_sequence,
        reference_accession=reference_accession,
    )

    assert isinstance(df, pd.DataFrame)

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

    row_full = df.loc[df["accession"] == "SAMPLE_FULL"].iloc[0]
    assert bool(row_full["scorable"]) is True
    assert row_full["skip_reason"] == ""
    assert row_full["num_mutations"] == 2
    assert row_full["risk_score"] == 3.0
    assert isinstance(row_full["genes_affected"], str) and "Spike" in row_full["genes_affected"]

    row_missing = df.loc[df["accession"] == "SAMPLE_MISSING_SEQ"].iloc[0]
    assert bool(row_missing["scorable"]) is False
    assert row_missing["skip_reason"] == "missing_sequence"
    assert row_missing["risk_level"] == "N/A"
    assert "Not scored:" in row_missing["risk_explanation"]

    row_short = df.loc[df["accession"] == "SAMPLE_TOO_SHORT"].iloc[0]
    assert bool(row_short["scorable"]) is False
    assert row_short["skip_reason"].startswith("too_short")
    assert row_short["risk_level"] == "N/A"


def test_summarize_genomes_requires_reference_sequence(monkeypatch):
    import src.ingest.analytics as analytics

    def fake_score_genome(rec):
        return {"accession": rec["accession"]}

    monkeypatch.setattr(analytics, "score_genome", fake_score_genome)

    records = [{"accession": "PX1", "source": "genbank", "sequence": "A" * 2000}]

    with pytest.raises(TypeError):
        analytics.summarize_genomes(records)


def test_summarize_genomes_handles_empty_input(monkeypatch):
    import src.ingest.analytics as analytics

    def fake_score_genome(_rec):
        raise AssertionError("score_genome should not be called for empty input")

    monkeypatch.setattr(analytics, "score_genome", fake_score_genome)

    df = analytics.summarize_genomes(
        [],
        reference_sequence="A" * 2000,
        reference_accession="NC_045512.2",
    )

    assert isinstance(df, pd.DataFrame)
    assert len(df) == 0


def test_summarize_genomes_logs_pre_scoring_skip(caplog):
    """Records that fail the pre-scoring gate must be logged at DEBUG."""
    import src.ingest.analytics as analytics

    with caplog.at_level(logging.DEBUG):
        analytics.summarize_genomes(
            [{"accession": "SHORT_REC", "sequence": "A" * 10}],
            reference_sequence="A" * 2000,
        )

    assert any("SHORT_REC" in r.message for r in caplog.records)


def test_summarize_genomes_logs_batch_summary(monkeypatch, caplog):
    """summarize_genomes must emit an INFO-level batch summary at the end."""
    import src.ingest.analytics as analytics

    def fake_score(rec):
        return {
            "accession": rec["accession"],
            "source": "genbank",
            "scorable": True,
            "qc_status": "PASS",
            "qc_reasons": [],
            "num_mutations": 0,
            "genes_affected": [],
            "risk_score": 0.0,
            "risk_level": "Low",
            "risk_explanation": "ok",
            "risk_by_gene": {},
        }

    monkeypatch.setattr(analytics, "score_genome", fake_score)

    with caplog.at_level(logging.INFO):
        analytics.summarize_genomes(
            [{"accession": "A1", "sequence": "A" * 2000}],
            reference_sequence="A" * 2000,
        )

    assert any(r.levelno == logging.INFO for r in caplog.records)

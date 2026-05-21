"""
REST API tests — TDD, written before implementation.
All tests expected to FAIL until src/ingest/api.py is implemented.

Uses FastAPI's TestClient (synchronous) with a dependency override to inject
a pre-built DataFrame so tests never touch disk or NCBI.
"""

from __future__ import annotations

import pandas as pd
import pytest
from fastapi.testclient import TestClient

from src.ingest.api import app, get_dataframe

# ── fixtures ──────────────────────────────────────────────────────────────────

_SAMPLE_DF = pd.DataFrame([
    {
        "accession": "MN908947",
        "lineage": "B.1.617.2",
        "lineage_display": "Delta",
        "who_label": "Delta",
        "who_class": "VOC",
        "lineage_confidence": 0.9,
        "country": "USA",
        "region": "California",
        "collection_date": "2021-10-01",
        "risk_level": "High",
        "risk_score": 7.0,
        "num_mutations": 12,
        "genes_affected": "S, ORF1ab",
        "escape_count": 2,
        "escape_antibodies": "Bamlanivimab, REGN10933",
        "has_critical_escape": True,
        "scorable": True,
        "sequence_length": 29903,
    },
    {
        "accession": "OX123456",
        "lineage": "BA.2",
        "lineage_display": "Omicron BA.2",
        "who_label": "Omicron",
        "who_class": "VOC",
        "lineage_confidence": 0.8,
        "country": "GBR",
        "region": "England",
        "collection_date": "2022-03-15",
        "risk_level": "Moderate",
        "risk_score": 4.0,
        "num_mutations": 8,
        "genes_affected": "S",
        "escape_count": 1,
        "escape_antibodies": "Bamlanivimab",
        "has_critical_escape": True,
        "scorable": True,
        "sequence_length": 29890,
    },
    {
        "accession": "PP000001",
        "lineage": "Unknown",
        "lineage_display": "Unknown",
        "who_label": "",
        "who_class": "",
        "lineage_confidence": 0.0,
        "country": "FRA",
        "region": None,
        "collection_date": None,
        "risk_level": "Low",
        "risk_score": 1.0,
        "num_mutations": 2,
        "genes_affected": "N",
        "escape_count": 0,
        "escape_antibodies": "",
        "has_critical_escape": False,
        "scorable": True,
        "sequence_length": 29500,
    },
])


@pytest.fixture()
def client():
    """TestClient with the sample DataFrame injected via DI override."""
    app.dependency_overrides[get_dataframe] = lambda: _SAMPLE_DF
    yield TestClient(app)
    app.dependency_overrides.clear()


# ── /health ───────────────────────────────────────────────────────────────────


def test_health_returns_200(client):
    r = client.get("/health")
    assert r.status_code == 200
    assert r.json()["status"] == "ok"


# ── GET /genomes ──────────────────────────────────────────────────────────────


def test_get_genomes_returns_list(client):
    r = client.get("/genomes")
    assert r.status_code == 200
    data = r.json()
    assert isinstance(data, list)
    assert len(data) == 3


def test_get_genomes_contains_expected_fields(client):
    r = client.get("/genomes")
    row = r.json()[0]
    for field in ("accession", "lineage", "country", "risk_level", "escape_count"):
        assert field in row, f"missing field: {field}"


def test_get_genomes_filter_by_country(client):
    r = client.get("/genomes?country=USA")
    assert r.status_code == 200
    data = r.json()
    assert len(data) == 1
    assert data[0]["accession"] == "MN908947"


def test_get_genomes_filter_by_lineage(client):
    r = client.get("/genomes?lineage=BA.2")
    assert r.status_code == 200
    data = r.json()
    assert len(data) == 1
    assert data[0]["accession"] == "OX123456"


def test_get_genomes_filter_country_case_insensitive(client):
    r = client.get("/genomes?country=usa")
    assert r.status_code == 200
    assert len(r.json()) == 1


def test_get_genomes_pagination_limit(client):
    r = client.get("/genomes?limit=2")
    assert r.status_code == 200
    assert len(r.json()) == 2


def test_get_genomes_pagination_offset(client):
    r = client.get("/genomes?limit=10&offset=2")
    assert r.status_code == 200
    assert len(r.json()) == 1


def test_get_genomes_empty_result(client):
    r = client.get("/genomes?country=ZZZZZ")
    assert r.status_code == 200
    assert r.json() == []


# ── GET /genomes/{accession} ──────────────────────────────────────────────────


def test_get_genome_found(client):
    r = client.get("/genomes/MN908947")
    assert r.status_code == 200
    data = r.json()
    assert data["accession"] == "MN908947"
    assert data["lineage"] == "B.1.617.2"


def test_get_genome_not_found(client):
    r = client.get("/genomes/NOTEXIST")
    assert r.status_code == 404


def test_get_genome_detail_has_extra_fields(client):
    r = client.get("/genomes/MN908947")
    data = r.json()
    for field in ("risk_score", "num_mutations", "escape_antibodies", "who_class"):
        assert field in data, f"missing detail field: {field}"


# ── GET /variants ─────────────────────────────────────────────────────────────


def test_get_variants_returns_list(client):
    r = client.get("/variants")
    assert r.status_code == 200
    data = r.json()
    assert isinstance(data, list)
    assert len(data) >= 1


def test_get_variants_has_count_field(client):
    r = client.get("/variants")
    for v in r.json():
        assert "lineage" in v
        assert "count" in v


def test_get_variants_counts_are_correct(client):
    r = client.get("/variants")
    counts = {v["lineage"]: v["count"] for v in r.json()}
    assert counts.get("B.1.617.2") == 1
    assert counts.get("BA.2") == 1


# ── GET /summary ──────────────────────────────────────────────────────────────


def test_get_summary_returns_totals(client):
    r = client.get("/summary")
    assert r.status_code == 200
    data = r.json()
    assert data["total_genomes"] == 3
    assert data["scorable_count"] == 3


def test_get_summary_risk_breakdown(client):
    r = client.get("/summary")
    rb = r.json()["risk_breakdown"]
    assert isinstance(rb, dict)
    assert rb.get("High") == 1
    assert rb.get("Moderate") == 1
    assert rb.get("Low") == 1


def test_get_summary_top_lineages(client):
    r = client.get("/summary")
    top = r.json()["top_lineages"]
    assert isinstance(top, list)
    lineage_names = [t["lineage"] for t in top]
    assert "B.1.617.2" in lineage_names


def test_get_summary_escape_stats(client):
    r = client.get("/summary")
    data = r.json()
    assert "escape_genome_count" in data
    # 2 out of 3 genomes have escape_count > 0
    assert data["escape_genome_count"] == 2

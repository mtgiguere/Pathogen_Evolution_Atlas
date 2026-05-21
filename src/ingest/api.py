"""
api.py — FastAPI REST interface for the Pathogen Evolution Atlas.

Start with:
    uvicorn src.ingest.api:app --reload

Endpoints
---------
GET /health                        — liveness check
GET /genomes                       — list genomes (filterable, paginated)
GET /genomes/{accession}           — single genome detail
GET /variants                      — lineage counts
GET /summary                       — aggregate dashboard summary
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Annotated, Any

import pandas as pd
from fastapi import Depends, FastAPI, HTTPException, Query

app = FastAPI(
    title="Pathogen Evolution Atlas",
    description="Genomic surveillance REST API for SARS-CoV-2.",
    version="1.0.0",
)

_GENOMES_PATH = Path(os.getenv("GENOMES_PATH", "data/raw/genomes.ndjson"))
_REF_PATH = Path(os.getenv("REF_PATH", "data/reference/genbank_reference.ndjson"))

# ── Data dependency ───────────────────────────────────────────────────────────
# Override `get_dataframe` in tests to inject a DataFrame without touching disk.


def get_dataframe() -> pd.DataFrame:
    """Load and score genomes from disk. Override via app.dependency_overrides in tests."""
    from .analytics import summarize_genomes
    from .io import read_ndjson

    if not _GENOMES_PATH.exists() or not _REF_PATH.exists():
        return pd.DataFrame()

    genomes = list(read_ndjson(_GENOMES_PATH))
    ref_records = list(read_ndjson(_REF_PATH))
    if not ref_records:
        return pd.DataFrame()

    ref_seq: str = ref_records[0].get("sequence", "")
    return summarize_genomes(genomes, reference_sequence=ref_seq)


DataFrameDep = Annotated[pd.DataFrame, Depends(get_dataframe)]


# ── helpers ───────────────────────────────────────────────────────────────────


def _row_to_summary(row: pd.Series) -> dict[str, Any]:
    return {
        "accession": row.get("accession"),
        "lineage": row.get("lineage", "Unknown"),
        "lineage_display": row.get("lineage_display", "Unknown"),
        "who_label": row.get("who_label", ""),
        "who_class": row.get("who_class", ""),
        "lineage_confidence": float(row.get("lineage_confidence") or 0.0),
        "country": row.get("country"),
        "region": row.get("region"),
        "collection_date": str(row.get("collection_date") or ""),
        "risk_level": row.get("risk_level", "N/A"),
        "risk_score": float(row.get("risk_score") or 0.0),
        "num_mutations": int(row.get("num_mutations") or 0),
        "genes_affected": row.get("genes_affected", ""),
        "escape_count": int(row.get("escape_count") or 0),
        "escape_antibodies": row.get("escape_antibodies", ""),
        "has_critical_escape": bool(row.get("has_critical_escape", False)),
        "scorable": bool(row.get("scorable", False)),
        "sequence_length": int(row.get("sequence_length") or 0),
    }


def _row_to_detail(row: pd.Series) -> dict[str, Any]:
    base = _row_to_summary(row)
    base["escape_mechanisms"] = row.get("escape_mechanisms", "")
    return base


# ── routes ────────────────────────────────────────────────────────────────────


@app.get("/health")
def health() -> dict[str, str]:
    return {"status": "ok"}


@app.get("/genomes")
def list_genomes(
    df: DataFrameDep,
    country: str | None = Query(default=None),
    lineage: str | None = Query(default=None),
    limit: int = Query(default=100, ge=1, le=10_000),
    offset: int = Query(default=0, ge=0),
) -> list[dict[str, Any]]:
    if df.empty:
        return []

    mask = pd.Series([True] * len(df), index=df.index)

    if country is not None:
        country_col = df["country"].fillna("").str.upper()
        mask &= country_col == country.upper()

    if lineage is not None:
        mask &= df["lineage"] == lineage

    filtered = df[mask]
    page = filtered.iloc[offset : offset + limit]

    return [_row_to_summary(row) for _, row in page.iterrows()]


@app.get("/genomes/{accession}")
def get_genome(accession: str, df: DataFrameDep) -> dict[str, Any]:
    if df.empty:
        raise HTTPException(status_code=404, detail="Genome not found")

    rows = df[df["accession"] == accession]
    if rows.empty:
        raise HTTPException(status_code=404, detail=f"Genome '{accession}' not found")

    return _row_to_detail(rows.iloc[0])


@app.get("/variants")
def list_variants(df: DataFrameDep) -> list[dict[str, Any]]:
    if df.empty:
        return []

    counts = (
        df.groupby("lineage")
        .size()
        .reset_index(name="count")
        .sort_values("count", ascending=False)
    )

    result: list[dict[str, Any]] = []
    for _, row in counts.iterrows():
        lineage = row["lineage"]
        # Pull display name and WHO info from first matching genome row
        meta = df[df["lineage"] == lineage].iloc[0]
        result.append({
            "lineage": lineage,
            "lineage_display": meta.get("lineage_display", lineage),
            "who_label": meta.get("who_label", ""),
            "who_class": meta.get("who_class", ""),
            "count": int(row["count"]),
        })

    return result


@app.get("/summary")
def summary(df: DataFrameDep) -> dict[str, Any]:
    if df.empty:
        return {
            "total_genomes": 0,
            "scorable_count": 0,
            "risk_breakdown": {},
            "top_lineages": [],
            "escape_genome_count": 0,
        }

    total = len(df)
    scorable = int(df["scorable"].sum()) if "scorable" in df.columns else 0

    risk_breakdown: dict[str, int] = {}
    if "risk_level" in df.columns:
        for level, cnt in df["risk_level"].value_counts().items():
            risk_breakdown[str(level)] = int(cnt)

    top_lineages: list[dict[str, Any]] = []
    if "lineage" in df.columns:
        top = df["lineage"].value_counts().head(10)
        for lin, cnt in top.items():
            top_lineages.append({"lineage": str(lin), "count": int(cnt)})

    escape_count = 0
    if "escape_count" in df.columns:
        escape_count = int((df["escape_count"] > 0).sum())

    return {
        "total_genomes": total,
        "scorable_count": scorable,
        "risk_breakdown": risk_breakdown,
        "top_lineages": top_lineages,
        "escape_genome_count": escape_count,
    }

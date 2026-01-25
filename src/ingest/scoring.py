"""
scoring.py
End-to-end scoring orchestration (dependency-injection friendly).
"""

from __future__ import annotations
from typing import Any
from collections.abc import Callable
from ingest.mutations import Mutation, diff_sequences
from ingest.risk import score_mutations

def _identify_mutations(record: Any):
    """
    Identify mutations for either:
      - dict-like records (record["sequence"], record["reference_sequence"])
      - CanonicalGenomeRecord dataclasses (record.sequence, record.reference_sequence)

    If sequence or reference is missing, return [].
    """
    # Dict-like path
def _identify_mutations(record: Any):
    # Dict-like path
    if isinstance(record, dict):
        ref = record.get("reference_sequence")
        sample = record.get("sequence")
        if not ref or not sample:
            return []
        try:
            return diff_sequences(ref, sample)
        except ValueError:
            # Length mismatch (partial/trimmed sequences) -> treat as unscorable for v1
            return []

    # Dataclass / object path
    ref = getattr(record, "reference_sequence", None)
    sample = getattr(record, "sequence", None)
    if not ref or not sample:
        return []

    try:
        return diff_sequences(ref, sample)
    except ValueError:
        return []


def _map_genes(mutations: list[Mutation]) -> list[str]:
    """
    Summarize affected genes from Mutation.gene annotations.
    Returns stable ordering for UI friendliness.
    """
    return sorted({m.gene for m in mutations if m.gene is not None})


def _compute_risk(mutations: list[Mutation]) -> dict[str, Any]:
    """
    Compute explainable risk summary from real Mutation objects.
    """
    if not mutations:
        return {
            "score": 0,
            "by_gene": {},
            "level": "Low",
            "explanation": "No mutations detected.",
        }

    return score_mutations(mutations)

def score_genome(
    record,
    *,
    reference_sequence: str | None = None,
    identify_mutations=None,
    map_genes=None,
    compute_risk=None,
) -> dict[str, Any]:
    identify_mutations = identify_mutations or _identify_mutations
    map_genes = map_genes or _map_genes
    compute_risk = compute_risk or score_mutations

    mutations = identify_mutations(record)
    genes_affected = map_genes(mutations)

    risk = compute_risk(mutations)
    risk_score = float(risk.get("score", 0.0))

    accession = record["accession"] if isinstance(record, dict) else record.accession
    source = record.get("source", "genbank") if isinstance(record, dict) else record.source

    return {
        "accession": accession,
        "source": source,
        "num_mutations": len(mutations),
        "genes_affected": genes_affected,
        "risk_score": risk_score,
        "risk_level": risk.get("level", ""),
        "risk_by_gene": risk.get("by_gene", {}),
        "risk_explanation": risk.get("explanation", ""),
    }

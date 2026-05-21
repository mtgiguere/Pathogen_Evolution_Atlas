"""
scoring.py
End-to-end scoring orchestration (dependency-injection friendly).
"""

from __future__ import annotations

import logging
from typing import Any

from .mutations import Mutation, diff_sequences, qc_compare_to_reference
from .risk import score_mutations

logger = logging.getLogger(__name__)


def _identify_mutations(record: Any):
    if isinstance(record, dict):
        ref = record.get("reference_sequence")
        sample = record.get("sequence")
        if not ref or not sample:
            return []
        try:
            return diff_sequences(ref, sample)
        except ValueError:
            return []

    ref = getattr(record, "reference_sequence", None)
    sample = getattr(record, "sequence", None)
    if not ref or not sample:
        return []

    try:
        return diff_sequences(ref, sample)
    except ValueError:
        return []


def _map_genes(mutations: list[Mutation]) -> list[str]:
    """Summarize affected genes; stable ordering for UI friendliness."""
    return sorted({m.gene for m in mutations if m.gene is not None})


def _compute_risk(mutations: list[Mutation]) -> dict[str, Any]:
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
    qc_fn=None,
    identify_mutations=None,
    map_genes=None,
    compute_risk=None,
) -> dict[str, Any]:
    identify_mutations = identify_mutations or _identify_mutations
    map_genes = map_genes or _map_genes
    compute_risk = compute_risk or _compute_risk
    qc_fn = qc_fn or qc_compare_to_reference

    sample_seq = (
        record.get("sequence") if isinstance(record, dict) else getattr(record, "sequence", None)
    )
    ref_seq = (
        record.get("reference_sequence")
        if isinstance(record, dict)
        else getattr(record, "reference_sequence", None)
    ) or reference_sequence

    accession = record["accession"] if isinstance(record, dict) else record.accession
    source = record.get("source", "genbank") if isinstance(record, dict) else record.source

    if not sample_seq:
        logger.warning("Missing sequence for %s", accession)
        return {
            "accession": accession,
            "source": source,
            "scorable": False,
            "qc_status": "FAIL",
            "qc_reasons": ["MISSING_SEQUENCE"],
            "num_mutations": 0,
            "genes_affected": [],
            "risk_score": None,
            "risk_level": "N/A",
            "risk_by_gene": {},
            "risk_explanation": "Not scored: MISSING_SEQUENCE",
        }

    if not ref_seq:
        logger.warning("Missing reference for %s", accession)
        return {
            "accession": accession,
            "source": source,
            "scorable": False,
            "qc_status": "FAIL",
            "qc_reasons": ["MISSING_REFERENCE"],
            "num_mutations": 0,
            "genes_affected": [],
            "risk_score": None,
            "risk_level": "N/A",
            "risk_by_gene": {},
            "risk_explanation": "Not scored: MISSING_REFERENCE",
        }

    # Require at least 80% overlap with the reference before scoring.
    min_overlap = int(0.8 * len(ref_seq))

    qc = qc_fn(ref_seq, sample_seq, min_overlap=min_overlap)

    if qc.status != "PASS":
        logger.warning("QC fail for %s: %s", accession, qc.reasons)
        return {
            "accession": accession,
            "source": source,
            "scorable": False,
            "qc_status": qc.status,
            "qc_reasons": qc.reasons,
            "num_mutations": 0,
            "genes_affected": [],
            "risk_score": None,
            "risk_level": "N/A",
            "risk_by_gene": {},
            "risk_explanation": f"Not scored: {', '.join(qc.reasons)}",
        }

    if isinstance(record, dict):
        record_for_scoring = dict(record)
        record_for_scoring.setdefault("reference_sequence", ref_seq)
    else:
        record_for_scoring = record
        if getattr(record_for_scoring, "reference_sequence", None) is None:
            record_for_scoring.reference_sequence = ref_seq

    mutations = identify_mutations(record_for_scoring)
    genes_affected = map_genes(mutations)

    risk = compute_risk(mutations)
    risk_score = float(risk.get("score", 0.0))

    logger.debug(
        "Scored %s: mutations=%d level=%s", accession, len(mutations), risk.get("level", "")
    )

    return {
        "accession": accession,
        "source": source,
        "scorable": True,
        "qc_status": "PASS",
        "qc_reasons": [],
        "num_mutations": len(mutations),
        "genes_affected": genes_affected,
        "risk_score": risk_score,
        "risk_level": risk.get("level", ""),
        "risk_by_gene": risk.get("by_gene", {}),
        "risk_explanation": risk.get("explanation", ""),
        "mutations_list": mutations,
    }

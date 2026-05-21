"""
analytics.py — Convert genome records into an analytics-ready DataFrame.
"""

from __future__ import annotations

import logging
from collections.abc import Iterable
from pathlib import Path

import pandas as pd

from .escape import EscapeMutation, escape_from_nt_mutations, escape_summary, load_escape_catalogue
from .growth import aggregate_by_week, estimate_growth_rates
from .lineage import LineageCall, LineageClassifier, load_signatures
from .scoring import score_genome

logger = logging.getLogger(__name__)

_DEFAULT_SIGNATURES_PATH = Path("data/lineages/signatures.json")
_DEFAULT_ESCAPE_PATH = Path("data/escape/catalogue.json")


def _escape_columns(matches: list) -> dict:
    summary = escape_summary(matches)
    antibodies = summary["antibodies_affected"]
    return {
        "escape_count": summary["total"],
        "escape_antibodies": ", ".join(antibodies),
        "has_critical_escape": summary["total"] > 0,
        "escape_mechanisms": ", ".join(sorted(summary["by_mechanism"])),
    }


def summarize_genomes(
    records: Iterable[dict],
    *,
    reference_sequence: str,
    reference_accession: str = "NC_045512.2",
    signatures_path: Path | None = None,
    escape_path: Path | None = None,
) -> pd.DataFrame:
    records = list(records)

    if not reference_sequence:
        raise ValueError("reference_sequence must be provided explicitly")

    ref_seq = reference_sequence

    sig_path = signatures_path or _DEFAULT_SIGNATURES_PATH
    try:
        sigs = load_signatures(sig_path)
        classifier: LineageClassifier | None = LineageClassifier(
            reference_sequence=ref_seq, signatures=sigs
        )
    except FileNotFoundError:
        logger.warning("Signatures file not found (%s); lineage classification disabled", sig_path)
        classifier = None

    esc_path = escape_path or _DEFAULT_ESCAPE_PATH
    try:
        escape_catalogue: list[EscapeMutation] = load_escape_catalogue(esc_path)
    except FileNotFoundError:
        logger.warning("Escape catalogue not found (%s); escape lookup disabled", esc_path)
        escape_catalogue = []

    def _get(rec, key, default=None):
        if isinstance(rec, dict):
            return rec.get(key, default)
        return getattr(rec, key, default)

    rows: list[dict] = []
    scorable_count = 0

    for r in records:
        seq = _get(r, "sequence") or ""

        if not seq:
            pre_scorable = False
            skip_reason = "missing_sequence"
        elif len(seq) < 1000:
            pre_scorable = False
            skip_reason = f"too_short ({len(seq)})"
        else:
            pre_scorable = True
            skip_reason = ""

        s = {
            "scorable": False,
            "qc_status": "FAIL",
            "qc_reasons": ["NOT_SCORED"],
            "num_mutations": 0,
            "genes_affected": [],
            "risk_score": None,
            "risk_level": "N/A",
            "risk_explanation": "Not scored",
            "risk_by_gene": {},
        }

        lineage_call: LineageCall | None = None
        esc_matches = []

        if pre_scorable:
            rec_for_scoring = {
                "accession": _get(r, "accession"),
                "source": _get(r, "source", "genbank"),
                "sequence": seq,
                "reference_sequence": ref_seq,
            }

            s = score_genome(rec_for_scoring)

            if bool(s.get("scorable", False)):
                scorable_count += 1
                mutations_list = s.get("mutations_list") or []
                if classifier is not None:
                    lineage_call = classifier.classify(mutations_list)
                if escape_catalogue:
                    esc_matches = escape_from_nt_mutations(
                        mutations_list, ref_seq, escape_catalogue
                    )
            else:
                reasons = list(s.get("qc_reasons", []) or [])
                skip_reason = "qc_fail: " + ", ".join(reasons) if reasons else "qc_fail"
        else:
            logger.debug(
                "Pre-scoring skip: accession=%s reason=%s", _get(r, "accession"), skip_reason
            )
            s["risk_explanation"] = "Not scored: " + skip_reason
            s["qc_reasons"] = [skip_reason]
            s["qc_status"] = "FAIL"

        num_mutations = int(s.get("num_mutations", 0) or 0)
        genes_list = list(s.get("genes_affected", []) or [])
        risk_score_raw = s.get("risk_score", None)
        risk_score = float(risk_score_raw) if risk_score_raw is not None else 0.0
        risk_level = s.get("risk_level", "N/A") or "N/A"
        risk_explanation = s.get("risk_explanation", "") or ""
        qc_status = s.get("qc_status", "FAIL") or "FAIL"
        qc_reasons_list = list(s.get("qc_reasons", []) or [])
        qc_reasons_str = ", ".join(qc_reasons_list)

        genes_str = ", ".join(genes_list)

        rows.append(
            {
                "accession": _get(r, "accession"),
                "source": _get(r, "source", "genbank"),
                "sequence_length": len(seq),
                "scorable": bool(s.get("scorable", False)) if pre_scorable else False,
                "skip_reason": skip_reason,
                "num_mutations": num_mutations,
                "genes_affected": genes_str,
                "risk_score": risk_score,
                "risk_level": risk_level,
                "risk_explanation": risk_explanation,
                "date": _get(r, "collection_date"),
                "lat": _get(r, "lat"),
                "lon": _get(r, "lon"),
                "qc_status": qc_status,
                "qc_reasons": qc_reasons_str,
                "qc_reasons_list": qc_reasons_list,
                "collection_date": _get(r, "collection_date"),
                "organism": _get(r, "organism"),
                "host": _get(r, "host"),
                "country": _get(r, "country"),
                "region": _get(r, "region"),
                "genes_affected_list": genes_list,
                "genes_affected_count": len(set(genes_list)),
                "mutations": list(s.get("mutations", []) or []),
                "mutations_str": ", ".join(list(s.get("mutations", []) or [])),
                "reference_accession": reference_accession,
                "reference_length": len(ref_seq),
                "lineage": lineage_call.lineage if lineage_call else "Unknown",
                "lineage_display": lineage_call.display_name if lineage_call else "Unknown",
                "who_label": lineage_call.who_label if lineage_call else "",
                "who_class": lineage_call.who_class if lineage_call else "",
                "lineage_confidence": lineage_call.confidence if lineage_call else 0.0,
                "lineage_supporting": (
                    len(lineage_call.supporting_mutations) if lineage_call else 0
                ),
                **_escape_columns(esc_matches),
            }
        )

    logger.info("summarize_genomes: %d / %d records scorable", scorable_count, len(records))
    return pd.DataFrame(rows)


def compute_variant_growth(
    df: pd.DataFrame,
    *,
    min_timepoints: int = 3,
    lineage_col: str = "lineage",
    date_col: str = "collection_date",
) -> pd.DataFrame:
    """
    Convenience wrapper: aggregate df by ISO week × lineage, then estimate growth rates.

    Returns a DataFrame with columns:
      lineage, growth_rate, doubling_time_days, r_squared, n_timepoints, trend
    """
    weekly = aggregate_by_week(df, lineage_col=lineage_col, date_col=date_col)
    rates = estimate_growth_rates(weekly, min_timepoints=min_timepoints)
    logger.info(
        "compute_variant_growth: %d lineages assessed (%d Growing, %d Declining)",
        len(rates),
        int((rates["trend"] == "Growing").sum()) if not rates.empty else 0,
        int((rates["trend"] == "Declining").sum()) if not rates.empty else 0,
    )
    return rates

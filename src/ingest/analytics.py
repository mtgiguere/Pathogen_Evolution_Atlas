"""
analytics.py
"""
from __future__ import annotations

from collections.abc import Iterable

import pandas as pd

from ingest.scoring import score_genome

"""
Convert genome records into an analytics-ready DataFrame.
"""
"""
Convert genome records into an analytics-ready DataFrame.
"""
def summarize_genomes(
    records: Iterable[dict],
    *,
    reference_sequence: str,
    reference_accession: str = "NC_045512.2",
) -> pd.DataFrame:
    records = list(records)

    if not reference_sequence:
        raise ValueError("reference_sequence must be provided explicitly")

    ref_seq = reference_sequence

    def _get(rec, key, default=None):
        if isinstance(rec, dict):
            return rec.get(key, default)
        return getattr(rec, key, default)

    rows: list[dict] = []
    for r in records:
        seq = _get(r, "sequence") or ""

        # --- Keep existing lightweight gates for backwards compatibility ---
        if not seq:
            pre_scorable = False
            skip_reason = "missing_sequence"
        elif len(seq) < 1000:
            pre_scorable = False
            skip_reason = f"too_short ({len(seq)})"
        else:
            pre_scorable = True
            skip_reason = ""

        # Defaults (not scored)
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

        if pre_scorable:
            rec_for_scoring = {
                "accession": _get(r, "accession"),
                "source": _get(r, "source", "genbank"),
                "sequence": seq,
                "reference_sequence": ref_seq,  # embed for compatibility with test stubs
            }

            s = score_genome(rec_for_scoring)

            # If scoring says "not scorable" (QC fail), reflect that in skip_reason
            if not bool(s.get("scorable", False)):
                # Keep older-style skip_reason string but tie it to QC
                reasons = list(s.get("qc_reasons", []) or [])
                skip_reason = "qc_fail: " + ", ".join(reasons) if reasons else "qc_fail"
        else:
            # keep your older explanation format
            s["risk_explanation"] = "Not scored: " + skip_reason
            s["qc_reasons"] = [skip_reason]
            s["qc_status"] = "FAIL"

        # Normalize outputs (backward compatible)
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
                # ---- columns your existing tests expect ----
                "accession": _get(r, "accession"),
                "source": _get(r, "source", "genbank"),
                "sequence_length": len(seq),
                "scorable": bool(s.get("scorable", False)) if pre_scorable else False,
                "skip_reason": skip_reason,
                "num_mutations": num_mutations,
                "genes_affected": genes_str,  # keep as string for backwards compat
                "risk_score": risk_score,
                "risk_level": risk_level,
                "risk_explanation": risk_explanation,
                "date": _get(r, "collection_date"),
                "lat": _get(r, "lat"),
                "lon": _get(r, "lon"),

                # ---- QC surfaced (new) ----
                "qc_status": qc_status,
                "qc_reasons": qc_reasons_str,
                "qc_reasons_list": qc_reasons_list,

                # ---- extra columns (existing optional) ----
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
            }
        )

    return pd.DataFrame(rows)

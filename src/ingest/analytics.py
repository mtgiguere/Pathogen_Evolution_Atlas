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

        # Decide if we can score this record
        if not ref_seq:
            scorable = False
            skip_reason = "missing_reference"
        elif not seq:
            scorable = False
            skip_reason = "missing_sequence"
        elif len(seq) < 1000:
            scorable = False
            skip_reason = f"too_short ({len(seq)})"
        else:
            scorable = True
            skip_reason = ""

        if scorable:
            rec_for_scoring = {
                "accession": _get(r, "accession"),
                "source": _get(r, "source", "genbank"),
                "sequence": seq,
                "reference_sequence": ref_seq,
            }
            s = score_genome(rec_for_scoring)

            num_mutations = int(s.get("num_mutations", 0))
            genes_list = list(s.get("genes_affected", []) or [])
            mutations_list = list(s.get("mutations", []) or [])

            risk_score = float(s.get("risk_score", 0.0))
            risk_level = s.get("risk_level", "N/A")
            risk_explanation = s.get("risk_explanation", "")
        else:
            num_mutations = 0
            genes_list = []
            mutations_list = []
            risk_score = 0.0
            risk_level = "N/A"
            risk_explanation = "Not scored: " + skip_reason

        genes_str = ", ".join(genes_list)
        muts_str = ", ".join(mutations_list)

        rows.append(
            {
                # ---- columns your existing tests expect ----
                "accession": _get(r, "accession"),
                "source": _get(r, "source", "genbank"),
                "sequence_length": len(seq),
                "scorable": scorable,
                "skip_reason": skip_reason,
                "num_mutations": num_mutations,
                "genes_affected": genes_str,  # NOTE: keep as string for backwards compat
                "risk_score": risk_score,
                "risk_level": risk_level,
                "risk_explanation": risk_explanation,
                "date": _get(r, "collection_date"),  # NOTE: keep name 'date' for tests
                "lat": _get(r, "lat"),
                "lon": _get(r, "lon"),

                # ---- extra “Tableau surface” columns (new, optional) ----
                "collection_date": _get(r, "collection_date"),
                "organism": _get(r, "organism"),
                "host": _get(r, "host"),
                "country": _get(r, "country"),
                "region": _get(r, "region"),
                "genes_affected_list": genes_list,
                "genes_affected_count": len(set(genes_list)),
                "mutations": mutations_list,
                "mutations_str": muts_str,
                "reference_accession": reference_accession,
                "reference_length": len(ref_seq),
            }
        )

    return pd.DataFrame(rows)

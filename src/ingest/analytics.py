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
def summarize_genomes(records: Iterable[dict]) -> pd.DataFrame:
    records = list(records)

    def _get(rec, key, default=None):
        if isinstance(rec, dict):
            return rec.get(key, default)
        return getattr(rec, key, default)

    # Select reference genome once:
    # Prefer NC_045512 if present; otherwise fall back to the longest sequence available.
    ref_rec = next((r for r in records if _get(r, "accession", "").startswith("NC_045512")), None)

    if ref_rec is None:
        ref_rec = max(records, key=lambda r: len(_get(r, "sequence") or ""), default=None)

    ref_seq = _get(ref_rec, "sequence") if ref_rec is not None else None

    # show first record sequence length too
    first = records[0] if records else None


    rows = []
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
            # Still skip obvious fragments (tweak threshold as desired)
            scorable = False
            skip_reason = f"too_short ({len(seq)})"
        # elif len(seq) != len(ref_seq):
        #     scorable = False
        #     skip_reason = f"length_mismatch ({len(seq)} != {len(ref_seq)})"
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
            num_mutations = s["num_mutations"]
            genes_affected = ", ".join(s["genes_affected"])
            risk_score = s["risk_score"]
            risk_level = s["risk_level"]
            risk_explanation = s["risk_explanation"]
        else:
            # Don't score; keep app running and make the data quality visible
            num_mutations = 0
            genes_affected = ""
            risk_score = 0.0
            risk_level = "N/A"
            risk_explanation = "Not scored: " + skip_reason

        rows.append(
            {
                "accession": _get(r, "accession"),
                "source": _get(r, "source", "genbank"),
                "sequence_length": len(seq),
                "scorable": scorable,
                "skip_reason": skip_reason,
                "num_mutations": num_mutations,
                "genes_affected": genes_affected,
                "risk_score": risk_score,
                "risk_level": risk_level,
                "risk_explanation": risk_explanation,
                "date": _get(r, "collection_date"),
                "lat": _get(r, "lat"),
                "lon": _get(r, "lon"),
            }
        )

    df = pd.DataFrame(rows)
    return df


    return pd.DataFrame(rows)

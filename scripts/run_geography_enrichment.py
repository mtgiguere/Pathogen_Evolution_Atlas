"""
One-shot runner for geography enrichment.

Reads canonical genome records, enriches missing geography
(country/region) via BioSample, and writes an enriched NDJSON.

This script is intentionally NOT imported anywhere.
It is meant to be run manually.
"""

import os
from pathlib import Path

from ingest.geography_enrichment import enrich_many_locations
from ingest.io import load_ndjson, write_ndjson

# --- Config ---
EMAIL = os.environ.get("NCBI_EMAIL")
IN_PATH = Path("data/raw/genomes.ndjson")
OUT_PATH = Path("data/raw/genomes.enriched.ndjson")


def main() -> None:
    if not EMAIL:
        raise RuntimeError("NCBI_EMAIL environment variable is required for Entrez access.")

    if not IN_PATH.exists():
        raise FileNotFoundError(f"Input file not found: {IN_PATH}")

    records = load_ndjson(IN_PATH)

    before_missing = sum(1 for r in records if not r.country)

    enriched = enrich_many_locations(records, email=EMAIL)

    after_missing = sum(1 for r in enriched if not r.country)
    filled = before_missing - after_missing

    write_ndjson(enriched, OUT_PATH)

    print(
        f"[geo] country missing: {before_missing}/{len(records)} "
        f"-> {after_missing}/{len(records)} (filled {filled})"
    )
    print(f"[geo] wrote enriched file: {OUT_PATH}")


if __name__ == "__main__":
    main()

"""
fetch_genbank_accessions.py
"""
from __future__ import annotations

import argparse
import os

from ingest.genbank import fetch_and_normalize_many
from ingest.io import write_ndjson


def main() -> None:
    p = argparse.ArgumentParser(description="Fetch GenBank accessions and write canonical NDJSON.")
    p.add_argument("--accessions", nargs="+", required=True, help="One or more GenBank accessions")
    p.add_argument("--out", required=True, help="Output path for NDJSON")
    args = p.parse_args()

    email = os.getenv("NCBI_EMAIL")
    if not email:
        raise SystemExit("NCBI_EMAIL environment variable not set (required by NCBI).")

    records = fetch_and_normalize_many(args.accessions, email=email)
    write_ndjson(records, args.out)

    print(f"Wrote {len(records)} records -> {args.out}")


if __name__ == "__main__":
    main()

"""
run_ingest.py — Automated incremental ingestion CLI.

Meant to be invoked periodically (cron, Task Scheduler, or manually):

    python scripts/run_ingest.py
    python scripts/run_ingest.py --interval-hours 12 --max-new 500 --force

On each run it:
  1. Checks whether enough time has elapsed since the last run (skip if not, unless --force).
  2. Searches NCBI for SARS-CoV-2 complete genomes submitted/updated since the last known date.
  3. Deduplicates against accessions already in the local data store.
  4. Fetches and normalises new records.
  5. Appends them to the genomes NDJSON file.
  6. Persists updated state for the next run.
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
from datetime import UTC, date
from pathlib import Path

from Bio import Entrez

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.ingest.genbank import fetch_many_genbank_minimal
from src.ingest.io import read_ndjson, write_ndjson
from src.ingest.scheduler import (
    IngestState,
    load_state,
    run_incremental_ingest,
    save_state,
    should_run,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s — %(message)s",
    datefmt="%Y-%m-%dT%H:%M:%S",
)
logger = logging.getLogger("run_ingest")

_STATE_PATH = Path("data/state/ingest_state.json")
_GENOMES_PATH = Path(os.getenv("GENOMES_PATH", "data/raw/genomes.ndjson"))
_DEFAULT_INTERVAL_HOURS = 24.0
_DEFAULT_MAX_NEW = 200


def _build_search_query(min_len: int = 29000, max_len: int = 31000) -> str:
    return (
        '"Severe acute respiratory syndrome coronavirus 2"[Organism] '
        'AND "complete genome"[Title] '
        f"AND {min_len}:{max_len}[SLEN]"
    )


def _search_ncbi(since: date | None, email: str, retmax: int) -> list[str]:
    Entrez.email = email
    query = _build_search_query()

    kwargs: dict = {"db": "nuccore", "term": query, "retmax": retmax}
    if since is not None:
        kwargs["mindate"] = since.strftime("%Y/%m/%d")
        kwargs["datetype"] = "pdat"

    with Entrez.esearch(**kwargs) as handle:
        res = Entrez.read(handle)

    ids = res.get("IdList", [])
    if not ids:
        return []

    with Entrez.efetch(db="nuccore", id=",".join(ids), rettype="acc", retmode="text") as handle:
        text = handle.read()

    accessions = [line.strip() for line in text.splitlines() if line.strip()]
    seen: set[str] = set()
    deduped: list[str] = []
    for a in accessions:
        if a not in seen:
            seen.add(a)
            deduped.append(a)
    return deduped


def _load_existing_accessions(path: Path) -> set[str]:
    if not path.exists():
        return set()
    records = read_ndjson(path)
    return {r.get("accession", r.get("accession", "")) for r in records if r.get("accession")}


def main() -> None:
    ap = argparse.ArgumentParser(description="Incremental SARS-CoV-2 genome ingestion from NCBI.")
    ap.add_argument("--interval-hours", type=float, default=_DEFAULT_INTERVAL_HOURS)
    ap.add_argument("--max-new", type=int, default=_DEFAULT_MAX_NEW)
    ap.add_argument("--force", action="store_true", help="Run even if interval has not elapsed")
    ap.add_argument("--state-path", type=Path, default=_STATE_PATH)
    ap.add_argument("--out", type=Path, default=_GENOMES_PATH)
    args = ap.parse_args()

    email = os.getenv("NCBI_EMAIL")
    if not email:
        raise SystemExit("NCBI_EMAIL environment variable not set.")

    state = load_state(args.state_path)

    if not args.force and not should_run(state, args.interval_hours):
        logger.info(
            "Last run was recent (< %.1f h); skipping. Use --force to override.",
            args.interval_hours,
        )
        return

    existing = _load_existing_accessions(args.out)
    logger.info("Existing accessions in store: %d", len(existing))

    def search_fn(since: date | None) -> list[str]:
        return _search_ncbi(since, email=email, retmax=args.max_new)

    def fetch_fn(accs: list[str]) -> list[dict]:
        return fetch_many_genbank_minimal(accs, email=email)

    result = run_incremental_ingest(
        state=state,
        search_fn=search_fn,
        fetch_fn=fetch_fn,
        existing_accessions=existing,
    )

    if result.error_count:
        logger.error("Fetch failed — %d error(s). State not updated.", result.error_count)
        sys.exit(1)

    if result.new_count > 0:
        # Append new raw records to the genomes file
        args.out.parent.mkdir(parents=True, exist_ok=True)
        existing_records = list(read_ndjson(args.out)) if args.out.exists() else []
        write_ndjson(existing_records + result.records, args.out)
        logger.info("Appended %d new records to %s", result.new_count, args.out)

    # Update state
    new_total = state.total_fetched + result.new_count
    # Advance last_accession_date to today so the next run won't re-scan old records
    new_state = IngestState(
        last_run=result.run_at,
        last_accession_date=result.run_at.astimezone(UTC).date(),
        total_fetched=new_total,
    )
    save_state(new_state, args.state_path)
    logger.info(
        "Done. new=%d skipped=%d total_fetched=%d",
        result.new_count,
        result.skipped_count,
        new_total,
    )


if __name__ == "__main__":
    main()

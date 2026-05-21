"""
scheduler.py — Incremental ingestion state management and pipeline runner.

Tracks when the pipeline last ran and which accessions have already been
fetched so that each run only pulls genuinely new sequences from NCBI.
"""

from __future__ import annotations

import json
import logging
from collections.abc import Callable
from dataclasses import dataclass, field
from datetime import UTC, date, datetime
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)

# ── State ─────────────────────────────────────────────────────────────────────


@dataclass
class IngestState:
    last_run: datetime | None = None
    last_accession_date: date | None = None
    total_fetched: int = 0


@dataclass
class IngestResult:
    new_count: int = 0
    skipped_count: int = 0
    error_count: int = 0
    run_at: datetime = field(default_factory=lambda: datetime.now(tz=UTC))
    records: list[dict[str, Any]] = field(default_factory=list)


def load_state(path: Path) -> IngestState:
    """Load persisted ingest state; returns a default state if the file is absent."""
    if not path.exists():
        return IngestState()

    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        logger.warning("Could not parse state file %s; using defaults", path)
        return IngestState()

    last_run = None
    if raw.get("last_run"):
        try:
            last_run = datetime.fromisoformat(raw["last_run"])
        except ValueError:
            pass

    last_acc_date = None
    if raw.get("last_accession_date"):
        try:
            last_acc_date = date.fromisoformat(raw["last_accession_date"])
        except ValueError:
            pass

    return IngestState(
        last_run=last_run,
        last_accession_date=last_acc_date,
        total_fetched=int(raw.get("total_fetched", 0)),
    )


def save_state(state: IngestState, path: Path) -> None:
    """Persist ingest state to JSON, creating parent directories as needed."""
    path.parent.mkdir(parents=True, exist_ok=True)
    payload: dict[str, Any] = {
        "last_run": state.last_run.isoformat() if state.last_run else None,
        "last_accession_date": (
            state.last_accession_date.isoformat() if state.last_accession_date else None
        ),
        "total_fetched": state.total_fetched,
    }
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


# ── Scheduling logic ──────────────────────────────────────────────────────────


def should_run(
    state: IngestState,
    interval_hours: float,
    now: datetime | None = None,
) -> bool:
    """Return True if enough time has passed since the last run (or it has never run)."""
    if state.last_run is None:
        return True
    _now = now or datetime.now(tz=UTC)
    last = state.last_run
    if last.tzinfo is None:
        last = last.replace(tzinfo=UTC)
    elapsed_hours = (_now - last).total_seconds() / 3600.0
    return elapsed_hours >= interval_hours


# ── Pipeline runner ───────────────────────────────────────────────────────────


def run_incremental_ingest(
    state: IngestState,
    *,
    search_fn: Callable[[date | None], list[str]],
    fetch_fn: Callable[[list[str]], list[dict[str, Any]]],
    existing_accessions: set[str],
) -> IngestResult:
    """
    Run one incremental ingestion cycle.

    Parameters
    ----------
    state:
        Current ingest state (read-only here; caller updates after writing records).
    search_fn:
        Callable(since: date | None) -> list[str] — returns accession IDs from NCBI.
        Receives the last known accession date so only newer records are returned.
    fetch_fn:
        Callable(accessions: list[str]) -> list[dict] — fetches raw GenBank dicts.
        Called once with the filtered (de-duplicated) accession list.
    existing_accessions:
        Set of accession strings already in the local data store; used for dedup.

    Returns
    -------
    IngestResult with counts and the newly fetched raw record dicts.
    """
    run_at = datetime.now(tz=UTC)
    result = IngestResult(run_at=run_at)

    found = search_fn(state.last_accession_date)

    new_accs = [a for a in found if a not in existing_accessions]
    result.skipped_count = len(found) - len(new_accs)

    if not new_accs:
        logger.info("run_incremental_ingest: no new accessions found")
        return result

    try:
        records = fetch_fn(new_accs)
    except Exception:
        logger.exception("run_incremental_ingest: fetch failed for %d accessions", len(new_accs))
        result.error_count = 1
        return result

    result.records = records
    result.new_count = len(records)
    logger.info(
        "run_incremental_ingest: fetched %d new, skipped %d existing",
        result.new_count,
        result.skipped_count,
    )
    return result

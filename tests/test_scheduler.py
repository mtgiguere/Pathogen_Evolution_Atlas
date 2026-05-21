"""
Periodic ingestion scheduler tests — TDD, written before implementation.
All tests expected to FAIL until src/ingest/scheduler.py is implemented.
"""

from __future__ import annotations

import json
import tempfile
from datetime import UTC, date, datetime, timedelta
from pathlib import Path

from src.ingest.scheduler import (
    IngestResult,
    IngestState,
    load_state,
    run_incremental_ingest,
    save_state,
    should_run,
)

# ── IngestState defaults ──────────────────────────────────────────────────────


def test_ingest_state_defaults():
    s = IngestState()
    assert s.last_run is None
    assert s.last_accession_date is None
    assert s.total_fetched == 0


# ── load_state / save_state ───────────────────────────────────────────────────


def test_load_state_missing_file_returns_default():
    state = load_state(Path("data/state/nonexistent.json"))
    assert state.last_run is None
    assert state.total_fetched == 0


def test_save_and_load_state_roundtrip():
    ts = datetime(2025, 6, 1, 12, 0, 0, tzinfo=UTC)
    acc_date = date(2025, 5, 31)
    state = IngestState(last_run=ts, last_accession_date=acc_date, total_fetched=42)

    with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
        tmp = Path(f.name)

    save_state(state, tmp)
    loaded = load_state(tmp)

    assert loaded.last_run == ts
    assert loaded.last_accession_date == acc_date
    assert loaded.total_fetched == 42


def test_save_state_creates_parent_dirs(tmp_path):
    state = IngestState(total_fetched=5)
    deep = tmp_path / "a" / "b" / "state.json"
    save_state(state, deep)
    assert deep.exists()
    loaded = load_state(deep)
    assert loaded.total_fetched == 5


def test_load_state_partial_json():
    """A file with only some fields set should not crash; missing fields fall back to defaults."""
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".json", delete=False, encoding="utf-8"
    ) as f:
        json.dump({"total_fetched": 7}, f)
        tmp = Path(f.name)

    loaded = load_state(tmp)
    assert loaded.total_fetched == 7
    assert loaded.last_run is None
    assert loaded.last_accession_date is None


# ── should_run ────────────────────────────────────────────────────────────────

_NOW = datetime(2025, 6, 1, 12, 0, 0, tzinfo=UTC)


def test_should_run_when_never_run():
    state = IngestState(last_run=None)
    assert should_run(state, interval_hours=24, now=_NOW) is True


def test_should_run_when_interval_elapsed():
    last = _NOW - timedelta(hours=25)
    state = IngestState(last_run=last)
    assert should_run(state, interval_hours=24, now=_NOW) is True


def test_should_run_false_when_recently_run():
    last = _NOW - timedelta(hours=1)
    state = IngestState(last_run=last)
    assert should_run(state, interval_hours=24, now=_NOW) is False


def test_should_run_exactly_at_interval_boundary():
    last = _NOW - timedelta(hours=24)
    state = IngestState(last_run=last)
    # Exactly at boundary → True (elapsed >= interval)
    assert should_run(state, interval_hours=24, now=_NOW) is True


# ── run_incremental_ingest ────────────────────────────────────────────────────

_ACC_DATE = date(2025, 5, 15)


def _make_raw(accession: str) -> dict:
    return {
        "accession": accession,
        "organism": "SARS-CoV-2",
        "collection_date": "2025-05-20",
        "location": "USA: California",
        "host": "Homo sapiens",
        "sequence": "ATCG" * 100,
        "lat_lon": None,
        "biosample": None,
    }


def test_run_incremental_ingest_no_new_accessions():
    result = run_incremental_ingest(
        state=IngestState(),
        search_fn=lambda since: [],
        fetch_fn=lambda accs: [],
        existing_accessions=set(),
    )
    assert isinstance(result, IngestResult)
    assert result.new_count == 0
    assert result.skipped_count == 0
    assert result.error_count == 0


def test_run_incremental_ingest_all_new():
    accs = ["ACC001", "ACC002", "ACC003"]
    result = run_incremental_ingest(
        state=IngestState(),
        search_fn=lambda since: accs,
        fetch_fn=lambda a: [_make_raw(x) for x in a],
        existing_accessions=set(),
    )
    assert result.new_count == 3
    assert result.skipped_count == 0
    assert result.error_count == 0


def test_run_incremental_ingest_deduplicates_existing():
    accs = ["ACC001", "ACC002", "ACC003"]
    existing = {"ACC001", "ACC002"}  # already in store
    result = run_incremental_ingest(
        state=IngestState(),
        search_fn=lambda since: accs,
        fetch_fn=lambda a: [_make_raw(x) for x in a],
        existing_accessions=existing,
    )
    assert result.new_count == 1    # only ACC003 is new
    assert result.skipped_count == 2


def test_run_incremental_ingest_passes_since_date_to_search():
    received: list = []

    def search(since):
        received.append(since)
        return []

    state = IngestState(last_accession_date=_ACC_DATE)
    run_incremental_ingest(
        state=state,
        search_fn=search,
        fetch_fn=lambda a: [],
        existing_accessions=set(),
    )
    assert received == [_ACC_DATE]


def test_run_incremental_ingest_search_gets_none_when_no_prior_date():
    received: list = []

    def search(since):
        received.append(since)
        return []

    run_incremental_ingest(
        state=IngestState(),
        search_fn=search,
        fetch_fn=lambda a: [],
        existing_accessions=set(),
    )
    assert received == [None]


def test_run_incremental_ingest_result_has_run_at_timestamp():
    result = run_incremental_ingest(
        state=IngestState(),
        search_fn=lambda since: [],
        fetch_fn=lambda a: [],
        existing_accessions=set(),
    )
    assert isinstance(result.run_at, datetime)
    assert result.run_at.tzinfo is not None  # timezone-aware


def test_run_incremental_ingest_fetch_error_counted():
    def bad_fetch(accs):
        raise RuntimeError("NCBI unavailable")

    result = run_incremental_ingest(
        state=IngestState(),
        search_fn=lambda since: ["ACC001"],
        fetch_fn=bad_fetch,
        existing_accessions=set(),
    )
    assert result.error_count == 1
    assert result.new_count == 0


def test_run_incremental_ingest_returns_fetched_records():
    accs = ["ACC001", "ACC002"]
    result = run_incremental_ingest(
        state=IngestState(),
        search_fn=lambda since: accs,
        fetch_fn=lambda a: [_make_raw(x) for x in a],
        existing_accessions=set(),
    )
    assert len(result.records) == 2
    assert result.records[0]["accession"] == "ACC001"


def test_run_incremental_ingest_empty_search_returns_empty_records():
    result = run_incremental_ingest(
        state=IngestState(),
        search_fn=lambda since: [],
        fetch_fn=lambda a: [],
        existing_accessions=set(),
    )
    assert result.records == []

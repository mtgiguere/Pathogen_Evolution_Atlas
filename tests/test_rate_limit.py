"""
Testing rate-limit behavior: timing, thread safety, and logging.
"""

import logging


def test_rate_limiter_sleeps_when_called_too_fast(monkeypatch):
    """
    Deterministic unit test: if two requests happen at the same 'time',
    we should call sleep() to enforce <= 3 req/sec.
    """
    from src.ingest import genbank

    slept = []
    monkeypatch.setattr(genbank.time, "sleep", lambda s: slept.append(s))

    t = {"now": 100.0}
    monkeypatch.setattr(genbank.time, "monotonic", lambda: t["now"])
    monkeypatch.setattr(genbank, "_LAST_REQUEST_TS", None)

    genbank._rate_limit()
    genbank._rate_limit()

    assert len(slept) == 1
    assert slept[0] >= (1.0 / 3.0)


def test_rate_limiter_thread_safe(monkeypatch):
    """Two threads hitting _rate_limit at the same instant: exactly one must sleep."""
    import threading

    from src.ingest import genbank

    slept = []
    mu = threading.Lock()

    monkeypatch.setattr(
        genbank.time,
        "sleep",
        lambda s: (mu.acquire(), slept.append(s), mu.release()),
    )
    monkeypatch.setattr(genbank.time, "monotonic", lambda: 100.0)
    monkeypatch.setattr(genbank, "_LAST_REQUEST_TS", None)

    barrier = threading.Barrier(2)

    def worker():
        barrier.wait()
        genbank._rate_limit()

    threads = [threading.Thread(target=worker) for _ in range(2)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()

    assert len(slept) == 1


def test_rate_limiter_logs_sleep(monkeypatch, caplog):
    """A rate-limit sleep must emit a DEBUG log."""
    from src.ingest import genbank

    monkeypatch.setattr(genbank.time, "sleep", lambda s: None)
    monkeypatch.setattr(genbank.time, "monotonic", lambda: 100.0)
    monkeypatch.setattr(genbank, "_LAST_REQUEST_TS", None)

    with caplog.at_level(logging.DEBUG):
        genbank._rate_limit()  # first call — no sleep
        genbank._rate_limit()  # second at same instant — sleeps and logs

    assert any("sleep" in r.message.lower() or "rate" in r.message.lower() for r in caplog.records)

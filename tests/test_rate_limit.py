"""
Testing DNS/HTTPS/NCBI response
"""
def test_rate_limiter_sleeps_when_called_too_fast(monkeypatch):
    """
    Deterministic unit test: if two requests happen at the same 'time',
    we should call sleep() to enforce <= 3 req/sec.
    """
    from src.ingest import genbank

    # Capture sleep durations instead of actually sleeping
    slept = []
    monkeypatch.setattr(genbank.time, "sleep", lambda s: slept.append(s))

    # Control time so both calls happen "instantly"
    t = {"now": 100.0}
    monkeypatch.setattr(genbank.time, "monotonic", lambda: t["now"])

    # Reset internal state for the test
    genbank._LAST_REQUEST_TS = None

    # Call twice without advancing time
    genbank._rate_limit()
    genbank._rate_limit()

    assert len(slept) == 1
    assert slept[0] >= (1.0 / 3.0)



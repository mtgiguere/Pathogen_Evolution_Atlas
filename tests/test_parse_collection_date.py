from datetime import date

import pytest

from ingest.genbank import parse_collection_date


@pytest.mark.parametrize(
    "raw, expected",
    [
        (None, None),
        ("", None),
        ("unknown", None),
        ("n/a", None),
        ("2020", date(2020, 1, 1)),
        ("2020-12", date(2020, 12, 1)),
        ("2020-12-31", date(2020, 12, 31)),
        ("Dec-2020", date(2020, 12, 1)),
        ("2020-Dec-01", date(2020, 12, 1)),
        # malformed real-world-ish cases -> None (don’t crash)
        ("2020-Feb-30", None),
        ("2020-02-30", None),
        ("2020-Apr-31", None),
    ],
)
def test_parse_collection_date_handles_common_and_invalid_formats(raw, expected):
    assert parse_collection_date(raw) == expected

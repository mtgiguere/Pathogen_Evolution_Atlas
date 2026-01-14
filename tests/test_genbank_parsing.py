from datetime import date

# We are importing the function we are testing.
# This is the heart of Test-Driven Development:
#   "Does this function behave the way we expect?"
from src.ingest.genbank import parse_collection_date


def test_parse_collection_date_full():
    """
    If we are given a full ISO-style date (YYYY-MM-DD),
    the function should return the exact same date
    as a Python `date` object.
    """
    assert parse_collection_date("2024-08-19") == date(2024, 8, 19)


def test_parse_collection_date_year_month():
    """
    If we are only given year and month (YYYY-MM),
    we assume the day is the first of that month.
    This keeps our data consistent and sortable.
    """
    assert parse_collection_date("2024-08") == date(2024, 8, 1)

def test_parse_collection_date_year_only():
    """
    If we are only given a year (YYYY),
    we assume January 1st of that year.
    This keeps ordering and timelines consistent.
    """
    assert parse_collection_date("2024") == date(2024, 1, 1)

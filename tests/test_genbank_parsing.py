"""
Unit Tests genbank.py.
"""

#Imports needed for this python page
from datetime import date

# We are importing the function we are testing.
# This is the heart of Test-Driven Development:
#   "Does this function behave the way we expect?"
from src.ingest.genbank import parse_collection_date, parse_location


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

def test_parse_collection_date_unknown_returns_none():
    """
    If GenBank provides 'unknown', an empty string,
    or None, we should return None instead of crashing.
    """
    assert parse_collection_date("unknown") is None
    assert parse_collection_date("") is None
    assert parse_collection_date(None) is None

def test_parse_location_country_region():
    """
    If location is formatted like 'Country: Region',
    we split into (country, region).
    """
    assert parse_location("USA: Illinois") == ("USA", "Illinois")

def test_parse_location_country_only():
    """
    If only a country is provided, region should be None.
    """
    assert parse_location("USA") == ("USA", None)

def test_parse_location_empty():
    """
    Missing or empty location values should return (None, None).
    """
    assert parse_location("") == (None, None)
    assert parse_location(None) == (None, None)

def test_parse_collection_date_month_name_year():
    """
    GenBank sometimes uses month abbreviations like 'Dec-2019'.
    We normalize that to the first day of that month.
    """
    assert parse_collection_date("Dec-2019") == date(2019, 12, 1)

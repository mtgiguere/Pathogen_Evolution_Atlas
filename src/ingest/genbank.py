from datetime import date
from typing import Optional


def parse_collection_date(raw: Optional[str]) -> Optional[date]:
    """
    Convert a GenBank-style collection date string into a Python date object.

    GenBank often stores dates in partially-known forms:
      - "YYYY-MM-DD"  (full date)
      - "YYYY-MM"     (month known, day unknown)
      - "YYYY"        (only the year is known)

    This function normalizes all of those into real `date` objects
    so that downstream code can sort, compare, and model them reliably.

    If the date is missing or unknown, we return None.
    """

    # If the input is None or an empty string, we cannot parse a date.
    # Returning None allows the rest of the pipeline to handle "unknown".
    if not raw:
        return None

    # Split the string on "-" so:
    #   "2024-08-19" -> ["2024", "08", "19"]
    #   "2024-08"    -> ["2024", "08"]
    #   "2024"       -> ["2024"]
    parts = raw.split("-")

    # Case 1: full date (year, month, day)
    if len(parts) == 3:
        y, m, d = parts
        return date(int(y), int(m), int(d))

    # Case 2: year and month only → assume day = 1
    # This lets us still place the sample on a timeline.
    if len(parts) == 2:
        y, m = parts
        return date(int(y), int(m), 1)

    # Case 3: year only → assume January 1st
    if len(parts) == 1:
        y = parts[0]
        return date(int(y), 1, 1)

    # Anything else is malformed → treat as unknown
    return None

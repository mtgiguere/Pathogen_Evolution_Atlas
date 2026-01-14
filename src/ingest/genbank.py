"""
genbank.py
"""
from datetime import date
from typing import Optional
from typing import Tuple
from typing import Any, Dict
from .models import CanonicalGenomeRecord
from typing import Any, Dict, Optional
from Bio import Entrez, SeqIO



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
    # Month mapping
    MONTHS = {
    "jan": 1, "feb": 2, "mar": 3, "apr": 4, "may": 5, "jun": 6,
    "jul": 7, "aug": 8, "sep": 9, "oct": 10, "nov": 11, "dec": 12,
    }


    # If the input is None or an empty string, we cannot parse a date.
    # Returning None allows the rest of the pipeline to handle "unknown".
    if not raw:
        return None

    # GenBank sometimes uses strings like "unknown" instead of a real date.
    # Treat those as missing data.
    s = str(raw).strip()
    if s.lower() in {"unknown", "na", "n/a", "none"}:
        return None
    
    # Split the string on "-" so:
    #   "2024-08-19" -> ["2024", "08", "19"]
    #   "2024-08"    -> ["2024", "08"]
    #   "2024"       -> ["2024"]
    parts = s.split("-")

    # Case 1: full date (year, month, day)
    if len(parts) == 3:
        y, m, d = parts
        return date(int(y), int(m), int(d))

    # Case 2: year and month only → assume day = 1
    # This lets us still place the sample on a timeline.
    if len(parts) == 2:
        a, b = parts[0].strip(), parts[1].strip()

        # GenBank sometimes uses "Dec-2019" (month name + year)
        if a[:3].lower() in MONTHS and b.isdigit():
            return date(int(b), MONTHS[a[:3].lower()], 1)

        # Otherwise assume numeric "YYYY-MM"
        y, m = a, b
        return date(int(y), int(m), 1)


    # Case 3: year only → assume January 1st
    if len(parts) == 1:
        y = parts[0]
        return date(int(y), 1, 1)

    # Anything else is malformed → treat as unknown
    return None

def parse_location(raw: Optional[str]) -> Tuple[Optional[str], Optional[str]]:
    """
    Normalize a location string into (country, region).

    Expected common format:
      - "USA: Illinois" -> ("USA", "Illinois")
    """
    if not raw:
        return (None, None)

    s = str(raw).strip()
    if not s:
        return (None, None)

    if ":" in s:
        country, region = s.split(":", 1)
        return (country.strip() or None, region.strip() or None)

    return (s, None)

def fetch_genbank_minimal(accession: str, email: str) -> Dict[str, Any]:
    """
    Fetch a single GenBank record from NCBI and extract only the fields we need
    for our MVP "minimal dict" intake format.

    We keep this function small and explicit so it's easy for a beginner
    to understand and so our downstream normalization stays stable.
    """
    # NCBI requires a real contact email for automated access.
    Entrez.email = email

    # Fetch the GenBank flatfile for this accession.
    with Entrez.efetch(db="nuccore", id=accession, rettype="gb", retmode="text") as handle:
        record = SeqIO.read(handle, "genbank")

    # GenBank stores most useful metadata on the "source" feature.
    # It is typically the first feature in the record.
    source_feature = record.features[0] if record.features else None
    qualifiers = source_feature.qualifiers if source_feature else {}

    # Common qualifier keys we care about:
    # - collection_date
    # - country (often formatted like "USA: California")
    # - host
    collection_date = (qualifiers.get("collection_date") or [None])[0]
    location = (qualifiers.get("country") or [None])[0]
    host = (qualifiers.get("host") or [None])[0]

    organism = record.annotations.get("organism", "") or ""
    sequence = str(record.seq)

    # Return a minimal dict that our normalizer already knows how to handle.
    return {
        "accession": accession,
        "organism": organism,
        "collection_date": collection_date,
        "location": location,
        "host": host,
        "sequence": sequence,
    }

def normalize_genbank_minimal(raw: Dict[str, Any]) -> CanonicalGenomeRecord:
    """
    Convert a minimal GenBank-like dict into our CanonicalGenomeRecord.

    Expected keys (MVP):
      - accession
      - organism
      - collection_date (ISO string or partial)
      - location ("Country: Region" or "Country")
      - host (optional)
      - sequence (string, optional)
    """
    accession = str(raw.get("accession", "")).strip()
    organism = str(raw.get("organism", "")).strip()

    # Required fields for our canonical contract.
    if not accession:
        raise ValueError("accession is required")
    if not organism:
        raise ValueError("organism is required")


    collection_date = parse_collection_date(raw.get("collection_date"))
    country, region = parse_location(raw.get("location"))

    host_raw = raw.get("host")
    host = str(host_raw).strip() if host_raw else None

    seq = raw.get("sequence") or ""
    seq_str = str(seq).strip()
    sequence_length = len(seq_str)

    return CanonicalGenomeRecord(
        accession=accession,
        organism=organism,
        collection_date=collection_date,
        country=country,
        region=region,
        host=host,
        sequence_length=sequence_length,
        source="genbank",
    )


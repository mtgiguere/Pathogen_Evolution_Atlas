"""
genbank.py
"""

import logging
import threading
import time
from collections.abc import Iterable
from datetime import date
from typing import Any

from Bio import Entrez, SeqIO

from .models import CanonicalGenomeRecord

# NCBI guideline: no more than ~3 requests/second
_MIN_SECONDS_BETWEEN_REQUESTS = 1.0 / 3.0
_LAST_REQUEST_TS = None
_rate_lock = threading.Lock()

logger = logging.getLogger(__name__)


def _first_qualifier(qualifiers: dict[str, Any], key: str) -> str | None:
    v = qualifiers.get(key)
    if not v:
        return None
    if isinstance(v, list):
        return v[0] if v else None
    return str(v)


def _find_source_feature(record):
    for f in getattr(record, "features", []) or []:
        if getattr(f, "type", None) == "source":
            return f
    return None


def _extract_biosample(qualifiers: dict[str, Any]) -> str | None:
    dbx = qualifiers.get("db_xref") or []
    if not isinstance(dbx, list):
        dbx = [dbx]
    for item in dbx:
        s = str(item)
        if s.startswith("BioSample:"):
            return s.split(":", 1)[1].strip() or None
    return None


def parse_collection_date(raw: str | None) -> date | None:
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
    MONTHS = {
        "jan": 1,
        "feb": 2,
        "mar": 3,
        "apr": 4,
        "may": 5,
        "jun": 6,
        "jul": 7,
        "aug": 8,
        "sep": 9,
        "oct": 10,
        "nov": 11,
        "dec": 12,
    }

    if not raw:
        return None

    s = str(raw).strip()
    if s.lower() in {"unknown", "na", "n/a", "none"}:
        return None

    parts = s.split("-")

    if len(parts) == 3:
        y, m, d = (p.strip() for p in parts)

        try:
            m_key = m[:3].lower()
            if m_key in MONTHS and not m.isdigit():
                return date(int(y), MONTHS[m_key], int(d))

            return date(int(y), int(m), int(d))
        except ValueError:
            return None

    if len(parts) == 2:
        a, b = parts[0].strip(), parts[1].strip()

        if a[:3].lower() in MONTHS and b.isdigit():
            return date(int(b), MONTHS[a[:3].lower()], 1)

        y, m = a, b
        return date(int(y), int(m), 1)

    if len(parts) == 1:
        y = parts[0]
        return date(int(y), 1, 1)

    return None


def parse_location(raw: str | None) -> tuple[str | None, str | None]:
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


def fetch_genbank_minimal(accession: str, email: str) -> dict[str, Any]:
    """
    Fetch a single GenBank record from NCBI and extract only the fields we need
    for our MVP "minimal dict" intake format.
    """
    Entrez.email = email

    logger.debug("Fetching GenBank record %s", accession)
    _rate_limit()
    with Entrez.efetch(db="nuccore", id=accession, rettype="gb", retmode="text") as handle:
        record = SeqIO.read(handle, "genbank")

    source_feature = _find_source_feature(record)
    qualifiers = source_feature.qualifiers if source_feature else {}

    collection_date = (qualifiers.get("collection_date") or [None])[0]
    location = _first_qualifier(qualifiers, "country") or _first_qualifier(
        qualifiers, "geo_loc_name"
    )
    host = (qualifiers.get("host") or [None])[0]
    lat_lon = _first_qualifier(qualifiers, "lat_lon")
    biosample = _extract_biosample(qualifiers)
    organism = record.annotations.get("organism", "") or ""
    sequence = str(record.seq)

    return {
        "accession": accession,
        "organism": organism,
        "collection_date": collection_date,
        "location": location,
        "host": host,
        "sequence": sequence,
        "lat_lon": lat_lon,
        "biosample": biosample,
    }


def normalize_genbank_minimal(raw: dict[str, Any]) -> CanonicalGenomeRecord:
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
        sequence=seq_str,
        source="genbank",
    )


def _rate_limit() -> None:
    """
    Enforce a minimum delay between NCBI requests so we stay under
    the recommended rate limit (~3 requests per second).

    A threading.Lock ensures this is safe under concurrent access.
    """
    global _LAST_REQUEST_TS

    with _rate_lock:
        now = time.monotonic()
        if _LAST_REQUEST_TS is not None:
            elapsed = now - _LAST_REQUEST_TS
            remaining = _MIN_SECONDS_BETWEEN_REQUESTS - elapsed
            if remaining > 0:
                logger.debug("Rate-limit sleep: %.3f s", remaining)
                time.sleep(remaining)
        _LAST_REQUEST_TS = time.monotonic()


def fetch_many_genbank_minimal(accessions: Iterable[str], email: str) -> list[dict]:
    """
    Fetch many GenBank records and return a list of minimal dicts.

    Note: This function intentionally keeps behavior simple:
    - preserves input order
    - uses the single-record fetch function internally
    """
    return [fetch_genbank_minimal(a, email=email) for a in accessions]


def normalize_many_genbank_minimal(
    raw_records: Iterable[dict[str, Any]],
) -> list[CanonicalGenomeRecord]:
    """
    Normalize many minimal GenBank-like dicts into CanonicalGenomeRecord objects.
    Preserves input order.
    """
    return [normalize_genbank_minimal(r) for r in raw_records]


def fetch_and_normalize_many(accessions: Iterable[str], email: str):
    """
    Convenience helper: fetch many GenBank records and immediately normalize them
    into CanonicalGenomeRecord objects.
    """
    raws = fetch_many_genbank_minimal(accessions=accessions, email=email)
    return normalize_many_genbank_minimal(raws)


"""
geography_enrichment.py
Geographic normalization and resolution logic.

This module is responsible for turning raw location metadata
into map-ready coordinates or regions.clear

"""
from dataclasses import replace
from typing import Optional, Iterable, List
from Bio import Entrez

from .models import CanonicalGenomeRecord
from .genbank import parse_location


def enrich_location_from_biosample(
    record: CanonicalGenomeRecord,
    email: str,
) -> CanonicalGenomeRecord:
    if record.country:
        return record

    geo = _fetch_biosample_geo_loc(record.accession, email)
    if not geo:
        return record

    country, region = parse_location(geo)
    return replace(record, country=country, region=region)



def _fetch_biosample_geo_loc(accession: str, email: str) -> Optional[str]:
    Entrez.email = email

    # nuccore → biosample
    with Entrez.elink(dbfrom="nuccore", db="biosample", id=accession) as h:
        link = Entrez.read(h)

    ids = link[0].get("LinkSetDb", [])
    if not ids:
        return None

    biosample_ids = ids[0]["Link"]
    if not biosample_ids:
        return None

    bs_id = biosample_ids[0]["Id"]

    with Entrez.efetch(db="biosample", id=bs_id, retmode="xml") as h:
        try:
            doc = Entrez.read(h, validate=False)
        except Exception:
            # BioSample XML is often missing a DTD
            # Enrichment is best-effort; never crash
            return None

    attrs = doc[0].get("Attributes", [])
    for a in attrs:
        if a.get("attribute_name") == "geo_loc_name":
            return a.get("content")

    return None

def enrich_many_locations(
    records: Iterable[CanonicalGenomeRecord],
    email: str,
) -> List[CanonicalGenomeRecord]:
    """
    Enrich a collection of CanonicalGenomeRecords with geographic metadata.

    Only records missing country information will be enriched.
    Preserves input order.
    """
    enriched: List[CanonicalGenomeRecord] = []

    for r in records:
        if r.country:
            enriched.append(r)
        else:
            enriched.append(enrich_location_from_biosample(r, email=email))

    return enriched
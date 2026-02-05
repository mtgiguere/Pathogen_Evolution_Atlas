
"""
geography_enrichment.py
Geographic normalization and resolution logic.

This module is responsible for turning raw location metadata
into map-ready coordinates or regions.clear

"""
from dataclasses import replace
from typing import Optional
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
        doc = Entrez.read(h)

    attrs = doc[0].get("Attributes", [])
    for a in attrs:
        if a.get("attribute_name") == "geo_loc_name":
            return a.get("content")

    return None

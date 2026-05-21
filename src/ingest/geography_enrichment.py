from __future__ import annotations

from collections.abc import Iterable
from dataclasses import replace
from xml.etree import ElementTree as ET

from Bio import Entrez

from .genbank import parse_location
from .models import CanonicalGenomeRecord


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


def _fetch_biosample_geo_loc(accession: str, email: str) -> str | None:
    """
    Try to enrich geography via:
      nuccore accession -> elink to biosample -> efetch biosample XML -> find geo_loc_name
    Uses ElementTree because some BioSample XML responses do not include a DTD/schema,
    which makes Bio.Entrez.read() raise.
    """
    Entrez.email = email

    # 1) nuccore → biosample ids
    with Entrez.elink(dbfrom="nuccore", db="biosample", id=accession) as h:
        link = Entrez.read(h)

    linksetdb = link[0].get("LinkSetDb", []) if link else []
    if not linksetdb:
        return None

    links = linksetdb[0].get("Link", [])
    if not links:
        return None

    bs_id = links[0].get("Id")
    if not bs_id:
        return None

    # 2) biosample efetch -> raw XML -> ElementTree parse
    with Entrez.efetch(db="biosample", id=bs_id, retmode="xml") as h:
        xml_bytes = h.read()

    if not xml_bytes:
        return None

    # Entrez handle may return bytes; ET wants bytes or str
    try:
        root = ET.fromstring(xml_bytes)
    except ET.ParseError:
        # Sometimes NCBI returns an HTML error page or truncated XML
        return None

    # 3) Find any attribute named geo_loc_name, return its text/content
    #
    # BioSample XML commonly looks like:
    # <BioSampleSet>
    #   <BioSample ...>
    #     <Attributes>
    #       <Attribute attribute_name="geo_loc_name">USA: California</Attribute>
    #
    for attr in root.findall(".//Attribute"):
        if attr.get("attribute_name") == "geo_loc_name":
            val = (attr.text or "").strip()
            return val or None

    return None


def enrich_many_locations(
    records: Iterable[CanonicalGenomeRecord],
    email: str,
) -> list[CanonicalGenomeRecord]:
    enriched: list[CanonicalGenomeRecord] = []
    for r in records:
        if r.country:
            enriched.append(r)
        else:
            enriched.append(enrich_location_from_biosample(r, email=email))
    return enriched

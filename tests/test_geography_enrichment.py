from ingest.geography_enrichment import enrich_location_from_biosample
from ingest.models import CanonicalGenomeRecord


def test_enrich_location_noop_if_country_present(mocker):
    rec = CanonicalGenomeRecord(
        accession="X",
        organism="Virus",
        collection_date=None,
        country="USA",
        region="Illinois",
        host=None,
        sequence_length=0,
        sequence="",
        source="genbank",
    )

    out = enrich_location_from_biosample(rec, email="x@y.z")
    assert out.country == "USA"
    assert out.region == "Illinois"


def test_enrich_location_uses_biosample_when_missing(mocker):
    mocker.patch(
        "ingest.geography_enrichment._fetch_biosample_geo_loc",
        return_value="USA: California",
    )

    rec = CanonicalGenomeRecord(
        accession="MN908947.3",
        organism="Virus",
        collection_date=None,
        country=None,
        region=None,
        host=None,
        sequence_length=0,
        sequence="",
        source="genbank",
    )

    out = enrich_location_from_biosample(rec, email="x@y.z")
    assert out.country == "USA"
    assert out.region == "California"


def test_enrich_location_no_biosample_data(mocker):
    mocker.patch(
        "ingest.geography_enrichment._fetch_biosample_geo_loc",
        return_value=None,
    )

    rec = CanonicalGenomeRecord(
        accession="X",
        organism="Virus",
        collection_date=None,
        country=None,
        region=None,
        host=None,
        sequence_length=0,
        sequence="",
        source="genbank",
    )

    out = enrich_location_from_biosample(rec, email="x@y.z")
    assert out.country is None
    assert out.region is None

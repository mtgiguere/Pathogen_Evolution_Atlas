from src.ingest.geography_enrichment import enrich_location_from_biosample, enrich_many_locations
from src.ingest.models import CanonicalGenomeRecord


def _rec(accession="X", country=None, region=None):
    return CanonicalGenomeRecord(
        accession=accession,
        organism="Virus",
        collection_date=None,
        country=country,
        region=region,
        host=None,
        sequence_length=0,
        sequence="",
        source="genbank",
    )


def test_enrich_location_noop_if_country_present(mocker):
    rec = _rec(country="USA", region="Illinois")
    out = enrich_location_from_biosample(rec, email="x@y.z")
    assert out.country == "USA"
    assert out.region == "Illinois"


def test_enrich_location_uses_biosample_when_missing(mocker):
    mocker.patch(
        "src.ingest.geography_enrichment._fetch_biosample_geo_loc",
        return_value="USA: California",
    )
    out = enrich_location_from_biosample(_rec(accession="MN908947.3"), email="x@y.z")
    assert out.country == "USA"
    assert out.region == "California"


def test_enrich_location_no_biosample_data(mocker):
    mocker.patch(
        "src.ingest.geography_enrichment._fetch_biosample_geo_loc",
        return_value=None,
    )
    out = enrich_location_from_biosample(_rec(), email="x@y.z")
    assert out.country is None
    assert out.region is None


def test_enrich_many_locations_mixed_records(mocker):
    mocker.patch(
        "src.ingest.geography_enrichment._fetch_biosample_geo_loc",
        return_value="USA: California",
    )
    out = enrich_many_locations(
        [_rec("A", country="USA", region="Illinois"), _rec("B")],
        email="x@y.z",
    )
    assert out[0].country == "USA"
    assert out[0].region == "Illinois"
    assert out[1].country == "USA"
    assert out[1].region == "California"


# ── rate limiting ─────────────────────────────────────────────────────────────


def test_enrich_many_sleeps_between_api_calls(mocker):
    mocker.patch(
        "src.ingest.geography_enrichment._fetch_biosample_geo_loc",
        return_value="Germany",
    )
    mock_sleep = mocker.patch("src.ingest.geography_enrichment.time.sleep")

    recs = [_rec(f"X{i}") for i in range(3)]
    enrich_many_locations(recs, email="x@y.z", delay_seconds=0.4)

    assert mock_sleep.call_count == 3
    mock_sleep.assert_called_with(0.4)


def test_enrich_many_no_sleep_for_already_enriched(mocker):
    mocker.patch(
        "src.ingest.geography_enrichment._fetch_biosample_geo_loc",
        return_value="Germany",
    )
    mock_sleep = mocker.patch("src.ingest.geography_enrichment.time.sleep")

    recs = [
        _rec("A", country="USA"),
        _rec("B", country="France"),
        _rec("C"),  # only this one needs API call
    ]
    enrich_many_locations(recs, email="x@y.z")

    assert mock_sleep.call_count == 1


def test_enrich_many_accepts_generator(mocker):
    mocker.patch(
        "src.ingest.geography_enrichment._fetch_biosample_geo_loc",
        return_value="Japan",
    )
    mocker.patch("src.ingest.geography_enrichment.time.sleep")

    result = enrich_many_locations((_rec(f"G{i}") for i in range(3)), email="x@y.z")
    assert len(result) == 3
    assert all(r.country == "Japan" for r in result)

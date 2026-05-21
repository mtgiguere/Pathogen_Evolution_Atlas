"""
Growth rate estimation tests — TDD, written before implementation.
All tests expected to FAIL until src/ingest/growth.py is implemented.
"""

from __future__ import annotations

import math

import pandas as pd
import pytest

from src.ingest.growth import aggregate_by_week, estimate_growth_rates

# ── aggregate_by_week ─────────────────────────────────────────────────────────


def _df(*rows: tuple[str, str]) -> pd.DataFrame:
    """Build a minimal DataFrame from (lineage, collection_date) tuples."""
    return pd.DataFrame(rows, columns=["lineage", "collection_date"])


def test_aggregate_by_week_empty():
    df = _df()
    result = aggregate_by_week(df)
    assert result.empty
    assert list(result.columns) == ["lineage", "week", "count", "total", "frequency"]


def test_aggregate_by_week_single_lineage():
    df = _df(
        ("Delta", "2021-10-04"),   # 2021-W40
        ("Delta", "2021-10-05"),   # 2021-W40
        ("Delta", "2021-10-11"),   # 2021-W41
    )
    result = aggregate_by_week(df)
    assert len(result) == 2
    w40 = result[result["week"] == "2021-W40"].iloc[0]
    assert w40["count"] == 2
    assert w40["total"] == 2
    assert w40["frequency"] == pytest.approx(1.0)
    w41 = result[result["week"] == "2021-W41"].iloc[0]
    assert w41["count"] == 1


def test_aggregate_by_week_multiple_lineages():
    df = _df(
        ("Delta", "2021-10-04"),   # 2021-W40
        ("BA.2",  "2021-10-04"),   # 2021-W40
        ("BA.2",  "2021-10-05"),   # 2021-W40
    )
    result = aggregate_by_week(df)
    assert len(result) == 2  # Delta×W40, BA.2×W40
    total_delta = result[result["lineage"] == "Delta"].iloc[0]["total"]
    total_ba2   = result[result["lineage"] == "BA.2"].iloc[0]["total"]
    assert total_delta == 3   # both share the same week total
    assert total_ba2   == 3
    freq_delta = result[result["lineage"] == "Delta"].iloc[0]["frequency"]
    assert freq_delta == pytest.approx(1 / 3)


def test_aggregate_by_week_drops_missing_dates():
    df = _df(
        ("Delta", "2021-10-04"),
        ("Delta", None),
        ("Delta", "bad-date"),
    )
    result = aggregate_by_week(df)
    # Only the one valid-date row should survive
    assert len(result) == 1
    assert result.iloc[0]["count"] == 1


def test_aggregate_by_week_week_format():
    df = _df(("Delta", "2022-01-03"))  # Monday of 2022-W01
    result = aggregate_by_week(df)
    assert result.iloc[0]["week"] == "2022-W01"


# ── estimate_growth_rates ─────────────────────────────────────────────────────


def _weekly(lineage: str, weeks: list[str], counts: list[int]) -> pd.DataFrame:
    """Build a weekly-count DataFrame for a single lineage."""
    rows = []
    for week, count in zip(weeks, counts, strict=True):
        rows.append({"lineage": lineage, "week": week, "count": count})
    df = pd.DataFrame(rows)
    df["total"] = df["count"]
    df["frequency"] = 1.0
    return df


WEEKS_4 = ["2022-W01", "2022-W02", "2022-W03", "2022-W04"]


def test_estimate_growth_rates_empty():
    result = estimate_growth_rates(pd.DataFrame(columns=["lineage", "week", "count"]))
    assert result.empty


def test_estimate_growth_rates_growing():
    # Perfect doubling each week: counts [1, 2, 4, 8] → slope ≈ ln(2) per week
    df = _weekly("Delta", WEEKS_4, [1, 2, 4, 8])
    result = estimate_growth_rates(df)
    assert len(result) == 1
    row = result.iloc[0]
    assert row["lineage"] == "Delta"
    assert row["growth_rate"] == pytest.approx(math.log(2), abs=0.05)
    assert row["trend"] == "Growing"
    assert row["n_timepoints"] == 4
    # r_squared close to 1 for perfect exponential
    assert row["r_squared"] > 0.99


def test_estimate_growth_rates_declining():
    df = _weekly("Delta", WEEKS_4, [8, 4, 2, 1])
    result = estimate_growth_rates(df)
    row = result.iloc[0]
    assert row["growth_rate"] == pytest.approx(-math.log(2), abs=0.05)
    assert row["trend"] == "Declining"


def test_estimate_growth_rates_stable():
    df = _weekly("BA.2", WEEKS_4, [10, 10, 10, 10])
    result = estimate_growth_rates(df)
    row = result.iloc[0]
    assert row["growth_rate"] == pytest.approx(0.0, abs=1e-10)
    assert row["trend"] == "Stable"


def test_estimate_growth_rates_insufficient_data():
    # Only 2 timepoints — below min_timepoints=3 default
    df = _weekly("BA.5", ["2022-W01", "2022-W02"], [3, 6])
    result = estimate_growth_rates(df, min_timepoints=3)
    row = result.iloc[0]
    assert row["trend"] == "Insufficient data"
    assert math.isnan(row["growth_rate"])


def test_estimate_growth_rates_multiple_variants():
    growing = _weekly("Delta", WEEKS_4, [1, 2, 4, 8])
    declining = _weekly("BA.2",  WEEKS_4, [8, 4, 2, 1])
    df = pd.concat([growing, declining], ignore_index=True)
    result = estimate_growth_rates(df)
    assert len(result) == 2
    trends = dict(zip(result["lineage"], result["trend"], strict=True))
    assert trends["Delta"] == "Growing"
    assert trends["BA.2"] == "Declining"


def test_doubling_time_days_growing():
    # growth_rate = ln(2)/week → doubling_time = 7 days
    df = _weekly("Delta", WEEKS_4, [1, 2, 4, 8])
    result = estimate_growth_rates(df)
    row = result.iloc[0]
    assert row["doubling_time_days"] == pytest.approx(7.0, abs=0.5)


def test_doubling_time_days_nan_for_declining():
    df = _weekly("Delta", WEEKS_4, [8, 4, 2, 1])
    result = estimate_growth_rates(df)
    row = result.iloc[0]
    # Declining variant: doubling time is meaningless, should be NaN
    assert math.isnan(row["doubling_time_days"])


def test_r_squared_perfect_fit():
    df = _weekly("Delta", WEEKS_4, [1, 2, 4, 8])
    result = estimate_growth_rates(df)
    assert result.iloc[0]["r_squared"] == pytest.approx(1.0, abs=1e-6)


def test_r_squared_noisy_data():
    # Noisy but generally growing — r_squared < 1 but positive growth
    df = _weekly("BA.5", WEEKS_4, [2, 5, 3, 9])
    result = estimate_growth_rates(df)
    row = result.iloc[0]
    assert 0.0 <= row["r_squared"] <= 1.0
    assert row["growth_rate"] > 0  # net positive trend


def test_estimate_growth_rates_output_columns():
    df = _weekly("Delta", WEEKS_4, [1, 2, 4, 8])
    result = estimate_growth_rates(df)
    assert set(result.columns) >= {
        "lineage", "growth_rate", "doubling_time_days", "r_squared", "n_timepoints", "trend"
    }

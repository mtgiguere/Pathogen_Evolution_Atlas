"""
forecast.py tests — TDD, written before implementation.
All tests expected to FAIL until src/ingest/forecast.py is implemented.
"""

from __future__ import annotations

import math

import pandas as pd

from src.ingest.forecast import forecast_variant_frequencies

# ── helpers ───────────────────────────────────────────────────────────────────


def _weekly(rows: list[tuple[str, str, int, int]]) -> pd.DataFrame:
    """Build a weekly_df from (lineage, week, count, total) tuples."""
    df = pd.DataFrame(rows, columns=["lineage", "week", "count", "total"])
    df["frequency"] = df["count"] / df["total"]
    return df


def _rates(rows: list[tuple]) -> pd.DataFrame:
    """Build a growth_rates_df from (lineage, growth_rate, trend) tuples."""
    return pd.DataFrame(rows, columns=["lineage", "growth_rate", "trend"])


# ── empty / degenerate inputs ─────────────────────────────────────────────────


def test_empty_weekly_returns_empty():
    weekly = pd.DataFrame(columns=["lineage", "week", "count", "total", "frequency"])
    rates = _rates([("Alpha", 0.1, "Growing")])
    result = forecast_variant_frequencies(weekly, rates, n_weeks=4)
    assert result.empty


def test_empty_rates_returns_empty():
    weekly = _weekly([("Alpha", "2021-W40", 10, 10)])
    rates = pd.DataFrame(columns=["lineage", "growth_rate", "trend"])
    result = forecast_variant_frequencies(weekly, rates, n_weeks=4)
    assert result.empty


def test_n_weeks_zero_returns_empty():
    weekly = _weekly([("Alpha", "2021-W40", 10, 10)])
    rates = _rates([("Alpha", 0.1, "Growing")])
    result = forecast_variant_frequencies(weekly, rates, n_weeks=0)
    assert result.empty


# ── output schema ─────────────────────────────────────────────────────────────


def test_output_columns():
    weekly = _weekly([("Alpha", "2021-W40", 10, 10)])
    rates = _rates([("Alpha", 0.1, "Growing")])
    result = forecast_variant_frequencies(weekly, rates, n_weeks=2)
    expected_cols = {"lineage", "week", "projected_count", "projected_frequency"}
    assert expected_cols.issubset(set(result.columns))


def test_output_weeks_are_future():
    weekly = _weekly([("Alpha", "2021-W40", 10, 10)])
    rates = _rates([("Alpha", 0.1, "Growing")])
    result = forecast_variant_frequencies(weekly, rates, n_weeks=3)
    assert set(result["week"].unique()) == {"2021-W41", "2021-W42", "2021-W43"}


def test_output_row_count():
    weekly = _weekly(
        [
            ("Alpha", "2021-W40", 8, 10),
            ("Delta", "2021-W40", 2, 10),
        ]
    )
    rates = _rates(
        [
            ("Alpha", 0.1, "Growing"),
            ("Delta", -0.1, "Declining"),
        ]
    )
    result = forecast_variant_frequencies(weekly, rates, n_weeks=4)
    # 2 lineages × 4 weeks
    assert len(result) == 8


# ── growth direction ──────────────────────────────────────────────────────────


def test_growing_variant_count_increases():
    weekly = _weekly([("Alpha", "2021-W40", 10, 10)])
    rates = _rates([("Alpha", 0.2, "Growing")])
    result = forecast_variant_frequencies(weekly, rates, n_weeks=4)
    counts = result.sort_values("week")["projected_count"].tolist()
    assert counts == sorted(counts), "Projected counts should increase for a growing variant"


def test_declining_variant_count_decreases():
    weekly = _weekly([("Alpha", "2021-W40", 10, 10)])
    rates = _rates([("Alpha", -0.2, "Declining")])
    result = forecast_variant_frequencies(weekly, rates, n_weeks=4)
    counts = result.sort_values("week")["projected_count"].tolist()
    assert counts == sorted(counts, reverse=True), (
        "Projected counts should decrease for a declining variant"
    )


def test_stable_variant_count_roughly_flat():
    weekly = _weekly([("Alpha", "2021-W40", 100, 100)])
    rates = _rates([("Alpha", 0.0, "Stable")])
    result = forecast_variant_frequencies(weekly, rates, n_weeks=4)
    counts = result["projected_count"].tolist()
    assert all(math.isclose(c, 100.0, rel_tol=1e-6) for c in counts)


# ── frequency properties ──────────────────────────────────────────────────────


def test_frequencies_sum_to_one_per_week():
    weekly = _weekly(
        [
            ("Alpha", "2021-W40", 7, 10),
            ("Delta", "2021-W40", 3, 10),
        ]
    )
    rates = _rates(
        [
            ("Alpha", 0.15, "Growing"),
            ("Delta", -0.1, "Declining"),
        ]
    )
    result = forecast_variant_frequencies(weekly, rates, n_weeks=4)
    for week, grp in result.groupby("week"):
        total_freq = grp["projected_frequency"].sum()
        assert math.isclose(total_freq, 1.0, abs_tol=1e-9), (
            f"Frequencies don't sum to 1.0 in week {week}: {total_freq}"
        )


def test_frequency_bounds():
    weekly = _weekly(
        [
            ("Alpha", "2021-W40", 9, 10),
            ("Delta", "2021-W40", 1, 10),
        ]
    )
    rates = _rates(
        [
            ("Alpha", 0.3, "Growing"),
            ("Delta", -0.3, "Declining"),
        ]
    )
    result = forecast_variant_frequencies(weekly, rates, n_weeks=6)
    assert (result["projected_frequency"] >= 0).all()
    assert (result["projected_frequency"] <= 1).all()


# ── growth formula ────────────────────────────────────────────────────────────


def test_projected_count_matches_exponential_formula():
    last_count = 50
    rate = 0.1
    weekly = _weekly([("Alpha", "2021-W40", last_count, last_count)])
    rates = _rates([("Alpha", rate, "Growing")])
    result = forecast_variant_frequencies(weekly, rates, n_weeks=3)
    result = result.sort_values("week").reset_index(drop=True)
    for t, row in result.iterrows():
        expected = last_count * math.exp(rate * (t + 1))
        assert math.isclose(row["projected_count"], expected, rel_tol=1e-9)


# ── lineage filtering ─────────────────────────────────────────────────────────


def test_lineage_missing_from_rates_is_excluded():
    weekly = _weekly(
        [
            ("Alpha", "2021-W40", 8, 10),
            ("Delta", "2021-W40", 2, 10),
        ]
    )
    # Only Alpha has a rate
    rates = _rates([("Alpha", 0.1, "Growing")])
    result = forecast_variant_frequencies(weekly, rates, n_weeks=2)
    assert set(result["lineage"].unique()) == {"Alpha"}


def test_lineage_with_nan_growth_rate_is_excluded():
    weekly = _weekly(
        [
            ("Alpha", "2021-W40", 8, 10),
            ("Delta", "2021-W40", 2, 10),
        ]
    )
    rates = pd.DataFrame(
        [
            {"lineage": "Alpha", "growth_rate": 0.1, "trend": "Growing"},
            {"lineage": "Delta", "growth_rate": float("nan"), "trend": "Insufficient data"},
        ]
    )
    result = forecast_variant_frequencies(weekly, rates, n_weeks=2)
    assert set(result["lineage"].unique()) == {"Alpha"}


# ── multi-week baseline ───────────────────────────────────────────────────────


def test_uses_last_observed_week_as_baseline():
    """Forecast should anchor to the most recent week, not an earlier one."""
    weekly = _weekly(
        [
            ("Alpha", "2021-W38", 5, 10),
            ("Alpha", "2021-W39", 8, 10),
            ("Alpha", "2021-W40", 12, 10),
        ]
    )
    rates = _rates([("Alpha", 0.1, "Growing")])
    result = forecast_variant_frequencies(weekly, rates, n_weeks=1)
    # Only one forecast week: 2021-W41
    assert list(result["week"]) == ["2021-W41"]
    expected_count = 12 * math.exp(0.1 * 1)
    assert math.isclose(result.iloc[0]["projected_count"], expected_count, rel_tol=1e-9)

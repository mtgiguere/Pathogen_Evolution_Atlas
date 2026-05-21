"""
forecast.py — Project variant frequencies forward from estimated growth rates.

Uses the log-linear growth model already fit by growth.py:
  projected_count(t) = last_count * exp(growth_rate * t)
where t is weeks from the last observed week.

Frequencies are computed by normalising projected counts across all forecastable
lineages within each projected week, so they always sum to 1.
"""

from __future__ import annotations

import math
from datetime import datetime, timedelta

import pandas as pd


def _advance_week(week_str: str, n: int) -> str:
    """Return the ISO-week string n weeks after week_str (e.g. '2021-W40' + 2 → '2021-W42')."""
    year, w = week_str.split("-W")
    monday = datetime.strptime(f"{year}-W{int(w):02d}-1", "%G-W%V-%u")
    future = monday + timedelta(weeks=n)
    iso = future.isocalendar()
    return f"{iso.year}-W{iso.week:02d}"


def forecast_variant_frequencies(
    weekly_df: pd.DataFrame,
    growth_rates_df: pd.DataFrame,
    *,
    n_weeks: int = 4,
) -> pd.DataFrame:
    """
    Project variant frequencies forward n_weeks from the last observed week.

    Parameters
    ----------
    weekly_df:
        Output of aggregate_by_week — columns: lineage, week, count, total, frequency.
    growth_rates_df:
        Output of estimate_growth_rates — columns: lineage, growth_rate, trend.
    n_weeks:
        Number of future weeks to project.

    Returns
    -------
    DataFrame with columns: lineage, week, projected_count, projected_frequency.
    One row per (lineage, future_week). Lineages with NaN or missing growth rates
    are excluded. Frequencies sum to 1.0 within each projected week.
    """
    empty = pd.DataFrame(columns=["lineage", "week", "projected_count", "projected_frequency"])

    if weekly_df.empty or growth_rates_df.empty or n_weeks <= 0:
        return empty

    # Drop lineages with no usable growth rate
    valid_rates = growth_rates_df.dropna(subset=["growth_rate"])
    if valid_rates.empty:
        return empty

    # Find last observed count per lineage (anchor for projection)
    last_counts: dict[str, float] = {}
    for lineage, grp in weekly_df.groupby("lineage"):
        last_row = grp.sort_values("week").iloc[-1]
        last_counts[lineage] = float(last_row["count"])

    # Global last week — project forward from here
    last_week = weekly_df["week"].max()

    rows: list[dict] = []
    for _, rate_row in valid_rates.iterrows():
        lineage = rate_row["lineage"]
        if lineage not in last_counts:
            continue
        rate = float(rate_row["growth_rate"])
        base = last_counts[lineage]

        for t in range(1, n_weeks + 1):
            rows.append(
                {
                    "lineage": lineage,
                    "week": _advance_week(last_week, t),
                    "projected_count": base * math.exp(rate * t),
                }
            )

    if not rows:
        return empty

    df = pd.DataFrame(rows)

    # Normalise counts to frequencies within each projected week
    week_totals = df.groupby("week")["projected_count"].transform("sum")
    df["projected_frequency"] = df["projected_count"] / week_totals

    return df.reset_index(drop=True)

"""
growth.py — Per-variant growth rate estimation.

Bins sequences by ISO week and lineage, then fits a log-linear model
(ln(count) ~ week_index) to estimate exponential growth rates.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
import pandas as pd

# Thresholds for labelling trend direction (natural-log scale, per week).
# ln(1.05) ≈ 0.049 → 5% weekly growth is the minimum "Growing" signal.
_GROWING_THRESHOLD = 0.05
_DECLINING_THRESHOLD = -0.05


@dataclass
class GrowthRate:
    lineage: str
    growth_rate: float  # ln-scale slope per week; NaN if insufficient data
    doubling_time_days: float  # positive only for growing variants; NaN otherwise
    r_squared: float  # goodness-of-fit; NaN if insufficient data
    n_timepoints: int  # number of weekly bins used in the fit
    trend: str  # "Growing" | "Declining" | "Stable" | "Insufficient data"


def aggregate_by_week(
    df: pd.DataFrame,
    *,
    lineage_col: str = "lineage",
    date_col: str = "collection_date",
) -> pd.DataFrame:
    """
    Group sequences by ISO week and lineage.

    Returns a DataFrame with columns:
      lineage, week (YYYY-WWW format), count, total (all lineages that week), frequency
    """
    empty = pd.DataFrame(columns=["lineage", "week", "count", "total", "frequency"])

    if df.empty or lineage_col not in df.columns or date_col not in df.columns:
        return empty

    work = df[[lineage_col, date_col]].copy()
    work.columns = ["lineage", "date"]

    # Parse dates, coercing bad/missing values to NaT
    work["date"] = pd.to_datetime(work["date"], errors="coerce")
    work = work.dropna(subset=["date"])

    if work.empty:
        return empty

    # ISO week string "YYYY-Www"
    work["week"] = work["date"].apply(
        lambda d: f"{d.isocalendar().year}-W{d.isocalendar().week:02d}"
    )

    counts = work.groupby(["lineage", "week"], sort=True).size().reset_index(name="count")

    week_totals = counts.groupby("week")["count"].sum().rename("total")
    counts = counts.join(week_totals, on="week")
    counts["frequency"] = counts["count"] / counts["total"]

    return counts.reset_index(drop=True)


def estimate_growth_rates(
    weekly_df: pd.DataFrame,
    *,
    min_timepoints: int = 3,
) -> pd.DataFrame:
    """
    Fit a log-linear growth model per lineage.

    Expects a DataFrame with at least columns: lineage, week, count.
    Returns a DataFrame with columns:
      lineage, growth_rate, doubling_time_days, r_squared, n_timepoints, trend
    """
    required = {"lineage", "week", "count"}
    if weekly_df.empty or not required.issubset(weekly_df.columns):
        return pd.DataFrame(
            columns=[
                "lineage",
                "growth_rate",
                "doubling_time_days",
                "r_squared",
                "n_timepoints",
                "trend",
            ]
        )

    rows: list[dict] = []

    for lineage, grp in weekly_df.groupby("lineage"):
        grp = grp.sort_values("week")
        valid = grp[grp["count"] > 0]
        n = len(valid)

        if n < min_timepoints:
            rows.append(
                {
                    "lineage": lineage,
                    "growth_rate": float("nan"),
                    "doubling_time_days": float("nan"),
                    "r_squared": float("nan"),
                    "n_timepoints": n,
                    "trend": "Insufficient data",
                }
            )
            continue

        x = np.arange(n, dtype=float)
        y = np.log(valid["count"].to_numpy(dtype=float))

        slope, intercept = np.polyfit(x, y, 1)

        # R-squared
        y_pred = slope * x + intercept
        ss_res = float(np.sum((y - y_pred) ** 2))
        ss_tot = float(np.sum((y - y.mean()) ** 2))
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
        r2 = max(0.0, min(1.0, r2))  # clamp to [0, 1]

        if slope > _GROWING_THRESHOLD:
            trend = "Growing"
            doubling = math.log(2) / slope * 7.0  # weeks → days
        elif slope < _DECLINING_THRESHOLD:
            trend = "Declining"
            doubling = float("nan")
        else:
            trend = "Stable"
            doubling = float("nan")

        rows.append(
            {
                "lineage": lineage,
                "growth_rate": float(slope),
                "doubling_time_days": doubling,
                "r_squared": r2,
                "n_timepoints": n,
                "trend": trend,
            }
        )

    return pd.DataFrame(rows)

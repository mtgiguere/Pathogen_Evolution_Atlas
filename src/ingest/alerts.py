"""
alerts.py — Rule-based alert / notification system.

An AlertEngine evaluates a set of AlertRules against the analytics DataFrame
and dispatches fired alerts to one or more AlertChannels (log, file, webhook).

Usage
-----
    engine = AlertEngine(
        rules=[HighRiskGenomeRule(), CriticalEscapeRule(), FastGrowingVariantRule()],
        channels=[LogChannel(), FileChannel(Path("data/alerts/alerts.ndjson"))],
    )
    alerts = engine.run(df, growth_df=growth_rates_df)
"""

from __future__ import annotations

import json
import logging
import math
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import pandas as pd
import requests

logger = logging.getLogger(__name__)


# ── Alert payload ─────────────────────────────────────────────────────────────


@dataclass
class Alert:
    level: str  # "INFO" | "WARNING" | "CRITICAL"
    rule_name: str
    message: str
    context: dict[str, Any] = field(default_factory=dict)
    triggered_at: datetime = field(default_factory=lambda: datetime.now(tz=UTC))

    def to_dict(self) -> dict[str, Any]:
        return {
            "level": self.level,
            "rule_name": self.rule_name,
            "message": self.message,
            "context": self.context,
            "triggered_at": self.triggered_at.isoformat(),
        }


# ── Rules ─────────────────────────────────────────────────────────────────────


class HighRiskGenomeRule:
    """Fire a CRITICAL alert for every genome scored as High risk."""

    name = "high_risk_genome"

    def evaluate(self, df: pd.DataFrame, *, growth_df: pd.DataFrame | None = None) -> list[Alert]:
        if df.empty or "risk_level" not in df.columns:
            return []
        alerts: list[Alert] = []
        for _, row in df[df["risk_level"] == "High"].iterrows():
            acc = row.get("accession", "unknown")
            lineage = row.get("lineage", "Unknown")
            score = row.get("risk_score", "?")
            alerts.append(
                Alert(
                    level="CRITICAL",
                    rule_name=self.name,
                    message=f"High-risk genome detected: {acc} ({lineage}, score={score})",
                    context={
                        "accession": acc,
                        "lineage": lineage,
                        "risk_score": score,
                        "country": row.get("country"),
                        "collection_date": str(row.get("collection_date") or ""),
                    },
                )
            )
        return alerts


class CriticalEscapeRule:
    """Fire when a genome carries ≥ min_escape_count known escape mutations."""

    name = "critical_escape"

    def __init__(self, min_escape_count: int = 2) -> None:
        self.min_escape_count = min_escape_count

    def evaluate(self, df: pd.DataFrame, *, growth_df: pd.DataFrame | None = None) -> list[Alert]:
        if df.empty or "escape_count" not in df.columns:
            return []
        mask = df["has_critical_escape"].astype(bool) & (
            df["escape_count"] >= self.min_escape_count
        )
        alerts: list[Alert] = []
        for _, row in df[mask].iterrows():
            acc = row.get("accession", "unknown")
            lineage = row.get("lineage", "Unknown")
            n = int(row.get("escape_count", 0))
            antibodies = row.get("escape_antibodies", "")
            alerts.append(
                Alert(
                    level="CRITICAL",
                    rule_name=self.name,
                    message=(
                        f"Immune escape detected: {acc} ({lineage}) "
                        f"carries {n} known escape mutation(s)"
                    ),
                    context={
                        "accession": acc,
                        "lineage": lineage,
                        "escape_count": n,
                        "escape_antibodies": antibodies,
                    },
                )
            )
        return alerts


class FastGrowingVariantRule:
    """Fire when a variant's doubling time is below the threshold."""

    name = "fast_growing_variant"

    def __init__(self, max_doubling_time_days: float = 14.0) -> None:
        self.max_doubling_time_days = max_doubling_time_days

    def evaluate(self, df: pd.DataFrame, *, growth_df: pd.DataFrame | None = None) -> list[Alert]:
        if growth_df is None or growth_df.empty or "doubling_time_days" not in growth_df.columns:
            return []
        alerts: list[Alert] = []
        for _, row in growth_df.iterrows():
            dt = row.get("doubling_time_days")
            if dt is None or (isinstance(dt, float) and math.isnan(dt)):
                continue
            if dt < self.max_doubling_time_days:
                lineage = row.get("lineage", "Unknown")
                alerts.append(
                    Alert(
                        level="WARNING",
                        rule_name=self.name,
                        message=(f"Fast-growing variant: {lineage} doubling every {dt:.1f} days"),
                        context={
                            "lineage": lineage,
                            "doubling_time_days": dt,
                            "growth_rate": row.get("growth_rate"),
                            "trend": row.get("trend"),
                        },
                    )
                )
        return alerts


class NewVOCDetectedRule:
    """Fire once per newly observed WHO Variant of Concern not in known_vocs."""

    name = "new_voc_detected"

    def __init__(self, known_vocs: set[str]) -> None:
        self.known_vocs = set(known_vocs)

    def evaluate(self, df: pd.DataFrame, *, growth_df: pd.DataFrame | None = None) -> list[Alert]:
        if df.empty or "who_class" not in df.columns:
            return []
        voc_df = df[df["who_class"] == "VOC"]
        seen_this_batch: set[str] = set()
        alerts: list[Alert] = []
        for _, row in voc_df.iterrows():
            lineage = str(row.get("lineage", ""))
            if not lineage or lineage in self.known_vocs or lineage in seen_this_batch:
                continue
            seen_this_batch.add(lineage)
            alerts.append(
                Alert(
                    level="WARNING",
                    rule_name=self.name,
                    message=f"New WHO Variant of Concern detected: {lineage}",
                    context={
                        "lineage": lineage,
                        "who_label": row.get("who_label", ""),
                        "country": row.get("country"),
                    },
                )
            )
        return alerts


# ── Channels ──────────────────────────────────────────────────────────────────


class LogChannel:
    """Write alerts to the Python logging system."""

    _LEVEL_MAP = {"INFO": logging.INFO, "WARNING": logging.WARNING, "CRITICAL": logging.CRITICAL}

    def send(self, alert: Alert) -> None:
        level = self._LEVEL_MAP.get(alert.level, logging.WARNING)
        logger.log(level, "[%s] %s", alert.rule_name, alert.message)


class FileChannel:
    """Append alerts as NDJSON to a file."""

    def __init__(self, path: Path) -> None:
        self.path = path

    def send(self, alert: Alert) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        with self.path.open("a", encoding="utf-8") as fh:
            fh.write(json.dumps(alert.to_dict()) + "\n")


class WebhookChannel:
    """POST alert JSON to an HTTP endpoint (Slack, Teams, generic webhook)."""

    def __init__(self, url: str, timeout: float = 5.0) -> None:
        self.url = url
        self.timeout = timeout

    def send(self, alert: Alert) -> None:
        try:
            resp = requests.post(self.url, json=alert.to_dict(), timeout=self.timeout)
            if resp.status_code >= 400:
                logger.warning(
                    "Webhook %s returned HTTP %d for alert %s",
                    self.url,
                    resp.status_code,
                    alert.rule_name,
                )
        except Exception as exc:
            logger.warning("Webhook delivery failed for %s: %s", self.url, exc)


# ── Engine ────────────────────────────────────────────────────────────────────


class AlertEngine:
    def __init__(self, rules: list, channels: list) -> None:
        self.rules = rules
        self.channels = channels

    def evaluate(self, df: pd.DataFrame, *, growth_df: pd.DataFrame | None = None) -> list[Alert]:
        alerts: list[Alert] = []
        for rule in self.rules:
            alerts.extend(rule.evaluate(df, growth_df=growth_df))
        return alerts

    def dispatch(self, alerts: list[Alert]) -> None:
        for alert in alerts:
            for channel in self.channels:
                try:
                    channel.send(alert)
                except Exception:
                    logger.exception("Channel %s failed to deliver alert", channel)

    def run(self, df: pd.DataFrame, *, growth_df: pd.DataFrame | None = None) -> list[Alert]:
        alerts = self.evaluate(df, growth_df=growth_df)
        if alerts:
            logger.info("AlertEngine: %d alert(s) fired", len(alerts))
            self.dispatch(alerts)
        return alerts

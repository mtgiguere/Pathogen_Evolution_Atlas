"""
Alert / notification system tests — TDD, written before implementation.
All tests expected to FAIL until src/ingest/alerts.py is implemented.
"""

from __future__ import annotations

import json
import math
import tempfile
from datetime import datetime
from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd

from src.ingest.alerts import (
    Alert,
    AlertEngine,
    CriticalEscapeRule,
    FastGrowingVariantRule,
    FileChannel,
    HighRiskGenomeRule,
    LogChannel,
    NewVOCDetectedRule,
    WebhookChannel,
)

# ── test data helpers ─────────────────────────────────────────────────────────


def _genome_df(**overrides) -> pd.DataFrame:
    base = {
        "accession": "ACC001",
        "lineage": "B.1.617.2",
        "who_class": "VOC",
        "risk_level": "High",
        "risk_score": 9.0,
        "escape_count": 3,
        "has_critical_escape": True,
        "country": "USA",
        "collection_date": "2024-01-15",
    }
    base.update(overrides)
    return pd.DataFrame([base])


def _growth_df(lineage: str, doubling_days: float, trend: str = "Growing") -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "lineage": lineage,
                "growth_rate": math.log(2) / (doubling_days / 7),
                "doubling_time_days": doubling_days,
                "r_squared": 0.95,
                "n_timepoints": 4,
                "trend": trend,
            }
        ]
    )


# ── Alert dataclass ───────────────────────────────────────────────────────────


def test_alert_has_required_fields():
    alert = Alert(
        level="WARNING",
        rule_name="test_rule",
        message="Test alert",
        context={"accession": "ACC001"},
    )
    assert alert.level == "WARNING"
    assert alert.rule_name == "test_rule"
    assert alert.message == "Test alert"
    assert isinstance(alert.triggered_at, datetime)
    assert alert.context["accession"] == "ACC001"


# ── HighRiskGenomeRule ────────────────────────────────────────────────────────


def test_high_risk_rule_fires_on_high_risk():
    rule = HighRiskGenomeRule()
    alerts = rule.evaluate(_genome_df(risk_level="High"))
    assert len(alerts) == 1
    assert alerts[0].level == "CRITICAL"
    assert "ACC001" in alerts[0].message


def test_high_risk_rule_no_alert_for_moderate():
    rule = HighRiskGenomeRule()
    assert rule.evaluate(_genome_df(risk_level="Moderate")) == []


def test_high_risk_rule_no_alert_for_low():
    rule = HighRiskGenomeRule()
    assert rule.evaluate(_genome_df(risk_level="Low")) == []


def test_high_risk_rule_fires_for_each_high_risk_genome():
    df = pd.DataFrame(
        [
            {**_genome_df(risk_level="High", accession="A1").iloc[0]},
            {**_genome_df(risk_level="High", accession="A2").iloc[0]},
            {**_genome_df(risk_level="Low", accession="A3").iloc[0]},
        ]
    )
    rule = HighRiskGenomeRule()
    alerts = rule.evaluate(df)
    assert len(alerts) == 2


def test_high_risk_rule_no_alert_for_empty_df():
    rule = HighRiskGenomeRule()
    assert rule.evaluate(pd.DataFrame()) == []


# ── CriticalEscapeRule ────────────────────────────────────────────────────────


def test_critical_escape_rule_fires():
    rule = CriticalEscapeRule(min_escape_count=2)
    alerts = rule.evaluate(_genome_df(escape_count=3, has_critical_escape=True))
    assert len(alerts) == 1
    assert alerts[0].level in ("WARNING", "CRITICAL")


def test_critical_escape_rule_respects_threshold():
    rule = CriticalEscapeRule(min_escape_count=4)
    # escape_count=3 < 4 → no alert
    assert rule.evaluate(_genome_df(escape_count=3, has_critical_escape=True)) == []


def test_critical_escape_rule_requires_has_critical_escape_true():
    rule = CriticalEscapeRule(min_escape_count=1)
    assert rule.evaluate(_genome_df(escape_count=5, has_critical_escape=False)) == []


def test_critical_escape_rule_context_includes_antibodies():
    rule = CriticalEscapeRule(min_escape_count=1)
    df = _genome_df(
        escape_count=2, has_critical_escape=True, escape_antibodies="Bamlanivimab, REGN10933"
    )
    alerts = rule.evaluate(df)
    assert len(alerts) == 1
    assert "antibodies" in alerts[0].context or "escape_antibodies" in alerts[0].context


# ── FastGrowingVariantRule ────────────────────────────────────────────────────


def test_fast_growing_rule_fires_when_doubling_time_short():
    rule = FastGrowingVariantRule(max_doubling_time_days=14)
    growth = _growth_df("Delta", doubling_days=7)  # 7 days < 14
    alerts = rule.evaluate(pd.DataFrame(), growth_df=growth)
    assert len(alerts) == 1
    assert "Delta" in alerts[0].message


def test_fast_growing_rule_no_alert_when_slow():
    rule = FastGrowingVariantRule(max_doubling_time_days=14)
    growth = _growth_df("Delta", doubling_days=30)  # 30 > 14
    assert rule.evaluate(pd.DataFrame(), growth_df=growth) == []


def test_fast_growing_rule_no_alert_for_declining():
    rule = FastGrowingVariantRule(max_doubling_time_days=14)
    growth = _growth_df("Delta", doubling_days=float("nan"), trend="Declining")
    assert rule.evaluate(pd.DataFrame(), growth_df=growth) == []


def test_fast_growing_rule_no_growth_df_returns_empty():
    rule = FastGrowingVariantRule(max_doubling_time_days=14)
    assert rule.evaluate(_genome_df(), growth_df=None) == []


# ── NewVOCDetectedRule ────────────────────────────────────────────────────────


def test_new_voc_rule_fires_for_unseen_voc():
    rule = NewVOCDetectedRule(known_vocs=set())
    alerts = rule.evaluate(_genome_df(who_class="VOC", lineage="XBB.1.5"))
    assert len(alerts) == 1
    assert "XBB.1.5" in alerts[0].message


def test_new_voc_rule_no_alert_for_known_voc():
    rule = NewVOCDetectedRule(known_vocs={"B.1.617.2"})
    assert rule.evaluate(_genome_df(who_class="VOC", lineage="B.1.617.2")) == []


def test_new_voc_rule_no_alert_for_non_voc():
    rule = NewVOCDetectedRule(known_vocs=set())
    assert rule.evaluate(_genome_df(who_class="", lineage="Unknown")) == []


def test_new_voc_rule_deduplicates_within_df():
    # Same new VOC appears twice — one alert, not two
    df = pd.DataFrame(
        [
            _genome_df(who_class="VOC", lineage="XBB.1.5").iloc[0],
            _genome_df(who_class="VOC", lineage="XBB.1.5").iloc[0],
        ]
    )
    rule = NewVOCDetectedRule(known_vocs=set())
    alerts = rule.evaluate(df)
    assert len(alerts) == 1


# ── AlertEngine ───────────────────────────────────────────────────────────────


def test_engine_evaluates_all_rules():
    high_rule = MagicMock()
    high_rule.evaluate.return_value = [Alert(level="CRITICAL", rule_name="r1", message="m1")]
    esc_rule = MagicMock()
    esc_rule.evaluate.return_value = [Alert(level="WARNING", rule_name="r2", message="m2")]

    engine = AlertEngine(rules=[high_rule, esc_rule], channels=[])
    alerts = engine.evaluate(pd.DataFrame())
    assert len(alerts) == 2
    high_rule.evaluate.assert_called_once()
    esc_rule.evaluate.assert_called_once()


def test_engine_dispatches_to_all_channels():
    ch1 = MagicMock()
    ch2 = MagicMock()
    engine = AlertEngine(rules=[], channels=[ch1, ch2])
    alert = Alert(level="INFO", rule_name="r", message="test")
    engine.dispatch([alert])
    ch1.send.assert_called_once_with(alert)
    ch2.send.assert_called_once_with(alert)


def test_engine_run_returns_alerts():
    rule = HighRiskGenomeRule()
    engine = AlertEngine(rules=[rule], channels=[])
    alerts = engine.run(_genome_df(risk_level="High"))
    assert len(alerts) == 1


def test_engine_run_no_alerts_when_conditions_not_met():
    rule = HighRiskGenomeRule()
    engine = AlertEngine(rules=[rule], channels=[])
    alerts = engine.run(_genome_df(risk_level="Low"))
    assert alerts == []


def test_engine_passes_growth_df_to_rules():
    rule = MagicMock()
    rule.evaluate.return_value = []
    growth = _growth_df("Delta", 7)
    engine = AlertEngine(rules=[rule], channels=[])
    engine.run(pd.DataFrame(), growth_df=growth)
    _, kwargs = rule.evaluate.call_args
    assert kwargs.get("growth_df") is not None or rule.evaluate.call_args[0][1] is not None


# ── LogChannel ────────────────────────────────────────────────────────────────


def test_log_channel_logs_alert(caplog):
    import logging

    ch = LogChannel()
    alert = Alert(level="WARNING", rule_name="test", message="Watch out")
    with caplog.at_level(logging.WARNING, logger="src.ingest.alerts"):
        ch.send(alert)
    assert "Watch out" in caplog.text


# ── FileChannel ───────────────────────────────────────────────────────────────


def test_file_channel_writes_alert():
    with tempfile.NamedTemporaryFile(suffix=".ndjson", delete=False) as f:
        tmp = Path(f.name)
    ch = FileChannel(path=tmp)
    alert = Alert(
        level="CRITICAL", rule_name="rule1", message="Critical finding", context={"acc": "X"}
    )
    ch.send(alert)
    lines = tmp.read_text(encoding="utf-8").strip().splitlines()
    assert len(lines) == 1
    data = json.loads(lines[0])
    assert data["level"] == "CRITICAL"
    assert data["message"] == "Critical finding"


def test_file_channel_appends():
    with tempfile.NamedTemporaryFile(suffix=".ndjson", delete=False) as f:
        tmp = Path(f.name)
    ch = FileChannel(path=tmp)
    for i in range(3):
        ch.send(Alert(level="INFO", rule_name="r", message=f"msg{i}"))
    lines = tmp.read_text(encoding="utf-8").strip().splitlines()
    assert len(lines) == 3


def test_file_channel_creates_parent_dirs(tmp_path):
    ch = FileChannel(path=tmp_path / "sub" / "alerts.ndjson")
    ch.send(Alert(level="INFO", rule_name="r", message="hi"))
    assert (tmp_path / "sub" / "alerts.ndjson").exists()


# ── WebhookChannel ────────────────────────────────────────────────────────────


def test_webhook_channel_posts_to_url():
    ch = WebhookChannel(url="https://hooks.example.com/test")
    alert = Alert(level="WARNING", rule_name="r", message="test webhook")
    with patch("src.ingest.alerts.requests.post") as mock_post:
        mock_post.return_value = MagicMock(status_code=200)
        ch.send(alert)
    mock_post.assert_called_once()
    call_kwargs = mock_post.call_args
    assert call_kwargs[0][0] == "https://hooks.example.com/test"


def test_webhook_channel_logs_on_failure(caplog):
    import logging

    ch = WebhookChannel(url="https://hooks.example.com/test")
    alert = Alert(level="CRITICAL", rule_name="r", message="important")
    with patch("src.ingest.alerts.requests.post", side_effect=Exception("timeout")):
        with caplog.at_level(logging.WARNING, logger="src.ingest.alerts"):
            ch.send(alert)
    assert "timeout" in caplog.text or "webhook" in caplog.text.lower()

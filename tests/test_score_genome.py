import logging

from src.ingest.mutations import QCResult, qc_compare_to_reference
from src.ingest.scoring import score_genome


def qc_lenient(ref, sample, **kw):
    return qc_compare_to_reference(ref, sample, min_overlap=1)


def test_score_genome_sets_num_mutations_from_identified_mutations():
    from src.ingest.mutations import Mutation

    canonical_record = {
        "accession": "TEST0002",
        "sequence": "ACGTACGTACGT",
        "source": "genbank",
    }

    def identify(_record):
        return [
            Mutation(pos=1, ref="A", alt="G", gene="S"),
            Mutation(pos=2, ref="C", alt="T", gene="N"),
            Mutation(pos=3, ref="G", alt="A", gene="ORF1ab"),
        ]

    summary = score_genome(
        canonical_record,
        reference_sequence="A",
        qc_fn=qc_lenient,
        identify_mutations=identify,
    )

    assert summary["num_mutations"] == 3


def test_score_genome_sets_genes_affected_from_gene_mapping():
    from src.ingest.mutations import Mutation

    canonical_record = {
        "accession": "TEST0003",
        "sequence": "ACGTACGTACGT",
        "source": "genbank",
    }

    def identify(_record):
        return [
            Mutation(pos=1, ref="A", alt="G", gene="S"),
            Mutation(pos=2, ref="C", alt="T", gene="N"),
        ]

    def map_genes(_mutations):
        return ["Spike", "N"]

    summary = score_genome(
        canonical_record,
        reference_sequence="A",
        qc_fn=qc_lenient,
        identify_mutations=identify,
        map_genes=map_genes,
    )

    assert summary["genes_affected"] == ["Spike", "N"]


def test_score_genome_sets_risk_score_from_risk_model():
    from src.ingest.mutations import Mutation

    canonical_record = {
        "accession": "TEST0004",
        "sequence": "ACGTACGTACGT",
        "source": "genbank",
    }

    def identify(_record):
        return [Mutation(pos=1, ref="A", alt="G", gene="S")]

    def compute_risk(_mutations):
        return {"score": 7.5, "by_gene": {}, "level": "Low", "explanation": ""}

    summary = score_genome(
        canonical_record,
        reference_sequence="A",
        qc_fn=qc_lenient,
        identify_mutations=identify,
        compute_risk=compute_risk,
    )

    assert summary["risk_score"] == 7.5


def test_score_genome_computes_real_risk_score_from_mutations():
    from src.ingest.mutations import Mutation
    from src.ingest.risk import score_mutations

    canonical_record = {
        "accession": "TEST2000",
        "sequence": "AG",
        "source": "genbank",
    }

    def identify(_record):
        return [Mutation(pos=1, ref="A", alt="G", gene="S")]

    summary = score_genome(
        canonical_record,
        reference_sequence="A",
        qc_fn=qc_lenient,
        identify_mutations=identify,
    )

    assert summary["risk_score"] == float(score_mutations(identify(canonical_record))["score"])


def test_score_genome_includes_risk_explainability_fields():
    from src.ingest.mutations import Mutation

    canonical_record = {
        "accession": "TEST3000",
        "sequence": "AA",
        "source": "genbank",
    }

    def identify(_record):
        return [
            Mutation(pos=1, ref="A", alt="G", gene="S"),
            Mutation(pos=2, ref="C", alt="T", gene="S"),
        ]

    summary = score_genome(
        canonical_record,
        reference_sequence="A",
        qc_fn=qc_lenient,
        identify_mutations=identify,
    )

    assert summary["risk_level"] in {"Low", "Moderate", "High"}
    assert isinstance(summary["risk_by_gene"], dict)
    assert "S" in summary["risk_by_gene"]
    assert isinstance(summary["risk_explanation"], str)
    assert len(summary["risk_explanation"]) > 0


def test_score_genome_uses_injected_reference_sequence_when_record_missing_reference():
    record = {
        "accession": "PX90_TEST",
        "source": "genbank",
        "sequence": "AGGT",
    }
    reference_sequence = "ACGT"

    out = score_genome(
        record,
        reference_sequence=reference_sequence,
        qc_fn=qc_lenient,
    )

    assert out["num_mutations"] == 1
    assert out["risk_score"] == 0.0


def test_score_genome_blocks_high_N_fraction():
    ref = "A" * 100
    sample = ("A" * 40) + ("N" * 60)

    record = {
        "accession": "PX90_BAD",
        "source": "genbank",
        "sequence": sample,
    }

    out = score_genome(
        record,
        reference_sequence=ref,
        qc_fn=lambda r, s, **kw: qc_compare_to_reference(r, s, min_overlap=1, max_n_fraction=0.10),
    )

    assert out["scorable"] is False
    assert out["qc_status"] == "FAIL"
    assert "HIGH_N_FRACTION" in out["qc_reasons"]
    assert out["risk_level"] == "N/A"
    assert out["risk_score"] is None


def test_score_genome_qc_fail_logs_warning(caplog):
    """A QC failure must emit a WARNING log containing the accession."""
    with caplog.at_level(logging.WARNING):
        result = score_genome(
            {"accession": "QC_FAIL_TEST", "source": "test", "sequence": "A" * 50},
            reference_sequence="A" * 100,
        )

    assert result["scorable"] is False
    assert any("QC_FAIL_TEST" in r.message for r in caplog.records)


def test_score_genome_success_logs_debug(caplog):
    """A successful score must emit a DEBUG log containing the accession."""

    def _qc_pass(ref, sample, *, min_overlap=0, **kw):
        return QCResult(
            status="PASS",
            reasons=[],
            overlap_len=4,
            n_fraction=0.0,
            non_acgt_fraction=0.0,
        )

    with caplog.at_level(logging.DEBUG):
        result = score_genome(
            {"accession": "SCORED_TEST", "source": "test", "sequence": "ACGT"},
            reference_sequence="ACGT",
            qc_fn=_qc_pass,
            identify_mutations=lambda _: [],
        )

    assert result["scorable"] is True
    assert any("SCORED_TEST" in r.message for r in caplog.records)

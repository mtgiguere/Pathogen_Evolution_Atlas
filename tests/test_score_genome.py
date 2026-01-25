"""
xxxxxxxxxx
"""

def test_score_genome_sets_num_mutations_from_identified_mutations():
    from ingest.mutations import Mutation
    from ingest.scoring import score_genome

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

    summary = score_genome(canonical_record, identify_mutations=identify)

    assert summary["num_mutations"] == 3


def test_score_genome_sets_genes_affected_from_gene_mapping():
    from ingest.mutations import Mutation
    from ingest.scoring import score_genome

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
        identify_mutations=identify,
        map_genes=map_genes,
    )

    assert summary["genes_affected"] == ["Spike", "N"]


def test_score_genome_sets_risk_score_from_risk_model():
    from ingest.mutations import Mutation
    from ingest.scoring import score_genome

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
        identify_mutations=identify,
        compute_risk=compute_risk,
    )

    assert summary["risk_score"] == 7.5


def test_score_genome_computes_real_risk_score_from_mutations():
    # This one uses the REAL risk model, but still injects mutations to avoid diffing.
    from ingest.mutations import Mutation
    from ingest.risk import score_mutations
    from ingest.scoring import score_genome

    canonical_record = {
        "accession": "TEST2000",
        "sequence": "AG",
        "source": "genbank",
    }

    def identify(_record):
        return [Mutation(pos=1, ref="A", alt="G", gene="S")]

    summary = score_genome(canonical_record, identify_mutations=identify)

    assert summary["risk_score"] == float(score_mutations(identify(canonical_record))["score"])


def test_score_genome_includes_risk_explainability_fields():
    from ingest.mutations import Mutation
    from ingest.scoring import score_genome

    canonical_record = {
        "accession": "TEST3000",
        "sequence": "XX",
        "source": "genbank",
    }

    def identify(_record):
        return [
            Mutation(pos=1, ref="A", alt="G", gene="S"),
            Mutation(pos=2, ref="C", alt="T", gene="S"),
        ]

    summary = score_genome(canonical_record, identify_mutations=identify)

    assert summary["risk_level"] in {"Low", "Moderate", "High"}
    assert isinstance(summary["risk_by_gene"], dict)
    assert "S" in summary["risk_by_gene"]
    assert isinstance(summary["risk_explanation"], str)
    assert len(summary["risk_explanation"]) > 0

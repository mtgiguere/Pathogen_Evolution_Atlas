from src.ingest.mutations import Mutation
from src.ingest.risk import _GENE_WEIGHTS, score_mutations


def test_score_mutations_basic():
    muts = [
        Mutation(pos=1000, ref="A", alt="G", gene="ORF1ab"),
        Mutation(pos=2000, ref="C", alt="T", gene="ORF1ab"),
        Mutation(pos=22000, ref="A", alt="G", gene="S"),
        Mutation(pos=29000, ref="T", alt="C", gene="N"),
        Mutation(pos=30000, ref="G", alt="A", gene=None),
    ]

    result = score_mutations(muts)

    assert result["score"] == 6  # 2*1 (ORF1ab) + 1*3 (S) + 1*1 (N) = 6
    assert result["by_gene"] == {
        "ORF1ab": 2,
        "S": 1,
        "N": 1,
    }


def test_score_mutations_level_and_explanation():
    muts = [
        Mutation(pos=22000, ref="A", alt="G", gene="S"),
        Mutation(pos=22010, ref="C", alt="T", gene="S"),
        Mutation(pos=1000, ref="A", alt="G", gene="ORF1ab"),
    ]
    # Score = 2*3 (S) + 1*1 (ORF1ab) = 7

    result = score_mutations(muts)

    assert result["score"] == 7
    assert result["level"] == "High"
    assert "Spike" in result["explanation"]


def test_score_mutations_all_structural_genes_have_weight():
    """Every structural gene must have a positive weight."""
    structural_genes = {
        "S", "ORF1ab", "N",
        "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "ORF10",
    }
    for gene in structural_genes:
        assert gene in _GENE_WEIGHTS, f"Missing weight for {gene}"
        assert _GENE_WEIGHTS[gene] > 0, f"Weight for {gene} must be positive"


def test_score_mutations_immune_evasion_genes_have_elevated_weight():
    """ORF6 (IFN antagonist) and ORF8 (MHC-I downregulation) must weigh >= 2."""
    assert _GENE_WEIGHTS["ORF6"] >= 2
    assert _GENE_WEIGHTS["ORF8"] >= 2


def test_score_mutations_new_structural_genes_produce_nonzero_score():
    """A single mutation in any new structural gene must produce a scorable result."""
    for gene in ("ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "ORF10"):
        muts = [Mutation(pos=1, ref="A", alt="G", gene=gene)]
        result = score_mutations(muts)
        assert result["score"] > 0, f"score should be > 0 for gene {gene}"

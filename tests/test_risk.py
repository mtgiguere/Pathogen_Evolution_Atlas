from src.ingest.mutations import Mutation
from src.ingest.risk import score_mutations


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

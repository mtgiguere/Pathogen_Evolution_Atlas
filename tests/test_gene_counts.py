from src.ingest.counts import count_mutations_by_gene
from src.ingest.mutations import Mutation


def test_count_mutations_by_gene_basic():
    muts = [
        Mutation(pos=1000, ref="A", alt="G", gene="ORF1ab"),
        Mutation(pos=2000, ref="C", alt="T", gene="ORF1ab"),
        Mutation(pos=22000, ref="A", alt="G", gene="S"),
        Mutation(pos=29000, ref="T", alt="C", gene="N"),
    ]

    counts = count_mutations_by_gene(muts)

    assert counts == {
        "ORF1ab": 2,
        "S": 1,
        "N": 1,
    }

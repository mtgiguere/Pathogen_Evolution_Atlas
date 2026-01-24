import pytest

from src.ingest.mutations import Mutation, diff_sequences


def test_diff_sequences_simple_substitution():
    ref = "ACGT"
    sample = "AGGT"

    muts = diff_sequences(ref, sample)

    assert muts == [Mutation(pos=2, ref="C", alt="G")]


def test_diff_sequences_no_mutations():
    assert diff_sequences("ACGT", "ACGT") == []


def test_diff_sequences_length_mismatch_raises():
    with pytest.raises(ValueError):
        diff_sequences("ACGT", "ACG")


def test_diff_sequences_ignores_ambiguous_n():
    # N means "unknown" base; we don't count it as a real mutation.
    assert diff_sequences("ACGT", "ANGT") == []


def test_diff_sequences_returns_mutation_objects():
    muts = diff_sequences("ACGT", "AGGT")

    assert len(muts) == 1
    m = muts[0]

    assert isinstance(m, Mutation)
    assert m.pos == 2
    assert m.ref == "C"
    assert m.alt == "G"


def test_diff_sequences_mutations_include_gene():
    # Position 21563 is the first base of Spike (S).
    ref = "A" * 21562 + "C"
    sample = "A" * 21562 + "G"

    muts = diff_sequences(ref, sample)

    assert muts == [Mutation(pos=21563, ref="C", alt="G", gene="S")]

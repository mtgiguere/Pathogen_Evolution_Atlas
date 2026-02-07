
from src.ingest.mutations import Mutation, diff_sequences
from src.ingest.mutations import qc_compare_to_reference


def test_qc_fails_low_overlap():
    ref = "A" * 100
    sample = "A" * 50

    qc = qc_compare_to_reference(ref, sample, min_overlap=80)

    assert qc.status == "FAIL"
    assert "LOW_OVERLAP" in qc.reasons
    assert qc.overlap_len == 50


def test_qc_fails_high_n_fraction():
    ref = "A" * 100
    sample = ("A" * 50) + ("N" * 50)

    qc = qc_compare_to_reference(ref, sample, min_overlap=80, max_n_fraction=0.10)

    assert qc.status == "FAIL"
    assert "HIGH_N_FRACTION" in qc.reasons
    # In overlap, N fraction should be 0.5
    assert abs(qc.n_fraction - 0.5) < 1e-9


def test_qc_fails_high_non_acgt_fraction():
    ref = "A" * 100
    sample = ("A" * 90) + ("R" * 10)  # IUPAC ambiguity

    qc = qc_compare_to_reference(ref, sample, min_overlap=80, max_non_acgt_fraction=0.05)

    assert qc.status == "FAIL"
    assert "HIGH_NON_ACGT_FRACTION" in qc.reasons


def test_qc_passes_clean_sample():
    ref = "ACGT" * 25  # 100 bases
    sample = "ACGT" * 25

    qc = qc_compare_to_reference(ref, sample, min_overlap=80)

    assert qc.status == "PASS"
    assert qc.reasons == []


def test_diff_sequences_simple_substitution():
    ref = "ACGT"
    sample = "AGGT"

    muts = diff_sequences(ref, sample)

    assert muts == [Mutation(pos=2, ref="C", alt="G")]


def test_diff_sequences_no_mutations():
    assert diff_sequences("ACGT", "ACGT") == []


# def test_diff_sequences_length_mismatch_raises():
#     with pytest.raises(ValueError):
#         diff_sequences("ACGT", "ACG")


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

def test_diff_sequences_ignores_ambiguous_and_gap_bases():
    ref = "ACGTACGT"
    sample = "ACNT-RYS"  # N, gap, and IUPAC ambiguity codes should be ignored
    muts = diff_sequences(ref, sample)
    assert muts == []
from src.ingest.genes import gene_for_position


def test_gene_for_position_spike():
    # Spike starts at 21563 in the SARS-CoV-2 reference genome (1-based)
    assert gene_for_position(21563) == "S"

def test_gene_for_position_orf1ab():
    # ORF1ab starts at position 266 in SARS-CoV-2 reference
    assert gene_for_position(266) == "ORF1ab"

def test_gene_for_position_n():
    # N gene: 28274–29533 (SARS-CoV-2 reference)
    assert gene_for_position(28274) == "N"

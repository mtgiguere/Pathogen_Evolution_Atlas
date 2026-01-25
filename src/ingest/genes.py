"""
genes.py
"""
from __future__ import annotations


def gene_for_position(pos: int) -> str | None:
    """
    Return the SARS-CoV-2 gene name for a 1-based nucleotide position.

    v0 scope:
    - ORF1ab
    - Spike (S)
    - Nucleocapsid (N)
    - Return None if position is outside known ranges
    """
    # SARS-CoV-2 reference (NC_045512.2), 1-based inclusive coordinates

    # ORF1ab: 266–21555
    if 266 <= pos <= 21555:
        return "ORF1ab"

    # Spike (S): 21563–25384
    if 21563 <= pos <= 25384:
        return "S"

    # Nucleocapsid (N): 28274–29533
    if 28274 <= pos <= 29533:
        return "N"

    return None

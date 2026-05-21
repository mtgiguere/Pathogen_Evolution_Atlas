"""
risk.py
"""

from __future__ import annotations

from collections.abc import Iterable
from typing import Any

from .counts import count_mutations_by_gene
from .mutations import Mutation

# Weights reflect biological importance of each gene as a surveillance signal.
# S=3: primary immune target, most selection pressure.
# ORF3a/E/ORF6/ORF8=2: virulence and immune-evasion factors.
# Remaining structural genes=1: important but lower direct immune pressure.
_GENE_WEIGHTS: dict[str, int] = {
    "S": 3,  # Spike — primary antibody target
    "ORF3a": 2,  # viroporin — virulence, apoptosis
    "E": 2,  # envelope — virulence, ion channel
    "ORF6": 2,  # IFN antagonist — interferon evasion
    "ORF8": 2,  # MHC-I downregulation — immune evasion
    "ORF1ab": 1,  # polyprotein (NSPs 1-16)
    "M": 1,  # membrane — structural assembly
    "N": 1,  # nucleocapsid — structural
    "ORF7a": 1,  # immune modulation
    "ORF7b": 1,  # minor accessory
    "ORF10": 1,  # minor / uncertain expression
}

# Simple label thresholds (v1)
_LOW_MAX = 2
_MODERATE_MAX = 6  # 3–6 inclusive is Moderate; 7+ is High


def _gene_label(gene: str) -> str:
    return {
        "S": "Spike",
        "ORF1ab": "ORF1ab (polyprotein)",
        "ORF3a": "ORF3a (viroporin)",
        "E": "Envelope",
        "M": "Membrane",
        "ORF6": "ORF6 (IFN antagonist)",
        "ORF7a": "ORF7a",
        "ORF7b": "ORF7b",
        "ORF8": "ORF8 (immune evasion)",
        "N": "Nucleocapsid",
        "ORF10": "ORF10",
    }.get(gene, gene)


def score_mutations(mutations: Iterable[Mutation]) -> dict[str, Any]:
    by_gene = count_mutations_by_gene(mutations)

    score = 0
    for gene, n in by_gene.items():
        score += n * _GENE_WEIGHTS.get(gene, 0)

    if score <= _LOW_MAX:
        level = "Low"
    elif score <= _MODERATE_MAX:
        level = "Moderate"
    else:
        level = "High"

    top_gene = None
    top_impact = -1
    for gene, n in by_gene.items():
        impact = n * _GENE_WEIGHTS.get(gene, 0)
        if impact > top_impact:
            top_gene = gene
            top_impact = impact

    if top_gene is None:
        explanation = "No gene-attributed mutations detected."
    else:
        explanation = (
            f"{level} risk driven mostly by {_gene_label(top_gene)} "
            f"({by_gene[top_gene]} mutations; weight {_GENE_WEIGHTS.get(top_gene, 0)})."
        )

    return {
        "score": score,
        "by_gene": by_gene,
        "level": level,
        "explanation": explanation,
    }

from __future__ import annotations

from collections.abc import Iterable
from typing import Any

from .counts import count_mutations_by_gene
from .mutations import Mutation

_GENE_WEIGHTS: dict[str, int] = {
    "S": 3,
    "ORF1ab": 1,
    "N": 1,
}

# Simple label thresholds (v1)
_LOW_MAX = 2
_MODERATE_MAX = 6  # 3–6 inclusive is Moderate; 7+ is High


def _gene_label(gene: str) -> str:
    """Human-friendly names for explanations."""
    return {
        "S": "Spike",
        "ORF1ab": "ORF1ab",
        "N": "Nucleocapsid",
    }.get(gene, gene)


def score_mutations(mutations: Iterable[Mutation]) -> dict[str, Any]:
    by_gene = count_mutations_by_gene(mutations)

    score = 0
    for gene, n in by_gene.items():
        score += n * _GENE_WEIGHTS.get(gene, 0)

    # Level based on simple thresholds
    if score <= _LOW_MAX:
        level = "Low"
    elif score <= _MODERATE_MAX:
        level = "Moderate"
    else:
        level = "High"

    # Explanation: mention the biggest contributor (by weighted impact)
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

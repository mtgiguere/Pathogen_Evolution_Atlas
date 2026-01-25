"""
counts.py
"""
from __future__ import annotations

from collections import Counter
from collections.abc import Iterable

from .mutations import Mutation


def count_mutations_by_gene(mutations: Iterable[Mutation]) -> dict[str, int]:
    """
    Count how many mutations fall in each gene.
    Mutations with gene=None are ignored.
    """
    return dict(
        Counter(m.gene for m in mutations if m.gene is not None)
    )

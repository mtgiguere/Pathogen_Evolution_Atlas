"""
mutations.py
"""
from __future__ import annotations

from dataclasses import dataclass

from .genes import gene_for_position


@dataclass(frozen=True)
class Mutation:
    pos: int
    ref: str
    alt: str
    gene: str | None = None


def diff_sequences(ref: str, sample: str) -> list[Mutation]:
    # Real-world sequences are often trimmed/partial.
    # For v1, diff only the overlapping region.
    L = min(len(ref), len(sample))
    ref = ref[:L]
    sample = sample[:L]

    # if len(ref) != len(sample):
    #     raise ValueError("Sequences must be the same length")

    mutations: list[Mutation] = []

    for i, (r, s) in enumerate(zip(ref, sample, strict=True), start=1):
        if r.upper() == "N" or s.upper() == "N":
            continue
        if r != s:
            mutations.append(Mutation(pos=i, ref=r, alt=s, gene=gene_for_position(i)))

    return mutations

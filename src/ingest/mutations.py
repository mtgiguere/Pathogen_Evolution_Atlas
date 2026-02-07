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
    ref = ref[:L].upper()
    sample = sample[:L].upper()

    # Only treat strict A/C/G/T as comparable; skip everything else (N, gaps, IUPAC ambiguity)
    VALID = {"A", "C", "G", "T"}

    mutations: list[Mutation] = []

    for i, (r, s) in enumerate(zip(ref, sample, strict=True), start=1):
        if r not in VALID or s not in VALID:
            continue
        if r != s:
            mutations.append(Mutation(pos=i, ref=r, alt=s, gene=gene_for_position(i)))

    return mutations

@dataclass(frozen=True)
class QCResult:
    status: str  # "PASS" or "FAIL"
    reasons: list[str]
    overlap_len: int
    n_fraction: float
    non_acgt_fraction: float


def qc_compare_to_reference(
    ref: str,
    sample: str,
    *,
    min_overlap: int = 20000,
    max_n_fraction: float = 0.10,
    max_non_acgt_fraction: float = 0.02,
) -> QCResult:
    """
    Quick sanity checks so we don't score partial/garbage sequences.

    - overlap_len: how many bases we can compare (min length of the two)
    - n_fraction: fraction of 'N' bases in the overlapping region
    - non_acgt_fraction: fraction of bases that are not A/C/G/T/N in the overlap
    """
    L = min(len(ref), len(sample))
    if L == 0:
        return QCResult(
            status="FAIL",
            reasons=["LOW_OVERLAP"],
            overlap_len=0,
            n_fraction=1.0,
            non_acgt_fraction=0.0,
        )

    ref_olap = ref[:L].upper()
    sample_olap = sample[:L].upper()

    n_count = sum(1 for b in sample_olap if b == "N")
    non_acgt_count = sum(1 for b in sample_olap if b not in {"A", "C", "G", "T", "N"})

    n_fraction = n_count / L
    non_acgt_fraction = non_acgt_count / L

    reasons: list[str] = []
    if L < min_overlap:
        reasons.append("LOW_OVERLAP")
    if n_fraction > max_n_fraction:
        reasons.append("HIGH_N_FRACTION")
    if non_acgt_fraction > max_non_acgt_fraction:
        reasons.append("HIGH_NON_ACGT_FRACTION")

    status = "FAIL" if reasons else "PASS"
    return QCResult(
        status=status,
        reasons=reasons,
        overlap_len=L,
        n_fraction=n_fraction,
        non_acgt_fraction=non_acgt_fraction,
    )

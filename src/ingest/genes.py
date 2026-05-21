"""
genes.py — Map nucleotide positions to SARS-CoV-2 gene and sub-gene annotations.

All coordinates are 1-based, inclusive, relative to NC_045512.2.
"""

from __future__ import annotations

# ── Structural genes (11 genes) ───────────────────────────────────────────────
# Each entry: (start, end, gene_name)
# Order matters for overlapping regions: first match wins.
_STRUCTURAL_GENES: list[tuple[int, int, str]] = [
    (266, 21555, "ORF1ab"),
    (21563, 25384, "S"),
    (25393, 26220, "ORF3a"),
    (26245, 26472, "E"),
    (26523, 27191, "M"),
    (27202, 27387, "ORF6"),
    (27394, 27759, "ORF7a"),
    (27756, 27887, "ORF7b"),  # overlaps ORF7a tail (27756-27759); ORF7a listed first
    (27894, 28259, "ORF8"),
    (28274, 29533, "N"),
    (29558, 29674, "ORF10"),
]

# ── NSPs within ORF1ab ────────────────────────────────────────────────────────
# Ribosomal frameshifting at ~13468 produces ORF1b; we use continuous genomic coords.
_NSPS: list[tuple[int, int, str]] = [
    (266, 805, "nsp1"),
    (806, 2719, "nsp2"),
    (2720, 8554, "nsp3"),    # papain-like protease
    (8555, 10054, "nsp4"),
    (10055, 10972, "nsp5"),  # 3CLpro/Mpro — Paxlovid/nirmatrelvir target
    (10973, 11842, "nsp6"),
    (11843, 12091, "nsp7"),
    (12092, 12685, "nsp8"),
    (12686, 13024, "nsp9"),
    (13025, 13441, "nsp10"),
    (13442, 16236, "nsp12"),  # RdRp — remdesivir target (nsp11 is 7 aa, merged here)
    (16237, 18039, "nsp13"),  # helicase
    (18040, 19620, "nsp14"),  # ExoN proofreading + methyltransferase
    (19621, 20658, "nsp15"),  # EndoU
    (20659, 21552, "nsp16"),  # 2'-O-methyltransferase
]

# ── Spike protein subdomains ──────────────────────────────────────────────────
# Derived from aa boundaries (Lan et al. 2020, Walls et al. 2020).
# aa → nt: start = 21563 + (aa_start - 1) * 3
_SPIKE_DOMAINS: list[tuple[int, int, str]] = [
    (21599, 22477, "NTD"),  # N-terminal domain (aa 13-305)
    (22517, 23185, "RBD"),  # receptor-binding domain (aa 319-541)
    (24008, 24061, "FP"),   # fusion peptide (aa 816-833)
    (24296, 24514, "HR1"),  # heptad repeat 1 (aa 912-984)
    (25049, 25201, "HR2"),  # heptad repeat 2 (aa 1163-1213)
]


def gene_for_position(pos: int) -> str | None:
    """Return the structural gene name for a 1-based nucleotide position, or None."""
    for start, end, name in _STRUCTURAL_GENES:
        if start <= pos <= end:
            return name
    return None


def subgene_for_position(pos: int) -> str | None:
    """
    Return the sub-gene annotation for a position, or None.

    Within ORF1ab: returns the NSP name (e.g. "nsp5", "nsp12").
    Within S:      returns the Spike domain (e.g. "RBD", "NTD").
    All other genes and intergenic positions return None.
    """
    gene = gene_for_position(pos)
    if gene == "ORF1ab":
        for start, end, name in _NSPS:
            if start <= pos <= end:
                return name
    elif gene == "S":
        for start, end, name in _SPIKE_DOMAINS:
            if start <= pos <= end:
                return name
    return None

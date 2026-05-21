"""
Gene position mapping tests — structural genes, NSPs, and Spike subdomains.
All coordinates based on NC_045512.2 (SARS-CoV-2 reference), 1-based.
"""

import pytest

from src.ingest.genes import gene_coordinates, gene_for_position, subgene_for_position

# ── Structural gene coverage ──────────────────────────────────────────────────


@pytest.mark.parametrize(
    "pos,expected",
    [
        # Before first gene
        (1, None),
        (265, None),
        # ORF1ab (266-21555)
        (266, "ORF1ab"),
        (10000, "ORF1ab"),
        (21555, "ORF1ab"),
        # Intergenic (21556-21562)
        (21556, None),
        (21562, None),
        # S — Spike (21563-25384)
        (21563, "S"),
        (23000, "S"),
        (25384, "S"),
        # Intergenic (25385-25392)
        (25392, None),
        # ORF3a (25393-26220)
        (25393, "ORF3a"),
        (26220, "ORF3a"),
        # Intergenic (26221-26244)
        (26221, None),
        # E — Envelope (26245-26472)
        (26245, "E"),
        (26472, "E"),
        # Intergenic (26473-26522)
        (26473, None),
        # M — Membrane (26523-27191)
        (26523, "M"),
        (27191, "M"),
        # Intergenic (27192-27201)
        (27192, None),
        # ORF6 (27202-27387)
        (27202, "ORF6"),
        (27387, "ORF6"),
        # Intergenic (27388-27393)
        (27393, None),
        # ORF7a (27394-27759)
        (27394, "ORF7a"),
        (27759, "ORF7a"),
        # ORF7a/ORF7b overlap (27756-27759): ORF7a is listed first, takes priority
        (27756, "ORF7a"),
        # ORF7b (27756-27887) — positions past the overlap
        (27760, "ORF7b"),
        (27887, "ORF7b"),
        # Intergenic (27888-27893)
        (27888, None),
        # ORF8 (27894-28259)
        (27894, "ORF8"),
        (28259, "ORF8"),
        # Intergenic (28260-28273)
        (28273, None),
        # N — Nucleocapsid (28274-29533)
        (28274, "N"),
        (29533, "N"),
        # Intergenic (29534-29557)
        (29534, None),
        # ORF10 (29558-29674)
        (29558, "ORF10"),
        (29674, "ORF10"),
        # After last gene
        (29675, None),
        (30000, None),
    ],
)
def test_gene_for_position(pos, expected):
    assert gene_for_position(pos) == expected


# ── NSP sub-annotation within ORF1ab ─────────────────────────────────────────


@pytest.mark.parametrize(
    "pos,expected_subgene",
    [
        # nsp1 (266-805)
        (266, "nsp1"),
        (805, "nsp1"),
        # nsp2 (806-2719)
        (806, "nsp2"),
        (2719, "nsp2"),
        # nsp3 (2720-8554) — papain-like protease
        (2720, "nsp3"),
        (8554, "nsp3"),
        # nsp4 (8555-10054)
        (8555, "nsp4"),
        (10054, "nsp4"),
        # nsp5 (10055-10972) — 3CLpro/Mpro, Paxlovid target
        (10055, "nsp5"),
        (10500, "nsp5"),
        (10972, "nsp5"),
        # nsp6 (10973-11842)
        (10973, "nsp6"),
        (11842, "nsp6"),
        # nsp7 (11843-12091)
        (11843, "nsp7"),
        (12091, "nsp7"),
        # nsp8 (12092-12685)
        (12092, "nsp8"),
        (12685, "nsp8"),
        # nsp9 (12686-13024)
        (12686, "nsp9"),
        (13024, "nsp9"),
        # nsp10 (13025-13441)
        (13025, "nsp10"),
        (13441, "nsp10"),
        # nsp12 (13442-16236) — RdRp, remdesivir target
        (13442, "nsp12"),
        (15000, "nsp12"),
        (16236, "nsp12"),
        # nsp13 (16237-18039) — helicase
        (16237, "nsp13"),
        (18039, "nsp13"),
        # nsp14 (18040-19620) — ExoN proofreading
        (18040, "nsp14"),
        (19620, "nsp14"),
        # nsp15 (19621-20658) — EndoU
        (19621, "nsp15"),
        (20658, "nsp15"),
        # nsp16 (20659-21552) — 2'-O-methyltransferase
        (20659, "nsp16"),
        (21552, "nsp16"),
        # ORF1ab tail beyond last NSP — no subgene
        (21553, None),
        (21555, None),
    ],
)
def test_subgene_for_position_nsps(pos, expected_subgene):
    assert subgene_for_position(pos) == expected_subgene


# ── Spike subdomain sub-annotation ───────────────────────────────────────────


@pytest.mark.parametrize(
    "pos,expected_subgene",
    [
        # Signal peptide region (21563-21598) — no specific domain yet
        (21563, None),
        (21598, None),
        # NTD — N-terminal domain (21599-22477)
        (21599, "NTD"),
        (22000, "NTD"),
        (22477, "NTD"),
        # Gap between NTD and RBD (22478-22516)
        (22478, None),
        (22516, None),
        # RBD — receptor-binding domain (22517-23185) — primary immune target
        (22517, "RBD"),
        (22800, "RBD"),
        (23185, "RBD"),
        # Between RBD and FP
        (23186, None),
        # FP — fusion peptide (24008-24061)
        (24008, "FP"),
        (24061, "FP"),
        # HR1 — heptad repeat 1 (24296-24514)
        (24296, "HR1"),
        (24514, "HR1"),
        # HR2 — heptad repeat 2 (25049-25201)
        (25049, "HR2"),
        (25201, "HR2"),
        # End of Spike — no specific domain
        (25202, None),
        (25384, None),
    ],
)
def test_subgene_for_position_spike_domains(pos, expected_subgene):
    assert subgene_for_position(pos) == expected_subgene


# ── Other structural genes have no subgene annotation ────────────────────────


@pytest.mark.parametrize(
    "pos,gene",
    [
        (25393, "ORF3a"),
        (26245, "E"),
        (26523, "M"),
        (27202, "ORF6"),
        (27394, "ORF7a"),
        (27760, "ORF7b"),
        (27894, "ORF8"),
        (28274, "N"),
        (29558, "ORF10"),
    ],
)
def test_subgene_none_for_other_structural_genes(pos, gene):
    assert gene_for_position(pos) == gene
    assert subgene_for_position(pos) is None


# ── Intergenic positions have no subgene ─────────────────────────────────────


def test_subgene_none_for_intergenic():
    assert subgene_for_position(265) is None
    assert subgene_for_position(21556) is None
    assert subgene_for_position(29675) is None


# ── gene_coordinates ──────────────────────────────────────────────────────────


def test_gene_coordinates_returns_eleven_genes():
    coords = gene_coordinates()
    assert len(coords) == 11


def test_gene_coordinates_contains_expected_genes():
    names = {name for _, _, name in gene_coordinates()}
    assert names == {"ORF1ab", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF10"}


def test_gene_coordinates_orf1ab_bounds():
    coords = {name: (s, e) for s, e, name in gene_coordinates()}
    assert coords["ORF1ab"] == (266, 21555)


def test_gene_coordinates_returns_copy():
    c1 = gene_coordinates()
    c1.append((0, 0, "fake"))
    assert len(gene_coordinates()) == 11

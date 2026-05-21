"""
Immune escape mutation lookup tests — TDD, written before implementation.
All tests expected to FAIL until src/ingest/escape.py is implemented.
"""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

import pytest

from src.ingest.escape import (
    EscapeMatch,
    EscapeMutation,
    escape_from_nt_mutations,
    escape_summary,
    find_escape_mutations,
    load_escape_catalogue,
)
from src.ingest.lineage import AminoAcidMutation
from src.ingest.mutations import Mutation

# ── helpers ───────────────────────────────────────────────────────────────────


def _entry(
    gene: str,
    aa_pos: int,
    ref_aa: str,
    alt_aa: str,
    mechanism: str = "antibody_escape",
    antibodies: list[str] | None = None,
    notes: str = "",
) -> EscapeMutation:
    return EscapeMutation(
        gene=gene,
        aa_pos=aa_pos,
        ref_aa=ref_aa,
        alt_aa=alt_aa,
        mechanism=mechanism,
        antibodies_affected=antibodies or [],
        notes=notes,
    )


E484K = _entry("S", 484, "E", "K", antibodies=["Bamlanivimab", "REGN10933"])
E484A = _entry("S", 484, "E", "A", antibodies=["LY-CoV016"])
K417N = _entry("S", 417, "K", "N", antibodies=["REGN10933", "LY-CoV016"])
L452R = _entry("S", 452, "L", "R", antibodies=["Bamlanivimab"])

CATALOGUE = [E484K, E484A, K417N, L452R]


# ── find_escape_mutations ─────────────────────────────────────────────────────


def test_find_escape_mutations_exact_match():
    observed = [AminoAcidMutation(gene="S", aa_pos=484, ref_aa="E", alt_aa="K")]
    matches = find_escape_mutations(observed, CATALOGUE)
    assert len(matches) == 1
    assert matches[0].escape_entry == E484K
    assert matches[0].observed_mutation == observed[0]


def test_find_escape_mutations_no_match():
    observed = [AminoAcidMutation(gene="S", aa_pos=999, ref_aa="A", alt_aa="T")]
    assert find_escape_mutations(observed, CATALOGUE) == []


def test_find_escape_mutations_multiple_matches():
    observed = [
        AminoAcidMutation(gene="S", aa_pos=484, ref_aa="E", alt_aa="K"),
        AminoAcidMutation(gene="S", aa_pos=417, ref_aa="K", alt_aa="N"),
    ]
    matches = find_escape_mutations(observed, CATALOGUE)
    assert len(matches) == 2
    matched_entries = [m.escape_entry for m in matches]
    assert E484K in matched_entries
    assert K417N in matched_entries


def test_find_escape_mutations_same_position_different_alt():
    # E484A should NOT match the E484K catalogue entry
    observed = [AminoAcidMutation(gene="S", aa_pos=484, ref_aa="E", alt_aa="A")]
    matches = find_escape_mutations(observed, CATALOGUE)
    assert len(matches) == 1
    assert matches[0].escape_entry == E484A  # matched E484A, not E484K


def test_find_escape_mutations_empty_observed():
    assert find_escape_mutations([], CATALOGUE) == []


def test_find_escape_mutations_empty_catalogue():
    observed = [AminoAcidMutation(gene="S", aa_pos=484, ref_aa="E", alt_aa="K")]
    assert find_escape_mutations(observed, []) == []


def test_find_escape_mutations_gene_mismatch():
    # Same position and aa change but different gene — should NOT match
    observed = [AminoAcidMutation(gene="N", aa_pos=484, ref_aa="E", alt_aa="K")]
    assert find_escape_mutations(observed, CATALOGUE) == []


# ── escape_from_nt_mutations ──────────────────────────────────────────────────


def _make_ref(length: int = 30000, overrides: dict[int, str] | None = None) -> str:
    bases = ["A"] * length
    if overrides:
        for pos_1based, base in overrides.items():
            bases[pos_1based - 1] = base
    return "".join(bases)


def test_escape_from_nt_mutations_translates_and_matches():
    # S aa_pos=1 codon: nt 21563-21565 = GCT (Ala)
    # Mutate nt 21563 G→T → codon TCT = Ser → nonsynonymous → gene=S,aa1,A→S
    # That won't match the catalogue, so let's make an entry for S:A1S
    cat = [_entry("S", 1, "A", "S")]
    ref = _make_ref(overrides={21563: "G", 21564: "C", 21565: "T"})
    nt_muts = [Mutation(pos=21563, ref="G", alt="T", gene="S")]
    matches = escape_from_nt_mutations(nt_muts, ref, cat)
    assert len(matches) == 1
    assert matches[0].escape_entry.aa_pos == 1


def test_escape_from_nt_mutations_synonymous_no_match():
    # GCT → GCC both Ala → synonymous → no AA change → no escape match
    cat = [_entry("S", 1, "A", "A")]  # degenerate entry, should never match
    ref = _make_ref(overrides={21563: "G", 21564: "C", 21565: "T"})
    nt_muts = [Mutation(pos=21565, ref="T", alt="C", gene="S")]
    assert escape_from_nt_mutations(nt_muts, ref, cat) == []


def test_escape_from_nt_mutations_empty():
    assert escape_from_nt_mutations([], _make_ref(), CATALOGUE) == []


# ── escape_summary ────────────────────────────────────────────────────────────


def test_escape_summary_empty():
    summary = escape_summary([])
    assert summary["total"] == 0
    assert summary["by_mechanism"] == {}
    assert summary["antibodies_affected"] == []


def test_escape_summary_counts_by_mechanism():
    matches = [
        EscapeMatch(
            observed_mutation=AminoAcidMutation(gene="S", aa_pos=484, ref_aa="E", alt_aa="K"),
            escape_entry=E484K,
        ),
        EscapeMatch(
            observed_mutation=AminoAcidMutation(gene="S", aa_pos=417, ref_aa="K", alt_aa="N"),
            escape_entry=K417N,
        ),
    ]
    summary = escape_summary(matches)
    assert summary["total"] == 2
    assert summary["by_mechanism"]["antibody_escape"] == 2


def test_escape_summary_unique_antibodies():
    # E484K affects [Bamlanivimab, REGN10933], K417N affects [REGN10933, LY-CoV016]
    # unique antibodies = {Bamlanivimab, REGN10933, LY-CoV016}
    matches = [
        EscapeMatch(
            observed_mutation=AminoAcidMutation(gene="S", aa_pos=484, ref_aa="E", alt_aa="K"),
            escape_entry=E484K,
        ),
        EscapeMatch(
            observed_mutation=AminoAcidMutation(gene="S", aa_pos=417, ref_aa="K", alt_aa="N"),
            escape_entry=K417N,
        ),
    ]
    summary = escape_summary(matches)
    assert set(summary["antibodies_affected"]) == {"Bamlanivimab", "REGN10933", "LY-CoV016"}


def test_escape_summary_mixed_mechanisms():
    ab_entry = _entry("S", 484, "E", "K", mechanism="antibody_escape")
    vax_entry = _entry("S", 501, "N", "Y", mechanism="vaccine_reduced_neutralization")
    matches = [
        EscapeMatch(
            observed_mutation=AminoAcidMutation(gene="S", aa_pos=484, ref_aa="E", alt_aa="K"),
            escape_entry=ab_entry,
        ),
        EscapeMatch(
            observed_mutation=AminoAcidMutation(gene="S", aa_pos=501, ref_aa="N", alt_aa="Y"),
            escape_entry=vax_entry,
        ),
    ]
    summary = escape_summary(matches)
    assert summary["by_mechanism"]["antibody_escape"] == 1
    assert summary["by_mechanism"]["vaccine_reduced_neutralization"] == 1


# ── load_escape_catalogue ─────────────────────────────────────────────────────


def test_load_escape_catalogue_from_json():
    path = Path("data/escape/catalogue.json")
    entries = load_escape_catalogue(path)
    assert len(entries) >= 5
    genes = {e.gene for e in entries}
    assert "S" in genes
    # All mechanisms must be known strings
    known = {"antibody_escape", "vaccine_reduced_neutralization", "antiviral_resistance"}
    for e in entries:
        assert e.mechanism in known


def test_load_escape_catalogue_missing_file_raises():
    with pytest.raises(FileNotFoundError):
        load_escape_catalogue(Path("data/escape/nonexistent.json"))


def test_load_escape_catalogue_roundtrip():
    entry = _entry("S", 484, "E", "K", mechanism="antibody_escape", antibodies=["TestMAb"])
    payload = [
        {
            "gene": entry.gene,
            "aa_pos": entry.aa_pos,
            "ref_aa": entry.ref_aa,
            "alt_aa": entry.alt_aa,
            "mechanism": entry.mechanism,
            "antibodies_affected": entry.antibodies_affected,
            "notes": entry.notes,
        }
    ]
    with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False, encoding="utf-8") as f:
        json.dump(payload, f)
        tmp = Path(f.name)

    loaded = load_escape_catalogue(tmp)
    assert len(loaded) == 1
    assert loaded[0].gene == "S"
    assert loaded[0].aa_pos == 484
    assert loaded[0].antibodies_affected == ["TestMAb"]

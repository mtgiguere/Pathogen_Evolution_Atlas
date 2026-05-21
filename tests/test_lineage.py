"""
Lineage classification tests — TDD, written before implementation.
All tests expected to FAIL until src/ingest/lineage.py is implemented.
"""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

import pytest

from src.ingest.lineage import (
    AminoAcidMutation,
    LineageClassifier,
    LineageSignature,
    load_signatures,
    translate_mutation,
)
from src.ingest.mutations import Mutation

# ── translate_mutation ────────────────────────────────────────────────────────

# Reference codon helpers
# S gene starts at NT pos 266+0 = 266 (ORF1ab=266), S=21563
# S aa_pos=1 → nt 21563,21564,21565 → ref codon = first 3 nt of S gene
# We build a minimal fake reference large enough to hold the codons we test.


def _make_ref(length: int = 30000, overrides: dict[int, str] | None = None) -> str:
    """Build a reference string of `length` A's with optional per-position overrides (1-based)."""
    bases = ["A"] * length
    if overrides:
        for pos_1based, base in overrides.items():
            bases[pos_1based - 1] = base
    return "".join(bases)


# S gene aa=1 codon occupies nt positions 21563, 21564, 21565 (1-based)
# We set them to G, C, T  → codon GCT = Ala (A)
# Mutate nt 21563 G→T → codon TCT = Ser (S) → nonsynonymous


def test_translate_mutation_nonsynonymous():
    ref = _make_ref(overrides={21563: "G", 21564: "C", 21565: "T"})
    mut = Mutation(pos=21563, ref="G", alt="T", gene="S")
    aa_mut = translate_mutation(mut, ref)
    assert aa_mut is not None
    assert aa_mut.gene == "S"
    assert aa_mut.aa_pos == 1
    assert aa_mut.ref_aa == "A"   # GCT = Ala
    assert aa_mut.alt_aa == "S"   # TCT = Ser


def test_translate_mutation_synonymous_returns_none():
    # GCT → GCC both encode Ala → synonymous → None
    ref = _make_ref(overrides={21563: "G", 21564: "C", 21565: "T"})
    mut = Mutation(pos=21565, ref="T", alt="C", gene="S")  # GCT → GCC
    assert translate_mutation(mut, ref) is None


def test_translate_mutation_intergenic_returns_none():
    ref = _make_ref()
    mut = Mutation(pos=100, ref="A", alt="T", gene=None)
    assert translate_mutation(mut, ref) is None


def test_translate_mutation_second_codon():
    # S aa_pos=2 → nt 21566, 21567, 21568
    # Set codon to TTT = Phe (F); mutate first base T→A → ATT = Ile (I)
    ref = _make_ref(overrides={21566: "T", 21567: "T", 21568: "T"})
    mut = Mutation(pos=21566, ref="T", alt="A", gene="S")
    aa_mut = translate_mutation(mut, ref)
    assert aa_mut is not None
    assert aa_mut.aa_pos == 2
    assert aa_mut.ref_aa == "F"   # TTT = Phe
    assert aa_mut.alt_aa == "I"   # ATT = Ile


def test_translate_mutation_unknown_gene_returns_none():
    ref = _make_ref()
    mut = Mutation(pos=100, ref="A", alt="T", gene="UNKNOWN_GENE")
    assert translate_mutation(mut, ref) is None


# ── LineageClassifier ─────────────────────────────────────────────────────────


def _make_sig(
    name: str,
    aa_muts: list[tuple[str, int, str, str]],
    min_hit: float = 0.6,
) -> LineageSignature:
    defining = [
        AminoAcidMutation(gene=g, aa_pos=p, ref_aa=r, alt_aa=a)
        for g, p, r, a in aa_muts
    ]
    return LineageSignature(
        name=name,
        display_name=name,
        who_label="",
        who_class="VOC",
        min_hit_fraction=min_hit,
        defining_mutations=defining,
    )


def test_classify_empty_mutations_returns_unknown():
    sig = _make_sig("Delta", [("S", 452, "L", "R"), ("S", 478, "T", "K")])
    classifier = LineageClassifier(reference_sequence="A" * 30000, signatures=[sig])
    result = classifier.classify([])
    assert result.lineage == "Unknown"


def test_classify_all_defining_mutations_match():
    sig = _make_sig(
        "Delta",
        [("S", 452, "L", "R"), ("S", 478, "T", "K"), ("S", 614, "D", "G")],
        min_hit=0.6,
    )
    # Provide all three AA mutations
    observed = [
        AminoAcidMutation(gene="S", aa_pos=452, ref_aa="L", alt_aa="R"),
        AminoAcidMutation(gene="S", aa_pos=478, ref_aa="T", alt_aa="K"),
        AminoAcidMutation(gene="S", aa_pos=614, ref_aa="D", alt_aa="G"),
    ]
    classifier = LineageClassifier(reference_sequence="A" * 30000, signatures=[sig])
    result = classifier.classify_aa(observed)
    assert result.lineage == "Delta"
    assert result.confidence == pytest.approx(1.0)


def test_classify_partial_match_above_threshold():
    sig = _make_sig(
        "Delta",
        [("S", 452, "L", "R"), ("S", 478, "T", "K"), ("S", 614, "D", "G")],
        min_hit=0.6,
    )
    # 2 out of 3 → confidence 0.667 > 0.6
    observed = [
        AminoAcidMutation(gene="S", aa_pos=452, ref_aa="L", alt_aa="R"),
        AminoAcidMutation(gene="S", aa_pos=614, ref_aa="D", alt_aa="G"),
    ]
    classifier = LineageClassifier(reference_sequence="A" * 30000, signatures=[sig])
    result = classifier.classify_aa(observed)
    assert result.lineage == "Delta"
    assert result.confidence == pytest.approx(2 / 3)


def test_classify_below_threshold_returns_unknown():
    sig = _make_sig(
        "Delta",
        [("S", 452, "L", "R"), ("S", 478, "T", "K"), ("S", 614, "D", "G")],
        min_hit=0.8,
    )
    # 1 out of 3 → 0.333 < 0.8
    observed = [AminoAcidMutation(gene="S", aa_pos=452, ref_aa="L", alt_aa="R")]
    classifier = LineageClassifier(reference_sequence="A" * 30000, signatures=[sig])
    result = classifier.classify_aa(observed)
    assert result.lineage == "Unknown"


def test_classify_picks_highest_confidence_lineage():
    sig_delta = _make_sig(
        "Delta",
        [("S", 452, "L", "R"), ("S", 478, "T", "K"), ("S", 614, "D", "G")],
        min_hit=0.5,
    )
    sig_ba2 = _make_sig(
        "BA.2",
        [("S", 477, "S", "N"), ("S", 478, "T", "K"), ("S", 614, "D", "G")],
        min_hit=0.5,
    )
    # Observed: Delta gets 2/3, BA.2 gets 3/3
    observed = [
        AminoAcidMutation(gene="S", aa_pos=477, ref_aa="S", alt_aa="N"),
        AminoAcidMutation(gene="S", aa_pos=478, ref_aa="T", alt_aa="K"),
        AminoAcidMutation(gene="S", aa_pos=614, ref_aa="D", alt_aa="G"),
    ]
    classifier = LineageClassifier(
        reference_sequence="A" * 30000, signatures=[sig_delta, sig_ba2]
    )
    result = classifier.classify_aa(observed)
    assert result.lineage == "BA.2"
    assert result.confidence == pytest.approx(1.0)


def test_classify_supporting_and_missing_mutations_populated():
    sig = _make_sig(
        "Delta",
        [("S", 452, "L", "R"), ("S", 478, "T", "K"), ("S", 614, "D", "G")],
        min_hit=0.5,
    )
    observed = [
        AminoAcidMutation(gene="S", aa_pos=452, ref_aa="L", alt_aa="R"),
        AminoAcidMutation(gene="S", aa_pos=614, ref_aa="D", alt_aa="G"),
    ]
    classifier = LineageClassifier(reference_sequence="A" * 30000, signatures=[sig])
    result = classifier.classify_aa(observed)
    assert result.lineage == "Delta"
    assert len(result.supporting_mutations) == 2
    assert len(result.missing_mutations) == 1
    missing = result.missing_mutations[0]
    assert missing.aa_pos == 478


# ── load_signatures ───────────────────────────────────────────────────────────


def test_load_signatures_from_json():
    path = Path("data/lineages/signatures.json")
    sigs = load_signatures(path)
    names = {s.name for s in sigs}
    assert len(sigs) >= 5
    assert "B.1.617.2" in names  # Delta
    assert "JN.1" in names       # recent Omicron descendant
    for sig in sigs:
        assert len(sig.defining_mutations) >= 2
        assert 0 < sig.min_hit_fraction <= 1.0


def test_load_signatures_missing_file_raises():
    with pytest.raises(FileNotFoundError):
        load_signatures(Path("data/lineages/nonexistent.json"))


def test_load_signatures_roundtrip():
    """Serialise a signature to JSON, reload it, check it matches."""
    sig = _make_sig(
        "TestVariant",
        [("S", 452, "L", "R"), ("S", 614, "D", "G")],
        min_hit=0.7,
    )
    payload = [
        {
            "name": sig.name,
            "display_name": sig.display_name,
            "who_label": sig.who_label,
            "who_class": sig.who_class,
            "min_hit_fraction": sig.min_hit_fraction,
            "defining_mutations": [
                {
                    "gene": m.gene,
                    "aa_pos": m.aa_pos,
                    "ref_aa": m.ref_aa,
                    "alt_aa": m.alt_aa,
                }
                for m in sig.defining_mutations
            ],
        }
    ]
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".json", delete=False, encoding="utf-8"
    ) as f:
        json.dump(payload, f)
        tmp_path = Path(f.name)

    loaded = load_signatures(tmp_path)
    assert len(loaded) == 1
    assert loaded[0].name == "TestVariant"
    assert loaded[0].min_hit_fraction == pytest.approx(0.7)
    assert len(loaded[0].defining_mutations) == 2
    assert loaded[0].defining_mutations[0].aa_pos == 452

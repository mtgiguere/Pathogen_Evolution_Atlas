"""
lineage.py — Pango lineage classification from nucleotide mutations.

Classification is protein-level: nucleotide mutations are translated to
amino-acid mutations, then scored against a catalogue of defining mutations
per lineage.  The lineage with the highest hit fraction (above its threshold)
is returned; ties are broken by hit count.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .mutations import Mutation

# ── Standard genetic code (NCBI codon table 1) ────────────────────────────────

_CODON_TABLE: dict[str, str] = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}

# 1-based start position of each gene's first nucleotide (NC_045512.2)
_GENE_STARTS: dict[str, int] = {
    "ORF1ab": 266,
    "S": 21563,
    "ORF3a": 25393,
    "E": 26245,
    "M": 26523,
    "ORF6": 27202,
    "ORF7a": 27394,
    "ORF7b": 27756,
    "ORF8": 27894,
    "N": 28274,
    "ORF10": 29558,
}


# ── Data classes ──────────────────────────────────────────────────────────────


@dataclass(frozen=True)
class AminoAcidMutation:
    gene: str
    aa_pos: int  # 1-based amino-acid position within the gene
    ref_aa: str  # single-letter reference amino acid
    alt_aa: str  # single-letter mutant amino acid
    nt_pos: int = 0  # originating nucleotide position (0 = unknown)
    subgene: str | None = None


@dataclass(frozen=True)
class LineageSignature:
    name: str  # Pango name, e.g. "B.1.617.2"
    display_name: str  # human-readable, e.g. "Delta"
    who_label: str  # WHO Greek letter or "" if none
    who_class: str  # "VOC", "VOI", "VBM", "former", ""
    min_hit_fraction: float  # minimum fraction to assign this lineage
    defining_mutations: list[AminoAcidMutation] = field(default_factory=list)


@dataclass(frozen=True)
class LineageCall:
    lineage: str  # Pango name or "Unknown"
    display_name: str
    who_label: str
    who_class: str
    confidence: float  # fraction of defining mutations observed
    supporting_mutations: list[AminoAcidMutation] = field(default_factory=list)
    missing_mutations: list[AminoAcidMutation] = field(default_factory=list)


_UNKNOWN_CALL = LineageCall(
    lineage="Unknown",
    display_name="Unknown",
    who_label="",
    who_class="",
    confidence=0.0,
)


# ── Translation ───────────────────────────────────────────────────────────────


def translate_mutation(mutation: Mutation, reference_sequence: str) -> AminoAcidMutation | None:
    """
    Convert a nucleotide Mutation to an AminoAcidMutation using the reference sequence.

    Returns None when:
    - the mutation's gene is None or not in _GENE_STARTS (intergenic / unknown)
    - the mutation is synonymous (ref_aa == alt_aa)
    - the affected codon contains a stop codon in the reference
    - the position is out of range
    """
    if mutation.gene is None or mutation.gene not in _GENE_STARTS:
        return None

    gene_start = _GENE_STARTS[mutation.gene]
    offset = mutation.pos - gene_start  # 0-based offset within the CDS
    if offset < 0:
        return None

    aa_pos = offset // 3 + 1  # 1-based amino-acid position
    codon_start_0 = gene_start - 1 + (aa_pos - 1) * 3  # 0-based index into ref string

    ref_str = reference_sequence.upper()
    if codon_start_0 + 3 > len(ref_str):
        return None

    ref_codon = ref_str[codon_start_0 : codon_start_0 + 3]
    if len(ref_codon) != 3 or not all(c in "ACGT" for c in ref_codon):
        return None

    ref_aa = _CODON_TABLE.get(ref_codon)
    if ref_aa is None or ref_aa == "*":
        return None

    # Build the alt codon by substituting the mutated base
    codon_offset = offset % 3
    alt_codon_list = list(ref_codon)
    alt_codon_list[codon_offset] = mutation.alt.upper()
    alt_codon = "".join(alt_codon_list)
    alt_aa = _CODON_TABLE.get(alt_codon)
    if alt_aa is None or alt_aa == "*":
        return None

    if ref_aa == alt_aa:
        return None  # synonymous

    return AminoAcidMutation(
        gene=mutation.gene,
        aa_pos=aa_pos,
        ref_aa=ref_aa,
        alt_aa=alt_aa,
        nt_pos=mutation.pos,
        subgene=mutation.subgene,
    )


# ── Classifier ────────────────────────────────────────────────────────────────


class LineageClassifier:
    def __init__(self, reference_sequence: str, signatures: list[LineageSignature]) -> None:
        self._ref = reference_sequence
        self._signatures = signatures

    def classify(self, mutations: list[Mutation]) -> LineageCall:
        """Translate NT mutations to AA, then classify."""
        aa_muts: list[AminoAcidMutation] = []
        for m in mutations:
            aa = translate_mutation(m, self._ref)
            if aa is not None:
                aa_muts.append(aa)
        return self.classify_aa(aa_muts)

    def classify_aa(self, aa_mutations: list[AminoAcidMutation]) -> LineageCall:
        """
        Classify using already-translated amino-acid mutations.

        Each signature's defining mutations are matched by (gene, aa_pos, ref_aa, alt_aa)
        key equality so that nt_pos and subgene fields do not affect matching.
        """
        observed_keys: set[tuple[str, int, str, str]] = {
            (m.gene, m.aa_pos, m.ref_aa, m.alt_aa) for m in aa_mutations
        }

        best: LineageCall | None = None

        for sig in self._signatures:
            if not sig.defining_mutations:
                continue

            supporting: list[AminoAcidMutation] = []
            missing: list[AminoAcidMutation] = []
            for dm in sig.defining_mutations:
                key = (dm.gene, dm.aa_pos, dm.ref_aa, dm.alt_aa)
                if key in observed_keys:
                    supporting.append(dm)
                else:
                    missing.append(dm)

            confidence = len(supporting) / len(sig.defining_mutations)
            if confidence < sig.min_hit_fraction:
                continue

            if (
                best is None
                or confidence > best.confidence
                or (
                    confidence == best.confidence
                    and len(supporting) > len(best.supporting_mutations)
                )
            ):
                best = LineageCall(
                    lineage=sig.name,
                    display_name=sig.display_name,
                    who_label=sig.who_label,
                    who_class=sig.who_class,
                    confidence=confidence,
                    supporting_mutations=supporting,
                    missing_mutations=missing,
                )

        return best if best is not None else _UNKNOWN_CALL


# ── Signature catalogue I/O ───────────────────────────────────────────────────


def load_signatures(path: Path) -> list[LineageSignature]:
    """Load lineage signatures from a JSON file."""
    if not path.exists():
        raise FileNotFoundError(f"Signatures file not found: {path}")

    raw: list[dict[str, Any]] = json.loads(path.read_text(encoding="utf-8"))
    signatures: list[LineageSignature] = []

    for entry in raw:
        defining = [
            AminoAcidMutation(
                gene=m["gene"],
                aa_pos=int(m["aa_pos"]),
                ref_aa=m["ref_aa"],
                alt_aa=m["alt_aa"],
            )
            for m in entry.get("defining_mutations", [])
        ]
        signatures.append(
            LineageSignature(
                name=entry["name"],
                display_name=entry.get("display_name", entry["name"]),
                who_label=entry.get("who_label", ""),
                who_class=entry.get("who_class", ""),
                min_hit_fraction=float(entry.get("min_hit_fraction", 0.6)),
                defining_mutations=defining,
            )
        )

    return signatures

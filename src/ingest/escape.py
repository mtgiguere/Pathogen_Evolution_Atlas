"""
escape.py — Known immune escape mutation lookup.

Maps amino-acid mutations against a curated catalogue of variants
known to reduce neutralization by monoclonal antibodies, vaccines,
or antiviral drugs.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .lineage import AminoAcidMutation, translate_mutation
from .mutations import Mutation


@dataclass(frozen=True)
class EscapeMutation:
    gene: str
    aa_pos: int
    ref_aa: str
    alt_aa: str
    mechanism: str  # "antibody_escape" | "vaccine_reduced_neutralization" | "antiviral_resistance"
    antibodies_affected: list[str] = field(default_factory=list)
    notes: str = ""


@dataclass(frozen=True)
class EscapeMatch:
    observed_mutation: AminoAcidMutation
    escape_entry: EscapeMutation


# Lookup key: (gene, aa_pos, ref_aa, alt_aa)
_EscapeKey = tuple[str, int, str, str]


def _key(e: EscapeMutation | AminoAcidMutation) -> _EscapeKey:
    return (e.gene, e.aa_pos, e.ref_aa, e.alt_aa)


def find_escape_mutations(
    aa_mutations: list[AminoAcidMutation],
    catalogue: list[EscapeMutation],
) -> list[EscapeMatch]:
    """Return all catalogue entries matched by the observed amino-acid mutations."""
    index: dict[_EscapeKey, EscapeMutation] = {_key(e): e for e in catalogue}
    matches: list[EscapeMatch] = []
    for obs in aa_mutations:
        entry = index.get(_key(obs))
        if entry is not None:
            matches.append(EscapeMatch(observed_mutation=obs, escape_entry=entry))
    return matches


def escape_from_nt_mutations(
    nt_mutations: list[Mutation],
    reference_sequence: str,
    catalogue: list[EscapeMutation],
) -> list[EscapeMatch]:
    """Translate nucleotide mutations to AA then look up escape catalogue."""
    aa_muts: list[AminoAcidMutation] = []
    for m in nt_mutations:
        aa = translate_mutation(m, reference_sequence)
        if aa is not None:
            aa_muts.append(aa)
    return find_escape_mutations(aa_muts, catalogue)


def escape_summary(matches: list[EscapeMatch]) -> dict[str, Any]:
    """Summarise escape matches: total count, breakdown by mechanism, unique antibodies."""
    by_mechanism: dict[str, int] = {}
    antibodies: set[str] = set()

    for m in matches:
        mech = m.escape_entry.mechanism
        by_mechanism[mech] = by_mechanism.get(mech, 0) + 1
        antibodies.update(m.escape_entry.antibodies_affected)

    return {
        "total": len(matches),
        "by_mechanism": by_mechanism,
        "antibodies_affected": sorted(antibodies),
    }


def load_escape_catalogue(path: Path) -> list[EscapeMutation]:
    """Load the escape mutation catalogue from a JSON file."""
    if not path.exists():
        raise FileNotFoundError(f"Escape catalogue not found: {path}")

    raw: list[dict[str, Any]] = json.loads(path.read_text(encoding="utf-8"))
    return [
        EscapeMutation(
            gene=e["gene"],
            aa_pos=int(e["aa_pos"]),
            ref_aa=e["ref_aa"],
            alt_aa=e["alt_aa"],
            mechanism=e.get("mechanism", "antibody_escape"),
            antibodies_affected=e.get("antibodies_affected", []),
            notes=e.get("notes", ""),
        )
        for e in raw
    ]

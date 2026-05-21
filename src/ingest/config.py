"""
config.py — Per-pathogen configuration.

A PathogenConfig holds everything that differs between pathogens:
gene coordinates, risk weights, NCBI search terms, paths to lineage
signatures and escape catalogues.  Loading from JSON makes it trivial
to add a new pathogen without touching Python code.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any


@dataclass
class GeneEntry:
    name: str
    start: int   # 1-based, inclusive
    end: int     # 1-based, inclusive


@dataclass
class PathogenConfig:
    pathogen_id: str              # slug, e.g. "sars-cov-2"
    display_name: str             # human-readable, e.g. "SARS-CoV-2"
    organism_query: str           # NCBI Entrez search term
    reference_accession: str      # e.g. "NC_045512.2"
    min_genome_length: int
    max_genome_length: int
    genes: list[GeneEntry] = field(default_factory=list)    # structural gene map; first match wins
    subgenes: list[GeneEntry] = field(default_factory=list) # NSPs, domains, etc.
    gene_weights: dict[str, int] = field(default_factory=dict)
    signatures_path: str = ""     # path to lineage signatures JSON (relative to project root)
    escape_path: str = ""         # path to escape catalogue JSON
    description: str = ""


# ── Config-driven gene lookup ─────────────────────────────────────────────────


def gene_for_position_config(pos: int, config: PathogenConfig) -> str | None:
    """Return the gene name for a position using config-defined coordinates."""
    for entry in config.genes:
        if entry.start <= pos <= entry.end:
            return entry.name
    return None


def subgene_for_position_config(pos: int, config: PathogenConfig) -> str | None:
    """Return the sub-gene name (NSP, domain, etc.) for a position using config-defined coordinates."""
    for entry in config.subgenes:
        if entry.start <= pos <= entry.end:
            return entry.name
    return None


# ── Serialisation ─────────────────────────────────────────────────────────────


def _config_to_dict(cfg: PathogenConfig) -> dict[str, Any]:
    d = asdict(cfg)
    return d


def _gene_entry_from_dict(g: dict[str, Any]) -> GeneEntry:
    return GeneEntry(name=g["name"], start=int(g["start"]), end=int(g["end"]))


def _dict_to_config(d: dict[str, Any]) -> PathogenConfig:
    genes = [_gene_entry_from_dict(g) for g in d.get("genes", [])]
    subgenes = [_gene_entry_from_dict(g) for g in d.get("subgenes", [])]
    return PathogenConfig(
        pathogen_id=d["pathogen_id"],
        display_name=d.get("display_name", d["pathogen_id"]),
        organism_query=d.get("organism_query", ""),
        reference_accession=d.get("reference_accession", ""),
        min_genome_length=int(d.get("min_genome_length", 0)),
        max_genome_length=int(d.get("max_genome_length", 999_999)),
        genes=genes,
        subgenes=subgenes,
        gene_weights=dict(d.get("gene_weights", {})),
        signatures_path=d.get("signatures_path", ""),
        escape_path=d.get("escape_path", ""),
        description=d.get("description", ""),
    )


def load_config(path: Path) -> PathogenConfig:
    """Load a PathogenConfig from a JSON file."""
    if not path.exists():
        raise FileNotFoundError(f"Pathogen config not found: {path}")
    raw = json.loads(path.read_text(encoding="utf-8"))
    return _dict_to_config(raw)


def save_config(config: PathogenConfig, path: Path) -> None:
    """Persist a PathogenConfig to JSON, creating parent directories as needed."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(_config_to_dict(config), indent=2), encoding="utf-8")


def list_pathogen_configs(directory: Path) -> list[PathogenConfig]:
    """Load all *.json configs from a directory. Returns [] if the directory is missing."""
    if not directory.exists():
        return []
    configs: list[PathogenConfig] = []
    for p in sorted(directory.glob("*.json")):
        try:
            configs.append(load_config(p))
        except (KeyError, ValueError):
            pass
    return configs


# ── Built-in defaults ─────────────────────────────────────────────────────────


def get_default_sars_cov2_config() -> PathogenConfig:
    """Return the canonical SARS-CoV-2 config (NC_045512.2 coordinates)."""
    return PathogenConfig(
        pathogen_id="sars-cov-2",
        display_name="SARS-CoV-2",
        organism_query=(
            '"Severe acute respiratory syndrome coronavirus 2"[Organism] '
            'AND "complete genome"[Title]'
        ),
        reference_accession="NC_045512.2",
        min_genome_length=29000,
        max_genome_length=31000,
        description="SARS-CoV-2 (COVID-19 causative agent). Reference: NC_045512.2 (Wuhan-Hu-1).",
        genes=[
            GeneEntry(name="ORF1ab", start=266,   end=21555),
            GeneEntry(name="S",      start=21563,  end=25384),
            GeneEntry(name="ORF3a",  start=25393,  end=26220),
            GeneEntry(name="E",      start=26245,  end=26472),
            GeneEntry(name="M",      start=26523,  end=27191),
            GeneEntry(name="ORF6",   start=27202,  end=27387),
            GeneEntry(name="ORF7a",  start=27394,  end=27759),
            GeneEntry(name="ORF7b",  start=27756,  end=27887),
            GeneEntry(name="ORF8",   start=27894,  end=28259),
            GeneEntry(name="N",      start=28274,  end=29533),
            GeneEntry(name="ORF10",  start=29558,  end=29674),
        ],
        subgenes=[
            # NSPs within ORF1ab
            GeneEntry(name="nsp1",  start=266,   end=805),
            GeneEntry(name="nsp2",  start=806,   end=2719),
            GeneEntry(name="nsp3",  start=2720,  end=8554),
            GeneEntry(name="nsp4",  start=8555,  end=10054),
            GeneEntry(name="nsp5",  start=10055, end=10972),   # 3CLpro / Paxlovid target
            GeneEntry(name="nsp6",  start=10973, end=11842),
            GeneEntry(name="nsp7",  start=11843, end=12091),
            GeneEntry(name="nsp8",  start=12092, end=12685),
            GeneEntry(name="nsp9",  start=12686, end=13024),
            GeneEntry(name="nsp10", start=13025, end=13441),
            GeneEntry(name="nsp12", start=13442, end=16236),   # RdRp / remdesivir target
            GeneEntry(name="nsp13", start=16237, end=18039),
            GeneEntry(name="nsp14", start=18040, end=19620),
            GeneEntry(name="nsp15", start=19621, end=20658),
            GeneEntry(name="nsp16", start=20659, end=21552),
            # Spike subdomains
            GeneEntry(name="NTD",   start=21599, end=22477),
            GeneEntry(name="RBD",   start=22517, end=23185),
            GeneEntry(name="FP",    start=24008, end=24061),
            GeneEntry(name="HR1",   start=24296, end=24514),
            GeneEntry(name="HR2",   start=25049, end=25201),
        ],
        gene_weights={
            "S":      3,
            "ORF3a":  2,
            "E":      2,
            "ORF6":   2,
            "ORF8":   2,
            "ORF1ab": 1,
            "M":      1,
            "N":      1,
            "ORF7a":  1,
            "ORF7b":  1,
            "ORF10":  1,
        },
        signatures_path="data/lineages/signatures.json",
        escape_path="data/escape/catalogue.json",
    )

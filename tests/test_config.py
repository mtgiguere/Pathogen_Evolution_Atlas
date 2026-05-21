"""
Multi-pathogen configuration tests — TDD, written before implementation.
All tests expected to FAIL until src/ingest/config.py is implemented.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import pytest

from src.ingest.config import (
    GeneEntry,
    PathogenConfig,
    gene_for_position_config,
    get_default_sars_cov2_config,
    list_pathogen_configs,
    load_config,
    save_config,
    subgene_for_position_config,
)

# ── PathogenConfig structure ──────────────────────────────────────────────────


def test_pathogen_config_has_required_fields():
    cfg = get_default_sars_cov2_config()
    assert cfg.pathogen_id == "sars-cov-2"
    assert cfg.display_name
    assert cfg.reference_accession == "NC_045512.2"
    assert cfg.organism_query
    assert cfg.min_genome_length > 0
    assert cfg.max_genome_length > cfg.min_genome_length


def test_default_sars_cov2_config_has_genes():
    cfg = get_default_sars_cov2_config()
    names = {g.name for g in cfg.genes}
    assert "S" in names
    assert "ORF1ab" in names
    assert "N" in names
    assert len(cfg.genes) == 11  # all structural genes


def test_default_sars_cov2_config_has_subgenes():
    cfg = get_default_sars_cov2_config()
    subnames = {g.name for g in cfg.subgenes}
    assert "nsp5" in subnames   # Paxlovid target
    assert "nsp12" in subnames  # RdRp
    assert "RBD" in subnames    # Spike receptor-binding domain


def test_default_sars_cov2_config_has_gene_weights():
    cfg = get_default_sars_cov2_config()
    assert cfg.gene_weights.get("S") == 3
    assert cfg.gene_weights.get("ORF1ab") == 1
    assert "E" in cfg.gene_weights


def test_default_sars_cov2_config_has_data_paths():
    cfg = get_default_sars_cov2_config()
    assert cfg.signatures_path
    assert cfg.escape_path


# ── config-driven gene lookup ─────────────────────────────────────────────────


def _minimal_config(genes: list[tuple[int, int, str]]) -> PathogenConfig:
    return PathogenConfig(
        pathogen_id="test",
        display_name="Test Pathogen",
        organism_query="test[Organism]",
        reference_accession="TEST001",
        min_genome_length=100,
        max_genome_length=10000,
        genes=[GeneEntry(start=s, end=e, name=n) for s, e, n in genes],
        subgenes=[],
        gene_weights={},
        signatures_path="",
        escape_path="",
        description="",
    )


def test_gene_for_position_config_basic():
    cfg = _minimal_config([(100, 500, "GeneA"), (600, 900, "GeneB")])
    assert gene_for_position_config(100, cfg) == "GeneA"
    assert gene_for_position_config(500, cfg) == "GeneA"
    assert gene_for_position_config(350, cfg) == "GeneA"
    assert gene_for_position_config(600, cfg) == "GeneB"
    assert gene_for_position_config(900, cfg) == "GeneB"


def test_gene_for_position_config_intergenic():
    cfg = _minimal_config([(100, 500, "GeneA"), (600, 900, "GeneB")])
    assert gene_for_position_config(50, cfg) is None    # before first
    assert gene_for_position_config(550, cfg) is None   # gap between genes
    assert gene_for_position_config(950, cfg) is None   # after last


def test_gene_for_position_config_empty_genes():
    cfg = _minimal_config([])
    assert gene_for_position_config(100, cfg) is None


def test_gene_for_position_config_overlap_first_wins():
    cfg = _minimal_config([(100, 400, "GeneA"), (300, 500, "GeneB")])
    assert gene_for_position_config(350, cfg) == "GeneA"  # first listed wins


def test_subgene_for_position_config_basic():
    cfg = PathogenConfig(
        pathogen_id="test",
        display_name="Test",
        organism_query="test",
        reference_accession="T1",
        min_genome_length=100,
        max_genome_length=5000,
        genes=[GeneEntry(start=100, end=2000, name="PolyProtein")],
        subgenes=[
            GeneEntry(start=100, end=600, name="nsp1"),
            GeneEntry(start=601, end=1200, name="nsp2"),
        ],
        gene_weights={},
        signatures_path="",
        escape_path="",
        description="",
    )
    assert subgene_for_position_config(100, cfg) == "nsp1"
    assert subgene_for_position_config(600, cfg) == "nsp1"
    assert subgene_for_position_config(601, cfg) == "nsp2"
    assert subgene_for_position_config(1201, cfg) is None  # beyond subgenes


def test_subgene_for_position_config_intergenic_returns_none():
    cfg = _minimal_config([(100, 500, "GeneA")])
    assert subgene_for_position_config(50, cfg) is None  # not in any gene


# ── load_config / save_config ─────────────────────────────────────────────────


def test_save_and_load_config_roundtrip():
    cfg = get_default_sars_cov2_config()
    with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
        tmp = Path(f.name)
    save_config(cfg, tmp)
    loaded = load_config(tmp)
    assert loaded.pathogen_id == cfg.pathogen_id
    assert loaded.reference_accession == cfg.reference_accession
    assert len(loaded.genes) == len(cfg.genes)
    assert len(loaded.subgenes) == len(cfg.subgenes)
    assert loaded.gene_weights == cfg.gene_weights


def test_load_config_missing_file_raises():
    with pytest.raises(FileNotFoundError):
        load_config(Path("data/pathogens/nonexistent.json"))


def test_save_config_creates_parent_dirs(tmp_path):
    cfg = get_default_sars_cov2_config()
    deep = tmp_path / "sub" / "dir" / "cfg.json"
    save_config(cfg, deep)
    assert deep.exists()


def test_load_sars_cov2_from_data_directory():
    path = Path("data/pathogens/sars-cov-2.json")
    cfg = load_config(path)
    assert cfg.pathogen_id == "sars-cov-2"
    assert cfg.reference_accession == "NC_045512.2"
    assert len(cfg.genes) >= 11


# ── list_pathogen_configs ─────────────────────────────────────────────────────


def test_list_pathogen_configs_finds_sars_cov2():
    cfgs = list_pathogen_configs(Path("data/pathogens"))
    ids = {c.pathogen_id for c in cfgs}
    assert "sars-cov-2" in ids


def test_list_pathogen_configs_empty_dir(tmp_path):
    result = list_pathogen_configs(tmp_path)
    assert result == []


def test_list_pathogen_configs_missing_dir_returns_empty():
    result = list_pathogen_configs(Path("data/pathogens/nonexistent_dir"))
    assert result == []


def test_list_pathogen_configs_finds_multiple(tmp_path):
    for pid in ("pathogen-a", "pathogen-b"):
        cfg = PathogenConfig(
            pathogen_id=pid,
            display_name=pid.title(),
            organism_query=f"{pid}[Organism]",
            reference_accession=f"REF_{pid}",
            min_genome_length=1000,
            max_genome_length=50000,
            genes=[],
            subgenes=[],
            gene_weights={},
            signatures_path="",
            escape_path="",
            description="",
        )
        save_config(cfg, tmp_path / f"{pid}.json")

    cfgs = list_pathogen_configs(tmp_path)
    assert len(cfgs) == 2
    ids = {c.pathogen_id for c in cfgs}
    assert ids == {"pathogen-a", "pathogen-b"}


# ── second pathogen sanity check ──────────────────────────────────────────────


def test_influenza_h3n2_config_loadable():
    path = Path("data/pathogens/influenza-h3n2.json")
    cfg = load_config(path)
    assert cfg.pathogen_id == "influenza-h3n2"
    names = {g.name for g in cfg.genes}
    assert "HA" in names   # haemagglutinin — primary immune target
    assert "NA" in names   # neuraminidase — Tamiflu target
    assert cfg.gene_weights.get("HA", 0) >= 2

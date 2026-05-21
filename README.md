# Pathogen Evolution Atlas

A genomic surveillance pipeline for tracking how pathogens evolve and spread —
built in Python, grounded in molecular biology, designed for real-world use.

---

## What This Is

When a new SARS-CoV-2 variant appears, public health teams need to answer three
questions quickly:

1. **What mutations does it carry?** — Are they in the Spike RBD, where antibodies bind?
   In nsp5, where Paxlovid acts? Do they match known immune escape patterns?
2. **Is it a known lineage?** — Delta, Omicron BA.2, JN.1, or something new?
3. **Is it growing?** — What is the doubling time this week compared to last?

This project automates that workflow end-to-end: from fetching raw sequences out of
NCBI GenBank, through mutation calling, lineage classification, and growth rate
estimation, to a REST API and alert system that fires when something dangerous appears.

---

## Capabilities

| Feature | Details |
|---|---|
| **Gene annotation** | All 11 SARS-CoV-2 structural genes; 15 NSPs within ORF1ab (nsp1–nsp16); 5 Spike subdomains (NTD, RBD, FP, HR1, HR2) — coordinates relative to NC\_045512.2 |
| **Lineage classification** | Translates nucleotide mutations to amino-acid changes; matches against a catalogue of Pango lineage signatures (Delta, BA.2, BA.5, XBB.1.5, JN.1) |
| **Immune escape lookup** | Checks each genome against 11 curated escape mutations — E484K, K417N, L452R, F486P, R346T, K444T and more — with the specific antibodies each evades |
| **Growth rate estimation** | Fits a log-linear model (ln count ~ week) per variant; reports growth rate (ln/week), doubling time (days), R², and trend label |
| **Risk scoring** | Weighted mutation burden per gene — Spike (×3), viroporins/IFN antagonists (×2), structural genes (×1) — yields Low / Moderate / High |
| **Automated ingestion** | Incremental NCBI search with date filtering; deduplicates against existing store; persists state across runs |
| **REST API** | FastAPI with endpoints for genomes, variants, growth summary, and escape stats; auto-generated OpenAPI docs at `/docs` |
| **Alert system** | Rule-based engine: high-risk genome, critical escape, fast-growing variant, new VOC — dispatches to log, file, or webhook (Slack/Teams) |
| **Multi-pathogen** | JSON-driven `PathogenConfig` — drop a new file in `data/pathogens/` to add Influenza H3N2, RSV, or any sequenced pathogen without code changes |

---

## Scientific Grounding

Gene coordinates, NSP boundaries, and Spike subdomain definitions follow
[NC\_045512.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2) (Wuhan-Hu-1 reference).

Key biological annotations built into the pipeline:

- **nsp5** (nt 10055–10972) — 3C-like protease (3CLpro/Mpro); target of nirmatrelvir (Paxlovid)
- **nsp12** (nt 13442–16236) — RNA-dependent RNA polymerase (RdRp); target of remdesivir
- **RBD** (aa 319–541 of Spike) — receptor-binding domain; primary target of neutralising antibodies
- **E484K / E484A** — class 2 antibody escape; reduces neutralisation by convalescent sera and several monoclonal antibodies
- **K417N / K417T** — class 1 antibody escape; disrupts binding of REGN10933, LY-CoV016
- **L452R** — class 2/3 escape; Delta and BA.4/BA.5 signature; increases ACE2 affinity
- **F486P** — XBB.1.5 signature; potent class 3 escape while restoring ACE2 binding lost by F486V

Lineage classification uses protein-level amino-acid mutation matching — the same
approach used by Nextclade and similar tools — rather than raw nucleotide similarity.

---

## Quick Start

```bash
# 1. Clone and create a virtual environment
git clone https://github.com/mtgiguere/Pathogen_Evolution_Atlas.git
cd Pathogen_Evolution_Atlas
python -m venv .venv
source .venv/bin/activate       # Windows: .venv\Scripts\activate

# 2. Install dependencies
pip install -r requirements.txt

# 3. Set your NCBI email (required by NCBI Entrez policy)
export NCBI_EMAIL="your@email.com"   # Windows: $env:NCBI_EMAIL="..."

# 4. Fetch a small set of SARS-CoV-2 genomes
python scripts/build_covid_accessions.py --n 20 --out data/accessions/covid_accessions.txt
python scripts/fetch_genbank_accessions.py \
    --accessions $(cat data/accessions/covid_accessions.txt | tr '\n' ' ') \
    --out data/raw/genomes.ndjson

# 5. Launch the dashboard
streamlit run scripts/app.py

# 6. Or start the REST API
uvicorn src.ingest.api:app --reload
# then open http://localhost:8000/docs
```

---

## Running Tests

```bash
pytest                                  # all 314 tests
pytest tests/test_lineage.py -v         # lineage classification
pytest tests/test_escape.py -v          # immune escape lookup
pytest tests/test_growth.py -v          # growth rate estimation
pytest tests/test_alerts.py -v          # alert system
```

---

## Automated Ingestion

Run once to fetch new sequences since the last run:

```bash
python scripts/run_ingest.py --max-new 100
```

State is persisted at `data/state/ingest_state.json`. Schedule via Task Scheduler
(Windows) or cron (Linux/macOS) to run daily:

```
0 6 * * * cd /path/to/project && python scripts/run_ingest.py
```

---

## Project Structure

```
src/ingest/
    genbank.py              — NCBI GenBank fetch + normalisation
    genes.py                — Gene/NSP/Spike-domain position lookup (NC_045512.2)
    mutations.py            — Nucleotide diff → Mutation objects
    lineage.py              — NT→AA translation + Pango lineage classification
    escape.py               — Immune escape mutation catalogue lookup
    growth.py               — Per-variant growth rate estimation
    scoring.py              — End-to-end genome scoring orchestration
    risk.py                 — Weighted mutation burden → risk level
    analytics.py            — Summarise a genome collection → DataFrame
    alerts.py               — Rule-based alert engine + channels
    scheduler.py            — Incremental ingestion state + runner
    config.py               — Multi-pathogen PathogenConfig (JSON-driven)
    api.py                  — FastAPI REST interface
    geography_enrichment.py — Lat/lon enrichment from BioSample

scripts/
    app.py                      — Streamlit dashboard
    run_ingest.py               — Automated ingestion CLI
    build_covid_accessions.py   — Build NCBI accession list
    fetch_genbank_accessions.py
    run_geography_enrichment.py

data/
    pathogens/      — Per-pathogen config JSON (sars-cov-2, influenza-h3n2)
    lineages/       — Pango lineage signature catalogues
    escape/         — Immune escape mutation catalogues
    accessions/     — Accession ID lists
    reference/      — Reference genome NDJSON

tests/              — 314 tests, strict TDD throughout
```

---

## Adding a New Pathogen

Create a JSON config in `data/pathogens/` — no code changes needed:

```json
{
  "pathogen_id": "rsv-b",
  "display_name": "RSV B",
  "organism_query": "\"Human respiratory syncytial virus B\"[Organism] AND \"complete genome\"[Title]",
  "reference_accession": "KX765876",
  "min_genome_length": 14000,
  "max_genome_length": 16000,
  "genes": [
    {"name": "F", "start": 5660, "end": 7390},
    {"name": "G", "start": 4300, "end": 5565},
    {"name": "N", "start": 1141, "end": 2316}
  ],
  "gene_weights": {"F": 3, "G": 2, "N": 1},
  "signatures_path": "data/lineages/rsv-b-signatures.json",
  "escape_path": "data/escape/rsv-b-catalogue.json"
}
```

---

## REST API Endpoints

| Method | Endpoint | Description |
|---|---|---|
| GET | `/health` | Liveness check |
| GET | `/genomes` | List genomes — filter by `?country=` / `?lineage=`, paginate with `?limit=&offset=` |
| GET | `/genomes/{accession}` | Single genome detail |
| GET | `/variants` | Lineage counts with WHO label/class |
| GET | `/summary` | Total genomes, risk breakdown, top lineages, escape genome count |

Interactive Swagger docs auto-generated at `http://localhost:8000/docs`.

---

## Roadmap

- [ ] Streamlit → Dash migration for production-grade map views (choropleth spread maps via Mapbox)
- [ ] Phylogenetic tree rendering
- [ ] Expanded lineage catalogue (all current WHO variants under monitoring)
- [ ] GISAID integration (requires data-use agreement)
- [ ] Docker + cloud deployment

---

## Data Sources

Genomic sequences are fetched from **NCBI GenBank** via the Entrez API.
A registered email address is required per
[NCBI policy](https://www.ncbi.nlm.nih.gov/home/about/policies/) — set via the
`NCBI_EMAIL` environment variable (never hardcoded).

Lineage signatures and immune escape catalogues are curated from published literature:
- Lan et al. 2020 — Spike protein structure
- Walls et al. 2020 — Spike subdomains and ACE2 interaction
- Starr et al. 2021 — Deep mutational scanning of ACE2 binding and antibody escape
- WHO variant classification updates

---

## License

MIT — see [LICENSE](LICENSE).

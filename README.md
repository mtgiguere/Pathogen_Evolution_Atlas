# Pathogen Evolution Atlas
⚠️ This project is under active development.


# 🧬 Pathogen Evolution Atlas

**Genomics • Geography • Machine Learning**

*A platform for tracking how pathogens evolve and spread across space and time — and predicting what comes next.*

---

## Overview

Pathogen Evolution Atlas is a genomic surveillance and forecasting system that integrates pathogen genome sequences, mutation data, collection time, and geographic metadata to model how variants emerge and spread.

The platform is designed to mirror the systems used by public health agencies, vaccine developers, and genomic surveillance teams to monitor evolutionary risk and anticipate future outbreaks.

---

## What This System Does

For each pathogen, the system produces:

- Lineage and variant tracking  
- Mutation frequency trends  
- Geographic spread maps  
- Variant growth rates  
- Short-term forecasts  
- Explainable risk scores  

Outputs can be explored through an interactive dashboard or queried via an API.

---

### ⚠️ Note on Geographic Metadata

Not all genomic datasets include geographic location information.

The current development dataset consists of recent (2025–2026) human
surveillance genomes that intentionally omit country and regional metadata.
This reflects modern privacy-preserving genomic surveillance practices.

As a result:
- Evolutionary and temporal analyses are fully supported
- Geographic analyses are conditionally supported when location metadata exists
- The system architecture supports geography, but availability depends on source data

---

## Scientific Focus

This project focuses on:

- Molecular evolution  
- Mutation dynamics  
- Protein-coding changes  
- Selection pressure  
- Phylogenetic structure  

The system is grounded in biological meaning — not just statistical pattern matching.

---

## Data Sources

The platform integrates data from:

- **GISAID**
- **NCBI GenBank**
- **Nextstrain**
- (optional) public epidemiological datasets

Each genome record includes:

- Nucleotide sequence  
- Collection date  
- Geographic location  (when available)  
- Lineage or variant annotation  

---

## Machine Learning

Machine learning is used to:

- Predict geographic spread of variants  
- Forecast mutation trends  
- Identify high-risk lineages  

Models include:

- Baseline linear models  
- Random forests  
- Gradient-boosted trees  

Predictions are explainable using feature importance and SHAP analysis.

---

## System Architecture

High-level pipeline:


The system is designed to be:

- Reproducible  
- Auditable  
- Scalable  

---

## Repository Structure

/data
/raw
/features
/src
/ingest
/features
/modeling
/notebooks
/reports
/api
/app
/docs


---

## Getting Started

Instructions for:

- Cloning the repository  
- Installing dependencies  
- Running notebooks  
- Launching the dashboard  

(added as the project develops)

---

## Roadmap

- Phase 1 — Genomic data ingestion  
- Phase 2 — Feature engineering  
- Phase 3 — Baseline evolutionary analysis  
- Phase 4 — Machine learning forecasting  
- Phase 5 — Interactive dashboard and API  
- Phase 6 — Cloud deployment  

---

## Why This Project Exists

This project is built to demonstrate:

- Computational biology  
- Machine learning  
- Geospatial analysis  
- Real-world pathogen surveillance  

It is intended as a professional, portfolio-grade research platform.

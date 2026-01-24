# Usage Guide

This document explains how to run the Pathogen Evolution Atlas code locally
and what outputs to expect.

No prior computer science background is assumed.

---

### 1. Setup

    ## Requirements
        - Python 3.12
        - Internet access (for GenBank data)

    ## Create and activate a virtual environment

        ```bash
        python -m venv .venv

    ## On Windows:

        .venv\Scripts\Activate.ps1

    ## On macOS/Linux:

        source .venv/bin/activate

    ## Install dependencies
        pip install -r requirements.txt

### 2. Required environment variables

        - NCBI requires a contact email when downloading genome data.

    ## Set this once:

        setx NCBI_EMAIL "your.email@example.com"

    ## Restart your terminal afterward.

### 3. Fetch genome data from GenBank

    ## Example: fetch two SARS-CoV-2 reference genomes.

    python -m scripts.fetch_genbank_accessions \
    --accessions NC_045512.2 MN908947.3 \
    --out data/raw/genbank.ndjson

    ## Output

        - data/raw/genbank.ndjson
        - Each line is one normalized genome record

### 4. What happens next (current state)

    ## At this stage, the pipeline can:

        - Load genome records
        - Compare sequences to a reference
        - Identify mutations
        - Map mutations to genes
        - Compute an explainable risk score

        ## End-to-end scoring scripts will be added next.

### 5. Common issues

    ## NCBI_EMAIL not set
    - The script will refuse to run
    - Fix by setting the environment variable

    ## Slow downloads
    - GenBank rate limits requests
    - This is intentional and respectful of NCBI resources

### 6. Where outputs go

    - data/raw/ → downloaded source data
    - data/derived/ → processed/scored outputs (future)
    - These directories are ignored by git.


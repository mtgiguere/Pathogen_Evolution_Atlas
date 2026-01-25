
---

# 🧠 `docs/ARCHITECTURE.md`

```md
# Architecture Overview

This document explains how the system is structured
and why it is designed this way.

---

## Design philosophy

The system is built around three ideas:

1. **Biology first** — code reflects biological concepts
2. **Canonical data** — all sources normalize to the same structure
3. **Explainability** — every score can be explained

---

## High-level data flow


---

# 🧠 `docs/ARCHITECTURE.md`

```md
# Architecture Overview

This document explains how the system is structured
and why it is designed this way.

---

## Design philosophy

The system is built around three ideas:

1. **Biology first** — code reflects biological concepts
2. **Canonical data** — all sources normalize to the same structure
3. **Explainability** — every score can be explained

---

## High-level data flow

Public genome data
→ Ingest + normalize
→ Canonical genome record
→ Compare to reference genome
→ Identify mutations
→ Map mutations to genes
→ Compute risk score


---

## Canonical genome records

Different data sources store genome data differently.
We convert all sources into a **canonical record** with fields like:

- accession
- organism
- collection date
- host
- sequence
- source

This ensures downstream analysis does **not** depend on the data source.

---

## Mutation analysis

Each genome sequence is compared to a reference genome.

For each position:
- If nucleotides differ → mutation
- If base is ambiguous (`N`) → ignored

Positions use **1-based indexing**, which matches biological conventions.

---

## Gene mapping

Each mutation position is mapped to a known gene
using reference genome coordinates.

Currently implemented:
- ORF1ab
- Spike (S)
- Nucleocapsid (N)

---

## Risk scoring

Risk is computed using a simple, explainable model:

- Mutations are counted by gene
- Each gene has a weight
- The total score is the weighted sum

This model is intentionally simple and interpretable.
More advanced models may be added later.

---

## Extending the system

To add a new data source:
1. Write a source-specific ingest module
2. Normalize into canonical records
3. Reuse all downstream logic unchanged

This is how Nextstrain will be integrated.

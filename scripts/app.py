import os
from pathlib import Path

import pandas as pd
import streamlit as st

from ingest.analytics import summarize_genomes
from ingest.io import load_ndjson

EMAIL = os.getenv("NCBI_EMAIL", "you@domain.com")

_REF_PATH = Path(os.getenv("REF_PATH", "data/reference/genbank_reference.ndjson"))
_GENOMES_PATH = Path(os.getenv("GENOMES_PATH", "data/raw/genomes.ndjson"))

st.set_page_config(page_title="Pathogen Evolution Atlas", layout="wide")
st.title("🧬 Pathogen Evolution Atlas")

# --- Load data ---
for p in (_REF_PATH, _GENOMES_PATH):
    if not p.exists():
        st.error(
            f"Data file not found: {p}. Set REF_PATH / GENOMES_PATH env vars or run the ingest scripts first."
        )
        st.stop()

@st.cache_data
def load_data(ref_path: Path, genomes_path: Path) -> pd.DataFrame:
    ref_rec = next(iter(load_ndjson(ref_path)))
    ref_seq = ref_rec["sequence"] if isinstance(ref_rec, dict) else ref_rec.sequence
    ref_acc = ref_rec["accession"] if isinstance(ref_rec, dict) else ref_rec.accession
    records = list(load_ndjson(genomes_path))
    return summarize_genomes(records, reference_sequence=ref_seq, reference_accession=ref_acc)


df = load_data(_REF_PATH, _GENOMES_PATH)

# Make sure "date" behaves like a date for charts
if "date" in df.columns:
    df["date"] = pd.to_datetime(df["date"], errors="coerce")

# --- Sidebar filters ---
st.sidebar.header("Filters")

include_unscored = st.sidebar.toggle("Include unscored", value=False)

risk_level_choices = sorted([x for x in df["risk_level"].dropna().unique()])
if not include_unscored:
    risk_level_choices = [x for x in risk_level_choices if x != "N/A"]

default_levels = list(risk_level_choices)

risk_levels = st.sidebar.multiselect(
    "Risk level",
    options=risk_level_choices,
    default=default_levels,
)

filtered = df[df["risk_level"].isin(risk_levels)].copy()

if not include_unscored and "scorable" in filtered.columns:
    filtered = filtered[filtered["scorable"] == True]  # noqa: E712

# --- KPIs ---
c1, c2, c3 = st.columns(3)
c1.metric("Genomes", int(len(filtered)))

avg_risk = float(filtered["risk_score"].mean()) if len(filtered) else 0.0
c2.metric("Avg Risk", round(avg_risk, 2))

unique_genes = int(filtered["genes_affected"].nunique()) if len(filtered) else 0
c3.metric("Unique Genes", unique_genes)

# --- Quick visuals ---
st.subheader("Overview")

colA, colB = st.columns(2)

with colA:
    st.caption("Risk level distribution")
    if len(filtered):
        counts = (
            filtered["risk_level"]
            .fillna("N/A")
            .value_counts()
            .rename_axis("risk_level")
            .to_frame("count")
        )
        st.bar_chart(counts)
    else:
        st.info("No rows match the current filters.")

with colB:
    st.caption("Average risk over time")
    if len(filtered) and filtered["date"].notna().any():
        trend = (
            filtered.dropna(subset=["date"])
            .groupby(pd.Grouper(key="date", freq="W"))["risk_score"]
            .mean()
            .to_frame("avg_risk_score")
            .sort_index()
        )
        st.line_chart(trend)
    else:
        st.info("No valid dates available for a time trend yet.")

# --- Table ---
st.subheader("Genome Summary")
st.dataframe(filtered, use_container_width=True)

# --- Map (only if lat/lon exist AND non-null) ---
if {"lat", "lon"}.issubset(filtered.columns):
    map_df = filtered.dropna(subset=["lat", "lon"])
    if not map_df.empty:
        st.subheader("Geographic Distribution")
        st.map(map_df[["lat", "lon"]])

# --- Details (the "drawer") ---
st.subheader("Details")
if len(filtered):
    selected = st.selectbox("Select accession", filtered["accession"].tolist())
    row = filtered.loc[filtered["accession"] == selected].iloc[0]

    cL, cR = st.columns([1, 2])

    with cL:
        st.markdown(f"**Accession:** {row['accession']}")
        st.markdown(f"**Risk Level:** {row['risk_level']}")
        st.markdown(f"**Risk Score:** {row['risk_score']}")
        st.markdown(f"**Mutations:** {row.get('num_mutations', 0)}")
        if pd.notna(row.get("date", None)):
            st.markdown(f"**Collection date:** {row['date'].date()}")
        st.markdown(f"**Scorable:** {bool(row.get('scorable', False))}")
        if row.get("skip_reason", ""):
            st.markdown(f"**Skip reason:** `{row['skip_reason']}`")
        if "qc_status" in row:
            st.markdown(f"**QC:** {row.get('qc_status', '')}")
        if row.get("qc_reasons", ""):
            st.markdown(f"**QC reasons:** `{row['qc_reasons']}`")

    with cR:
        st.markdown("**Genes affected**")
        st.code(row.get("genes_affected", ""), language="text")
        st.markdown("**Risk explanation**")
        st.write(row.get("risk_explanation", ""))
else:
    st.info("Nothing to show. Try widening filters or enabling 'Include unscored'.")

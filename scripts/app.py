import streamlit as st
import pandas as pd

from ingest.analytics import summarize_genomes

from ingest.io import load_ndjson

records = load_ndjson("data/raw/genomes.ndjson")
df = summarize_genomes(records)

st.set_page_config(page_title="Pathogen Evolution Atlas", layout="wide")
st.title("🧬 Pathogen Evolution Atlas")

# --- Sidebar filters ---
st.sidebar.header("Filters")
risk_levels = st.sidebar.multiselect(
    "Risk level",
    sorted(df["risk_level"].unique()),
    default=list(df["risk_level"].unique()),
)

filtered = df[df["risk_level"].isin(risk_levels)]

# --- KPIs ---
c1, c2, c3 = st.columns(3)
c1.metric("Genomes", len(filtered))
c2.metric("Avg Risk", round(filtered["risk_score"].mean(), 2))
c3.metric("Unique Genes", filtered["genes_affected"].nunique())

# --- Table ---
st.subheader("Genome Summary")
st.dataframe(filtered, use_container_width=True)

# --- Map ---
if {"lat", "lon"}.issubset(filtered.columns):
    st.subheader("Geographic Distribution")
    map_df = filtered.dropna(subset=["lat", "lon"])
    if not map_df.empty:
        st.map(map_df[["lat", "lon"]])

# --- Explainability ---
st.subheader("Risk Explanation")
selected = st.selectbox("Select accession", filtered["accession"])
row = filtered[filtered["accession"] == selected].iloc[0]
st.markdown(f"**Risk Level:** {row['risk_level']}")
st.markdown(row["risk_explanation"])

import streamlit as st

from ingest.analytics import summarize_genomes
from ingest.io import load_ndjson
from ingest.geography_enrichment import enrich_many_locations  # <-- add

EMAIL = "you@domain.com"  # <-- put your Entrez email here

records = load_ndjson("data/raw/genomes.ndjson")

# --- GEO enrichment impact summary ---
before_country_missing = sum(1 for r in records if not getattr(r, "country", None))
before_region_missing = sum(1 for r in records if not getattr(r, "region", None))

#records = enrich_many_locations(records, email=EMAIL)

after_country_missing = sum(1 for r in records if not getattr(r, "country", None))
after_region_missing = sum(1 for r in records if not getattr(r, "region", None))

total = len(records)
print(
    f"[geo] country missing: {before_country_missing}/{total} -> {after_country_missing}/{total} "
    f"(filled {before_country_missing - after_country_missing})"
)
print(
    f"[geo] region missing: {before_region_missing}/{total} -> {after_region_missing}/{total} "
    f"(filled {before_region_missing - after_region_missing})"
)

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

import json

# Load the mutation calls we saved from call_all.py (the pipeline handoff).
with open("mutation_calls.json") as f:
    genome_mutations = json.load(f)

# We also need each genome's collection date. That lives in the original
# ndjson, so read it again and build a simple accession -> date lookup.
dates = {}
with open("data/raw/genomes.ndjson") as f:
    for line in f:
        record = json.loads(line)
        dates[record["accession"]] = record["collection_date"]

# Turn a date like "2025-04-17" into a quarter label like "2025-Q2".
def to_quarter(date_str):
    year = date_str[:4]
    month = int(date_str[5:7])
    quarter = (month - 1) // 3 + 1     # months 1-3 -> Q1, 4-6 -> Q2, etc.
    return f"{year}-Q{quarter}"

# Bin every genome into its quarter, and count how many genomes per quarter.
from collections import Counter
quarter_counts = Counter()
for accession in genome_mutations:
    q = to_quarter(dates[accession])
    quarter_counts[q] += 1

print("=== Genomes per quarter ===")
for q in sorted(quarter_counts):
    print(f" {q}: {quarter_counts[q]} genomes")
# The three quarters with enough genomes to trust.
GOOD_QUARTERS = ["2025-Q3", "2025-Q4", "2026-Q1"]

# Group genome accessions by quarter (only the good ones).
genomes_by_quarter = {q: [] for q in GOOD_QUARTERS}
for accession in genome_mutations:
    q = to_quarter(dates[accession])
    if q in genomes_by_quarter:
        genomes_by_quarter[q].append(accession)

# Helper: within one quarter, what fraction of genomes carry a given mutation?
# Each mutation was saved as a [pos, ref, alt] list, so we compare as lists.
def frequency_in_quarter(mutation, accessions):
    carriers = 0
    for acc in accessions:
        if list(mutation) in genome_mutations[acc]:
            carriers += 1
    return 100 * carriers / len(accessions)

# The middle-band mutations we want to track through time.
watch_list = [
    [22995, "C", "A"], [23277, "C", "T"], [23021, "A", "G"],
    [20178, "C", "T"], [13608, "C", "T"], [10615, "C", "T"], [7113, "C", "T"],
]

print()
print("=== Mutation frequency by quarter (Q3'25 -> Q4'25 -> Q1'26) ===")
print(f"{'mutation':>12}    {'2025-Q3':>8} {'2025-Q4':>8} {'2026-Q1':>8}    trend")
for mut in watch_list:
    pos, ref, alt = mut
    freqs = [frequency_in_quarter(mut, genomes_by_quarter[q]) for q in GOOD_QUARTERS]
    arrow = "rising" if freqs[-1] > freqs[0] + 10 else "falling" if freqs[-1] < freqs[0] - 10 else "flat"
    name = f"{ref}{pos}{alt}"
    print(f"{name:>12}  {freqs[0]:7.0f}% {freqs[1]:7.0f}% {freqs[2]:7.0f}%  {arrow}")

## Load each genome's lineage from the Nextclade results (a tab-separated file).
# We collapse sub-lineages to their family: "XFG.3.16.1" -> "XFG".
lineage = {}
with open("nextclade_output/nextclade.tsv") as f:
    header = f.readline().split("\t")
    name_col = header.index("seqName")
    pango_col = header.index("Nextclade_pango")
    for line in f:
        cols = line.split("\t")
        full_lineage = cols[pango_col].strip()
        family = full_lineage.split(".")[0] # part before the first dot
        lineage[cols[name_col].strip()] = family

# Quick check: what families do we have, and how many genomes each?
# We collapse sub-lineages to their family: "XFG.3.16.1" -> "XFG".
lineage = {}
with open("nextclade_output/nextclade.tsv") as f:
    header = f.readline().split("\t")
    name_col = header.index("seqName")
    pango_col = header.index("Nextclade_pango")
    for line in f:
        cols = line.split("\t")
        full_lineage = cols[pango_col].strip()
        family = full_lineage.split(".")[0] # part before the first dot
        lineage[cols[name_col].strip()] = family
# Quick check: what families do we have, and how many genomes each?
family_counts = Counter(lineage.values())
print()
print("=== Lineage families in the dataset ===")
for fam, n in family_counts.most_common():
    print(f" {fam}: {n} genomes")

# Track each lineage family's share within each trustworthy quarter.
print()
print("=== Lineage family share by quarter (Q3'25 -> Q4'25 -> Q1'26) ===")
print(f"{'family':>8}   {'2025-Q3':>8}  {'2025-Q4':>8}  {'2026-Q1':>8}  trend")

for fam in ["XFG", "NB", "PQ"]:
    shares = []
    for q in GOOD_QUARTERS:
        genomes_this_q = genomes_by_quarter[q]
        carriers = sum(1 for acc in genomes_this_q if lineage.get(acc) == fam)
        shares.append(100 * carriers / len(genomes_this_q))
    arrow = "rising" if shares[-1] > shares[0] + 10 else "falling" if shares[-1] < shares[0] - 10 else "flat"
    print(f"{fam:>8}    {shares[0]:7.0f}% {shares[1]:7.0f}% {shares[2]:7.0f}%   {arrow}")
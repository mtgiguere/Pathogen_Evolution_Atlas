import csv
from collections import Counter

# Read Nextclade's results. It's a TSV (tab-separated), and we look columns
# up by name rather than position (the file has 60+ columns).
xfg_aa_changes = []    # amino-acid change lists, one per XFG genome
xfg_count = 0

with open("nextclade_output/nextclade.tsv") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        family = row["Nextclade_pango"].split(".")[0]
        if family != "XFG":
            continue
        xfg_count += 1
        # aaSubstitutions is a comma-separated string like "S:F456L, N:R203K".
        # Split it into individual changes; skip genomes with none listed.
        aa_field = row["aaSubstitutions"].strip()
        if aa_field:
            xfg_aa_changes.append(aa_field.split(","))

print(f"XFG genomes: {xfg_count}")
print(f"Example amino-acid changes form one genome:")
print(" ", xfg_aa_changes[0][:12], "...")

# Count how often each SPIKE amino-acid change appears across XFG genomes.
# Spike changes are the ones whose gene prefix is "S".
spike_counter = Counter()
for changes in xfg_aa_changes:
    for change in changes:
        if change.startswith("S:"):
            spike_counter[change] += 1

print()
print(f"=== Spike amino-acid changes across {len(xfg_aa_changes)} XFG genomes ===")
print(f"(total distinct spike changes seen: {len(spike_counter)})")
print()
print("Near-universal in XFG (the defining spike signature):")
for change, n in spike_counter.most_common():
    pct = 100 * n / len(xfg_aa_changes)
    if pct >= 90:
        print(f"    {change:<10} {pct:.0f}%")
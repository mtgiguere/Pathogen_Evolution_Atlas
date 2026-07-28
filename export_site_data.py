import csv
import json
from collections import defaultdict

def changes_at_positions(family_name, positions):
    """For one family, collect every spike change at the given position.
    Returns (genome_count, {position: {change: count}})."""
    n = 0
    found = {pos: defaultdict(int) for pos in positions}
    with open("nextclade_output/nextclade.tsv") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row["Nextclade_pango"].split(".")[0] != family_name:
                continue
            n += 1
            aa = row["aaSubstitutions"].strip()
            if not aa:
                continue
            for change in aa.split(","):
                if not change.startswith("S:"):
                    continue
                pos = int(change[2:][1:-1])
                if pos in positions:
                    found[pos][change] += 1
    return n, found

# --- Build the data structure for the website ---
positions = {445, 478, 487, 346, 444}
families = ["XFG", "PQ", "NB"]

# We want, per positoin, per family: which changes and what % of that family
export = {}
for pos in sorted(positions):
    export[pos] = {}
    for fam in families:
        n, found = changes_at_positions(fam, {pos})
        changes = found[pos]
        # convert raw counts to a list of {change, count, pct}
        export[pos][fam] = {
            "n": n,
            "changes": [
                {"change": chg, "count": cnt, "pct": round(100 * cnt / n, 1)}
                for chg, cnt in changes.items()
            ],
        }

# --- Write it to JSON ---
with open("data/site/residues.json", "w") as out:
    json.dump(export, out, indent=2)

print("Wrote data/site/residues.json")
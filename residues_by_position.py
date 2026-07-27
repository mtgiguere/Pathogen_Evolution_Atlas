import csv 
from collections import defaultdict 

def changes_at_positions(family_name, positions):
    """For one family, collect every spike change whose position is in 'positions'.
    Returns: {position: Counter of changes} and the genome count."""
    n = 0
    # position -> {change_string: how many genomes have it}
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
                pos = int(change[2:][1:-1])  # "S:T478K" -> 478
                if pos in positions:
                    found[pos][change] += 1
    return n, found

positions = {445, 478, 487}

for fam in ["XFG", "PQ", "NB"]:
    n, found = changes_at_positions(fam, positions)
    print(f"\n=== {fam} (n={n}) ===")
    for pos in sorted(positions):
        changes = found[pos]
        if not changes:
            print(f" position {pos}: (no substitution called - matches Wuhan reference)")
        else:
            parts = [f"{chg} in {cnt}/{n}" for chg, cnt in changes.items()]
            print(f" position {pos}: " + "; ".join(parts))
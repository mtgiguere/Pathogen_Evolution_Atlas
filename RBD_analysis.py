import csv
from collections import Counter

# --- Rebuild the spike counter (same logic as functional_analysis.py) ---
xfg_aa_changes = []
xfg_count = 0
with open("nextclade_output/nextclade.tsv") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        family = row["Nextclade_pango"].split(".")[0]
        if family != "XFG":
            continue
        xfg_count += 1
        aa_field = row["aaSubstitutions"].strip()
        if aa_field:
            xfg_aa_changes.append(aa_field.split(","))

spike_counter = Counter()
for changes in xfg_aa_changes:
    for change in changes:
        if change.startswith("S:"):
            spike_counter[change] += 1

print(f"XFG genomes: {xfg_count}")
print(f"Distinct spike changes: {len(spike_counter)}")

# --- Filter to RBD positions and sort ---
RBD_START = 319
RBD_END = 541

rbd_hits = []
for change, n in spike_counter.items():
    core = change[2:]       # drop "S:", e.g. "L455S"
    pos = int(core[1:-1])   # pull the number, e.g. 455
    if RBD_START <= pos <+ RBD_END:
        pct = 100 * n / len(xfg_aa_changes)
        rbd_hits.append((pos, change, pct))

# Sort by position so the cluster reads left-to-right along the spike
rbd_hits.sort()

# --- Mark which positions fall in the receptor-binding motif (RBM) ---
# The RBM (~437-508) is the part of the RBD that physically cotacts the 
# ACE2 receptor. It's where most antibody-escape pressure concentrates,
# so it's worth flagging separately from the wider RBD.
RBM_START = 437
RBM_END = 508

print()
print(f"=== RBD spike changes in XFG (positions {RBD_START}-{RBD_END}) ===")
print(f"{len(rbd_hits)} changes found   [*] = in receptor-binding motif ({RBM_START}-{RBM_END})")
print()
for pos, change, pct in rbd_hits:
    marker = "[*]" if RBM_START <= pos <= RBM_END else " "
    print(f" {marker} {change:<10} {pct:.0f}%")
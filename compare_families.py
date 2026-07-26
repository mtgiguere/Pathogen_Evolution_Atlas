import csv
from collections import Counter

def spike_counter_for(family_name):
    aa_changes = []
    genome_count = 0
    with open("nextclade_output/nextclade.tsv") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            family = row["Nextclade_pango"].split(".")[0]
            if family != family_name:       # the blank, not hardcoded "XFG"
                continue
            genome_count += 1
            aa_field = row["aaSubstitutions"].strip()
            if aa_field:
                aa_changes.append(aa_field.split(","))

    spike_counter = Counter()
    for changes in aa_changes:
        for change in changes:
            if change.startswith("S:"):
                spike_counter[change] += 1
    return genome_count, spike_counter
def rbd_signature(genome_count, spike_counter, threshold=90):
    """Return the set of RBD spike changes present in >= threshold% of genomes."""
    signature = set()
    for change, n in spike_counter.items():
        pos = int(change[2:][1:-1])         # "S:L455S" -> 455
        if 319 <= pos <= 541:               # in the RBD
            pct = 100 * n / genome_count
            if pct >= threshold:
                signature.add(change)
    return signature    
xfg_n, xfg_spikes = spike_counter_for("XFG")
pq_n, pq_spikes = spike_counter_for("PQ")

xfg_sig = rbd_signature(xfg_n, xfg_spikes)
pq_sig = rbd_signature(pq_n, pq_spikes)

shared = xfg_sig & pq_sig           # in BOTH families
xfg_only = xfg_sig - pq_sig         # in XFG, NOT in PQ
pq_only = pq_sig - xfg_sig          # in PQ, NOT in XFG

print(f"XFG RBD signature: {len(xfg_sig)} changes")
print(f"PQ  RBD signature: {len(pq_sig)} changes")
print()
print(f"Shared (inherited background): {len(shared)}")
print(f"XFG-only (distinguishes XFG):  {len(xfg_only)}")
print(f"PQ-only:                       {len(pq_only)}")
def show(label, sig):
    print(f"\n{label} ({len(sig)})")
    for change in sorted(sig, key=lambda c: int(c[2:][1:-1])):
        pos = int(change[2:][1:-1])
        marker = "[*]" if 437 <= pos <= 508 else "  "
        print(f" {marker} {change}")

show("XFG-only (distinguishes XFG from PQ)", xfg_only)
show("PQ-only (distinguishes PQ from XFG)", pq_only)
show("Shared background", shared)

print("\n=== Threshold robustness: XFG-only RBD changes ===")
for t in [80, 85, 90, 95]:
    xfg_sig_t = rbd_signature(xfg_n, xfg_spikes, threshold=t)
    pq_sig_t = rbd_signature(pq_n, pq_spikes, threshold=t)
    xfg_only_t = sorted(xfg_sig_t - pq_sig_t, key=lambda c: int(c[2:][1:-1]))
    print(f"\n threshold {t}% -> {len(xfg_only_t)} XFG-only:")
    for change in xfg_only_t:
        print(f"    {change}")
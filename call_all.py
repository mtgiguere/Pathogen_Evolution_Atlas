import json

# Read every genome record and write them all into one multi-FASTA file.
# A multi-FASTA is jst many ">header / sequence" blocks stacked in one file.
count = 0
with open("data/raw/genomes.ndjson") as infile, open("all_samples.fasta", "w") as out:
    for line in infile:
        record = json.loads(line)
        out.write(">" + record["accession"] + "\n")
        out.write(record["sequence"] + "\n")
        count += 1

print(f"Wrote {count} sequences to all_samples.fasta")

import pysam
from collections import Counter

# Load the reference sequence (same as before).
with open("data/reference/NC_045512.2.fasta") as f:
    reference = "".join(line for line in f.read().splitlines() if not line.startswith(">"))

# We'll tally two things:
recurrence = Counter() # how many genomes carry each specific mutation
per_genome = {}        # how many mutations each genome has
genome_mutations = {}  # the actual mutation list for each genome, to save for reuse

alignment = pysam.AlignmentFile("all_samples.sam", "r")
for read in alignment.fetch(until_eof=True):
    # Skip anything that isn't a clean primary alignment.
    if read.is_unmapped or read.is_secondary or read.is_supplementary:
        continue

    mutations = []
    for sample_pos, ref_pos in read.get_aligned_pairs():
        if sample_pos is None or ref_pos is None:
            continue # gap (indel), skip
        ref_base = reference[ref_pos]
        sample_base = read.query_sequence[sample_pos]
        if sample_base == "N":
            continue # missing data, not a mutation
        if ref_base != sample_base:
            mut = (ref_pos + 1, ref_base, sample_base) # 1-based
            mutations.append(mut)
            recurrence[mut] += 1

    per_genome[read.query_name] = len(mutations)
    genome_mutations[read.query_name] = mutations

print("Genomes processed:", len(per_genome))
print("Distinct mutations seen:", len(recurrence))
print()
print("=== Top 20 most recurrent mutations (position, ref->, # of genomes) ===")
for (pos, ref, alt), n in recurrence.most_common(20):
    print(f" {ref}{pos}{alt}  in {n} of {len(per_genome)} genomes")
print()
print("=== Mutations in a MIDDLE band (present in 15% to 85% of genomes) ===")
print("    These are the ones actively dividing the population - the interesting ones.")
low = int(0.15 * len(per_genome))
high = int(0.85 * len(per_genome))
middle = [(m,n) for m, n in recurrence.items() if low <= n <= high]
middle.sort(key=lambda x: x[1], reverse=True)
print(f"    ({len(middle)} mutations fall in this band)")
print()
for (pos, ref, alt), n in middle[:25]:
    pct = 100 * n / len(per_genome)
    print(f" {ref}{pos}{alt}    in {n} of {len(per_genome)} genomes ({pct:.0f}%)")
# Save each genome's mutation calls so the temporal script can reuse them
# without re-aligning. (Pipeline pattern: one stage writes, the enxt reads.)
with open("mutation_calls.json", "w") as f:
    json.dump(genome_mutations, f)
print("Saved per-genome mutation calls to mutation_calls.json")
"""Export mutations from the detection system to JSON for web visualization."""

import json

from src.ingest.mutations import Mutation
from src.ingest.risk import score_mutations


def create_sample_data():
    """Create sample genome sequences to demonstrate mutation detection."""
    mutations_data = [
        {"pos": 1000, "ref": "A", "alt": "G", "gene": "ORF1ab"},
        {"pos": 22000, "ref": "C", "alt": "T", "gene": "S"},
        {"pos": 29000, "ref": "T", "alt": "C", "gene": "N"},
    ]
    return [Mutation(pos=m["pos"], ref=m["ref"], alt=m["alt"], gene=m["gene"]) for m in mutations_data]


def export_mutation_data():
    """Run mutation detection and export results as JSON."""
    mutations = create_sample_data()
    print(f"Created {len(mutations)} sample mutations")

    risk_analysis = score_mutations(mutations)
    print(f"Risk analysis complete: {risk_analysis['level']}")

    mutation_list = [
        {
            "position": mutation.pos,
            "gene": mutation.gene or "Unknown",
            "type": f"{mutation.ref}->{mutation.alt}",
            "ref": mutation.ref,
            "alt": mutation.alt,
        }
        for mutation in mutations
    ]

    export_data = {
        "mutations": mutation_list,
        "risk_analysis": risk_analysis,
        "summary": {
            "total_mutations": len(mutations),
            "genes_affected": list(risk_analysis["by_gene"].keys()),
        },
    }

    with open("mutation_data.json", "w") as f:
        json.dump(export_data, f, indent=2)

    print("Mutation data exported to mutation_data.json")
    print(f"Found {len(mutations)} mutations")
    print(f"Risk level: {risk_analysis['level']}")


if __name__ == "__main__":
    export_mutation_data()

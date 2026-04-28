""" Export mutations from the detection system to JSON for web visualization """


import json
import sys
sys.path.append('../src')

from ingest.mutations import Mutation, diff_sequences
from ingest.risk import score_mutations
from ingest.genes import gene_for_position

def create_sample_data():
    """Create sample genome sequences to demonstrate mutation detection"""
    # Create sample mutations directly
    mutations_data = [
        {'pos': 1000, 'ref': 'A', 'alt': 'G', 'gene': 'ORF1ab'},
        {'pos': 22000, 'ref': 'C', 'alt': 'T', 'gene': 'S'},
        {'pos': 29000, 'ref': 'T', 'alt': 'C', 'gene': 'N'}
    ]

    # Convert to Mutation objects
    mutations = []
    for m in mutations_data:
        mutation = Mutation(pos=m['pos'], ref=m['ref'], alt=m['alt'], gene=m['gene'])
        mutations.append(mutation)

    return mutations    

def export_mutation_data():
    """Run mutation detection and export results as JSON"""
    # Get sample data
    mutations = create_sample_data()
    print(f"Created {len(mutations)} sample mutations")
    # Calculate risk scores using your system
    risk_analysis = score_mutations(mutations)
    print(f"Risk analysis complete: {risk_analysis['level']}")
    # Convert mutations to web-friendly format
    mutation_list = []
    for mutation in mutations:
        mutation_data = {
            'position': mutation.pos,
            'gene': mutation.gene or 'Unknown',
            'type': f'{mutation.ref}->{mutation.alt}',
            'ref': mutation.ref,
            'alt': mutation.alt
        }
        mutation_list.append(mutation_data)


    # Create complete data structure for web visualization
    export_data = {
        'mutations': mutation_list,
        'risk_analysis': risk_analysis,
        'summary': {
            'total_mutations': len(mutations),
            'genes_affected': list(risk_analysis['by_gene'].keys())
    }
}

    # Save to JSON file
    with open('mutation_data.json', 'w') as f:
        json.dump(export_data, f, indent=2)

        print("Mutation data exported to mutation_data.json")
        print(f"Found {len(mutations)} mutations")
        print(f"Risk level: {risk_analysis['level']}")

if __name__ == "__main__":
    export_mutation_data()            
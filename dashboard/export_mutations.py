""" Export mutations from the detection system to JSON for web visualization """


import json
import sys
sys.path.append('../src')

from ingest.mutations import Mutation, diff_sequences
from ingest.risk import score_mutations
from ingest.genes import gene_for_position

def create_sample_data():
    """Create sample genome sequences to demonstrate mutation detection"""
    # ... mutations_data stays the same ...

    
    # Create sample mutations directly
    mutations_data = [
    # ORF1ab mutations (replication complex) - positions 266-21555
    {'pos': 3037, 'ref': 'C', 'alt': 'T', 'gene': 'ORF1ab', 'location': {'country': 'China', 'city': 'Wuhan', 'lat': 30.5928, 'lng': 114.3055, 'collection_date': '2020-01-05'}},
    {'pos': 8782, 'ref': 'C', 'alt': 'T', 'gene': 'ORF1ab', 'location': {'country': 'Italy', 'city': 'Milan', 'lat': 45.4642, 'lng': 9.1900, 'collection_date': '2020-02-20'}},
    {'pos': 11083, 'ref': 'G', 'alt': 'T', 'gene': 'ORF1ab', 'location': {'country': 'Spain', 'city': 'Madrid', 'lat': 40.4168, 'lng': -3.7038, 'collection_date': '2020-03-15'}},
    {'pos': 14408, 'ref': 'C', 'alt': 'T', 'gene': 'ORF1ab', 'location': {'country': 'United States', 'city': 'Seattle', 'lat': 47.6062, 'lng': -122.3321, 'collection_date': '2020-02-29'}},
    {'pos': 15324, 'ref': 'C', 'alt': 'T', 'gene': 'ORF1ab', 'location': {'country': 'France', 'city': 'Paris', 'lat': 48.8566, 'lng': 2.3522, 'collection_date': '2020-03-01'}},

    # Spike protein mutations (positions 21563-25384) - vaccine target
    {'pos': 21762, 'ref': 'C', 'alt': 'T', 'gene': 'S', 'location': {'country': 'South Africa', 'city': 'Pretoria', 'lat': -25.7479, 'lng': 28.2293, 'collection_date': '2021-11-15'}},
    {'pos': 21765, 'ref': 'A', 'alt': 'T', 'gene': 'S', 'location': {'country': 'United Kingdom', 'city': 'London', 'lat': 51.5074, 'lng': -0.1278, 'collection_date': '2020-12-18'}},
    {'pos': 22578, 'ref': 'G', 'alt': 'A', 'gene': 'S', 'location': {'country': 'Botswana', 'city': 'Gaborone', 'lat': -24.6282, 'lng': 25.9231, 'collection_date': '2021-11-18'}},
    {'pos': 22995, 'ref': 'C', 'alt': 'A', 'gene': 'S', 'location': {'country': 'India', 'city': 'Mumbai', 'lat': 19.0760, 'lng': 72.8777, 'collection_date': '2021-02-10'}},
    {'pos': 23063, 'ref': 'A', 'alt': 'T', 'gene': 'S', 'location': {'country': 'United Kingdom', 'city': 'London', 'lat': 51.5074, 'lng': -0.1278, 'collection_date': '2020-12-15'}},
    {'pos': 23271, 'ref': 'A', 'alt': 'T', 'gene': 'S', 'location': {'country': 'Netherlands', 'city': 'Amsterdam', 'lat': 52.3676, 'lng': 4.9041, 'collection_date': '2021-12-01'}},

    # ORF3a mutations (positions 25393-26220) - accessory protein
    {'pos': 25563, 'ref': 'G', 'alt': 'T', 'gene': 'ORF3a', 'location': {'country': 'Brazil', 'city': 'São Paulo', 'lat': -23.5558, 'lng': -46.6396, 'collection_date': '2020-04-15'}},
    {'pos': 25904, 'ref': 'C', 'alt': 'T', 'gene': 'ORF3a', 'location': {'country': 'Australia', 'city': 'Sydney', 'lat': -33.8688, 'lng': 151.2093, 'collection_date': '2020-03-20'}},

    # Envelope protein mutations (positions 26245-26472)
    {'pos': 26340, 'ref': 'C', 'alt': 'T', 'gene': 'E', 'location': {'country': 'South Korea', 'city': 'Seoul', 'lat': 37.5665, 'lng': 126.9780, 'collection_date': '2020-02-18'}},

    # Membrane protein mutations (positions 26523-27191)
    {'pos': 26735, 'ref': 'C', 'alt': 'T', 'gene': 'M', 'location': {'country': 'Canada', 'city': 'Toronto', 'lat': 43.6532, 'lng': -79.3832, 'collection_date': '2020-03-12'}},

    # ORF6 mutations (positions 27202-27387)
    {'pos': 27259, 'ref': 'A', 'alt': 'T', 'gene': 'ORF6', 'location': {'country': 'Japan', 'city': 'Tokyo', 'lat': 35.6762, 'lng': 139.6503, 'collection_date': '2020-02-25'}},

    # ORF7a mutations (positions 27394-27759)
    {'pos': 27638, 'ref': 'T', 'alt': 'C', 'gene': 'ORF7a', 'location': {'country': 'Singapore', 'city': 'Singapore', 'lat': 1.3521, 'lng': 103.8198, 'collection_date': '2020-02-08'}},

    # ORF8 mutations (positions 27894-28259)
    {'pos': 28144, 'ref': 'T', 'alt': 'C', 'gene': 'ORF8', 'location': {'country': 'Thailand', 'city': 'Bangkok', 'lat': 13.7563, 'lng': 100.5018, 'collection_date': '2020-01-25'}},

    # Nucleocapsid mutations (positions 28274-29533)
    {'pos': 28311, 'ref': 'C', 'alt': 'T', 'gene': 'N', 'location': {'country': 'Germany', 'city': 'Berlin', 'lat': 52.5200, 'lng': 13.4050, 'collection_date': '2020-03-08'}},
    {'pos': 28881, 'ref': 'G', 'alt': 'A', 'gene': 'N', 'location': {'country': 'Iran', 'city': 'Tehran', 'lat': 35.6892, 'lng': 51.3890, 'collection_date': '2020-02-28'}},
    {'pos': 29000, 'ref': 'T', 'alt': 'C', 'gene': 'N', 'location': {'country': 'Netherlands', 'city': 'Amsterdam', 'lat': 52.3676, 'lng': 4.9041, 'collection_date': '2020-03-18'}},

    # ORF10 mutations (positions 29558-29674)  
    {'pos': 29645, 'ref': 'G', 'alt': 'T', 'gene': 'ORF10', 'location': {'country': 'Mexico', 'city': 'Mexico City', 'lat': 19.4326, 'lng': -99.1332, 'collection_date': '2020-03-25'}},
]

    # Convert to Mutation objects and store location data separately
    mutations = []
    locations = []
    for m in mutations_data:
        mutation = Mutation(pos=m['pos'], ref=m['ref'], alt=m['alt'], gene=m['gene'])
        mutations.append(mutation)


        # Store location data with mutation position as key
        location_data = m['location'].copy()
        location_data['position'] = m['pos']  # Link location to mutation position
        locations.append(location_data)

    return mutations, locations    

def export_mutation_data():
    """Run mutation detection and export results as JSON"""
    # Get sample data
    mutations, locations = create_sample_data()
    print(f"Debug - locations: {locations}") # ADD THIS LINE
    print(f"Created {len(mutations)} sample mutations")
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
        'locations': locations, # ADD THIS LINE
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
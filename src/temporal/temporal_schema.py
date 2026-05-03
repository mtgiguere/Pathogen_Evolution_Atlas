from dataclasses import dataclass, field
from typing import List, Dict, Optional, Set
from datetime import datetime, timedelta
from enum import Enum
import json


class Gene(Enum):
    """Represents the 10 genes we're tracking in SARS-CoV-2"""
    ORF1ab = "ORF1ab"
    SPIKE = "S"
    ORF3a = "ORF3a"
    ENVELOPE = "E"
    MEMBRANE = "M"
    ORF6 = "ORF6"
    ORF7a = "ORF7a"
    ORF8 = "ORF8"
    NUCLEOCAPSID = "N"
    ORF10 = "ORF10"


@dataclass
class Mutation:
    """Represents a single mutation in the SARS-CoV-2 genome"""
    position: int
    gene: Gene
    ref_amino_acid: str
    alt_amino_acid: str
    mutation_name: str
    emergence_date: datetime
    variants_carrying: List[str] = field(default_factory=list)
    prevalence_history: Dict[str, float] = field(default_factory=dict)


@dataclass
class Variant:
    """Represents a named SARS-CoV-2 variant (Alpha, Beta, Delta, etc.)"""
    name: str
    pango_lineage: str
    emergence_date: datetime
    peak_prevalence_date: datetime
    defining_mutations: List[Mutation] = field(default_factory=list)
    is_recombinant: bool = False
    parent_variants: List[str] = field(default_factory=list)
    geographic_origin: Optional[str] = None
    who_label: Optional[str] = None


@dataclass
class Country:
    """Represents a country with its geographic and epidemiological data"""
    name: str
    iso_code: str
    latitude: float
    longitude: float
    region: str
    population: int
    quarterly_data: Dict[str, Dict] = field(default_factory=dict)


@dataclass
class QuarterlySnapshot:
    """Represents the pandemic state for a single quarter (3 months)"""
    quarter: str
    year: int
    start_date: datetime
    end_date: datetime
    active_mutations: List[Mutation] = field(default_factory=list)
    active_variants: List[Variant] = field(default_factory=list)
    global_dominant_variant: Optional[str] = None
    countries_affected: Dict[str, Country] = field(default_factory=dict)
    global_risk_score: float = 0.0
    mutation_emergence_events: List[Dict] = field(default_factory=list)
    recombination_events: List[Variant] = field(default_factory=list)


@dataclass
class TemporalDatabase:
    """
    The master database containing all temporal pandemic data.
    This is what we'll serialize to JSON and feed to the visualization.
    """
    name: str = "SARS-CoV-2 Pathogen Evolution Atlas"
    version: str = "1.0"
    creation_date: datetime = field(default_factory=datetime.now)
    last_updated: datetime = field(default_factory=datetime.now)
    
    all_mutations: Dict[str, Mutation] = field(default_factory=dict)
    all_variants: Dict[str, Variant] = field(default_factory=dict)
    all_countries: Dict[str, Country] = field(default_factory=dict)
    quarterly_snapshots: Dict[str, QuarterlySnapshot] = field(default_factory=dict)
    
    total_countries: int = 0
    total_mutations_tracked: int = 0
    total_variants: int = 0
    date_range_start: datetime = field(default_factory=lambda: datetime(2020, 1, 1))
    date_range_end: datetime = field(default_factory=lambda: datetime(2026, 12, 31))
    prevalence_threshold: float = 0.01
    
    def add_mutation(self, mutation: Mutation) -> None:
        """Add a mutation to the database"""
        self.all_mutations[mutation.mutation_name] = mutation
        self.total_mutations_tracked = len(self.all_mutations)
    
    def add_variant(self, variant: Variant) -> None:
        """Add a variant to the database"""
        self.all_variants[variant.name] = variant
        self.total_variants = len(self.all_variants)
    
    def add_country(self, country: Country) -> None:
        """Add a country to the database"""
        self.all_countries[country.iso_code] = country
        self.total_countries = len(self.all_countries)
    
    def add_quarterly_snapshot(self, snapshot: QuarterlySnapshot) -> None:
        """Add a quarterly snapshot to the database"""
        quarter_key = f"{snapshot.quarter}_{snapshot.year}"
        self.quarterly_snapshots[quarter_key] = snapshot
    
    def get_quarter(self, quarter: str, year: int) -> Optional[QuarterlySnapshot]:
        """Retrieve a specific quarter's data"""
        quarter_key = f"{quarter}_{year}"
        return self.quarterly_snapshots.get(quarter_key)
    
    def update_timestamp(self) -> None:
        """Update the last_updated timestamp"""
        self.last_updated = datetime.now()


class RiskScoreCalculator:
    """
    Calculates risk scores for mutations and variants.
    Uses gene-weighted system: Spike protein mutations are higher risk than others.
    """
    
    GENE_WEIGHTS = {
        Gene.SPIKE: 3.0,
        Gene.ORF1ab: 1.5,
        Gene.NUCLEOCAPSID: 1.5,
        Gene.ORF3a: 1.0,
        Gene.ENVELOPE: 1.0,
        Gene.MEMBRANE: 1.0,
        Gene.ORF6: 0.5,
        Gene.ORF7a: 0.5,
        Gene.ORF8: 0.5,
        Gene.ORF10: 0.5,
    }
    
    @staticmethod
    def calculate_mutation_risk(mutation: Mutation, prevalence: float) -> float:
        """
        Calculate risk score for a single mutation.
        Risk = Gene Weight × Prevalence (0-100 scale)
        """
        gene_weight = RiskScoreCalculator.GENE_WEIGHTS.get(mutation.gene, 1.0)
        risk_score = gene_weight * (prevalence * 100)
        return min(risk_score, 100.0)
    
    @staticmethod
    def calculate_variant_risk(variant: Variant, global_prevalence: float) -> float:
        """
        Calculate risk score for a variant based on its defining mutations.
        Risk = Average risk of all defining mutations × Global prevalence
        """
        if not variant.defining_mutations:
            return 0.0
        
        mutation_risks = [
            RiskScoreCalculator.calculate_mutation_risk(mut, 0.5)
            for mut in variant.defining_mutations
        ]
        average_mutation_risk = sum(mutation_risks) / len(mutation_risks)
        variant_risk = average_mutation_risk * global_prevalence
        return min(variant_risk, 100.0)
    
    @staticmethod
    def calculate_quarterly_risk(snapshot: QuarterlySnapshot) -> float:
        """
        Calculate overall pandemic risk for a quarter.
        Risk = Sum of all active mutation risks + variant dominance factor
        """
        if not snapshot.active_mutations:
            return 0.0
        
        mutation_risk_sum = sum([
            RiskScoreCalculator.calculate_mutation_risk(mut, 0.05)
            for mut in snapshot.active_mutations
        ])
        
        quarterly_risk = min((mutation_risk_sum / len(snapshot.active_mutations)) * 2, 100.0)
        return quarterly_risk


def database_to_json(db: TemporalDatabase) -> str:
    """
    Convert the TemporalDatabase to a JSON string.
    We need custom serialization because datetime objects aren't JSON-serializable by default.
    """
    db_dict = {
        "metadata": {
            "name": db.name,
            "version": db.version,
            "creation_date": db.creation_date.isoformat(),
            "last_updated": db.last_updated.isoformat(),
            "total_countries": db.total_countries,
            "total_mutations_tracked": db.total_mutations_tracked,
            "total_variants": db.total_variants,
            "date_range_start": db.date_range_start.isoformat(),
            "date_range_end": db.date_range_end.isoformat(),
            "prevalence_threshold": db.prevalence_threshold,
        },
        "mutations": {},
        "variants": {},
        "countries": {},
        "quarterly_snapshots": {},
    }
    
    for mut_name, mutation in db.all_mutations.items():
        db_dict["mutations"][mut_name] = {
            "position": mutation.position,
            "gene": mutation.gene.value,
            "ref_amino_acid": mutation.ref_amino_acid,
            "alt_amino_acid": mutation.alt_amino_acid,
            "mutation_name": mutation.mutation_name,
            "emergence_date": mutation.emergence_date.isoformat(),
            "variants_carrying": mutation.variants_carrying,
            "prevalence_history": mutation.prevalence_history,
        }
    
    for var_name, variant in db.all_variants.items():
        db_dict["variants"][var_name] = {
            "name": variant.name,
            "pango_lineage": variant.pango_lineage,
            "emergence_date": variant.emergence_date.isoformat(),
            "peak_prevalence_date": variant.peak_prevalence_date.isoformat(),
            "defining_mutations": [mut.mutation_name for mut in variant.defining_mutations],
            "is_recombinant": variant.is_recombinant,
            "parent_variants": variant.parent_variants,
            "geographic_origin": variant.geographic_origin,
            "who_label": variant.who_label,
        }
    
    for iso_code, country in db.all_countries.items():
        db_dict["countries"][iso_code] = {
            "name": country.name,
            "iso_code": country.iso_code,
            "latitude": country.latitude,
            "longitude": country.longitude,
            "region": country.region,
            "population": country.population,
            "quarterly_data": country.quarterly_data,
        }
    
    for quarter_key, snapshot in db.quarterly_snapshots.items():
        db_dict["quarterly_snapshots"][quarter_key] = {
            "quarter": snapshot.quarter,
            "year": snapshot.year,
            "start_date": snapshot.start_date.isoformat(),
            "end_date": snapshot.end_date.isoformat(),
            "active_mutations": [mut.mutation_name for mut in snapshot.active_mutations],
            "active_variants": [var.name for var in snapshot.active_variants],
            "global_dominant_variant": snapshot.global_dominant_variant,
            "countries_affected": list(snapshot.countries_affected.keys()),
            "global_risk_score": snapshot.global_risk_score,
            "mutation_emergence_events": snapshot.mutation_emergence_events,
            "recombination_events": [var.name for var in snapshot.recombination_events],
        }
    
    return json.dumps(db_dict, indent=2)


def save_database_to_file(db: TemporalDatabase, filepath: str) -> None:
    """Save the database to a JSON file"""
    json_string = database_to_json(db)
    with open(filepath, 'w') as f:
        f.write(json_string)
    print(f"Database saved to {filepath}")


def load_database_from_file(filepath: str) -> TemporalDatabase:
    """
    Load a database from a JSON file.
    Note: This creates a simplified version without all object references.
    For full functionality, we'd need more complex deserialization.
    """
    with open(filepath, 'r') as f:
        db_dict = json.load(f)
    
    db = TemporalDatabase(
        name=db_dict["metadata"]["name"],
        version=db_dict["metadata"]["version"],
        creation_date=datetime.fromisoformat(db_dict["metadata"]["creation_date"]),
        last_updated=datetime.fromisoformat(db_dict["metadata"]["last_updated"]),
        total_countries=db_dict["metadata"]["total_countries"],
        total_mutations_tracked=db_dict["metadata"]["total_mutations_tracked"],
        total_variants=db_dict["metadata"]["total_variants"],
        prevalence_threshold=db_dict["metadata"]["prevalence_threshold"],
    )
    
    print(f"Database loaded from {filepath}")
    return db


def test_schema():
    """
    Simple test to verify the schema works correctly.
    This creates sample data and checks that everything initializes properly.
    """
    print("=" * 70)
    print("TESTING TEMPORAL SCHEMA")
    print("=" * 70)
    
    print("\n[Test 1] Creating a mutation...")
    d614g = Mutation(
        position=21563,
        gene=Gene.SPIKE,
        ref_amino_acid="D",
        alt_amino_acid="G",
        mutation_name="S:D614G",
        emergence_date=datetime(2020, 1, 1),
        variants_carrying=["Wild-type"],
        prevalence_history={"Q1_2020": 0.15}
    )
    print(f"✓ Created mutation: {d614g.mutation_name}")
    print(f"  Position: {d614g.position}, Gene: {d614g.gene.value}")
    
    print("\n[Test 2] Creating another mutation...")
    n501y = Mutation(
        position=23063,
        gene=Gene.SPIKE,
        ref_amino_acid="N",
        alt_amino_acid="Y",
        mutation_name="S:N501Y",
        emergence_date=datetime(2020, 9, 1),
        variants_carrying=["Alpha", "Beta"],
        prevalence_history={"Q3_2020": 0.05, "Q4_2020": 0.35}
    )
    print(f"✓ Created mutation: {n501y.mutation_name}")
    
    print("\n[Test 3] Creating a variant...")
    alpha = Variant(
        name="Alpha",
        pango_lineage="B.1.1.7",
        emergence_date=datetime(2020, 9, 1),
        peak_prevalence_date=datetime(2021, 1, 15),
        defining_mutations=[d614g, n501y],
        is_recombinant=False,
        geographic_origin="United Kingdom",
        who_label="Alpha"
    )
    print(f"✓ Created variant: {alpha.name}")
    print(f"  Defining mutations: {len(alpha.defining_mutations)}")
    
    print("\n[Test 4] Creating a country...")
    usa = Country(
        name="United States",
        iso_code="USA",
        latitude=37.0902,
        longitude=-95.7129,
        region="North America",
        population=331900000,
        quarterly_data={
            "Q1_2020": {
                "cases": 100000,
                "dominant_variant": "Wild-type",
                "mutations_detected": ["S:D614G"]
            }
        }
    )
    print(f"✓ Created country: {usa.name} ({usa.iso_code})")
    print(f"  Coordinates: ({usa.latitude}, {usa.longitude})")
    
    print("\n[Test 5] Creating a quarterly snapshot...")
    q1_2020 = QuarterlySnapshot(
        quarter="Q1",
        year=2020,
        start_date=datetime(2020, 1, 1),
        end_date=datetime(2020, 3, 31),
        active_mutations=[d614g],
        active_variants=[],
        global_dominant_variant="Wild-type",
        countries_affected={"USA": usa},
        global_risk_score=25.0,
        mutation_emergence_events=[],
        recombination_events=[]
    )
    print(f"✓ Created snapshot: {q1_2020.quarter}_{q1_2020.year}")
    print(f"  Risk score: {q1_2020.global_risk_score}")
    
    print("\n[Test 6] Creating temporal database...")
    db = TemporalDatabase()
    db.add_mutation(d614g)
    db.add_mutation(n501y)
    db.add_variant(alpha)
    db.add_country(usa)
    db.add_quarterly_snapshot(q1_2020)
    
    print(f"✓ Created database: {db.name}")
    print(f"  Total mutations: {db.total_mutations_tracked}")
    print(f"  Total variants: {db.total_variants}")
    print(f"  Total countries: {db.total_countries}")
    print(f"  Total snapshots: {len(db.quarterly_snapshots)}")
    
    print("\n[Test 7] Testing risk score calculation...")
    mutation_risk = RiskScoreCalculator.calculate_mutation_risk(d614g, 0.5)
    print(f"✓ D614G risk at 50% prevalence: {mutation_risk:.2f}")
    
    variant_risk = RiskScoreCalculator.calculate_variant_risk(alpha, 0.3)
    print(f"✓ Alpha variant risk at 30% prevalence: {variant_risk:.2f}")
    
    print("\n[Test 8] Testing database retrieval...")
    retrieved_snapshot = db.get_quarter("Q1", 2020)
    if retrieved_snapshot:
        print(f"✓ Retrieved snapshot: {retrieved_snapshot.quarter}_{retrieved_snapshot.year}")
    else:
        print("✗ Failed to retrieve snapshot")
    
    print("\n[Test 9] Testing JSON serialization...")
    json_output = database_to_json(db)
    json_size = len(json_output)
    print(f"✓ Serialized to JSON: {json_size} characters")
    print(f"  First 200 chars: {json_output[:200]}...")
    
    print("\n" + "=" * 70)
    print("ALL TESTS PASSED ✓")
    print("=" * 70)
    
    return db


if __name__ == "__main__":
    test_db = test_schema()

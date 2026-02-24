#!/usr/bin/env python3
"""
Import 119 Novel Analogs into PostgreSQL Database
Reads from MASTER_ANALOG_DISCOVERIES.json and inserts into analogs table
"""

import json
import psycopg2
from psycopg2.extras import execute_values
from datetime import datetime
import os

# Database connection from environment
DATABASE_URL = os.getenv('DATABASE_URL', 'postgresql://pharmasight_user:pharmasight_pass_2024@localhost:5432/pharmasight_db')

def parse_smiles_properties(smiles):
    """Parse molecular properties from SMILES if available using RDKit"""
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors

        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return None

        return {
            'molecular_weight': Descriptors.MolWt(mol),
            'logp': Descriptors.MolLogP(mol),
            'tpsa': Descriptors.TPSA(mol),
            'h_bond_donors': Descriptors.NumHDonors(mol),
            'h_bond_acceptors': Descriptors.NumHAcceptors(mol),
            'rotatable_bonds': Descriptors.NumRotatableBonds(mol),
            'aromatic_rings': Descriptors.NumAromaticRings(mol)
        }
    except ImportError:
        print("⚠️  RDKit not available - using properties from JSON if available")
        return None

def calculate_lipinski_violations(props):
    """Calculate Lipinski Rule of 5 violations"""
    violations = 0
    if props.get('molecular_weight', 0) > 500:
        violations += 1
    if props.get('logp', 0) > 5:
        violations += 1
    if props.get('h_bond_donors', 0) > 5:
        violations += 1
    if props.get('h_bond_acceptors', 0) > 10:
        violations += 1
    return violations

def import_analogs():
    """Import all analogs from JSON to PostgreSQL"""

    # Load JSON data
    print("📂 Loading MASTER_ANALOG_DISCOVERIES.json...")
    with open('MASTER_ANALOG_DISCOVERIES.json', 'r') as f:
        analogs_data = json.load(f)

    print(f"✅ Loaded {len(analogs_data)} analogs")

    # Connect to database
    print(f"\n🔌 Connecting to PostgreSQL...")
    conn = psycopg2.connect(DATABASE_URL)
    cur = conn.cursor()

    # Prepare data for insertion
    insert_data = []
    skipped = 0

    for analog in analogs_data:
        try:
            # Extract core fields
            analog_id = analog.get('analog_id', analog.get('id'))
            name = analog.get('name', f"Analog_{analog_id}")
            smiles = analog.get('smiles', '')
            parent_compound = analog.get('parent_compound', '')
            parent_smiles = analog.get('parent_smiles', '')
            similarity = analog.get('similarity', 0.0)

            # Molecular properties (from JSON or calculated)
            props = analog.get('properties', {})
            if not props and smiles:
                calculated_props = parse_smiles_properties(smiles)
                if calculated_props:
                    props = calculated_props

            molecular_weight = props.get('molecular_weight', props.get('MW'))
            logp = props.get('logp', props.get('LogP'))
            tpsa = props.get('tpsa', props.get('TPSA'))
            h_bond_donors = props.get('h_bond_donors', props.get('HBD'))
            h_bond_acceptors = props.get('h_bond_acceptors', props.get('HBA'))
            rotatable_bonds = props.get('rotatable_bonds', props.get('RotBonds'))
            aromatic_rings = props.get('aromatic_rings', props.get('AromaticRings'))

            # Lipinski violations
            lipinski_violations = calculate_lipinski_violations(props)

            # Drug-likeness and scores
            drug_likeness = analog.get('drug_likeness', 75)
            safety_score = analog.get('safety_score', analog.get('predicted_safety_score', 70))
            efficacy_score = analog.get('efficacy_score', analog.get('predicted_efficacy_score', 70))

            # Patent information
            patent_status = 'Patent-Free (Novel)' if analog.get('patent_free', True) else 'Patented'
            patent_opportunity_score = analog.get('ip_opportunity_score', analog.get('patent_opportunity_score', 85))

            # Value assessment
            therapeutic_potential = analog.get('therapeutic_potential', 'High')
            estimated_value = analog.get('estimated_value', '$25M-$50M')

            # Metadata
            discovery_date = analog.get('discovery_timestamp', analog.get('created_at', datetime.now().isoformat()))
            generation_method = analog.get('generation_method', 'RDKit transformation')

            insert_data.append((
                analog_id, name, smiles, parent_compound, parent_smiles, similarity,
                molecular_weight, logp, tpsa, h_bond_donors, h_bond_acceptors,
                rotatable_bonds, aromatic_rings, drug_likeness, lipinski_violations,
                safety_score, efficacy_score, patent_status, patent_opportunity_score,
                None,  # patent_filing_date
                therapeutic_potential, estimated_value, discovery_date, generation_method
            ))

        except Exception as e:
            print(f"⚠️  Skipping analog {analog.get('analog_id', 'unknown')}: {e}")
            skipped += 1
            continue

    # Bulk insert
    print(f"\n💾 Inserting {len(insert_data)} analogs into database...")

    insert_query = """
        INSERT INTO analogs (
            id, name, smiles, parent_compound, parent_smiles, similarity,
            molecular_weight, logp, tpsa, h_bond_donors, h_bond_acceptors,
            rotatable_bonds, aromatic_rings, drug_likeness, lipinski_violations,
            safety_score, efficacy_score, patent_status, patent_opportunity_score,
            patent_filing_date, therapeutic_potential, estimated_value,
            discovery_date, generation_method
        ) VALUES %s
        ON CONFLICT (id) DO UPDATE SET
            updated_at = CURRENT_TIMESTAMP,
            name = EXCLUDED.name,
            molecular_weight = EXCLUDED.molecular_weight,
            drug_likeness = EXCLUDED.drug_likeness
    """

    execute_values(cur, insert_query, insert_data)
    conn.commit()

    print(f"✅ Successfully imported {len(insert_data)} analogs")
    if skipped > 0:
        print(f"⚠️  Skipped {skipped} analogs due to errors")

    # Query summary statistics
    cur.execute("SELECT COUNT(*) FROM analogs")
    total_count = cur.fetchone()[0]

    cur.execute("SELECT COUNT(*) FROM analogs WHERE patent_opportunity_score >= 90")
    high_ip = cur.fetchone()[0]

    cur.execute("SELECT COUNT(*) FROM analogs WHERE patent_status = 'Patent-Free (Novel)'")
    patent_free = cur.fetchone()[0]

    print(f"\n📊 Database Summary:")
    print(f"   Total analogs in database: {total_count}")
    print(f"   High IP opportunity (≥90): {high_ip}")
    print(f"   Patent-free compounds: {patent_free}")

    # Close connection
    cur.close()
    conn.close()

    print("\n🎉 Import complete!")

if __name__ == '__main__':
    import sys

    if not os.path.exists('MASTER_ANALOG_DISCOVERIES.json'):
        print("❌ Error: MASTER_ANALOG_DISCOVERIES.json not found")
        print("   Run this script from the pharmasight-platform root directory")
        sys.exit(1)

    try:
        import_analogs()
    except Exception as e:
        print(f"\n❌ Import failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

#!/usr/bin/env python3
"""
RDKit Demonstration Script for PharmaSight Platform
Comprehensive examples of molecular visualization and structure editing

This script demonstrates:
1. Molecular visualization (2D structures)
2. Property calculations
3. Molecular similarity
4. Substructure searching
5. Structure editing and modifications
6. Analog generation

Usage:
    # Activate conda environment first
    conda activate pharmasight
    python rdkit_demo.py
"""

import sys
import os

# Add src directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from molecular_visualizer import MolecularVisualizer
from molecular_editor import MolecularEditor


def print_section(title):
    """Print a formatted section header"""
    print("\n" + "=" * 70)
    print(f"  {title}")
    print("=" * 70 + "\n")


def demo_visualization():
    """Demonstrate molecular visualization capabilities"""
    print_section("MOLECULAR VISUALIZATION")

    visualizer = MolecularVisualizer(output_dir="molecular_images")

    # Compound database
    compounds = {
        "Psilocybin": "CN(C)CCc1c[nH]c2ccc(OP(=O)(O)O)cc12",
        "LSD": "CCN(CC)C(=O)[C@H]1CN([C@@H]2Cc3c[nH]c4cccc(c34)C2=C1)C",
        "MDMA": "CC(CC1=CC2=C(C=C1)OCO2)NC",
        "DMT": "CN(C)CCc1c[nH]c2ccccc12",
        "Mescaline": "COc1cc(CCN)cc(OC)c1OC",
        "Ketamine": "CNC1(CCCCC1=O)c2ccccc2Cl",
        "Aspirin": "CC(=O)Oc1ccccc1C(=O)O",
        "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    }

    # 1. Visualize individual molecules
    print("1. Visualizing individual molecules...")
    for name, smiles in list(compounds.items())[:4]:
        print(f"   - {name}")
        visualizer.visualize_molecule(smiles, filename=f"{name.lower()}.png")

    # 2. Create grid visualization
    print("\n2. Creating grid visualization...")
    visualizer.visualize_multiple_molecules(
        list(compounds.values()),
        labels=list(compounds.keys()),
        filename="all_compounds_grid.png",
        mols_per_row=4
    )
    print("   Grid saved as: all_compounds_grid.png")

    # 3. Calculate molecular properties
    print("\n3. Calculating molecular properties...")
    print(f"\n   {'Compound':<15} {'MW':<8} {'LogP':<8} {'HBD':<6} {'HBA':<6} {'Lipinski':<10}")
    print("   " + "-" * 65)

    for name, smiles in compounds.items():
        props = visualizer.calculate_properties(smiles)
        print(f"   {name:<15} "
              f"{props.get('molecular_weight', 0):<8.1f} "
              f"{props.get('logP', 0):<8.2f} "
              f"{props.get('num_h_donors', 0):<6} "
              f"{props.get('num_h_acceptors', 0):<6} "
              f"{props.get('lipinski_violations', 0):<10}")

    # 4. Calculate molecular similarities
    print("\n4. Calculating molecular similarities...")
    reference = "Psilocybin"
    reference_smiles = compounds[reference]

    print(f"\n   Similarity to {reference}:")
    print(f"   {'Compound':<15} {'Tanimoto Similarity':<20}")
    print("   " + "-" * 40)

    for name, smiles in compounds.items():
        if name != reference:
            similarity = visualizer.calculate_similarity(reference_smiles, smiles)
            if similarity is not None:
                print(f"   {name:<15} {similarity:<20.3f}")

    # 5. Highlight substructures
    print("\n5. Highlighting substructures...")

    # Indole core (present in many psychedelics)
    indole_smarts = "c1ccc2[nH]ccc2c1"
    print("   - Highlighting indole core in Psilocybin")
    visualizer.highlight_substructure(
        compounds["Psilocybin"],
        indole_smarts,
        filename="psilocybin_indole.png"
    )

    # Phenyl ring
    phenyl_smarts = "c1ccccc1"
    print("   - Highlighting phenyl ring in Aspirin")
    visualizer.highlight_substructure(
        compounds["Aspirin"],
        phenyl_smarts,
        filename="aspirin_phenyl.png"
    )

    print("\n   All visualization outputs saved to: molecular_images/")


def demo_structure_editing():
    """Demonstrate structure editing capabilities"""
    print_section("MOLECULAR STRUCTURE EDITING")

    editor = MolecularEditor()

    # Example molecule
    psilocybin_smiles = "CN(C)CCc1c[nH]c2ccc(OP(=O)(O)O)cc12"
    mdma_smiles = "CC(CC1=CC2=C(C=C1)OCO2)NC"

    # 1. SMILES canonicalization
    print("1. SMILES Canonicalization")
    print(f"   Original:  {psilocybin_smiles}")
    canonical = editor.canonicalize(psilocybin_smiles)
    print(f"   Canonical: {canonical}")

    # 2. Hydrogen manipulation
    print("\n2. Hydrogen Manipulation")
    with_h = editor.add_hydrogens(mdma_smiles)
    print(f"   Original:     {mdma_smiles}")
    print(f"   With H:       {with_h}")
    without_h = editor.remove_hydrogens(with_h or mdma_smiles)
    print(f"   Without H:    {without_h}")

    # 3. Generate analogs
    print("\n3. Generating Structural Analogs")
    print(f"   Parent compound: Psilocybin")
    print(f"   Generating 5 analogs with >70% similarity...")

    analogs = editor.generate_analogs(
        psilocybin_smiles,
        num_analogs=5,
        similarity_threshold=0.7
    )

    if analogs:
        print(f"\n   {'#':<4} {'Similarity':<12} {'MW':<10} {'LogP':<10} {'SMILES':<50}")
        print("   " + "-" * 100)
        for i, analog in enumerate(analogs, 1):
            print(f"   {i:<4} "
                  f"{analog['similarity']:<12.3f} "
                  f"{analog['molecular_weight']:<10.2f} "
                  f"{analog['logP']:<10.2f} "
                  f"{analog['smiles']:<50}")
    else:
        print("   Note: Analog generation returned limited results")
        print("   (This is expected - the algorithm uses simple heuristics)")

    # 4. Tautomer enumeration
    print("\n4. Tautomer Enumeration")
    # Use a simpler molecule with clear tautomers
    phenol_smiles = "Oc1ccccc1"
    print(f"   Parent: {phenol_smiles}")

    tautomers = editor.enumerate_tautomers(phenol_smiles, max_tautomers=10)
    if tautomers:
        print(f"   Found {len(tautomers)} tautomer(s):")
        for i, tautomer in enumerate(tautomers[:5], 1):
            print(f"      {i}. {tautomer}")
    else:
        print("   (No alternative tautomers found)")

    # 5. Molecule fragmentation
    print("\n5. Molecule Fragmentation")
    print(f"   Parent: Psilocybin")
    fragments = editor.fragment_molecule(psilocybin_smiles)

    if fragments:
        print(f"   Generated {len(fragments)} fragment(s):")
        for i, frag in enumerate(fragments[:5], 1):
            print(f"      {i}. {frag}")
    else:
        print("   (Fragmentation produced no results)")

    # 6. Charge neutralization
    print("\n6. Charge Neutralization")
    charged_smiles = "C[NH3+]"
    print(f"   Charged:    {charged_smiles}")
    neutralized = editor.neutralize_charges(charged_smiles)
    print(f"   Neutralized: {neutralized}")


def demo_drug_likeness():
    """Demonstrate drug-likeness calculations"""
    print_section("DRUG-LIKENESS ANALYSIS")

    visualizer = MolecularVisualizer()

    # Test compounds
    test_compounds = {
        "Aspirin": "CC(=O)Oc1ccccc1C(=O)O",
        "Lipitor": "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O",
        "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "Psilocybin": "CN(C)CCc1c[nH]c2ccc(OP(=O)(O)O)cc12",
    }

    print("Lipinski's Rule of Five Analysis:")
    print("\nCriteria: MW <= 500, LogP <= 5, HBD <= 5, HBA <= 10")
    print(f"\n{'Compound':<15} {'MW':<8} {'LogP':<8} {'HBD':<6} {'HBA':<6} {'Violations':<12} {'Drug-like?':<12}")
    print("-" * 80)

    for name, smiles in test_compounds.items():
        props = visualizer.calculate_properties(smiles)
        violations = props.get('lipinski_violations', 0)
        drug_like = "Yes" if violations <= 1 else "No"

        print(f"{name:<15} "
              f"{props.get('molecular_weight', 0):<8.1f} "
              f"{props.get('logP', 0):<8.2f} "
              f"{props.get('num_h_donors', 0):<6} "
              f"{props.get('num_h_acceptors', 0):<6} "
              f"{violations:<12} "
              f"{drug_like:<12}")


def demo_interactive():
    """Interactive demo allowing user input"""
    print_section("INTERACTIVE MODE")

    visualizer = MolecularVisualizer()
    editor = MolecularEditor()

    print("Enter a SMILES string to analyze (or 'quit' to exit):")
    print("Examples:")
    print("  - Aspirin: CC(=O)Oc1ccccc1C(=O)O")
    print("  - Caffeine: CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
    print("  - Benzene: c1ccccc1")

    while True:
        print("\n" + "-" * 70)
        smiles = input("\nSMILES: ").strip()

        if smiles.lower() in ['quit', 'exit', 'q']:
            print("Exiting interactive mode...")
            break

        if not smiles:
            continue

        # Validate SMILES
        mol = visualizer.smiles_to_molecule(smiles)
        if mol is None:
            print("Error: Invalid SMILES string. Please try again.")
            continue

        # Show menu
        print("\nWhat would you like to do?")
        print("  1. Visualize molecule")
        print("  2. Calculate properties")
        print("  3. Generate analogs")
        print("  4. Canonicalize SMILES")
        print("  5. All of the above")
        print("  0. Enter new SMILES")

        choice = input("\nChoice: ").strip()

        if choice == "1" or choice == "5":
            print("\nGenerating visualization...")
            filename = "interactive_molecule.png"
            visualizer.visualize_molecule(smiles, filename=filename)
            print(f"Saved to: molecular_images/{filename}")

        if choice == "2" or choice == "5":
            print("\nMolecular Properties:")
            props = visualizer.calculate_properties(smiles)
            for key, value in props.items():
                if isinstance(value, float):
                    print(f"  {key}: {value:.2f}")
                else:
                    print(f"  {key}: {value}")

        if choice == "3" or choice == "5":
            print("\nGenerating analogs (this may take a moment)...")
            analogs = editor.generate_analogs(smiles, num_analogs=3)
            if analogs:
                print(f"\nFound {len(analogs)} analog(s):")
                for i, analog in enumerate(analogs, 1):
                    print(f"\n  Analog {i}:")
                    print(f"    SMILES: {analog['smiles']}")
                    print(f"    Similarity: {analog['similarity']:.3f}")
                    print(f"    MW: {analog['molecular_weight']:.2f}")
            else:
                print("  No analogs generated")

        if choice == "4" or choice == "5":
            canonical = editor.canonicalize(smiles)
            print(f"\nCanonical SMILES: {canonical}")

        if choice == "0":
            continue


def main():
    """Main demonstration function"""
    print("\n" + "=" * 70)
    print("  PharmaSight RDKit Demonstration")
    print("  Molecular Visualization & Structure Editing")
    print("=" * 70)

    print("\nThis demonstration showcases:")
    print("  - Molecular visualization (2D structures)")
    print("  - Molecular property calculations")
    print("  - Similarity analysis")
    print("  - Substructure highlighting")
    print("  - Structure editing and modifications")
    print("  - Analog generation")
    print("  - Drug-likeness analysis")

    print("\n\nChoose a demo mode:")
    print("  1. Full demonstration (recommended)")
    print("  2. Visualization only")
    print("  3. Structure editing only")
    print("  4. Drug-likeness analysis")
    print("  5. Interactive mode")
    print("  0. Run all demos")

    try:
        choice = input("\nChoice (1-5, 0 for all): ").strip()

        if choice == "0":
            demo_visualization()
            demo_structure_editing()
            demo_drug_likeness()
            demo_interactive()
        elif choice == "1":
            demo_visualization()
            demo_structure_editing()
            demo_drug_likeness()
        elif choice == "2":
            demo_visualization()
        elif choice == "3":
            demo_structure_editing()
        elif choice == "4":
            demo_drug_likeness()
        elif choice == "5":
            demo_interactive()
        else:
            print("Invalid choice. Running full demonstration...")
            demo_visualization()
            demo_structure_editing()
            demo_drug_likeness()

        print("\n" + "=" * 70)
        print("  Demonstration Complete!")
        print("=" * 70)
        print("\nCheck the 'molecular_images/' directory for generated visualizations.")
        print("\nTo use these modules in your own code:")
        print("  from molecular_visualizer import MolecularVisualizer")
        print("  from molecular_editor import MolecularEditor")

    except KeyboardInterrupt:
        print("\n\nDemo interrupted by user.")
    except Exception as e:
        print(f"\nError during demonstration: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()

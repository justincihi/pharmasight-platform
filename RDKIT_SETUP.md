# RDKit Setup Guide for PharmaSight Platform

This guide explains how to set up and use RDKit for molecular visualization and structure manipulation in the PharmaSight platform.

## Table of Contents

1. [Installation](#installation)
2. [Quick Start](#quick-start)
3. [Module Overview](#module-overview)
4. [Usage Examples](#usage-examples)
5. [API Reference](#api-reference)
6. [Troubleshooting](#troubleshooting)

## Installation

### Method 1: Using Conda (Recommended)

RDKit is already installed in the `pharmasight` conda environment. To activate it:

```bash
# Activate the conda environment
conda activate pharmasight

# Verify installation
python -c "from rdkit import Chem; print('RDKit version:', Chem.rdBase.rdkitVersion)"
```

### Method 2: Creating a New Environment

If you need to create a fresh environment:

```bash
# Create environment from file
conda env create -f environment.yml

# Or create manually
conda create -n pharmasight python=3.11 rdkit -c conda-forge

# Activate environment
conda activate pharmasight
```

## Quick Start

### Running the Demo

The easiest way to get started is to run the interactive demo:

```bash
# Make sure conda environment is activated
conda activate pharmasight

# Run the demo
python rdkit_demo.py
```

This will present you with several demo options:
1. Full demonstration (visualization + editing)
2. Visualization only
3. Structure editing only
4. Drug-likeness analysis
5. Interactive mode (enter your own SMILES)

### Basic Example

```python
from src.molecular_visualizer import MolecularVisualizer
from src.molecular_editor import MolecularEditor

# Initialize
visualizer = MolecularVisualizer()
editor = MolecularEditor()

# Visualize a molecule
smiles = "CC(=O)Oc1ccccc1C(=O)O"  # Aspirin
visualizer.visualize_molecule(smiles, filename="aspirin.png")

# Calculate properties
props = visualizer.calculate_properties(smiles)
print(f"Molecular Weight: {props['molecular_weight']:.2f}")
print(f"LogP: {props['logP']:.2f}")

# Generate analogs
analogs = editor.generate_analogs(smiles, num_analogs=5)
for analog in analogs:
    print(f"Analog: {analog['smiles']} (similarity: {analog['similarity']:.3f})")
```

## Module Overview

### molecular_visualizer.py

The `MolecularVisualizer` class provides tools for visualizing molecules and calculating properties.

**Key Features:**
- Convert SMILES to molecule objects
- Generate 2D molecular images
- Create grid visualizations of multiple molecules
- Calculate molecular properties (MW, LogP, H-donors/acceptors, etc.)
- Highlight substructures
- Calculate molecular similarity
- Generate 3D structures

### molecular_editor.py

The `MolecularEditor` class provides tools for modifying molecular structures.

**Key Features:**
- SMILES canonicalization
- Add/remove hydrogen atoms
- Substitute atoms
- Add functional groups
- Generate structural analogs
- Enumerate tautomers
- Fragment molecules
- Neutralize charges

## Usage Examples

### 1. Visualizing Molecules

```python
from src.molecular_visualizer import MolecularVisualizer

viz = MolecularVisualizer(output_dir="my_molecules")

# Single molecule
viz.visualize_molecule(
    "CN(C)CCc1c[nH]c2ccc(OP(=O)(O)O)cc12",
    filename="psilocybin.png",
    size=(600, 600)
)

# Multiple molecules in grid
compounds = {
    "Aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "Ibuprofen": "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
}

viz.visualize_multiple_molecules(
    list(compounds.values()),
    labels=list(compounds.keys()),
    filename="drugs_grid.png",
    mols_per_row=3
)
```

### 2. Calculating Properties

```python
from src.molecular_visualizer import MolecularVisualizer

viz = MolecularVisualizer()

smiles = "CC(=O)Oc1ccccc1C(=O)O"  # Aspirin
props = viz.calculate_properties(smiles)

print(f"Molecular Formula: {props['molecular_formula']}")
print(f"Molecular Weight: {props['molecular_weight']:.2f} Da")
print(f"LogP: {props['logP']:.2f}")
print(f"H-Bond Donors: {props['num_h_donors']}")
print(f"H-Bond Acceptors: {props['num_h_acceptors']}")
print(f"Rotatable Bonds: {props['num_rotatable_bonds']}")
print(f"Aromatic Rings: {props['num_aromatic_rings']}")
print(f"TPSA: {props['tpsa']:.2f} Ų")
print(f"Lipinski Violations: {props['lipinski_violations']}")
```

### 3. Similarity Analysis

```python
from src.molecular_visualizer import MolecularVisualizer

viz = MolecularVisualizer()

# Compare two molecules
smiles1 = "CN(C)CCc1c[nH]c2ccc(OP(=O)(O)O)cc12"  # Psilocybin
smiles2 = "CN(C)CCc1c[nH]c2ccccc12"  # DMT

similarity = viz.calculate_similarity(smiles1, smiles2, method="morgan")
print(f"Tanimoto Similarity: {similarity:.3f}")
```

### 4. Highlighting Substructures

```python
from src.molecular_visualizer import MolecularVisualizer

viz = MolecularVisualizer()

# Highlight indole core in psilocybin
viz.highlight_substructure(
    smiles="CN(C)CCc1c[nH]c2ccc(OP(=O)(O)O)cc12",
    substructure_smarts="c1ccc2[nH]ccc2c1",  # Indole pattern
    filename="psilocybin_indole.png"
)
```

### 5. Generating Analogs

```python
from src.molecular_editor import MolecularEditor

editor = MolecularEditor()

parent_smiles = "CN(C)CCc1c[nH]c2ccc(OP(=O)(O)O)cc12"  # Psilocybin

# Generate 10 analogs with at least 70% similarity
analogs = editor.generate_analogs(
    parent_smiles,
    num_analogs=10,
    similarity_threshold=0.7
)

for i, analog in enumerate(analogs, 1):
    print(f"\nAnalog {i}:")
    print(f"  SMILES: {analog['smiles']}")
    print(f"  Similarity: {analog['similarity']:.3f}")
    print(f"  MW: {analog['molecular_weight']:.2f}")
    print(f"  LogP: {analog['logP']:.2f}")
```

### 6. Structure Modifications

```python
from src.molecular_editor import MolecularEditor

editor = MolecularEditor()

smiles = "CC(=O)Oc1ccccc1C(=O)O"  # Aspirin

# Canonicalize SMILES
canonical = editor.canonicalize(smiles)
print(f"Canonical: {canonical}")

# Add hydrogens
with_h = editor.add_hydrogens(smiles)
print(f"With H: {with_h}")

# Enumerate tautomers
tautomers = editor.enumerate_tautomers(smiles, max_tautomers=10)
for i, tautomer in enumerate(tautomers, 1):
    print(f"Tautomer {i}: {tautomer}")

# Fragment molecule
fragments = editor.fragment_molecule(smiles)
for i, fragment in enumerate(fragments, 1):
    print(f"Fragment {i}: {fragment}")
```

### 7. Drug-likeness Analysis (Lipinski's Rule of Five)

```python
from src.molecular_visualizer import MolecularVisualizer

viz = MolecularVisualizer()

compounds = {
    "Aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "Lipitor": "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O",
}

print("Lipinski's Rule of Five: MW ≤ 500, LogP ≤ 5, HBD ≤ 5, HBA ≤ 10")
print()

for name, smiles in compounds.items():
    props = viz.calculate_properties(smiles)
    violations = props['lipinski_violations']
    drug_like = "Yes" if violations <= 1 else "No"

    print(f"{name}:")
    print(f"  MW: {props['molecular_weight']:.1f} Da")
    print(f"  LogP: {props['logP']:.2f}")
    print(f"  H-Bond Donors: {props['num_h_donors']}")
    print(f"  H-Bond Acceptors: {props['num_h_acceptors']}")
    print(f"  Violations: {violations}")
    print(f"  Drug-like: {drug_like}")
    print()
```

## API Reference

### MolecularVisualizer

#### Methods

- `smiles_to_molecule(smiles: str) -> Optional[Chem.Mol]`
  - Convert SMILES string to RDKit molecule object

- `visualize_molecule(smiles: str, filename: str, size: Tuple[int, int]) -> bool`
  - Generate 2D visualization of a molecule

- `visualize_multiple_molecules(smiles_list: List[str], labels: List[str], filename: str) -> bool`
  - Create a grid visualization of multiple molecules

- `calculate_properties(smiles: str) -> dict`
  - Calculate molecular properties including MW, LogP, Lipinski violations

- `calculate_similarity(smiles1: str, smiles2: str, method: str) -> float`
  - Calculate Tanimoto similarity between two molecules

- `highlight_substructure(smiles: str, substructure_smarts: str, filename: str) -> bool`
  - Highlight a substructure in a molecule

- `generate_3d_structure(smiles: str, optimize: bool) -> Optional[Chem.Mol]`
  - Generate 3D coordinates for a molecule

### MolecularEditor

#### Methods

- `smiles_to_molecule(smiles: str) -> Optional[Chem.Mol]`
  - Convert SMILES string to RDKit molecule object

- `molecule_to_smiles(mol: Chem.Mol) -> Optional[str]`
  - Convert RDKit molecule to SMILES string

- `add_hydrogens(smiles: str, explicit: bool) -> Optional[str]`
  - Add hydrogen atoms to molecule

- `remove_hydrogens(smiles: str) -> Optional[str]`
  - Remove explicit hydrogen atoms

- `canonicalize(smiles: str) -> Optional[str]`
  - Return canonical SMILES representation

- `generate_analogs(smiles: str, num_analogs: int, similarity_threshold: float) -> List[Dict]`
  - Generate structural analogs

- `enumerate_tautomers(smiles: str, max_tautomers: int) -> List[str]`
  - Enumerate possible tautomers

- `fragment_molecule(smiles: str) -> List[str]`
  - Fragment molecule at rotatable bonds

- `neutralize_charges(smiles: str) -> Optional[str]`
  - Neutralize charged molecules

## Troubleshooting

### ImportError: No module named 'rdkit'

Make sure you've activated the conda environment:
```bash
conda activate pharmasight
```

### Visualizations not appearing

Check that the output directory exists. The default is `molecular_images/` in the current working directory.

### Invalid SMILES errors

Verify your SMILES strings using a SMILES validator or RDKit directly:
```python
from rdkit import Chem
mol = Chem.MolFromSmiles("YOUR_SMILES_HERE")
if mol is None:
    print("Invalid SMILES")
```

### Conda environment issues

To recreate the environment:
```bash
conda deactivate
conda env remove -n pharmasight
conda env create -f environment.yml
conda activate pharmasight
```

## Resources

- [RDKit Documentation](https://www.rdkit.org/docs/)
- [RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html)
- [SMILES Tutorial](https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html)
- [SMARTS Tutorial](https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html)

## Support

For issues or questions:
1. Check this documentation
2. Review the example code in `rdkit_demo.py`
3. Consult the RDKit documentation
4. Check the module source code with inline comments

## License

This module is part of the PharmaSight Platform.

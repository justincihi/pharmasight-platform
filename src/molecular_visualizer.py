#!/usr/bin/env python3
"""
Molecular Visualization Module for PharmaSight Platform
Uses RDKit to visualize and display molecular structures
"""

from rdkit import Chem
from rdkit.Chem import Draw, AllChem, Descriptors, Lipinski
from rdkit.Chem import rdMolDescriptors, rdFingerprintGenerator
import os
from typing import List, Optional, Union, Tuple


class MolecularVisualizer:
    """
    Comprehensive molecular visualization toolkit using RDKit
    """

    def __init__(self, output_dir: str = "molecular_images"):
        """
        Initialize the molecular visualizer

        Args:
            output_dir: Directory to save generated molecular images
        """
        self.output_dir = output_dir
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

    def smiles_to_molecule(self, smiles: str) -> Optional[Chem.Mol]:
        """
        Convert SMILES string to RDKit molecule object

        Args:
            smiles: SMILES string representation of molecule

        Returns:
            RDKit molecule object or None if invalid
        """
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                print(f"Error: Invalid SMILES string: {smiles}")
                return None
            return mol
        except Exception as e:
            print(f"Error parsing SMILES: {e}")
            return None

    def visualize_molecule(
        self,
        smiles: str,
        filename: str = "molecule.png",
        size: Tuple[int, int] = (400, 400),
        highlight_atoms: Optional[List[int]] = None,
        highlight_bonds: Optional[List[int]] = None
    ) -> bool:
        """
        Generate 2D visualization of a molecule from SMILES

        Args:
            smiles: SMILES string
            filename: Output filename for the image
            size: Tuple of (width, height) for image
            highlight_atoms: List of atom indices to highlight
            highlight_bonds: List of bond indices to highlight

        Returns:
            True if successful, False otherwise
        """
        mol = self.smiles_to_molecule(smiles)
        if mol is None:
            return False

        try:
            # Generate 2D coordinates
            AllChem.Compute2DCoords(mol)

            # Create drawing
            if highlight_atoms or highlight_bonds:
                img = Draw.MolToImage(
                    mol,
                    size=size,
                    highlightAtoms=highlight_atoms or [],
                    highlightBonds=highlight_bonds or []
                )
            else:
                img = Draw.MolToImage(mol, size=size)

            # Save image
            output_path = os.path.join(self.output_dir, filename)
            img.save(output_path)
            print(f"Molecule visualization saved to: {output_path}")
            return True

        except Exception as e:
            print(f"Error visualizing molecule: {e}")
            return False

    def visualize_multiple_molecules(
        self,
        smiles_list: List[str],
        labels: Optional[List[str]] = None,
        filename: str = "molecules_grid.png",
        mols_per_row: int = 3,
        sub_img_size: Tuple[int, int] = (300, 300)
    ) -> bool:
        """
        Create a grid visualization of multiple molecules

        Args:
            smiles_list: List of SMILES strings
            labels: Optional labels for each molecule
            filename: Output filename
            mols_per_row: Number of molecules per row in grid
            sub_img_size: Size of each sub-image

        Returns:
            True if successful, False otherwise
        """
        mols = []
        for smiles in smiles_list:
            mol = self.smiles_to_molecule(smiles)
            if mol:
                AllChem.Compute2DCoords(mol)
                mols.append(mol)
            else:
                mols.append(None)

        if not any(mols):
            print("Error: No valid molecules to visualize")
            return False

        try:
            img = Draw.MolsToGridImage(
                mols,
                molsPerRow=mols_per_row,
                subImgSize=sub_img_size,
                legends=labels
            )

            output_path = os.path.join(self.output_dir, filename)
            img.save(output_path)
            print(f"Grid visualization saved to: {output_path}")
            return True

        except Exception as e:
            print(f"Error creating grid visualization: {e}")
            return False

    def calculate_properties(self, smiles: str) -> dict:
        """
        Calculate various molecular properties

        Args:
            smiles: SMILES string

        Returns:
            Dictionary of molecular properties
        """
        mol = self.smiles_to_molecule(smiles)
        if mol is None:
            return {}

        properties = {
            "molecular_formula": Chem.rdMolDescriptors.CalcMolFormula(mol),
            "molecular_weight": Descriptors.MolWt(mol),
            "exact_molecular_weight": Descriptors.ExactMolWt(mol),
            "logP": Descriptors.MolLogP(mol),
            "num_h_donors": Descriptors.NumHDonors(mol),
            "num_h_acceptors": Descriptors.NumHAcceptors(mol),
            "num_rotatable_bonds": Descriptors.NumRotatableBonds(mol),
            "num_aromatic_rings": Descriptors.NumAromaticRings(mol),
            "num_aliphatic_rings": Descriptors.NumAliphaticRings(mol),
            "tpsa": Descriptors.TPSA(mol),  # Topological Polar Surface Area
            "num_atoms": mol.GetNumAtoms(),
            "num_heavy_atoms": mol.GetNumHeavyAtoms(),
            "num_bonds": mol.GetNumBonds(),
            "fraction_csp3": Descriptors.FractionCSP3(mol),
            "num_saturated_rings": Descriptors.NumSaturatedRings(mol),
        }

        # Lipinski's Rule of Five
        properties["lipinski_violations"] = self._check_lipinski(mol)

        return properties

    def _check_lipinski(self, mol: Chem.Mol) -> int:
        """
        Check Lipinski's Rule of Five violations

        Returns:
            Number of violations (0-4)
        """
        violations = 0
        if Descriptors.MolWt(mol) > 500:
            violations += 1
        if Descriptors.MolLogP(mol) > 5:
            violations += 1
        if Descriptors.NumHDonors(mol) > 5:
            violations += 1
        if Descriptors.NumHAcceptors(mol) > 10:
            violations += 1
        return violations

    def generate_3d_structure(self, smiles: str, optimize: bool = True) -> Optional[Chem.Mol]:
        """
        Generate 3D coordinates for a molecule

        Args:
            smiles: SMILES string
            optimize: Whether to optimize geometry with force field

        Returns:
            Molecule with 3D coordinates or None
        """
        mol = self.smiles_to_molecule(smiles)
        if mol is None:
            return None

        try:
            # Add hydrogens
            mol = Chem.AddHs(mol)

            # Generate 3D coordinates
            success = AllChem.EmbedMolecule(mol, randomSeed=42)
            if success != 0:
                print("Warning: 3D embedding failed, using 2D coordinates")
                return mol

            # Optimize with force field if requested
            if optimize:
                AllChem.MMFFOptimizeMolecule(mol)

            return mol

        except Exception as e:
            print(f"Error generating 3D structure: {e}")
            return None

    def highlight_substructure(
        self,
        smiles: str,
        substructure_smarts: str,
        filename: str = "highlighted_substructure.png",
        size: Tuple[int, int] = (400, 400)
    ) -> bool:
        """
        Visualize molecule with substructure highlighted

        Args:
            smiles: SMILES string of the molecule
            substructure_smarts: SMARTS pattern for substructure
            filename: Output filename
            size: Image size

        Returns:
            True if successful, False otherwise
        """
        mol = self.smiles_to_molecule(smiles)
        if mol is None:
            return False

        try:
            # Parse substructure
            substructure = Chem.MolFromSmarts(substructure_smarts)
            if substructure is None:
                print(f"Error: Invalid SMARTS pattern: {substructure_smarts}")
                return False

            # Find matching atoms
            matches = mol.GetSubstructMatches(substructure)
            if not matches:
                print("Warning: No substructure matches found")
                return False

            # Highlight first match
            highlight_atoms = list(matches[0])

            # Get bonds to highlight
            highlight_bonds = []
            for bond in mol.GetBonds():
                if bond.GetBeginAtomIdx() in highlight_atoms and \
                   bond.GetEndAtomIdx() in highlight_atoms:
                    highlight_bonds.append(bond.GetIdx())

            return self.visualize_molecule(
                smiles,
                filename=filename,
                size=size,
                highlight_atoms=highlight_atoms,
                highlight_bonds=highlight_bonds
            )

        except Exception as e:
            print(f"Error highlighting substructure: {e}")
            return False

    def calculate_similarity(
        self,
        smiles1: str,
        smiles2: str,
        method: str = "morgan"
    ) -> Optional[float]:
        """
        Calculate molecular similarity between two molecules

        Args:
            smiles1: First SMILES string
            smiles2: Second SMILES string
            method: Fingerprint method ('morgan', 'rdkit', 'topological')

        Returns:
            Tanimoto similarity score (0-1) or None if error
        """
        mol1 = self.smiles_to_molecule(smiles1)
        mol2 = self.smiles_to_molecule(smiles2)

        if mol1 is None or mol2 is None:
            return None

        try:
            if method == "morgan":
                fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=2048)
                fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=2048)
            elif method == "rdkit":
                fp1 = Chem.RDKFingerprint(mol1)
                fp2 = Chem.RDKFingerprint(mol2)
            elif method == "topological":
                fp1 = Chem.RDKFingerprint(mol1)
                fp2 = Chem.RDKFingerprint(mol2)
            else:
                print(f"Unknown method: {method}")
                return None

            from rdkit import DataStructs
            similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
            return similarity

        except Exception as e:
            print(f"Error calculating similarity: {e}")
            return None


def main():
    """Example usage of the MolecularVisualizer"""

    visualizer = MolecularVisualizer()

    # Example molecules (from the compound database)
    examples = {
        "Psilocybin": "CN(C)CCc1c[nH]c2ccc(OP(=O)(O)O)cc12",
        "LSD": "CCN(CC)C(=O)[C@H]1CN([C@@H]2Cc3c[nH]c4cccc(c34)C2=C1)C",
        "MDMA": "CC(CC1=CC2=C(C=C1)OCO2)NC",
        "DMT": "CN(C)CCc1c[nH]c2ccccc12",
        "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    }

    print("=== Molecular Visualization Examples ===\n")

    # 1. Visualize single molecules
    for name, smiles in examples.items():
        print(f"Visualizing {name}...")
        visualizer.visualize_molecule(smiles, filename=f"{name.lower()}.png")

        # Calculate properties
        props = visualizer.calculate_properties(smiles)
        print(f"  MW: {props.get('molecular_weight', 0):.2f}")
        print(f"  LogP: {props.get('logP', 0):.2f}")
        print(f"  Lipinski violations: {props.get('lipinski_violations', 0)}")
        print()

    # 2. Create grid visualization
    print("Creating grid visualization...")
    visualizer.visualize_multiple_molecules(
        list(examples.values()),
        labels=list(examples.keys()),
        filename="psychedelics_grid.png"
    )

    # 3. Calculate similarity
    print("\nCalculating molecular similarities:")
    psilocybin_smiles = examples["Psilocybin"]
    for name, smiles in examples.items():
        if name != "Psilocybin":
            similarity = visualizer.calculate_similarity(psilocybin_smiles, smiles)
            print(f"  Psilocybin vs {name}: {similarity:.3f}")

    # 4. Highlight substructure (indole core)
    print("\nHighlighting indole substructure in Psilocybin...")
    indole_smarts = "c1ccc2[nH]ccc2c1"
    visualizer.highlight_substructure(
        examples["Psilocybin"],
        indole_smarts,
        filename="psilocybin_indole_highlighted.png"
    )

    print("\n=== All visualizations complete! ===")
    print(f"Images saved to: {visualizer.output_dir}/")


if __name__ == "__main__":
    main()

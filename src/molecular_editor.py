#!/usr/bin/env python3
"""
Molecular Structure Editor Module for PharmaSight Platform
Uses RDKit to alter and modify molecular structures
"""

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit.Chem import Fragments, Lipinski
from rdkit.Chem import rdFingerprintGenerator, DataStructs
from typing import List, Optional, Tuple, Dict
import random


class MolecularEditor:
    """
    Comprehensive molecular structure editing toolkit using RDKit
    """

    def __init__(self):
        """Initialize the molecular editor"""
        pass

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

    def molecule_to_smiles(self, mol: Chem.Mol) -> Optional[str]:
        """
        Convert RDKit molecule to SMILES string

        Args:
            mol: RDKit molecule object

        Returns:
            SMILES string or None if error
        """
        try:
            return Chem.MolToSmiles(mol)
        except Exception as e:
            print(f"Error converting to SMILES: {e}")
            return None

    def add_hydrogens(self, smiles: str, explicit: bool = True) -> Optional[str]:
        """
        Add hydrogen atoms to molecule

        Args:
            smiles: Input SMILES string
            explicit: If True, add explicit hydrogens

        Returns:
            Modified SMILES string or None
        """
        mol = self.smiles_to_molecule(smiles)
        if mol is None:
            return None

        try:
            if explicit:
                mol_h = Chem.AddHs(mol)
            else:
                mol_h = mol

            return self.molecule_to_smiles(mol_h)
        except Exception as e:
            print(f"Error adding hydrogens: {e}")
            return None

    def remove_hydrogens(self, smiles: str) -> Optional[str]:
        """
        Remove explicit hydrogen atoms from molecule

        Args:
            smiles: Input SMILES string

        Returns:
            Modified SMILES string or None
        """
        mol = self.smiles_to_molecule(smiles)
        if mol is None:
            return None

        try:
            mol_no_h = Chem.RemoveHs(mol)
            return self.molecule_to_smiles(mol_no_h)
        except Exception as e:
            print(f"Error removing hydrogens: {e}")
            return None

    def substitute_atom(
        self,
        smiles: str,
        atom_idx: int,
        new_atom: str
    ) -> Optional[str]:
        """
        Replace an atom at specific index with a new atom

        Args:
            smiles: Input SMILES string
            atom_idx: Index of atom to replace
            new_atom: Atomic symbol of new atom (e.g., 'N', 'O', 'S')

        Returns:
            Modified SMILES string or None
        """
        mol = self.smiles_to_molecule(smiles)
        if mol is None:
            return None

        try:
            # Create editable molecule
            em = Chem.EditableMol(mol)

            # Get the atom
            if atom_idx >= mol.GetNumAtoms():
                print(f"Error: Atom index {atom_idx} out of range")
                return None

            atom = mol.GetAtomWithIdx(atom_idx)

            # Create new atom
            new_atom_obj = Chem.Atom(new_atom)

            # Replace by removing and adding
            # Note: This is a simplified approach
            # For more complex substitutions, consider using reaction templates
            rwmol = Chem.RWMol(mol)
            rwmol.ReplaceAtom(atom_idx, new_atom_obj)

            return self.molecule_to_smiles(rwmol)

        except Exception as e:
            print(f"Error substituting atom: {e}")
            return None

    def add_functional_group(
        self,
        smiles: str,
        attachment_idx: int,
        functional_group: str
    ) -> Optional[str]:
        """
        Add a functional group to a molecule at specified atom

        Args:
            smiles: Input SMILES string
            attachment_idx: Atom index where to attach group
            functional_group: SMILES of functional group to add

        Returns:
            Modified SMILES string or None
        """
        mol = self.smiles_to_molecule(smiles)
        if mol is None:
            return None

        fg_mol = self.smiles_to_molecule(functional_group)
        if fg_mol is None:
            return None

        try:
            # Combine molecules (simplified approach)
            combo = Chem.CombineMols(mol, fg_mol)
            return self.molecule_to_smiles(combo)

        except Exception as e:
            print(f"Error adding functional group: {e}")
            return None

    def methylate(self, smiles: str, position: Optional[int] = None) -> List[str]:
        """
        Add methyl group(s) to molecule

        Args:
            smiles: Input SMILES string
            position: Specific atom index to methylate (None for all possibilities)

        Returns:
            List of methylated SMILES strings
        """
        mol = self.smiles_to_molecule(smiles)
        if mol is None:
            return []

        results = []

        try:
            if position is not None:
                # Methylate specific position
                combo = Chem.CombineMols(mol, self.smiles_to_molecule("C"))
                results.append(self.molecule_to_smiles(combo))
            else:
                # Generate multiple methylated versions
                # This is a simplified approach
                # In practice, use reaction SMARTS for proper methylation
                for atom in mol.GetAtoms():
                    if atom.GetSymbol() in ['N', 'O']:  # Methylate N and O
                        combo = Chem.CombineMols(mol, self.smiles_to_molecule("C"))
                        results.append(self.molecule_to_smiles(combo))

            return results

        except Exception as e:
            print(f"Error methylating molecule: {e}")
            return []

    def generate_analogs(
        self,
        smiles: str,
        num_analogs: int = 5,
        similarity_threshold: float = 0.7
    ) -> List[Dict[str, any]]:
        """
        Generate structural analogs of a molecule

        Args:
            smiles: Input SMILES string
            num_analogs: Number of analogs to generate
            similarity_threshold: Minimum Tanimoto similarity to parent

        Returns:
            List of dictionaries with analog information
        """
        mol = self.smiles_to_molecule(smiles)
        if mol is None:
            return []

        analogs = []

        try:
            # Strategy 1: Replace aromatic rings
            analogs.extend(self._replace_aromatic_rings(smiles))

            # Strategy 2: Modify substituents
            analogs.extend(self._modify_substituents(smiles))

            # Strategy 3: Add/remove methyl groups
            analogs.extend(self._add_remove_methyl(smiles))

            # Strategy 4: Halogenation
            analogs.extend(self._halogenate(smiles))

            # Calculate properties and filter
            result_analogs = []
            fpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
            parent_fp = fpgen.GetFingerprint(mol)

            for analog_smiles in analogs[:num_analogs * 2]:  # Generate extra for filtering
                analog_mol = self.smiles_to_molecule(analog_smiles)
                if analog_mol is None:
                    continue

                # Calculate similarity using modern generator API
                analog_fp = fpgen.GetFingerprint(analog_mol)
                similarity = DataStructs.TanimotoSimilarity(parent_fp, analog_fp)

                if similarity >= similarity_threshold:
                    result_analogs.append({
                        "smiles": analog_smiles,
                        "similarity": similarity,
                        "molecular_weight": Descriptors.MolWt(analog_mol),
                        "logP": Descriptors.MolLogP(analog_mol)
                    })

            # Sort by similarity and return top results
            result_analogs.sort(key=lambda x: x["similarity"], reverse=True)
            return result_analogs[:num_analogs]

        except Exception as e:
            print(f"Error generating analogs: {e}")
            return []

    def _replace_aromatic_rings(self, smiles: str) -> List[str]:
        """Replace aromatic rings with alternatives"""
        # Simplified: just return original for now
        # In practice, use reaction SMARTS to replace rings
        return []

    def _modify_substituents(self, smiles: str) -> List[str]:
        """Modify substituents on the molecule"""
        mol = self.smiles_to_molecule(smiles)
        if mol is None:
            return []

        results = []
        # Simple modifications - add variations
        modifications = [
            ("C", "CC"),  # Methyl to ethyl
            ("CC", "C"),  # Ethyl to methyl
            ("O", "S"),   # O to S
            ("N", "O"),   # N to O
        ]

        original_smiles = self.molecule_to_smiles(mol)
        for old, new in modifications:
            if old in original_smiles:
                modified = original_smiles.replace(old, new, 1)
                results.append(modified)

        return results

    def _add_remove_methyl(self, smiles: str) -> List[str]:
        """Add or remove methyl groups"""
        mol = self.smiles_to_molecule(smiles)
        if mol is None:
            return []

        results = []

        # Add methyl
        for i in range(min(3, mol.GetNumAtoms())):
            combo = Chem.CombineMols(mol, self.smiles_to_molecule("C"))
            results.append(self.molecule_to_smiles(combo))

        return results

    def _halogenate(self, smiles: str) -> List[str]:
        """Add halogen substitutions"""
        results = []
        halogens = ["F", "Cl", "Br"]

        for halogen in halogens:
            # Simple approach: combine with halogen
            mol = self.smiles_to_molecule(smiles)
            if mol:
                combo = Chem.CombineMols(mol, self.smiles_to_molecule(halogen))
                results.append(self.molecule_to_smiles(combo))

        return results

    def neutralize_charges(self, smiles: str) -> Optional[str]:
        """
        Neutralize charged molecules

        Args:
            smiles: Input SMILES string

        Returns:
            Neutralized SMILES string or None
        """
        mol = self.smiles_to_molecule(smiles)
        if mol is None:
            return None

        try:
            # Define neutralization reactions
            pattern = Chem.MolFromSmarts("[+1!h0]")
            if mol.HasSubstructMatch(pattern):
                # Remove proton from positively charged atoms
                for match in mol.GetSubstructMatches(pattern):
                    atom = mol.GetAtomWithIdx(match[0])
                    atom.SetFormalCharge(0)
                    atom.SetNumExplicitHs(atom.GetTotalNumHs() - 1)

            pattern = Chem.MolFromSmarts("[-1]")
            if mol.HasSubstructMatch(pattern):
                # Add proton to negatively charged atoms
                for match in mol.GetSubstructMatches(pattern):
                    atom = mol.GetAtomWithIdx(match[0])
                    atom.SetFormalCharge(0)
                    atom.SetNumExplicitHs(atom.GetTotalNumHs() + 1)

            return self.molecule_to_smiles(mol)

        except Exception as e:
            print(f"Error neutralizing charges: {e}")
            return None

    def canonicalize(self, smiles: str) -> Optional[str]:
        """
        Return canonical SMILES representation

        Args:
            smiles: Input SMILES string

        Returns:
            Canonical SMILES string or None
        """
        mol = self.smiles_to_molecule(smiles)
        if mol is None:
            return None

        try:
            return Chem.MolToSmiles(mol, canonical=True)
        except Exception as e:
            print(f"Error canonicalizing SMILES: {e}")
            return None

    def enumerate_tautomers(self, smiles: str, max_tautomers: int = 10) -> List[str]:
        """
        Enumerate possible tautomers of a molecule

        Args:
            smiles: Input SMILES string
            max_tautomers: Maximum number of tautomers to generate

        Returns:
            List of tautomer SMILES strings
        """
        mol = self.smiles_to_molecule(smiles)
        if mol is None:
            return []

        try:
            from rdkit.Chem.MolStandardize import rdMolStandardize

            enumerator = rdMolStandardize.TautomerEnumerator()
            enumerator.SetMaxTautomers(max_tautomers)

            tautomers = enumerator.Enumerate(mol)
            tautomer_smiles = [self.molecule_to_smiles(t) for t in tautomers]

            return [s for s in tautomer_smiles if s is not None]

        except Exception as e:
            print(f"Error enumerating tautomers: {e}")
            return []

    def fragment_molecule(self, smiles: str) -> List[str]:
        """
        Fragment molecule at rotatable bonds

        Args:
            smiles: Input SMILES string

        Returns:
            List of fragment SMILES strings
        """
        mol = self.smiles_to_molecule(smiles)
        if mol is None:
            return []

        try:
            fragments = []

            # Get rotatable bonds
            rotatable_bonds = Chem.Fragments.fr_Al_OH(mol)

            # Fragment at each bond (simplified)
            from rdkit.Chem import FragmentOnBonds

            # Get all bonds
            bonds_to_break = [bond.GetIdx() for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.SINGLE]

            if bonds_to_break:
                # Break one bond at a time
                for bond_idx in bonds_to_break[:5]:  # Limit to first 5
                    fragmented = Chem.FragmentOnBonds(mol, [bond_idx], addDummies=False)
                    frags = Chem.GetMolFrags(fragmented, asMols=True)
                    for frag in frags:
                        frag_smiles = self.molecule_to_smiles(frag)
                        if frag_smiles:
                            fragments.append(frag_smiles)

            return list(set(fragments))  # Remove duplicates

        except Exception as e:
            print(f"Error fragmenting molecule: {e}")
            return []


def main():
    """Example usage of the MolecularEditor"""

    editor = MolecularEditor()

    print("=== Molecular Structure Editing Examples ===\n")

    # Example molecule: Psilocybin
    psilocybin_smiles = "CN(C)CCc1c[nH]c2ccc(OP(=O)(O)O)cc12"
    print(f"Original molecule (Psilocybin): {psilocybin_smiles}\n")

    # 1. Canonicalize
    print("1. Canonicalizing SMILES:")
    canonical = editor.canonicalize(psilocybin_smiles)
    print(f"   Canonical: {canonical}\n")

    # 2. Add/Remove hydrogens
    print("2. Hydrogen manipulation:")
    with_h = editor.add_hydrogens(psilocybin_smiles)
    print(f"   With explicit H: {with_h}")
    without_h = editor.remove_hydrogens(with_h or psilocybin_smiles)
    print(f"   Without explicit H: {without_h}\n")

    # 3. Generate analogs
    print("3. Generating structural analogs:")
    analogs = editor.generate_analogs(psilocybin_smiles, num_analogs=5)
    for i, analog in enumerate(analogs, 1):
        print(f"   Analog {i}:")
        print(f"      SMILES: {analog['smiles']}")
        print(f"      Similarity: {analog['similarity']:.3f}")
        print(f"      MW: {analog['molecular_weight']:.2f}")
        print(f"      LogP: {analog['logP']:.2f}")
        print()

    # 4. Enumerate tautomers
    print("4. Enumerating tautomers:")
    # Use simpler molecule for tautomer example
    simple_mol = "Oc1ccccc1"
    tautomers = editor.enumerate_tautomers(simple_mol, max_tautomers=5)
    print(f"   Parent: {simple_mol}")
    for i, tautomer in enumerate(tautomers, 1):
        print(f"   Tautomer {i}: {tautomer}")
    print()

    # 5. Fragment molecule
    print("5. Fragmenting molecule:")
    fragments = editor.fragment_molecule(psilocybin_smiles)
    print(f"   Found {len(fragments)} fragments:")
    for i, frag in enumerate(fragments[:5], 1):  # Show first 5
        print(f"   Fragment {i}: {frag}")
    print()

    print("=== All structure editing examples complete! ===")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Unit tests for RDKit integration modules

Run with: pytest tests/test_rdkit_integration.py -v
"""

import pytest
import sys
import os

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from molecular_visualizer import MolecularVisualizer
from molecular_editor import MolecularEditor
from admet_predictor import ADMETPredictor


class TestMolecularVisualizer:
    """Tests for MolecularVisualizer class"""

    def setup_method(self):
        """Setup test fixtures"""
        self.visualizer = MolecularVisualizer(output_dir="test_output")
        self.valid_smiles = "CC(=O)O"  # Acetic acid
        self.invalid_smiles = "Invalid SMILES"

    def test_smiles_to_molecule_valid(self):
        """Test valid SMILES conversion"""
        mol = self.visualizer.smiles_to_molecule(self.valid_smiles)
        assert mol is not None
        assert mol.GetNumAtoms() > 0

    def test_smiles_to_molecule_invalid(self):
        """Test invalid SMILES conversion"""
        mol = self.visualizer.smiles_to_molecule(self.invalid_smiles)
        assert mol is None

    def test_calculate_properties(self):
        """Test property calculation"""
        props = self.visualizer.calculate_properties(self.valid_smiles)

        assert 'molecular_weight' in props
        assert 'logP' in props
        assert 'num_h_donors' in props
        assert 'num_h_acceptors' in props
        assert 'lipinski_violations' in props

        # Check specific values for acetic acid
        assert abs(props['molecular_weight'] - 60.05) < 0.1
        assert props['num_h_donors'] == 1
        assert props['num_h_acceptors'] >= 2

    def test_calculate_properties_invalid(self):
        """Test property calculation with invalid SMILES"""
        props = self.visualizer.calculate_properties(self.invalid_smiles)
        assert props == {}

    def test_calculate_similarity(self):
        """Test similarity calculation"""
        smiles1 = "CC(=O)O"  # Acetic acid
        smiles2 = "CCC(=O)O"  # Propionic acid

        similarity = self.visualizer.calculate_similarity(smiles1, smiles2)

        assert similarity is not None
        assert 0 <= similarity <= 1
        assert similarity > 0.5  # Should be quite similar

    def test_calculate_similarity_identical(self):
        """Test similarity of identical molecules"""
        similarity = self.visualizer.calculate_similarity(
            self.valid_smiles,
            self.valid_smiles
        )

        assert similarity == 1.0

    def test_calculate_similarity_invalid(self):
        """Test similarity with invalid SMILES"""
        similarity = self.visualizer.calculate_similarity(
            self.valid_smiles,
            self.invalid_smiles
        )

        assert similarity is None


class TestMolecularEditor:
    """Tests for MolecularEditor class"""

    def setup_method(self):
        """Setup test fixtures"""
        self.editor = MolecularEditor()
        self.valid_smiles = "CC(=O)O"
        self.benzene = "c1ccccc1"

    def test_canonicalize(self):
        """Test SMILES canonicalization"""
        # Different representations of the same molecule
        smiles1 = "C(C)O"
        smiles2 = "OCC"

        canonical1 = self.editor.canonicalize(smiles1)
        canonical2 = self.editor.canonicalize(smiles2)

        assert canonical1 == canonical2
        assert canonical1 is not None

    def test_canonicalize_invalid(self):
        """Test canonicalization with invalid SMILES"""
        result = self.editor.canonicalize("Invalid")
        assert result is None

    def test_add_hydrogens(self):
        """Test hydrogen addition"""
        with_h = self.editor.add_hydrogens(self.benzene, explicit=True)

        assert with_h is not None
        assert "[H]" in with_h or "H" in with_h  # Contains explicit hydrogens

    def test_remove_hydrogens(self):
        """Test hydrogen removal"""
        with_h = self.editor.add_hydrogens(self.benzene, explicit=True)
        without_h = self.editor.remove_hydrogens(with_h)

        assert without_h is not None
        # Should not contain explicit H notation
        assert "[H]" not in without_h

    def test_enumerate_tautomers(self):
        """Test tautomer enumeration"""
        phenol = "Oc1ccccc1"
        tautomers = self.editor.enumerate_tautomers(phenol, max_tautomers=5)

        assert isinstance(tautomers, list)
        assert len(tautomers) > 0
        # Original should be in tautomers
        assert phenol in tautomers or any(phenol in t for t in tautomers)

    def test_fragment_molecule(self):
        """Test molecule fragmentation"""
        ethanol = "CCO"
        fragments = self.editor.fragment_molecule(ethanol)

        assert isinstance(fragments, list)
        # May or may not have fragments depending on implementation
        # Just verify it returns a list

    def test_neutralize_charges(self):
        """Test charge neutralization"""
        charged = "C[NH3+]"
        neutralized = self.editor.neutralize_charges(charged)

        # Should succeed (though exact result depends on neutralization logic)
        assert neutralized is not None


class TestADMETPredictor:
    """Tests for ADMETPredictor class"""

    def setup_method(self):
        """Setup test fixtures"""
        self.predictor = ADMETPredictor()
        self.aspirin = "CC(=O)Oc1ccccc1C(=O)O"
        self.caffeine = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"

    def test_predict_all(self):
        """Test complete ADMET prediction"""
        predictions = self.predictor.predict_all(self.aspirin)

        assert 'absorption' in predictions
        assert 'distribution' in predictions
        assert 'metabolism' in predictions
        assert 'excretion' in predictions
        assert 'toxicity' in predictions
        assert 'drug_likeness' in predictions
        assert 'blood_brain_barrier' in predictions

    def test_predict_absorption(self):
        """Test absorption prediction"""
        absorption = self.predictor.predict_absorption(self.aspirin)

        assert 'class' in absorption
        assert 'score' in absorption
        assert 'lipinski_violations' in absorption
        assert absorption['class'] in ['Low', 'Moderate', 'High']
        assert 0 <= absorption['score'] <= 1

    def test_predict_distribution(self):
        """Test distribution prediction"""
        distribution = self.predictor.predict_distribution(self.aspirin)

        assert 'volume_of_distribution' in distribution
        assert 'plasma_protein_binding' in distribution
        assert 'logP' in distribution

    def test_predict_metabolism(self):
        """Test metabolism prediction"""
        metabolism = self.predictor.predict_metabolism(self.aspirin)

        assert 'metabolic_stability' in metabolism
        assert 'cyp_substrate_prediction' in metabolism
        assert metabolism['metabolic_stability']['class'] in ['Low', 'Moderate', 'High']

    def test_predict_excretion(self):
        """Test excretion prediction"""
        excretion = self.predictor.predict_excretion(self.aspirin)

        assert 'renal_excretion' in excretion
        assert 'biliary_excretion' in excretion
        assert 'clearance' in excretion

    def test_predict_toxicity(self):
        """Test toxicity prediction"""
        toxicity = self.predictor.predict_toxicity(self.aspirin)

        assert 'overall_risk' in toxicity
        assert 'structural_alerts' in toxicity
        assert 'cardiotoxicity_hERG_risk' in toxicity
        assert toxicity['overall_risk']['class'] in ['Low', 'Moderate', 'High']

    def test_predict_drug_likeness(self):
        """Test drug-likeness prediction"""
        drug_likeness = self.predictor.predict_drug_likeness(self.aspirin)

        assert 'drug_likeness' in drug_likeness
        assert 'score' in drug_likeness
        assert 'lipinski_rule_of_five' in drug_likeness
        assert 0 <= drug_likeness['score'] <= 100

    def test_predict_bbb(self):
        """Test BBB penetration prediction"""
        bbb = self.predictor.predict_bbb_penetration(self.caffeine)

        assert 'penetration' in bbb
        assert 'score' in bbb
        assert 'cns_mpo_score' in bbb
        assert bbb['penetration'] in ['Low', 'Moderate', 'High']

    def test_invalid_smiles(self):
        """Test predictions with invalid SMILES"""
        predictions = self.predictor.predict_all("Invalid")

        assert 'error' in predictions


class TestIntegration:
    """Integration tests combining multiple modules"""

    def setup_method(self):
        """Setup test fixtures"""
        self.visualizer = MolecularVisualizer()
        self.editor = MolecularEditor()
        self.predictor = ADMETPredictor()

    def test_full_workflow(self):
        """Test complete drug discovery workflow"""
        # Start with a SMILES
        parent_smiles = "CC(=O)O"

        # 1. Validate and canonicalize
        canonical = self.editor.canonicalize(parent_smiles)
        assert canonical is not None

        # 2. Calculate properties
        props = self.visualizer.calculate_properties(canonical)
        assert 'molecular_weight' in props

        # 3. Predict ADMET
        admet = self.predictor.predict_all(canonical)
        assert 'absorption' in admet

        # 4. Check drug-likeness
        assert admet['drug_likeness']['score'] > 0

    def test_analog_workflow(self):
        """Test analog generation and comparison"""
        parent = "c1ccccc1"  # Benzene

        # Generate analogs
        analogs = self.editor.generate_analogs(parent, num_analogs=3)

        # Should generate some analogs
        assert isinstance(analogs, list)

        # Each analog should have required fields
        for analog in analogs:
            assert 'smiles' in analog
            assert 'similarity' in analog
            assert 0 <= analog['similarity'] <= 1


if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])

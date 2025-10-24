#!/usr/bin/env python3
"""
ADMET Prediction Module for PharmaSight Platform

Predicts Absorption, Distribution, Metabolism, Excretion, and Toxicity
properties using RDKit molecular descriptors and rule-based approaches.

Note: These are simplified heuristic predictions. For production use,
consider integrating machine learning models or commercial ADMET tools.
"""

from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, Crippen, MolSurf
from rdkit.Chem import rdMolDescriptors
from typing import Dict, Optional, List, Tuple
import math


class ADMETPredictor:
    """
    ADMET prediction using molecular descriptors and empirical rules
    """

    def __init__(self):
        """Initialize ADMET predictor"""
        pass

    def predict_all(self, smiles: str) -> Dict:
        """
        Predict all ADMET properties

        Args:
            smiles: SMILES string

        Returns:
            Dictionary with all ADMET predictions
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {'error': 'Invalid SMILES string'}

        return {
            'absorption': self.predict_absorption(smiles),
            'distribution': self.predict_distribution(smiles),
            'metabolism': self.predict_metabolism(smiles),
            'excretion': self.predict_excretion(smiles),
            'toxicity': self.predict_toxicity(smiles),
            'drug_likeness': self.predict_drug_likeness(smiles),
            'blood_brain_barrier': self.predict_bbb_penetration(smiles)
        }

    def predict_absorption(self, smiles: str) -> Dict:
        """
        Predict oral absorption properties

        Based on:
        - Lipinski's Rule of Five
        - Veber's rules (rotatable bonds, TPSA)
        - Intestinal absorption models
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {'error': 'Invalid SMILES'}

        mw = Descriptors.MolWt(mol)
        logp = Crippen.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        tpsa = MolSurf.TPSA(mol)
        rotatable = Descriptors.NumRotatableBonds(mol)

        # Lipinski violations
        lipinski_violations = sum([
            mw > 500,
            logp > 5,
            hbd > 5,
            hba > 10
        ])

        # Veber violations
        veber_violations = sum([
            rotatable > 10,
            tpsa > 140
        ])

        # Predict oral bioavailability (simplified)
        if lipinski_violations == 0 and veber_violations == 0:
            absorption_class = "High"
            absorption_score = 0.85 + (0.15 * (1 - tpsa / 200))
        elif lipinski_violations <= 1 and veber_violations <= 1:
            absorption_class = "Moderate"
            absorption_score = 0.50 + (0.35 * (1 - tpsa / 200))
        else:
            absorption_class = "Low"
            absorption_score = 0.30 * (1 - tpsa / 200)

        # Intestinal absorption prediction (0-100%)
        # Based on TPSA (Polar Surface Area)
        if tpsa <= 60:
            intestinal_absorption = 95
        elif tpsa <= 140:
            intestinal_absorption = 95 - ((tpsa - 60) * 0.5)
        else:
            intestinal_absorption = max(0, 55 - ((tpsa - 140) * 0.4))

        return {
            'class': absorption_class,
            'score': round(absorption_score, 3),
            'intestinal_absorption_percent': round(intestinal_absorption, 1),
            'lipinski_violations': lipinski_violations,
            'veber_violations': veber_violations,
            'tpsa': round(tpsa, 2),
            'rotatable_bonds': rotatable,
            'recommendation': self._get_absorption_recommendation(
                absorption_class, lipinski_violations, tpsa
            )
        }

    def predict_distribution(self, smiles: str) -> Dict:
        """
        Predict distribution properties including volume of distribution
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {'error': 'Invalid SMILES'}

        logp = Crippen.MolLogP(mol)
        mw = Descriptors.MolWt(mol)
        hba = Descriptors.NumHAcceptors(mol)
        hbd = Descriptors.NumHDonors(mol)

        # Plasma protein binding prediction
        # High LogP and MW correlate with high protein binding
        if logp > 3 and mw > 400:
            ppb_class = "High"
            ppb_percent = 90 + min(10, (logp - 3) * 2)
        elif logp > 1.5:
            ppb_class = "Moderate"
            ppb_percent = 70 + (logp - 1.5) * 10
        else:
            ppb_class = "Low"
            ppb_percent = 50 + logp * 10

        # Volume of distribution prediction (simplified)
        # Based on logP and molecular size
        if logp > 3:
            vd_class = "High"  # Extensive tissue distribution
            vd_l_kg = 3 + logp
        elif logp > 0:
            vd_class = "Moderate"
            vd_l_kg = 0.5 + logp
        else:
            vd_class = "Low"  # Mainly in plasma
            vd_l_kg = 0.2 + abs(logp) * 0.1

        return {
            'volume_of_distribution': {
                'class': vd_class,
                'value_L_per_kg': round(vd_l_kg, 2)
            },
            'plasma_protein_binding': {
                'class': ppb_class,
                'percent': round(min(99, max(10, ppb_percent)), 1)
            },
            'tissue_distribution': 'Extensive' if logp > 2 else 'Limited',
            'logP': round(logp, 2)
        }

    def predict_metabolism(self, smiles: str) -> Dict:
        """
        Predict metabolic liability and CYP interaction

        Identifies potential sites of metabolism and CYP substrates
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {'error': 'Invalid SMILES'}

        # Common CYP substrate patterns (SMARTS)
        cyp_patterns = {
            'CYP3A4': [
                '[CH3][CH2]',  # N-dealkylation
                'c1ccccc1',    # Aromatic hydroxylation
                '[NX3;H2,H1;!$(NC=O)]'  # N-oxidation
            ],
            'CYP2D6': [
                '[NX3;H2,H1;!$(NC=O)]',  # Basic nitrogen
                'c1ccccc1[CH2][NH]'  # Aromatic amine
            ],
            'CYP2C9': [
                'c1ccccc1',    # Aromatic rings
                'CC(=O)O',     # Carboxylic acids
            ],
            'CYP1A2': [
                'c1ccc2ccccc2c1',  # Polycyclic aromatics
                '[nH]1cccc1'   # Indoles
            ]
        }

        cyp_substrate = {}
        for cyp, patterns in cyp_patterns.items():
            matches = 0
            for pattern in patterns:
                patt = Chem.MolFromSmarts(pattern)
                if patt and mol.HasSubstructMatch(patt):
                    matches += 1

            if matches >= 2:
                cyp_substrate[cyp] = "High probability"
            elif matches == 1:
                cyp_substrate[cyp] = "Moderate probability"
            else:
                cyp_substrate[cyp] = "Low probability"

        # Estimate metabolic stability
        num_aromatic = Descriptors.NumAromaticRings(mol)
        num_aliphatic = Descriptors.NumAliphaticRings(mol)
        num_heteroatoms = Descriptors.NumHeteroatoms(mol)

        stability_score = 50
        stability_score += num_aromatic * 10  # Aromatic rings more stable
        stability_score -= num_heteroatoms * 3  # Heteroatoms = metabolism sites
        stability_score -= Descriptors.NumRotatableBonds(mol) * 2

        stability_score = max(0, min(100, stability_score))

        if stability_score > 70:
            stability_class = "High"
        elif stability_score > 40:
            stability_class = "Moderate"
        else:
            stability_class = "Low"

        return {
            'metabolic_stability': {
                'class': stability_class,
                'score': round(stability_score, 1)
            },
            'cyp_substrate_prediction': cyp_substrate,
            'phase_1_sites': num_heteroatoms + Descriptors.NumAromaticRings(mol),
            'phase_2_sites': Descriptors.NumHDonors(mol)
        }

    def predict_excretion(self, smiles: str) -> Dict:
        """
        Predict excretion properties
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {'error': 'Invalid SMILES'}

        mw = Descriptors.MolWt(mol)
        logp = Crippen.MolLogP(mol)
        tpsa = MolSurf.TPSA(mol)

        # Renal excretion prediction
        # Low MW, hydrophilic compounds favor renal excretion
        if mw < 300 and logp < 0:
            renal_class = "High"
            renal_percent = 70
        elif mw < 500 and logp < 2:
            renal_class = "Moderate"
            renal_percent = 40
        else:
            renal_class = "Low"
            renal_percent = 15

        # Biliary excretion (high MW, lipophilic)
        if mw > 500 and logp > 2:
            biliary_class = "High"
        elif mw > 400:
            biliary_class = "Moderate"
        else:
            biliary_class = "Low"

        # Clearance estimation (simplified)
        clearance_class = "High" if logp < 1 else "Moderate" if logp < 3 else "Low"

        return {
            'renal_excretion': {
                'class': renal_class,
                'percent': renal_percent
            },
            'biliary_excretion': biliary_class,
            'clearance': clearance_class,
            'half_life_prediction': 'Short' if mw < 300 else 'Moderate' if mw < 500 else 'Long'
        }

    def predict_toxicity(self, smiles: str) -> Dict:
        """
        Predict potential toxicity based on structural alerts

        Uses common toxicophores and PAINS filters
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {'error': 'Invalid SMILES'}

        # Toxicophore patterns (structural alerts)
        toxicophores = {
            'Epoxide': 'C1OC1',
            'Aldehyde': '[CH]=O',
            'Quinone': 'O=C1C=CC(=O)C=C1',
            'Nitro aromatic': 'c[N+](=O)[O-]',
            'Aromatic amine': 'c[NH2]',
            'Halogenated aromatic': 'c[F,Cl,Br,I]',
            'Michael acceptor': 'C=CC(=O)',
            'Aziridine': 'C1CN1',
            'Peroxide': 'OO',
            'Hydrazine': 'N[NH]'
        }

        alerts = []
        for name, smarts in toxicophores.items():
            patt = Chem.MolFromSmarts(smarts)
            if patt and mol.HasSubstructMatch(patt):
                alerts.append(name)

        # Mutagenicity concerns
        mutagenicity_patterns = ['c[N+](=O)[O-]', 'c[NH2]', 'C1OC1']
        mutagenicity_risk = sum(
            1 for p in mutagenicity_patterns
            if Chem.MolFromSmarts(p) and mol.HasSubstructMatch(Chem.MolFromSmarts(p))
        )

        # Cardiotoxicity (hERG liability)
        # Typically basic amines + aromatic rings + certain LogP range
        num_basic_nitrogens = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[NX3;H2,H1;!$(NC=O)]')))
        num_aromatic_rings = Descriptors.NumAromaticRings(mol)
        logp = Crippen.MolLogP(mol)

        if num_basic_nitrogens > 0 and num_aromatic_rings >= 2 and 2 < logp < 5:
            herg_risk = "High"
        elif num_basic_nitrogens > 0 and num_aromatic_rings >= 1:
            herg_risk = "Moderate"
        else:
            herg_risk = "Low"

        # Overall toxicity score
        toxicity_score = len(alerts) * 15 + mutagenicity_risk * 10
        if herg_risk == "High":
            toxicity_score += 20
        elif herg_risk == "Moderate":
            toxicity_score += 10

        toxicity_score = min(100, toxicity_score)

        if toxicity_score < 20:
            risk_class = "Low"
        elif toxicity_score < 50:
            risk_class = "Moderate"
        else:
            risk_class = "High"

        return {
            'overall_risk': {
                'class': risk_class,
                'score': toxicity_score
            },
            'structural_alerts': alerts,
            'mutagenicity_risk': 'High' if mutagenicity_risk >= 2 else 'Moderate' if mutagenicity_risk == 1 else 'Low',
            'cardiotoxicity_hERG_risk': herg_risk,
            'num_alerts': len(alerts),
            'recommendation': 'Requires toxicity testing' if len(alerts) > 2 else 'Standard safety assessment recommended'
        }

    def predict_drug_likeness(self, smiles: str) -> Dict:
        """
        Predict drug-likeness using multiple rules
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {'error': 'Invalid SMILES'}

        mw = Descriptors.MolWt(mol)
        logp = Crippen.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        rotatable = Descriptors.NumRotatableBonds(mol)
        tpsa = MolSurf.TPSA(mol)

        # Lipinski's Rule of Five
        lipinski = {
            'mw_ok': mw <= 500,
            'logp_ok': logp <= 5,
            'hbd_ok': hbd <= 5,
            'hba_ok': hba <= 10,
            'violations': sum([mw > 500, logp > 5, hbd > 5, hba > 10])
        }

        # Veber's Rule
        veber = {
            'rotatable_ok': rotatable <= 10,
            'tpsa_ok': tpsa <= 140,
            'violations': sum([rotatable > 10, tpsa > 140])
        }

        # Ghose Filter
        ghose = {
            'mw_ok': 160 <= mw <= 480,
            'logp_ok': -0.4 <= logp <= 5.6,
            'atoms_ok': 20 <= mol.GetNumAtoms() <= 70,
            'violations': sum([
                not (160 <= mw <= 480),
                not (-0.4 <= logp <= 5.6),
                not (20 <= mol.GetNumAtoms() <= 70)
            ])
        }

        # Overall drug-likeness score
        score = 100
        score -= lipinski['violations'] * 15
        score -= veber['violations'] * 10
        score -= ghose['violations'] * 5

        if score >= 80:
            drug_likeness = "High"
        elif score >= 50:
            drug_likeness = "Moderate"
        else:
            drug_likeness = "Low"

        return {
            'drug_likeness': drug_likeness,
            'score': max(0, score),
            'lipinski_rule_of_five': lipinski,
            'veber_rule': veber,
            'ghose_filter': ghose,
            'passes_all_filters': (
                lipinski['violations'] == 0 and
                veber['violations'] == 0 and
                ghose['violations'] == 0
            )
        }

    def predict_bbb_penetration(self, smiles: str) -> Dict:
        """
        Predict blood-brain barrier penetration
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {'error': 'Invalid SMILES'}

        mw = Descriptors.MolWt(mol)
        logp = Crippen.MolLogP(mol)
        tpsa = MolSurf.TPSA(mol)
        hbd = Descriptors.NumHDonors(mol)

        # BBB penetration rules
        # TPSA < 90, MW < 450, LogP 1-3 optimal
        if tpsa < 60 and mw < 400 and 1 < logp < 3:
            bbb_class = "High"
            bbb_score = 0.9
        elif tpsa < 90 and mw < 450 and 0 < logp < 5:
            bbb_class = "Moderate"
            bbb_score = 0.6
        else:
            bbb_class = "Low"
            bbb_score = 0.2

        # CNS MPO Score (simplified version)
        # Considers LogP, LogD, TPSA, HBD, pKa, MW
        cns_score = 0
        if logp <= 3:
            cns_score += 1
        if tpsa <= 90:
            cns_score += 1
        if hbd <= 3:
            cns_score += 1
        if mw <= 360:
            cns_score += 1

        return {
            'penetration': bbb_class,
            'score': bbb_score,
            'cns_mpo_score': cns_score,
            'tpsa': round(tpsa, 2),
            'suitable_for_cns': cns_score >= 3,
            'recommendation': 'Good CNS candidate' if cns_score >= 3 else 'May have limited CNS exposure'
        }

    def _get_absorption_recommendation(self, absorption_class: str,
                                      lipinski_violations: int,
                                      tpsa: float) -> str:
        """Generate absorption recommendation"""
        if absorption_class == "High":
            return "Expected good oral bioavailability"
        elif absorption_class == "Moderate":
            if lipinski_violations > 0:
                return "Consider reducing molecular weight or lipophilicity"
            elif tpsa > 120:
                return "Consider reducing polar surface area"
            else:
                return "May have acceptable oral bioavailability"
        else:
            recommendations = []
            if lipinski_violations > 1:
                recommendations.append("significantly reduce MW and LogP")
            if tpsa > 140:
                recommendations.append("reduce polar surface area")

            if recommendations:
                return f"Poor oral absorption likely. Consider: {', '.join(recommendations)}"
            else:
                return "Poor oral absorption predicted"


def main():
    """Example usage"""
    predictor = ADMETPredictor()

    print("=" * 70)
    print("  ADMET Prediction Examples")
    print("=" * 70 + "\n")

    test_compounds = {
        "Aspirin": "CC(=O)Oc1ccccc1C(=O)O",
        "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "Psilocybin": "CN(C)CCc1c[nH]c2ccc(OP(=O)(O)O)cc12",
        "Ibuprofen": "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
    }

    for name, smiles in test_compounds.items():
        print(f"\n{name}: {smiles}")
        print("-" * 70)

        predictions = predictor.predict_all(smiles)

        print(f"\nAbsorption: {predictions['absorption']['class']}")
        print(f"  Intestinal Absorption: {predictions['absorption']['intestinal_absorption_percent']}%")
        print(f"  Lipinski Violations: {predictions['absorption']['lipinski_violations']}")

        print(f"\nDistribution:")
        print(f"  Volume of Distribution: {predictions['distribution']['volume_of_distribution']['class']}")
        print(f"  Plasma Protein Binding: {predictions['distribution']['plasma_protein_binding']['percent']}%")

        print(f"\nMetabolism:")
        print(f"  Metabolic Stability: {predictions['metabolism']['metabolic_stability']['class']}")

        print(f"\nToxicity:")
        print(f"  Overall Risk: {predictions['toxicity']['overall_risk']['class']}")
        if predictions['toxicity']['structural_alerts']:
            print(f"  Alerts: {', '.join(predictions['toxicity']['structural_alerts'])}")

        print(f"\nDrug-likeness: {predictions['drug_likeness']['drug_likeness']}")
        print(f"  Score: {predictions['drug_likeness']['score']}/100")

        print(f"\nBBB Penetration: {predictions['blood_brain_barrier']['penetration']}")
        print(f"  CNS MPO Score: {predictions['blood_brain_barrier']['cns_mpo_score']}/4")

        print()


if __name__ == "__main__":
    main()

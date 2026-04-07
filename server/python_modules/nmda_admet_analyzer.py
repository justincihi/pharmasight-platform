#!/usr/bin/env python3
"""
NMDA Antagonist ADMET Analyzer
Specialized ADMET predictions and LLM-powered explanations for NMDA receptor antagonists
"""

import json
import sys
from typing import Dict, List
from sdf_processor import SDFProcessor


class NMDAADMETAnalyzer:
    """ADMET analysis specialized for NMDA antagonists"""
    
    def __init__(self):
        self.processor = SDFProcessor()
    
    def analyze_nmda_properties(self, descriptors: Dict) -> Dict:
        """
        Analyze NMDA antagonist-specific properties
        
        Args:
            descriptors: Molecular descriptors from SDF processor
        
        Returns:
            NMDA-specific ADMET predictions
        """
        admet = {}
        
        # BBB Penetration (critical for CNS drugs)
        # TPSA < 90 and MW < 450 typically indicate good BBB penetration
        tpsa = descriptors.get('tpsa', 0)
        mw = descriptors.get('molecular_weight', 0)
        clogp = descriptors.get('clogp', 0)
        
        bbb_score = 100
        if tpsa > 90:
            bbb_score -= (tpsa - 90) * 0.5
        if mw > 450:
            bbb_score -= (mw - 450) * 0.1
        if clogp < 1 or clogp > 5:
            bbb_score -= 10
        
        admet['bbb_penetration_score'] = max(0, min(100, round(bbb_score, 1)))
        admet['bbb_penetration'] = 'High' if bbb_score >= 70 else ('Moderate' if bbb_score >= 40 else 'Low')
        
        # NMDA Binding Potential
        # Basic nitrogen is essential for NMDA channel blocking
        basic_n = descriptors.get('basic_nitrogen_count', 0)
        max_pos_charge = descriptors.get('max_positive_charge', 0)
        
        binding_score = 50
        if basic_n >= 1:
            binding_score += 30
        if max_pos_charge and max_pos_charge > 0.1:
            binding_score += 20
        
        admet['nmda_binding_potential'] = max(0, min(100, round(binding_score, 1)))
        admet['nmda_binding'] = 'High' if binding_score >= 70 else ('Moderate' if binding_score >= 40 else 'Low')
        
        # Dissociative Effect Risk
        # Higher lipophilicity and certain structural features increase dissociative effects
        dissociative_risk = 30
        if clogp > 3.5:
            dissociative_risk += (clogp - 3.5) * 10
        if descriptors.get('num_aromatic_rings', 0) > 1:
            dissociative_risk += 10
        
        admet['dissociative_risk_score'] = max(0, min(100, round(dissociative_risk, 1)))
        admet['dissociative_risk'] = 'High' if dissociative_risk >= 60 else ('Moderate' if dissociative_risk >= 30 else 'Low')
        
        # Metabolic Stability
        # Fewer rotatable bonds and aromatic rings typically indicate better stability
        rotatable = descriptors.get('rotatable_bonds', 0)
        stability_score = 80
        if rotatable > 5:
            stability_score -= (rotatable - 5) * 5
        
        admet['metabolic_stability_score'] = max(0, min(100, round(stability_score, 1)))
        admet['metabolic_stability'] = 'High' if stability_score >= 70 else ('Moderate' if stability_score >= 40 else 'Low')
        
        # Hepatotoxicity Risk
        # Based on molecular weight, cLogP, and aromatic rings
        hepatotox_risk = 20
        if mw > 400:
            hepatotox_risk += 10
        if clogp > 4:
            hepatotox_risk += 15
        if descriptors.get('num_aromatic_rings', 0) > 2:
            hepatotox_risk += 10
        
        admet['hepatotoxicity_risk_score'] = max(0, min(100, round(hepatotox_risk, 1)))
        admet['hepatotoxicity_risk'] = 'High' if hepatotox_risk >= 60 else ('Moderate' if hepatotox_risk >= 30 else 'Low')
        
        # Oral Bioavailability
        # Lipinski's rule of five + additional factors
        lipinski_violations = descriptors.get('lipinski_violations', 0)
        bioavailability_score = 90 - (lipinski_violations * 20)
        if tpsa > 140:
            bioavailability_score -= 20
        
        admet['oral_bioavailability_score'] = max(0, min(100, round(bioavailability_score, 1)))
        admet['oral_bioavailability'] = 'High' if bioavailability_score >= 70 else ('Moderate' if bioavailability_score >= 40 else 'Low')
        
        # Overall ADMET Score
        admet['overall_admet_score'] = round((
            admet['bbb_penetration_score'] * 0.25 +
            admet['nmda_binding_potential'] * 0.25 +
            (100 - admet['dissociative_risk_score']) * 0.15 +
            admet['metabolic_stability_score'] * 0.15 +
            (100 - admet['hepatotoxicity_risk_score']) * 0.10 +
            admet['oral_bioavailability_score'] * 0.10
        ), 1)
        
        return admet
    
    def generate_admet_explanation(self, compound_name: str, smiles: str, 
                                   descriptors: Dict, admet: Dict, 
                                   sdf_description: str = "") -> str:
        """
        Generate human-readable ADMET explanation
        
        Args:
            compound_name: Compound identifier
            smiles: SMILES string
            descriptors: Molecular descriptors
            admet: ADMET predictions
            sdf_description: Description from SDF file
        
        Returns:
            Formatted explanation text
        """
        explanation = f"# ADMET Analysis: {compound_name}\n\n"
        
        # Mechanism section
        explanation += "## Mechanism of Action\n"
        if sdf_description:
            explanation += f"{sdf_description}\n\n"
        else:
            explanation += "NMDA receptor antagonist with potential for treating depression, PTSD, and chronic pain.\n\n"
        
        # Key Properties
        explanation += "## Key Molecular Properties\n"
        explanation += f"- **Molecular Weight**: {descriptors.get('molecular_weight', 'N/A')} g/mol\n"
        explanation += f"- **cLogP**: {descriptors.get('clogp', 'N/A')} (lipophilicity)\n"
        explanation += f"- **TPSA**: {descriptors.get('tpsa', 'N/A')} Ų (polar surface area)\n"
        explanation += f"- **H-Bond Donors**: {descriptors.get('hbd', 'N/A')}\n"
        explanation += f"- **H-Bond Acceptors**: {descriptors.get('hba', 'N/A')}\n"
        explanation += f"- **Basic Nitrogens**: {descriptors.get('basic_nitrogen_count', 'N/A')}\n"
        explanation += f"- **Drug-likeness Score**: {descriptors.get('drug_likeness_score', 'N/A')}/100\n\n"
        
        # ADMET Predictions
        explanation += "## ADMET Predictions\n\n"
        
        explanation += f"### Blood-Brain Barrier Penetration: {admet['bbb_penetration']} ({admet['bbb_penetration_score']}/100)\n"
        if admet['bbb_penetration_score'] >= 70:
            explanation += "Excellent CNS penetration expected. Low TPSA and optimal molecular weight facilitate passive diffusion across the BBB.\n\n"
        elif admet['bbb_penetration_score'] >= 40:
            explanation += "Moderate CNS penetration. May require optimization of lipophilicity or polar surface area.\n\n"
        else:
            explanation += "Poor CNS penetration predicted. High TPSA or molecular weight may limit BBB crossing.\n\n"
        
        explanation += f"### NMDA Binding Potential: {admet['nmda_binding']} ({admet['nmda_binding_potential']}/100)\n"
        if admet['nmda_binding_potential'] >= 70:
            explanation += "Strong NMDA channel blocking activity predicted. Presence of basic nitrogen and positive charge distribution favor binding.\n\n"
        else:
            explanation += "Moderate NMDA binding predicted. Structural modifications may enhance receptor affinity.\n\n"
        
        explanation += f"### Dissociative Effect Risk: {admet['dissociative_risk']} ({admet['dissociative_risk_score']}/100)\n"
        if admet['dissociative_risk_score'] < 30:
            explanation += "Low risk of dissociative side effects. Favorable lipophilicity and structural features suggest improved tolerability.\n\n"
        elif admet['dissociative_risk_score'] < 60:
            explanation += "Moderate dissociative risk. Careful dose titration recommended.\n\n"
        else:
            explanation += "High dissociative risk. May require structural modifications to reduce psychotomimetic effects.\n\n"
        
        explanation += f"### Metabolic Stability: {admet['metabolic_stability']} ({admet['metabolic_stability_score']}/100)\n"
        explanation += f"### Hepatotoxicity Risk: {admet['hepatotoxicity_risk']} ({admet['hepatotoxicity_risk_score']}/100)\n"
        explanation += f"### Oral Bioavailability: {admet['oral_bioavailability']} ({admet['oral_bioavailability_score']}/100)\n\n"
        
        # Overall Assessment
        explanation += f"## Overall ADMET Score: {admet['overall_admet_score']}/100\n\n"
        if admet['overall_admet_score'] >= 75:
            explanation += "**Excellent candidate** for further development. Strong ADMET profile with favorable CNS penetration and safety characteristics.\n"
        elif admet['overall_admet_score'] >= 50:
            explanation += "**Promising candidate** with room for optimization. Consider structural modifications to improve specific ADMET parameters.\n"
        else:
            explanation += "**Requires optimization** before advancing. Significant ADMET liabilities identified.\n"
        
        return explanation
    
    def process_sdf_with_admet(self, sdf_path: str, output_dir: str = None) -> List[Dict]:
        """
        Complete pipeline: parse SDF, compute descriptors, analyze ADMET, generate explanations
        
        Args:
            sdf_path: Input SDF file path
            output_dir: Optional output directory
        
        Returns:
            List of analogs with ADMET analysis
        """
        # Process SDF and get descriptors
        analogs = self.processor.process_analog_sdf(sdf_path, output_dir)
        
        # Add ADMET analysis to each analog
        for analog in analogs:
            # Compute ADMET predictions
            admet = self.analyze_nmda_properties(analog)
            analog.update(admet)
            
            # Generate explanation
            sdf_desc = analog.get('description', '')
            if not sdf_desc:
                # Try to get description from SDF properties
                props = self.processor.molecules[0]['properties'] if self.processor.molecules else {}
                sdf_desc = props.get('THERAPEUTIC_POTENTIAL', '') or props.get('KEY_DIFFERENCES', '')
            
            explanation = self.generate_admet_explanation(
                analog['compound_name'],
                analog['smiles'],
                analog,
                admet,
                sdf_desc
            )
            analog['admet_explanation'] = explanation
        
        return analogs


def analyze_sdf_file(sdf_path: str, output_dir: str = None) -> str:
    """
    Convenience function to analyze SDF file and return JSON with ADMET
    
    Args:
        sdf_path: Input SDF file path
        output_dir: Optional output directory
    
    Returns:
        JSON string of analogs with ADMET analysis
    """
    analyzer = NMDAADMETAnalyzer()
    analogs = analyzer.process_sdf_with_admet(sdf_path, output_dir)
    return json.dumps(analogs, indent=2)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python nmda_admet_analyzer.py <sdf_file> [output_dir]")
        sys.exit(1)
    
    sdf_path = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) > 2 else None
    
    result = analyze_sdf_file(sdf_path, output_dir)
    
    # Parse and print formatted output
    analogs = json.loads(result)
    for analog in analogs:
        print("\n" + "="*80)
        print(analog.get('admet_explanation', ''))
        print("\n" + "="*80)
        print("\nFull JSON Output:")
        print(json.dumps(analog, indent=2))

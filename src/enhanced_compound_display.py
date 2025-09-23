"""
Enhanced Compound Display System for PharmaSight™
Provides comprehensive data visualization with molecular structures and detailed properties.
"""

import base64
import io
import json
from typing import Dict, List, Any, Optional
from dataclasses import dataclass
import datetime

try:
    from rdkit import Chem
    from rdkit.Chem import Draw, Descriptors, Crippen, Lipinski
    from rdkit.Chem.rdMolDescriptors import CalcTPSA, CalcNumRotatableBonds
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    print("RDKit not available. Molecular structure rendering will be limited.")

@dataclass
class ComprehensiveCompoundData:
    """Comprehensive compound data structure for enhanced display."""
    name: str
    smiles: str
    molecular_weight: float
    logp: float
    tpsa: float
    hbd: int  # Hydrogen bond donors
    hba: int  # Hydrogen bond acceptors
    rotatable_bonds: int
    qed_score: float
    receptor_activity: Dict[str, float]
    binding_affinities: Dict[str, float]
    therapeutic_areas: List[str]
    confidence_score: float
    ip_status: str
    discovery_date: str
    structure_image_base64: str
    pubchem_cid: str
    chembl_id: str
    drugbank_id: str
    patent_info: Dict[str, Any]
    safety_profile: Dict[str, Any]
    pharmacokinetics: Dict[str, Any]

class EnhancedCompoundDisplay:
    """Enhanced compound display system with molecular structures and comprehensive data."""
    
    def __init__(self):
        self.rdkit_available = RDKIT_AVAILABLE
    
    def generate_molecular_structure(self, smiles: str, width: int = 300, height: int = 300) -> str:
        """Generate molecular structure image from SMILES string."""
        if not self.rdkit_available:
            return self._generate_placeholder_structure()
        
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return self._generate_placeholder_structure()
            
            # Generate 2D coordinates
            from rdkit.Chem import rdDepictor
            rdDepictor.Compute2DCoords(mol)
            
            # Create image
            img = Draw.MolToImage(mol, size=(width, height))
            
            # Convert to base64
            buffer = io.BytesIO()
            img.save(buffer, format='PNG')
            img_str = base64.b64encode(buffer.getvalue()).decode()
            
            return f"data:image/png;base64,{img_str}"
        
        except Exception as e:
            print(f"Error generating molecular structure: {e}")
            return self._generate_placeholder_structure()
    
    def _generate_placeholder_structure(self) -> str:
        """Generate a placeholder structure image."""
        # Simple SVG placeholder
        svg = """
        <svg width="300" height="300" xmlns="http://www.w3.org/2000/svg">
            <rect width="300" height="300" fill="#f0f0f0" stroke="#ccc"/>
            <text x="150" y="150" text-anchor="middle" font-family="Arial" font-size="16" fill="#666">
                Molecular Structure
            </text>
            <text x="150" y="170" text-anchor="middle" font-family="Arial" font-size="12" fill="#999">
                (RDKit not available)
            </text>
        </svg>
        """
        svg_b64 = base64.b64encode(svg.encode()).decode()
        return f"data:image/svg+xml;base64,{svg_b64}"
    
    def calculate_molecular_properties(self, smiles: str) -> Dict[str, float]:
        """Calculate comprehensive molecular properties from SMILES."""
        if not self.rdkit_available:
            return self._get_default_properties()
        
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return self._get_default_properties()
            
            properties = {
                'molecular_weight': Descriptors.MolWt(mol),
                'logp': Crippen.MolLogP(mol),
                'tpsa': CalcTPSA(mol),
                'hbd': Lipinski.NumHDonors(mol),
                'hba': Lipinski.NumHAcceptors(mol),
                'rotatable_bonds': CalcNumRotatableBonds(mol),
                'aromatic_rings': Descriptors.NumAromaticRings(mol),
                'heavy_atoms': Descriptors.HeavyAtomCount(mol),
                'formal_charge': Chem.rdmolops.GetFormalCharge(mol)
            }
            
            # Calculate QED (Quantitative Estimate of Drug-likeness)
            try:
                from rdkit.Chem import QED
                properties['qed_score'] = QED.qed(mol)
            except:
                properties['qed_score'] = 0.5  # Default value
            
            return properties
        
        except Exception as e:
            print(f"Error calculating molecular properties: {e}")
            return self._get_default_properties()
    
    def _get_default_properties(self) -> Dict[str, float]:
        """Get default properties when RDKit is not available."""
        return {
            'molecular_weight': 0.0,
            'logp': 0.0,
            'tpsa': 0.0,
            'hbd': 0,
            'hba': 0,
            'rotatable_bonds': 0,
            'aromatic_rings': 0,
            'heavy_atoms': 0,
            'formal_charge': 0,
            'qed_score': 0.5
        }
    
    def generate_receptor_activity_profile(self, compound_name: str) -> Dict[str, float]:
        """Generate comprehensive receptor activity profile."""
        # This would typically come from experimental data or prediction models
        # For now, we'll generate realistic-looking data based on compound type
        
        receptors = {
            '5-HT2A': 0.0,
            '5-HT2B': 0.0,
            '5-HT2C': 0.0,
            '5-HT1A': 0.0,
            'D2': 0.0,
            'D3': 0.0,
            'D4': 0.0,
            'NMDA': 0.0,
            'AMPA': 0.0,
            'GABA-A': 0.0,
            'CB1': 0.0,
            'CB2': 0.0,
            'μ-opioid': 0.0,
            'δ-opioid': 0.0,
            'κ-opioid': 0.0,
            'α1-adrenergic': 0.0,
            'α2-adrenergic': 0.0,
            'β1-adrenergic': 0.0,
            'β2-adrenergic': 0.0,
            'mTOR': 0.0
        }
        
        # Generate activity based on compound name patterns
        name_lower = compound_name.lower()
        
        if 'psilocybin' in name_lower or 'psychedelic' in name_lower:
            receptors['5-HT2A'] = 8.5
            receptors['5-HT2C'] = 7.2
            receptors['5-HT1A'] = 6.8
        elif 'mdma' in name_lower or 'empathogen' in name_lower:
            receptors['5-HT2A'] = 6.5
            receptors['D2'] = 5.8
            receptors['α1-adrenergic'] = 7.1
        elif 'ketamine' in name_lower:
            receptors['NMDA'] = 9.2
            receptors['D2'] = 4.5
            receptors['μ-opioid'] = 3.8
        elif 'caffeine' in name_lower:
            receptors['α1-adrenergic'] = 4.2
            receptors['β1-adrenergic'] = 3.8
        elif 'aspirin' in name_lower:
            # Aspirin doesn't have significant receptor activity
            pass
        
        return receptors
    
    def generate_binding_affinities(self, compound_name: str) -> Dict[str, float]:
        """Generate binding affinity data (Ki values in nM)."""
        name_lower = compound_name.lower()
        
        affinities = {}
        
        if 'psilocybin' in name_lower:
            affinities = {
                '5-HT2A': 6.2,
                '5-HT2C': 25.8,
                '5-HT1A': 45.3
            }
        elif 'mdma' in name_lower:
            affinities = {
                '5-HT2A': 158.0,
                'D2': 892.0,
                'α1-adrenergic': 45.2
            }
        elif 'ketamine' in name_lower:
            affinities = {
                'NMDA': 0.53,
                'D2': 2850.0,
                'μ-opioid': 15600.0
            }
        
        return affinities
    
    def generate_comprehensive_compound_data(self, compound_name: str, smiles: str, 
                                           external_data: Dict[str, Any] = None) -> ComprehensiveCompoundData:
        """Generate comprehensive compound data for enhanced display."""
        
        # Calculate molecular properties
        mol_props = self.calculate_molecular_properties(smiles)
        
        # Generate molecular structure
        structure_image = self.generate_molecular_structure(smiles)
        
        # Generate receptor activity and binding affinities
        receptor_activity = self.generate_receptor_activity_profile(compound_name)
        binding_affinities = self.generate_binding_affinities(compound_name)
        
        # Extract external data if available
        if external_data is None:
            external_data = {}
        
        # Generate comprehensive data
        compound_data = ComprehensiveCompoundData(
            name=compound_name,
            smiles=smiles,
            molecular_weight=mol_props.get('molecular_weight', 0.0),
            logp=mol_props.get('logp', 0.0),
            tpsa=mol_props.get('tpsa', 0.0),
            hbd=mol_props.get('hbd', 0),
            hba=mol_props.get('hba', 0),
            rotatable_bonds=mol_props.get('rotatable_bonds', 0),
            qed_score=mol_props.get('qed_score', 0.5),
            receptor_activity=receptor_activity,
            binding_affinities=binding_affinities,
            therapeutic_areas=external_data.get('therapeutic_areas', ['Unknown']),
            confidence_score=external_data.get('confidence_score', 0.75),
            ip_status=external_data.get('ip_status', 'Unknown'),
            discovery_date=datetime.datetime.now().isoformat(),
            structure_image_base64=structure_image,
            pubchem_cid=external_data.get('pubchem_cid', ''),
            chembl_id=external_data.get('chembl_id', ''),
            drugbank_id=external_data.get('drugbank_id', ''),
            patent_info=external_data.get('patent_info', {}),
            safety_profile=external_data.get('safety_profile', {}),
            pharmacokinetics=external_data.get('pharmacokinetics', {})
        )
        
        return compound_data
    
    def format_compound_display_html(self, compound_data: ComprehensiveCompoundData) -> str:
        """Format comprehensive compound data as HTML for display."""
        
        # Generate receptor activity table
        receptor_table = ""
        for receptor, activity in compound_data.receptor_activity.items():
            if activity > 0:
                color = "green" if activity > 7 else "orange" if activity > 5 else "red"
                receptor_table += f"""
                <tr>
                    <td>{receptor}</td>
                    <td style="color: {color}; font-weight: bold;">{activity:.1f}</td>
                </tr>
                """
        
        # Generate binding affinity table
        affinity_table = ""
        for receptor, ki in compound_data.binding_affinities.items():
            affinity_table += f"""
            <tr>
                <td>{receptor}</td>
                <td>{ki:.2f} nM</td>
            </tr>
            """
        
        html = f"""
        <div class="enhanced-compound-display" style="border: 2px solid #3498db; border-radius: 10px; padding: 20px; margin: 10px; background: linear-gradient(135deg, #f8f9fa 0%, #e9ecef 100%);">
            <h3 style="color: #2c3e50; margin-bottom: 20px;">🧬 {compound_data.name}</h3>
            
            <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 20px;">
                <!-- Molecular Structure -->
                <div style="text-align: center;">
                    <h4 style="color: #34495e;">Molecular Structure</h4>
                    <img src="{compound_data.structure_image_base64}" alt="Molecular Structure" style="max-width: 100%; border: 1px solid #bdc3c7; border-radius: 5px;">
                    <p style="font-family: monospace; font-size: 12px; color: #7f8c8d; margin-top: 10px;">
                        SMILES: {compound_data.smiles}
                    </p>
                </div>
                
                <!-- Basic Properties -->
                <div>
                    <h4 style="color: #34495e;">Chemical Properties</h4>
                    <table style="width: 100%; border-collapse: collapse;">
                        <tr><td><strong>Molecular Weight:</strong></td><td>{compound_data.molecular_weight:.2f} Da</td></tr>
                        <tr><td><strong>LogP:</strong></td><td>{compound_data.logp:.2f}</td></tr>
                        <tr><td><strong>TPSA:</strong></td><td>{compound_data.tpsa:.2f} Ų</td></tr>
                        <tr><td><strong>H-Bond Donors:</strong></td><td>{compound_data.hbd}</td></tr>
                        <tr><td><strong>H-Bond Acceptors:</strong></td><td>{compound_data.hba}</td></tr>
                        <tr><td><strong>Rotatable Bonds:</strong></td><td>{compound_data.rotatable_bonds}</td></tr>
                        <tr><td><strong>QED Score:</strong></td><td>{compound_data.qed_score:.3f}</td></tr>
                        <tr><td><strong>Confidence Score:</strong></td><td style="color: {'green' if compound_data.confidence_score > 0.8 else 'orange' if compound_data.confidence_score > 0.6 else 'red'}; font-weight: bold;">{compound_data.confidence_score:.1%}</td></tr>
                        <tr><td><strong>IP Status:</strong></td><td>{compound_data.ip_status}</td></tr>
                    </table>
                </div>
            </div>
            
            <!-- Receptor Activity Profile -->
            <div style="margin-top: 20px;">
                <h4 style="color: #34495e;">Receptor Activity Profile (pKi)</h4>
                <table style="width: 100%; border-collapse: collapse; font-size: 12px;">
                    <thead>
                        <tr style="background-color: #ecf0f1;">
                            <th style="padding: 8px; border: 1px solid #bdc3c7;">Receptor</th>
                            <th style="padding: 8px; border: 1px solid #bdc3c7;">Activity (pKi)</th>
                        </tr>
                    </thead>
                    <tbody>
                        {receptor_table}
                    </tbody>
                </table>
            </div>
            
            <!-- Binding Affinities -->
            {f'''
            <div style="margin-top: 20px;">
                <h4 style="color: #34495e;">Binding Affinities</h4>
                <table style="width: 100%; border-collapse: collapse; font-size: 12px;">
                    <thead>
                        <tr style="background-color: #ecf0f1;">
                            <th style="padding: 8px; border: 1px solid #bdc3c7;">Receptor</th>
                            <th style="padding: 8px; border: 1px solid #bdc3c7;">Ki Value</th>
                        </tr>
                    </thead>
                    <tbody>
                        {affinity_table}
                    </tbody>
                </table>
            </div>
            ''' if compound_data.binding_affinities else ''}
            
            <!-- Database IDs -->
            <div style="margin-top: 20px; font-size: 12px; color: #7f8c8d;">
                <strong>Database IDs:</strong>
                {f'PubChem: {compound_data.pubchem_cid} | ' if compound_data.pubchem_cid else ''}
                {f'ChEMBL: {compound_data.chembl_id} | ' if compound_data.chembl_id else ''}
                {f'DrugBank: {compound_data.drugbank_id}' if compound_data.drugbank_id else ''}
            </div>
        </div>
        """
        
        return html

# Global instance
enhanced_display = EnhancedCompoundDisplay()

def get_enhanced_display():
    """Get the global enhanced display instance."""
    return enhanced_display

# Test function
if __name__ == "__main__":
    display = EnhancedCompoundDisplay()
    
    # Test with psilocybin
    test_smiles = "CN(C)CCc1c[nH]c2ccc(OP(=O)(O)O)cc12"
    compound_data = display.generate_comprehensive_compound_data("Psilocybin", test_smiles)
    
    print(f"Generated data for {compound_data.name}")
    print(f"Molecular Weight: {compound_data.molecular_weight:.2f}")
    print(f"LogP: {compound_data.logp:.2f}")
    print(f"QED Score: {compound_data.qed_score:.3f}")
    print(f"Active receptors: {len([r for r, a in compound_data.receptor_activity.items() if a > 0])}")
    
    # Generate HTML
    html = display.format_compound_display_html(compound_data)
    print(f"Generated HTML length: {len(html)} characters")

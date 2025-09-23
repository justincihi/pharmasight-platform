"""
API endpoints for Enhanced Compound Display System
Provides REST API for generating comprehensive compound visualizations
"""

from flask import Blueprint, request, jsonify, session
from enhanced_compound_display import get_enhanced_display
from unified_search_engine import unified_search
import json

# Create blueprint for enhanced display API
enhanced_display_bp = Blueprint('enhanced_display', __name__, url_prefix='/api/display')

@enhanced_display_bp.route('/compound/<compound_name>', methods=['GET'])
def get_enhanced_compound_display(compound_name):
    """Get enhanced compound display with molecular structures and comprehensive data."""
    try:
        display = get_enhanced_display()
        
        # Get comprehensive data from unified search
        search_results = unified_search.search_compound(compound_name)
        
        # Extract SMILES from search results
        smiles = None
        external_data = {}
        
        if search_results and 'molecular_data' in search_results:
            mol_data = search_results['molecular_data']
            smiles = mol_data.get('smiles')
            
            # Extract additional data
            external_data = {
                'pubchem_cid': mol_data.get('pubchem_cid'),
                'chembl_id': mol_data.get('chembl_id'),
                'drugbank_id': mol_data.get('drugbank_id'),
                'confidence_score': search_results.get('confidence_score', 0.75),
                'therapeutic_areas': mol_data.get('therapeutic_areas', ['Unknown']),
                'ip_status': mol_data.get('patent_status', 'Unknown'),
                'patent_info': mol_data.get('patent_info', {}),
                'safety_profile': mol_data.get('safety_profile', {}),
                'pharmacokinetics': mol_data.get('pharmacokinetics', {})
            }
        
        # Fallback SMILES for known compounds
        if not smiles:
            known_smiles = {
                'psilocybin': 'CN(C)CCc1c[nH]c2ccc(OP(=O)(O)O)cc12',
                'caffeine': 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
                'aspirin': 'CC(=O)OC1=CC=CC=C1C(=O)O',
                'mdma': 'CC(CC1=CC2=C(C=C1)OCO2)NC',
                'ketamine': 'CNC1(CCCCC1=O)C2=CC=CC=C2Cl'
            }
            smiles = known_smiles.get(compound_name.lower())
        
        if not smiles:
            return jsonify({'error': 'SMILES string not found for compound'}), 404
        
        # Generate comprehensive compound data
        compound_data = display.generate_comprehensive_compound_data(
            compound_name, smiles, external_data
        )
        
        # Format as HTML
        html_display = display.format_compound_display_html(compound_data)
        
        # Convert to dict for JSON response
        compound_dict = {
            'name': compound_data.name,
            'smiles': compound_data.smiles,
            'molecular_weight': compound_data.molecular_weight,
            'logp': compound_data.logp,
            'tpsa': compound_data.tpsa,
            'hbd': compound_data.hbd,
            'hba': compound_data.hba,
            'rotatable_bonds': compound_data.rotatable_bonds,
            'qed_score': compound_data.qed_score,
            'receptor_activity': compound_data.receptor_activity,
            'binding_affinities': compound_data.binding_affinities,
            'therapeutic_areas': compound_data.therapeutic_areas,
            'confidence_score': compound_data.confidence_score,
            'ip_status': compound_data.ip_status,
            'discovery_date': compound_data.discovery_date,
            'structure_image_base64': compound_data.structure_image_base64,
            'pubchem_cid': compound_data.pubchem_cid,
            'chembl_id': compound_data.chembl_id,
            'drugbank_id': compound_data.drugbank_id,
            'html_display': html_display
        }
        
        return jsonify({
            'compound_data': compound_dict,
            'status': 'success'
        })
    
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@enhanced_display_bp.route('/batch', methods=['POST'])
def get_batch_enhanced_display():
    """Get enhanced display for multiple compounds."""
    try:
        data = request.get_json()
        compound_names = data.get('compounds', [])
        
        if not compound_names:
            return jsonify({'error': 'No compounds provided'}), 400
        
        display = get_enhanced_display()
        results = []
        
        for compound_name in compound_names:
            try:
                # Get comprehensive data from unified search
                search_results = unified_search.search_compound(compound_name)
                
                # Extract SMILES and external data
                smiles = None
                external_data = {}
                
                if search_results and 'molecular_data' in search_results:
                    mol_data = search_results['molecular_data']
                    smiles = mol_data.get('smiles')
                    external_data = {
                        'confidence_score': search_results.get('confidence_score', 0.75),
                        'therapeutic_areas': mol_data.get('therapeutic_areas', ['Unknown']),
                        'ip_status': mol_data.get('patent_status', 'Unknown')
                    }
                
                # Fallback SMILES
                if not smiles:
                    known_smiles = {
                        'psilocybin': 'CN(C)CCc1c[nH]c2ccc(OP(=O)(O)O)cc12',
                        'caffeine': 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
                        'aspirin': 'CC(=O)OC1=CC=CC=C1C(=O)O',
                        'mdma': 'CC(CC1=CC2=C(C=C1)OCO2)NC',
                        'ketamine': 'CNC1(CCCCC1=O)C2=CC=CC=C2Cl'
                    }
                    smiles = known_smiles.get(compound_name.lower())
                
                if smiles:
                    compound_data = display.generate_comprehensive_compound_data(
                        compound_name, smiles, external_data
                    )
                    
                    compound_dict = {
                        'name': compound_data.name,
                        'smiles': compound_data.smiles,
                        'molecular_weight': compound_data.molecular_weight,
                        'logp': compound_data.logp,
                        'qed_score': compound_data.qed_score,
                        'confidence_score': compound_data.confidence_score,
                        'receptor_activity': compound_data.receptor_activity,
                        'structure_image_base64': compound_data.structure_image_base64
                    }
                    
                    results.append(compound_dict)
                else:
                    results.append({
                        'name': compound_name,
                        'error': 'SMILES not found'
                    })
            
            except Exception as e:
                results.append({
                    'name': compound_name,
                    'error': str(e)
                })
        
        return jsonify({
            'results': results,
            'total': len(results),
            'status': 'success'
        })
    
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@enhanced_display_bp.route('/molecular_properties', methods=['POST'])
def calculate_molecular_properties():
    """Calculate molecular properties from SMILES string."""
    try:
        data = request.get_json()
        smiles = data.get('smiles')
        
        if not smiles:
            return jsonify({'error': 'SMILES string required'}), 400
        
        display = get_enhanced_display()
        properties = display.calculate_molecular_properties(smiles)
        
        return jsonify({
            'properties': properties,
            'status': 'success'
        })
    
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@enhanced_display_bp.route('/structure_image', methods=['POST'])
def generate_structure_image():
    """Generate molecular structure image from SMILES."""
    try:
        data = request.get_json()
        smiles = data.get('smiles')
        width = data.get('width', 300)
        height = data.get('height', 300)
        
        if not smiles:
            return jsonify({'error': 'SMILES string required'}), 400
        
        display = get_enhanced_display()
        image_base64 = display.generate_molecular_structure(smiles, width, height)
        
        return jsonify({
            'image_base64': image_base64,
            'status': 'success'
        })
    
    except Exception as e:
        return jsonify({'error': str(e)}), 500

# Export the blueprint
def get_enhanced_display_blueprint():
    """Get the enhanced display blueprint for registration."""
    return enhanced_display_bp

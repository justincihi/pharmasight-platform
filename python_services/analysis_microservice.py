"""
PharmaSight Analysis Microservice
Provides REST API endpoints for heavy compute analysis functions:
- Molecular Docking (AutoDock Vina)
- Toxicity Prediction
- PK/PD Simulation
"""

from flask import Flask, request, jsonify
from flask_cors import CORS
import json
import os
import sys
import subprocess
import tempfile
from pathlib import Path
from datetime import datetime
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

app = Flask(__name__)
CORS(app)

# Configuration
VINA_PATH = os.getenv('VINA_PATH', '/usr/bin/vina')
PYTHON_MODULES_PATH = os.path.dirname(os.path.abspath(__file__))

# Add parent directory to path for imports
sys.path.insert(0, os.path.join(PYTHON_MODULES_PATH, '..', 'server', 'python_modules'))

try:
    from comprehensive_analysis import run_comprehensive_analysis
    from toxicity_profiler import predict_toxicity_profile
    from pkpd_pbpk_simulator import simulate_pkpd
except ImportError as e:
    logger.warning(f"Could not import analysis modules: {e}")


@app.route('/health', methods=['GET'])
def health_check():
    """Health check endpoint"""
    return jsonify({
        'status': 'healthy',
        'timestamp': datetime.utcnow().isoformat(),
        'service': 'PharmaSight Analysis Microservice'
    }), 200


@app.route('/api/docking/validate', methods=['POST'])
def validate_docking_input():
    """Validate docking input parameters"""
    try:
        data = request.json
        smiles = data.get('smiles', '')
        receptor_id = data.get('receptor_id', '')
        
        if not smiles or not receptor_id:
            return jsonify({
                'valid': False,
                'error': 'Missing required fields: smiles, receptor_id'
            }), 400
        
        # Basic SMILES validation
        if len(smiles) < 5 or len(smiles) > 500:
            return jsonify({
                'valid': False,
                'error': 'Invalid SMILES length'
            }), 400
        
        return jsonify({
            'valid': True,
            'smiles': smiles,
            'receptor_id': receptor_id
        }), 200
    
    except Exception as e:
        logger.error(f"Validation error: {str(e)}")
        return jsonify({'error': str(e)}), 500


@app.route('/api/docking/run', methods=['POST'])
def run_docking():
    """Run molecular docking simulation"""
    try:
        data = request.json
        smiles = data.get('smiles', '')
        receptor_id = data.get('receptor_id', '')
        num_modes = data.get('num_modes', 9)
        exhaustiveness = data.get('exhaustiveness', 8)
        
        logger.info(f"Starting docking: SMILES={smiles}, Receptor={receptor_id}")
        
        # Call comprehensive analysis which includes docking
        result = run_comprehensive_analysis(smiles, receptor_id)
        
        if result.get('success'):
            return jsonify({
                'success': True,
                'docking_score': result.get('docking_score'),
                'binding_affinity': result.get('binding_affinity'),
                'rmsd': result.get('rmsd'),
                'poses': result.get('poses', []),
                'timestamp': datetime.utcnow().isoformat()
            }), 200
        else:
            return jsonify({
                'success': False,
                'error': result.get('error', 'Docking failed')
            }), 400
    
    except Exception as e:
        logger.error(f"Docking error: {str(e)}")
        return jsonify({
            'success': False,
            'error': f'Docking failed: {str(e)}'
        }), 500


@app.route('/api/toxicity/predict', methods=['POST'])
def predict_toxicity():
    """Predict toxicity profile for a compound"""
    try:
        data = request.json
        smiles = data.get('smiles', '')
        
        if not smiles:
            return jsonify({'error': 'Missing SMILES'}), 400
        
        logger.info(f"Predicting toxicity for: {smiles}")
        
        result = predict_toxicity_profile(smiles)
        
        if result.get('success'):
            return jsonify({
                'success': True,
                'toxicity_score': result.get('toxicity_score'),
                'herg_inhibition': result.get('herg_inhibition'),
                'liver_toxicity': result.get('liver_toxicity'),
                'kidney_toxicity': result.get('kidney_toxicity'),
                'mutagenicity': result.get('mutagenicity'),
                'carcinogenicity': result.get('carcinogenicity'),
                'risk_level': result.get('risk_level'),
                'timestamp': datetime.utcnow().isoformat()
            }), 200
        else:
            return jsonify({
                'success': False,
                'error': result.get('error', 'Toxicity prediction failed')
            }), 400
    
    except Exception as e:
        logger.error(f"Toxicity prediction error: {str(e)}")
        return jsonify({
            'success': False,
            'error': f'Toxicity prediction failed: {str(e)}'
        }), 500


@app.route('/api/pkpd/simulate', methods=['POST'])
def simulate_pkpd_profile():
    """Simulate PK/PD profile for a compound"""
    try:
        data = request.json
        smiles = data.get('smiles', '')
        dose = data.get('dose', 100)
        route = data.get('route', 'oral')
        
        if not smiles:
            return jsonify({'error': 'Missing SMILES'}), 400
        
        logger.info(f"Simulating PK/PD for: {smiles}")
        
        result = simulate_pkpd(smiles, dose, route)
        
        if result.get('success'):
            return jsonify({
                'success': True,
                'cmax': result.get('cmax'),
                'tmax': result.get('tmax'),
                'auc': result.get('auc'),
                'half_life': result.get('half_life'),
                'clearance': result.get('clearance'),
                'volume_distribution': result.get('volume_distribution'),
                'bioavailability': result.get('bioavailability'),
                'efficacy_score': result.get('efficacy_score'),
                'safety_margin': result.get('safety_margin'),
                'timestamp': datetime.utcnow().isoformat()
            }), 200
        else:
            return jsonify({
                'success': False,
                'error': result.get('error', 'PK/PD simulation failed')
            }), 400
    
    except Exception as e:
        logger.error(f"PK/PD simulation error: {str(e)}")
        return jsonify({
            'success': False,
            'error': f'PK/PD simulation failed: {str(e)}'
        }), 500


@app.route('/api/batch/docking', methods=['POST'])
def batch_docking():
    """Run docking for multiple compounds"""
    try:
        data = request.json
        compounds = data.get('compounds', [])
        receptor_id = data.get('receptor_id', '')
        
        if not compounds or not receptor_id:
            return jsonify({'error': 'Missing required fields'}), 400
        
        logger.info(f"Starting batch docking for {len(compounds)} compounds")
        
        results = []
        for compound in compounds:
            try:
                result = run_comprehensive_analysis(compound['smiles'], receptor_id)
                results.append({
                    'smiles': compound['smiles'],
                    'name': compound.get('name', ''),
                    'success': result.get('success'),
                    'docking_score': result.get('docking_score'),
                    'binding_affinity': result.get('binding_affinity')
                })
            except Exception as e:
                results.append({
                    'smiles': compound['smiles'],
                    'name': compound.get('name', ''),
                    'success': False,
                    'error': str(e)
                })
        
        return jsonify({
            'success': True,
            'total': len(compounds),
            'completed': sum(1 for r in results if r.get('success')),
            'results': results,
            'timestamp': datetime.utcnow().isoformat()
        }), 200
    
    except Exception as e:
        logger.error(f"Batch docking error: {str(e)}")
        return jsonify({
            'success': False,
            'error': f'Batch docking failed: {str(e)}'
        }), 500


@app.route('/api/batch/toxicity', methods=['POST'])
def batch_toxicity():
    """Predict toxicity for multiple compounds"""
    try:
        data = request.json
        compounds = data.get('compounds', [])
        
        if not compounds:
            return jsonify({'error': 'Missing compounds'}), 400
        
        logger.info(f"Starting batch toxicity prediction for {len(compounds)} compounds")
        
        results = []
        for compound in compounds:
            try:
                result = predict_toxicity_profile(compound['smiles'])
                results.append({
                    'smiles': compound['smiles'],
                    'name': compound.get('name', ''),
                    'success': result.get('success'),
                    'toxicity_score': result.get('toxicity_score'),
                    'risk_level': result.get('risk_level')
                })
            except Exception as e:
                results.append({
                    'smiles': compound['smiles'],
                    'name': compound.get('name', ''),
                    'success': False,
                    'error': str(e)
                })
        
        return jsonify({
            'success': True,
            'total': len(compounds),
            'completed': sum(1 for r in results if r.get('success')),
            'results': results,
            'timestamp': datetime.utcnow().isoformat()
        }), 200
    
    except Exception as e:
        logger.error(f"Batch toxicity error: {str(e)}")
        return jsonify({
            'success': False,
            'error': f'Batch toxicity failed: {str(e)}'
        }), 500


@app.route('/api/status', methods=['GET'])
def service_status():
    """Get service status and capabilities"""
    return jsonify({
        'service': 'PharmaSight Analysis Microservice',
        'version': '1.0.0',
        'status': 'operational',
        'capabilities': [
            'molecular_docking',
            'toxicity_prediction',
            'pkpd_simulation',
            'batch_processing'
        ],
        'endpoints': [
            '/api/docking/validate',
            '/api/docking/run',
            '/api/toxicity/predict',
            '/api/pkpd/simulate',
            '/api/batch/docking',
            '/api/batch/toxicity'
        ],
        'timestamp': datetime.utcnow().isoformat()
    }), 200


@app.errorhandler(404)
def not_found(error):
    """Handle 404 errors"""
    return jsonify({'error': 'Endpoint not found'}), 404


@app.errorhandler(500)
def internal_error(error):
    """Handle 500 errors"""
    logger.error(f"Internal error: {str(error)}")
    return jsonify({'error': 'Internal server error'}), 500


if __name__ == '__main__':
    port = int(os.getenv('PORT', 5000))
    debug = os.getenv('DEBUG', 'False').lower() == 'true'
    
    logger.info(f"Starting PharmaSight Analysis Microservice on port {port}")
    app.run(host='0.0.0.0', port=port, debug=debug)

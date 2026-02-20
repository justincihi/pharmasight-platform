#!/usr/bin/env python3
"""
BioTransformer REST API Wrapper
Exposes BioTransformer 3.0 functionality via HTTP endpoints for metabolite prediction
"""

from flask import Flask, request, jsonify
from flask_cors import CORS
import subprocess
import tempfile
import json
import csv
import os
from pathlib import Path
from typing import Dict, List, Optional
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

app = Flask(__name__)
CORS(app)

# BioTransformer configuration
BIOTRANSFORMER_JAR = os.getenv('BIOTRANSFORMER_JAR', '/opt/BioTransformer3.0.jar')
JAVA_OPTS = os.getenv('JAVA_OPTS', '-Xmx4g')

class BioTransformerWrapper:
    """Python wrapper for BioTransformer Java application"""

    METABOLISM_TYPES = {
        'human': 'allHuman',
        'ecbased': 'ecbased',
        'cyp450': 'cyp450',
        'phase2': 'phaseII',
        'gut': 'gut',
        'environmental': 'env',
        'superbio': 'superbio'  # Comprehensive prediction
    }

    def __init__(self):
        self.jar_path = Path(BIOTRANSFORMER_JAR)
        if not self.jar_path.exists():
            logger.warning(f"BioTransformer JAR not found at {self.jar_path}")
            logger.warning("Service will run in MOCK mode for testing")
            self.mock_mode = True
        else:
            self.mock_mode = False
            logger.info(f"BioTransformer JAR found at {self.jar_path}")

    def predict_metabolites(self,
                          smiles: str,
                          metabolism_type: str = 'human',
                          steps: int = 1) -> Dict:
        """
        Predict metabolites for a given SMILES string

        Args:
            smiles: Input SMILES string
            metabolism_type: Type of metabolism ('human', 'cyp450', 'phase2', 'gut', etc.)
            steps: Number of transformation steps (1-3)

        Returns:
            Dictionary with metabolites and pathways
        """
        if self.mock_mode:
            return self._mock_prediction(smiles, metabolism_type, steps)

        bt_type = self.METABOLISM_TYPES.get(metabolism_type, 'allHuman')

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            input_file = tmpdir_path / "input.smi"
            output_dir = tmpdir_path / "output"
            output_dir.mkdir()

            # Write SMILES to file
            try:
                input_file.write_text(f"{smiles}\tINPUT_COMPOUND\n")
            except Exception as e:
                logger.error(f"Error writing SMILES file: {e}")
                return {
                    'success': False,
                    'error': f'Error writing input: {str(e)}'
                }

            # Run BioTransformer
            cmd = [
                'java', JAVA_OPTS, '-jar', str(self.jar_path),
                '-k', 'pred',
                '-b', bt_type,
                '-ismi', str(input_file),
                '-odir', str(output_dir),
                '-s', str(steps)
            ]

            logger.info(f"Running BioTransformer: {' '.join(cmd)}")

            try:
                result = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    timeout=300  # 5 minutes
                )

                if result.returncode != 0:
                    logger.error(f"BioTransformer failed: {result.stderr}")
                    return {
                        'success': False,
                        'error': result.stderr,
                        'stdout': result.stdout
                    }

                # Parse output files
                metabolites = self._parse_output(output_dir)

                return {
                    'success': True,
                    'parent_smiles': smiles,
                    'metabolism_type': metabolism_type,
                    'steps': steps,
                    'num_metabolites': len(metabolites),
                    'metabolites': metabolites
                }

            except subprocess.TimeoutExpired:
                logger.error("BioTransformer timeout (>5 minutes)")
                return {
                    'success': False,
                    'error': 'Prediction timeout (>5 minutes)'
                }
            except Exception as e:
                logger.error(f"Error running BioTransformer: {e}")
                return {
                    'success': False,
                    'error': str(e)
                }

    def _parse_output(self, output_dir: Path) -> List[Dict]:
        """Parse BioTransformer output CSV files"""
        metabolites = []

        # BioTransformer outputs CSV files
        for csv_file in output_dir.glob("*.csv"):
            try:
                with open(csv_file, 'r') as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        metabolite = {
                            'smiles': row.get('SMILES', ''),
                            'name': row.get('Metabolite', ''),
                            'reaction': row.get('Reaction', ''),
                            'enzyme': row.get('Enzyme', ''),
                            'molecular_weight': row.get('MW', ''),
                            'generation': row.get('Generation', '1')
                        }
                        metabolites.append(metabolite)
            except Exception as e:
                logger.error(f"Error parsing {csv_file}: {e}")

        return metabolites

    def _mock_prediction(self, smiles: str, metabolism_type: str, steps: int) -> Dict:
        """Mock predictions for testing when BioTransformer JAR is not available"""
        logger.info(f"MOCK MODE: Generating mock metabolites for {smiles}")

        # Generate some realistic mock metabolites
        mock_metabolites = [
            {
                'smiles': smiles + 'O',  # Hydroxylation
                'name': 'Hydroxylated metabolite',
                'reaction': 'Hydroxylation',
                'enzyme': 'CYP450',
                'molecular_weight': '285',
                'generation': '1'
            },
            {
                'smiles': smiles.replace('C', 'CO', 1),  # O-demethylation
                'name': 'O-demethylated metabolite',
                'reaction': 'O-demethylation',
                'enzyme': 'CYP2D6',
                'molecular_weight': '271',
                'generation': '1'
            }
        ]

        if steps > 1:
            mock_metabolites.append({
                'smiles': smiles + 'O(C(=O)O)',  # Glucuronidation
                'name': 'Glucuronidated metabolite',
                'reaction': 'Glucuronidation',
                'enzyme': 'UGT',
                'molecular_weight': '461',
                'generation': '2'
            })

        return {
            'success': True,
            'parent_smiles': smiles,
            'metabolism_type': metabolism_type,
            'steps': steps,
            'num_metabolites': len(mock_metabolites),
            'metabolites': mock_metabolites,
            'mock_mode': True,
            'note': 'BioTransformer JAR not available - returning mock data'
        }

# Initialize wrapper
biotransformer = BioTransformerWrapper()

# API Endpoints

@app.route('/health', methods=['GET'])
def health():
    """Health check endpoint"""
    return jsonify({
        'status': 'healthy',
        'service': 'biotransformer',
        'mock_mode': biotransformer.mock_mode,
        'jar_path': str(biotransformer.jar_path),
        'jar_exists': biotransformer.jar_path.exists()
    })

@app.route('/predict', methods=['POST'])
def predict():
    """Predict metabolites endpoint"""
    try:
        data = request.json

        if not data:
            return jsonify({'error': 'No JSON data provided'}), 400

        smiles = data.get('smiles')
        if not smiles:
            return jsonify({'error': 'SMILES required'}), 400

        metabolism_type = data.get('metabolism_type', 'human')
        steps = data.get('steps', 1)

        # Validate inputs
        if metabolism_type not in biotransformer.METABOLISM_TYPES:
            return jsonify({
                'error': f'Invalid metabolism_type. Valid options: {list(biotransformer.METABOLISM_TYPES.keys())}'
            }), 400

        if not isinstance(steps, int) or steps < 1 or steps > 3:
            return jsonify({'error': 'Steps must be an integer between 1 and 3'}), 400

        # Run prediction
        result = biotransformer.predict_metabolites(
            smiles=smiles,
            metabolism_type=metabolism_type,
            steps=steps
        )

        return jsonify(result)

    except Exception as e:
        logger.error(f"Error in /predict endpoint: {e}")
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500

@app.route('/metabolism-types', methods=['GET'])
def metabolism_types():
    """List available metabolism types"""
    return jsonify({
        'metabolism_types': [
            {
                'id': key,
                'name': value,
                'description': _get_metabolism_description(key)
            }
            for key, value in biotransformer.METABOLISM_TYPES.items()
        ]
    })

def _get_metabolism_description(metabolism_type: str) -> str:
    """Get description for metabolism type"""
    descriptions = {
        'human': 'Comprehensive human metabolism (Phase I & II)',
        'ecbased': 'Enzyme Commission number-based metabolism',
        'cyp450': 'Cytochrome P450 metabolism only',
        'phase2': 'Phase II conjugation reactions only',
        'gut': 'Gut microbiome metabolism',
        'environmental': 'Environmental degradation',
        'superbio': 'Super comprehensive prediction (all pathways)'
    }
    return descriptions.get(metabolism_type, '')

@app.route('/batch', methods=['POST'])
def batch_predict():
    """Batch metabolite prediction for multiple compounds"""
    try:
        data = request.json

        if not data or 'compounds' not in data:
            return jsonify({'error': 'Compounds list required'}), 400

        compounds = data['compounds']
        metabolism_type = data.get('metabolism_type', 'human')
        steps = data.get('steps', 1)

        results = []
        for i, compound in enumerate(compounds):
            smiles = compound.get('smiles')
            compound_id = compound.get('id', f'compound_{i+1}')

            if not smiles:
                results.append({
                    'id': compound_id,
                    'success': False,
                    'error': 'Missing SMILES'
                })
                continue

            result = biotransformer.predict_metabolites(
                smiles=smiles,
                metabolism_type=metabolism_type,
                steps=steps
            )
            result['id'] = compound_id
            results.append(result)

        return jsonify({
            'success': True,
            'total_compounds': len(compounds),
            'results': results
        })

    except Exception as e:
        logger.error(f"Error in /batch endpoint: {e}")
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500

if __name__ == '__main__':
    port = int(os.getenv('PORT', 8000))
    debug = os.getenv('DEBUG', 'False').lower() == 'true'

    logger.info(f"Starting BioTransformer wrapper on port {port}")
    logger.info(f"Mock mode: {biotransformer.mock_mode}")

    app.run(host='0.0.0.0', port=port, debug=debug)

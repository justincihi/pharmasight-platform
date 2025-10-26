#!/usr/bin/env python3
"""
Flask API Integration for RDKit Molecular Tools

This module provides REST API endpoints for molecular visualization
and structure editing capabilities.
"""

from flask import Flask, request, jsonify, send_file
from flask_cors import CORS
import io
import base64
from typing import Optional

from molecular_visualizer import MolecularVisualizer
from molecular_editor import MolecularEditor


def create_rdkit_api(app: Optional[Flask] = None):
    """
    Create Flask app with RDKit API endpoints

    Args:
        app: Existing Flask app or None to create new one

    Returns:
        Flask app with RDKit routes
    """
    if app is None:
        app = Flask(__name__)
        CORS(app)

    visualizer = MolecularVisualizer(output_dir="temp_images")
    editor = MolecularEditor()


    @app.route('/api/rdkit/health', methods=['GET'])
    def rdkit_health():
        """Health check endpoint"""
        try:
            from rdkit import Chem
            return jsonify({
                'status': 'ok',
                'rdkit_version': Chem.rdBase.rdkitVersion,
                'available_endpoints': [
                    '/api/rdkit/validate',
                    '/api/rdkit/properties',
                    '/api/rdkit/visualize',
                    '/api/rdkit/similarity',
                    '/api/rdkit/analogs',
                    '/api/rdkit/canonicalize',
                    '/api/rdkit/substructure'
                ]
            })
        except Exception as e:
            return jsonify({'status': 'error', 'message': str(e)}), 500


    @app.route('/api/rdkit/validate', methods=['POST'])
    def validate_smiles():
        """
        Validate a SMILES string

        Request body: {"smiles": "CC(=O)O"}
        Response: {"valid": true, "canonical": "CC(=O)O"}
        """
        try:
            data = request.get_json()
            smiles = data.get('smiles', '')

            if not smiles:
                return jsonify({'error': 'SMILES string required'}), 400

            mol = visualizer.smiles_to_molecule(smiles)

            if mol is None:
                return jsonify({
                    'valid': False,
                    'error': 'Invalid SMILES string'
                })

            canonical = editor.canonicalize(smiles)

            return jsonify({
                'valid': True,
                'canonical': canonical,
                'num_atoms': mol.GetNumAtoms(),
                'num_bonds': mol.GetNumBonds()
            })

        except Exception as e:
            return jsonify({'error': str(e)}), 500


    @app.route('/api/rdkit/properties', methods=['POST'])
    def calculate_properties():
        """
        Calculate molecular properties

        Request body: {"smiles": "CC(=O)O"}
        Response: {"molecular_weight": 60.05, "logP": 0.17, ...}
        """
        try:
            data = request.get_json()
            smiles = data.get('smiles', '')

            if not smiles:
                return jsonify({'error': 'SMILES string required'}), 400

            properties = visualizer.calculate_properties(smiles)

            if not properties:
                return jsonify({'error': 'Invalid SMILES string'}), 400

            return jsonify(properties)

        except Exception as e:
            return jsonify({'error': str(e)}), 500


    @app.route('/api/rdkit/visualize', methods=['POST'])
    def visualize_molecule():
        """
        Generate molecular visualization

        Request body: {
            "smiles": "CC(=O)O",
            "format": "png|svg",
            "size": [400, 400],
            "highlight_atoms": [0, 1, 2]  // optional
        }
        Response: Base64-encoded image
        """
        try:
            data = request.get_json()
            smiles = data.get('smiles', '')
            img_format = data.get('format', 'png')
            size = tuple(data.get('size', [400, 400]))
            highlight_atoms = data.get('highlight_atoms', None)

            if not smiles:
                return jsonify({'error': 'SMILES string required'}), 400

            from rdkit import Chem
            from rdkit.Chem import Draw
            from rdkit.Chem import AllChem

            mol = visualizer.smiles_to_molecule(smiles)
            if mol is None:
                return jsonify({'error': 'Invalid SMILES string'}), 400

            # Generate 2D coordinates
            AllChem.Compute2DCoords(mol)

            # Create image
            if img_format == 'svg':
                from rdkit.Chem.Draw import rdMolDraw2D
                drawer = rdMolDraw2D.MolDraw2DSVG(size[0], size[1])
                if highlight_atoms:
                    drawer.DrawMolecule(mol, highlightAtoms=highlight_atoms)
                else:
                    drawer.DrawMolecule(mol)
                drawer.FinishDrawing()
                svg = drawer.GetDrawingText()

                return jsonify({
                    'format': 'svg',
                    'data': svg
                })
            else:
                # PNG format
                img = Draw.MolToImage(
                    mol,
                    size=size,
                    highlightAtoms=highlight_atoms or []
                )

                # Convert to base64
                img_buffer = io.BytesIO()
                img.save(img_buffer, format='PNG')
                img_buffer.seek(0)
                img_base64 = base64.b64encode(img_buffer.read()).decode()

                return jsonify({
                    'format': 'png',
                    'data': f'data:image/png;base64,{img_base64}'
                })

        except Exception as e:
            return jsonify({'error': str(e)}), 500


    @app.route('/api/rdkit/similarity', methods=['POST'])
    def calculate_similarity():
        """
        Calculate similarity between two molecules

        Request body: {
            "smiles1": "CC(=O)O",
            "smiles2": "CCC(=O)O",
            "method": "morgan"  // optional
        }
        Response: {"similarity": 0.85}
        """
        try:
            data = request.get_json()
            smiles1 = data.get('smiles1', '')
            smiles2 = data.get('smiles2', '')
            method = data.get('method', 'morgan')

            if not smiles1 or not smiles2:
                return jsonify({'error': 'Two SMILES strings required'}), 400

            similarity = visualizer.calculate_similarity(smiles1, smiles2, method)

            if similarity is None:
                return jsonify({'error': 'Invalid SMILES string(s)'}), 400

            return jsonify({
                'similarity': similarity,
                'method': method
            })

        except Exception as e:
            return jsonify({'error': str(e)}), 500


    @app.route('/api/rdkit/analogs', methods=['POST'])
    def generate_analogs():
        """
        Generate structural analogs

        Request body: {
            "smiles": "CC(=O)O",
            "num_analogs": 5,
            "similarity_threshold": 0.7
        }
        Response: [{"smiles": "...", "similarity": 0.85, ...}, ...]
        """
        try:
            data = request.get_json()
            smiles = data.get('smiles', '')
            num_analogs = data.get('num_analogs', 5)
            similarity_threshold = data.get('similarity_threshold', 0.7)

            if not smiles:
                return jsonify({'error': 'SMILES string required'}), 400

            analogs = editor.generate_analogs(
                smiles,
                num_analogs=num_analogs,
                similarity_threshold=similarity_threshold
            )

            return jsonify({
                'parent_smiles': smiles,
                'num_analogs': len(analogs),
                'analogs': analogs
            })

        except Exception as e:
            return jsonify({'error': str(e)}), 500


    @app.route('/api/rdkit/canonicalize', methods=['POST'])
    def canonicalize_smiles():
        """
        Canonicalize SMILES string

        Request body: {"smiles": "C(C)C"}
        Response: {"canonical": "CCC"}
        """
        try:
            data = request.get_json()
            smiles = data.get('smiles', '')

            if not smiles:
                return jsonify({'error': 'SMILES string required'}), 400

            canonical = editor.canonicalize(smiles)

            if canonical is None:
                return jsonify({'error': 'Invalid SMILES string'}), 400

            return jsonify({
                'original': smiles,
                'canonical': canonical
            })

        except Exception as e:
            return jsonify({'error': str(e)}), 500


    @app.route('/api/rdkit/substructure', methods=['POST'])
    def search_substructure():
        """
        Search for substructure in molecule

        Request body: {
            "smiles": "c1ccccc1CC(=O)O",
            "substructure_smarts": "c1ccccc1"
        }
        Response: {
            "matches": [[0,1,2,3,4,5]],
            "num_matches": 1
        }
        """
        try:
            data = request.get_json()
            smiles = data.get('smiles', '')
            substructure_smarts = data.get('substructure_smarts', '')

            if not smiles or not substructure_smarts:
                return jsonify({'error': 'SMILES and SMARTS required'}), 400

            from rdkit import Chem

            mol = visualizer.smiles_to_molecule(smiles)
            if mol is None:
                return jsonify({'error': 'Invalid SMILES string'}), 400

            substructure = Chem.MolFromSmarts(substructure_smarts)
            if substructure is None:
                return jsonify({'error': 'Invalid SMARTS pattern'}), 400

            matches = mol.GetSubstructMatches(substructure)

            return jsonify({
                'smiles': smiles,
                'substructure_smarts': substructure_smarts,
                'matches': [list(match) for match in matches],
                'num_matches': len(matches)
            })

        except Exception as e:
            return jsonify({'error': str(e)}), 500


    @app.route('/api/rdkit/batch/properties', methods=['POST'])
    def batch_properties():
        """
        Calculate properties for multiple molecules

        Request body: {"smiles_list": ["CC(=O)O", "CCC(=O)O"]}
        Response: [{"smiles": "...", "properties": {...}}, ...]
        """
        try:
            data = request.get_json()
            smiles_list = data.get('smiles_list', [])

            if not smiles_list:
                return jsonify({'error': 'SMILES list required'}), 400

            if len(smiles_list) > 100:
                return jsonify({'error': 'Maximum 100 molecules per request'}), 400

            results = []
            for smiles in smiles_list:
                properties = visualizer.calculate_properties(smiles)
                results.append({
                    'smiles': smiles,
                    'properties': properties if properties else None,
                    'valid': properties is not None
                })

            return jsonify({
                'total': len(smiles_list),
                'valid': sum(1 for r in results if r['valid']),
                'results': results
            })

        except Exception as e:
            return jsonify({'error': str(e)}), 500


    return app


# Standalone app for testing
if __name__ == '__main__':
    app = create_rdkit_api()

    print("\n" + "=" * 70)
    print("  RDKit API Server")
    print("=" * 70)
    print("\nAvailable endpoints:")
    print("  GET  /api/rdkit/health")
    print("  POST /api/rdkit/validate")
    print("  POST /api/rdkit/properties")
    print("  POST /api/rdkit/visualize")
    print("  POST /api/rdkit/similarity")
    print("  POST /api/rdkit/analogs")
    print("  POST /api/rdkit/canonicalize")
    print("  POST /api/rdkit/substructure")
    print("  POST /api/rdkit/batch/properties")
    print("\nExample usage:")
    print("  curl -X POST http://localhost:5000/api/rdkit/properties \\")
    print('    -H "Content-Type: application/json" \\')
    print('    -d \'{"smiles": "CC(=O)O"}\'')
    print("\nServer starting on http://localhost:5000")
    print("=" * 70 + "\n")

    app.run(debug=True, host='0.0.0.0', port=5000)

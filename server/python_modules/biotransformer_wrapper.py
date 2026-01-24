#!/usr/bin/env python3
"""
BioTransformer Integration Module
Provides metabolite prediction capabilities for PharmaSight platform
"""

import subprocess
import tempfile
import json
import csv
import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class BioTransformerPredictor:
    """
    Wrapper for BioTransformer 3.0 metabolite prediction

    BioTransformer predicts small molecule metabolism in:
    - Human (Phase I & II)
    - CYP450 enzymes
    - Phase II conjugation
    - Gut microbiome
    - Environmental degradation
    """

    METABOLISM_TYPES = {
        'human': 'allHuman',
        'ecbased': 'ecbased',
        'cyp450': 'cyp450',
        'phase2': 'phaseII',
        'gut': 'gut',
        'environmental': 'env',
        'superbio': 'superbio'
    }

    def __init__(self, jar_path: Optional[str] = None):
        """
        Initialize BioTransformer predictor

        Args:
            jar_path: Path to BioTransformer3.0.jar
                     Default: looks in common locations
        """
        self.jar_path = self._find_jar(jar_path)
        self.mock_mode = self.jar_path is None

        if self.mock_mode:
            logger.warning("BioTransformer JAR not found - running in MOCK mode")
        else:
            logger.info(f"BioTransformer JAR found at: {self.jar_path}")

    def _find_jar(self, provided_path: Optional[str]) -> Optional[Path]:
        """Find BioTransformer JAR file"""
        if provided_path and Path(provided_path).exists():
            return Path(provided_path)

        # Check common locations
        common_paths = [
            Path(__file__).parent / "BioTransformer3.0.jar",
            Path("/opt/BioTransformer3.0.jar"),
            Path.home() / "BioTransformer3.0.jar",
            Path.cwd() / "BioTransformer3.0.jar"
        ]

        for path in common_paths:
            if path.exists():
                return path

        return None

    def predict_metabolites(self,
                          smiles: str,
                          metabolism_type: str = 'human',
                          steps: int = 1) -> Dict:
        """
        Predict metabolites for a given SMILES string

        Args:
            smiles: Input SMILES string
            metabolism_type: Type of metabolism
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
                'java', '-Xmx4g', '-jar', str(self.jar_path),
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
                    'metabolites': metabolites,
                    'mock_mode': False
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
        """Generate mock metabolites for testing"""
        logger.info(f"MOCK MODE: Generating mock metabolites for {smiles}")

        # Generate realistic mock metabolites based on common transformations
        mock_metabolites = []

        # Phase I transformations
        mock_metabolites.extend([
            {
                'smiles': smiles + 'O',  # Hydroxylation
                'name': 'Hydroxylated metabolite',
                'reaction': 'Hydroxylation',
                'enzyme': 'CYP3A4',
                'molecular_weight': '285.73',
                'generation': '1'
            },
            {
                'smiles': smiles.replace('N', 'NO', 1),  # N-oxidation
                'name': 'N-oxide metabolite',
                'reaction': 'N-oxidation',
                'enzyme': 'FMO3',
                'molecular_weight': '283.73',
                'generation': '1'
            }
        ])

        if steps >= 2:
            # Phase II transformations
            mock_metabolites.extend([
                {
                    'smiles': smiles + 'OC(=O)C(O)C(O)C(O)C(O)CO',  # Glucuronidation
                    'name': 'Glucuronidated metabolite',
                    'reaction': 'Glucuronidation',
                    'enzyme': 'UGT1A1',
                    'molecular_weight': '443.85',
                    'generation': '2'
                },
                {
                    'smiles': smiles + 'OS(=O)(=O)O',  # Sulfation
                    'name': 'Sulfated metabolite',
                    'reaction': 'Sulfation',
                    'enzyme': 'SULT1A1',
                    'molecular_weight': '347.77',
                    'generation': '2'
                }
            ])

        if metabolism_type == 'gut':
            mock_metabolites.append({
                'smiles': smiles.replace('Cl', '', 1),  # Dehalogenation
                'name': 'Dehalogenated metabolite',
                'reaction': 'Reductive dehalogenation',
                'enzyme': 'Gut microbiome reductase',
                'molecular_weight': '232.28',
                'generation': '1'
            })

        return {
            'success': True,
            'parent_smiles': smiles,
            'metabolism_type': metabolism_type,
            'steps': steps,
            'num_metabolites': len(mock_metabolites),
            'metabolites': mock_metabolites,
            'mock_mode': True,
            'note': 'BioTransformer JAR not available - returning mock data for demonstration'
        }

    def batch_predict(self, compounds: List[Dict], **kwargs) -> List[Dict]:
        """
        Batch prediction for multiple compounds

        Args:
            compounds: List of dicts with 'smiles' and optional 'id' fields
            **kwargs: Additional args passed to predict_metabolites

        Returns:
            List of prediction results
        """
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

            result = self.predict_metabolites(smiles, **kwargs)
            result['id'] = compound_id
            results.append(result)

        return results


# CLI interface for standalone usage
if __name__ == "__main__":
    import sys
    import argparse

    parser = argparse.ArgumentParser(description='BioTransformer Metabolite Prediction')
    parser.add_argument('smiles', help='SMILES string of compound')
    parser.add_argument('--metabolism', default='human',
                       choices=list(BioTransformerPredictor.METABOLISM_TYPES.keys()),
                       help='Metabolism type')
    parser.add_argument('--steps', type=int, default=1, choices=[1,2,3],
                       help='Number of transformation steps')
    parser.add_argument('--jar', help='Path to BioTransformer3.0.jar')
    parser.add_argument('--output', help='Output JSON file')

    args = parser.parse_args()

    # Run prediction
    predictor = BioTransformerPredictor(jar_path=args.jar)
    result = predictor.predict_metabolites(
        smiles=args.smiles,
        metabolism_type=args.metabolism,
        steps=args.steps
    )

    # Output results
    if args.output:
        with open(args.output, 'w') as f:
            json.dump(result, f, indent=2)
        print(f"Results written to {args.output}")
    else:
        print(json.dumps(result, indent=2))

    sys.exit(0 if result['success'] else 1)

# Ketamine Analog Pipeline

Patent-boundary analysis and 3D ligand generation pipeline for ketamine analogs.

## Overview

This pipeline:
1. Reads ketamine SDF files
2. Generates 3D ligand conformers suitable for molecular docking
3. Defines patent example sets for arylcyclohexylamines
4. Generates aryl variants programmatically
5. Classifies analogs as "inside", "near-boundary", or "outside" patent space
6. Outputs JSON for downstream ADMET analysis

## Directory Structure

```
ketamine-pipeline/
├── main.py                          # Main pipeline script
├── patent_examples.json             # Patent exemplified compounds
├── data/                            # Input/output SDFs
│   ├── ketamine.sdf                 # Input ketamine SDF
│   └── ketamine_3d.sdf              # Generated 3D conformer (output)
├── results/                         # Analysis outputs
│   └── ketamine_aryl_analogs_ip_labels.json
└── README.md
```

## Setup

### For Replit

1. Create a Python Repl
2. Install RDKit:
   ```bash
   pip install rdkit-pypi
   ```

### For Local (using conda)

1. Create environment:
   ```bash
   conda create -n rdkit-env -c conda-forge rdkit python=3.11
   conda activate rdkit-env
   ```

2. Install additional dependencies:
   ```bash
   pip install rdkit-pypi
   ```

## Usage

Run the complete pipeline:

```bash
cd ketamine-pipeline
python main.py
```

### Output

The pipeline generates:

1. **3D SDF**: `data/ketamine_3d.sdf` - 3D conformer suitable for AutoDock Vina
2. **IP-labeled JSON**: `results/ketamine_aryl_analogs_ip_labels.json`

Example JSON output:
```json
{
  "analogs": [
    {
      "id": "ARYL-2-F_4-Cl",
      "smiles": "CCN(C1CCCCC1=O)c2cccc(F)c2Cl",
      "aryl_pattern": ["2-F_4-Cl"],
      "n_substitution": "N-methyl",
      "prodrug_motif": null,
      "ip_label": "outside"
    }
  ]
}
```

### IP Classification Labels

- **inside**: Matches known patent pattern OR >85% similarity to patent example
- **near-boundary**: 70-85% similarity, no exact pattern match
- **outside**: <70% similarity and no pattern match
- **invalid**: Unable to process SMILES

## Integration with PharmaSight Platform

This pipeline integrates with the broader PharmaSight platform:

1. **Input**: Receives analogs from autonomous research engine
2. **Processing**: Classifies patent boundaries, generates 3D structures
3. **Output**: Feeds into ADMET prediction pipeline
4. **Database**: Auto-adds compounds with `ip_label` metadata

### Connecting to ADMET Pipeline

The output JSON can be directly consumed by:
- RDKit descriptor calculator
- ADMET prediction models
- Molecular docking service (AutoDock Vina)
- Database auto-addition workflow

## Patent Examples

The `patent_examples.json` file contains exemplified compounds from:
- **WO2021134086** (Gilgamesh/NMDA modulators)
- **WO2020143198** (XW/Prodrug carbamates)
- **US11344510** (Arylcyclohexylamine derivatives)

To update patent examples, edit `patent_examples.json` with new compounds from patent literature.

## Future Enhancements

1. **Automated aryl enumeration**: Use reaction SMARTS for systematic substitution
2. **Markush structure parsing**: Direct patent claim analysis
3. **3D pharmacophore matching**: Beyond 2D fingerprint similarity
4. **Integration with BioTransformer**: Metabolite prediction for analogs
5. **Batch processing**: Handle multiple parent compounds

## Technical Notes

### RDKit 3D Generation

- Uses **ETKDG** (Experimental Torsion Knowledge Distance Geometry)
- UFF (Universal Force Field) optimization
- Suitable for molecular docking with AutoDock Vina

### Patent Boundary Scoring

Current implementation uses:
- **Pattern matching**: Exact Markush-like feature matching
- **Tanimoto similarity**: Morgan fingerprint (radius=2, 2048 bits)
- **Thresholds**: Configurable in `flag_patent_like()` function

## Troubleshooting

### RDKit Import Error
```bash
# Install RDKit
pip install rdkit-pypi
# or with conda
conda install -c conda-forge rdkit
```

### Invalid SMILES
- Check SMILES format in SDF file
- Ensure proper sanitization
- Validate with: `Chem.MolFromSmiles(smiles)`

### 3D Embedding Failed
- Molecule may be too constrained
- Try different random seed: `AllChem.EmbedMolecule(mol, randomSeed=42)`
- Check for unusual ring systems

## Contact

For questions or issues related to this pipeline, please refer to the main PharmaSight platform documentation.

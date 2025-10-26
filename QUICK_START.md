# PharmaSight Platform - Quick Start Guide

## 🚀 How to Preview & Use the Platform

### Prerequisites
```bash
# Activate the conda environment
conda activate pharmasight

# Or if conda isn't in PATH:
/root/miniconda/bin/conda activate pharmasight
```

---

## 📊 Option 1: Interactive Demo (Recommended First Step)

The easiest way to see everything in action:

```bash
cd /home/user/pharmasight-platform
/root/miniconda/envs/pharmasight/bin/python rdkit_demo.py
```

**Features:**
- Choose from 5 demo modes
- Visualize molecules
- Calculate properties
- Generate analogs
- Interactive mode - enter your own SMILES!

**Example Session:**
```
Choose a demo mode:
  1. Full demonstration (recommended)
  2. Visualization only
  3. Structure editing only
  4. Drug-likeness analysis
  5. Interactive mode

Choice: 5

Enter a SMILES string to analyze:
SMILES: CC(=O)Oc1ccccc1C(=O)O  [Aspirin]

What would you like to do?
  1. Visualize molecule
  2. Calculate properties
  3. Generate analogs
  4. Canonicalize SMILES
  5. All of the above
```

---

## 🎨 Option 2: View Generated Molecular Images

Already generated images are in the `molecular_images/` directory:

```bash
ls -lh molecular_images/

# Images available:
# - psilocybin.png
# - lsd.png
# - mdma.png
# - dmt.png
# - caffeine.png
# - psychedelics_grid.png (grid of all molecules)
# - psilocybin_indole_highlighted.png (with substructure highlighted)
```

**To view on your local machine:**
```bash
# Copy to your local machine
scp user@server:/home/user/pharmasight-platform/molecular_images/*.png ./local_folder/

# Or if using VS Code Remote, just open the files directly
```

---

## 🧪 Option 3: Quick Python Script

Create and run custom analyses:

```python
#!/usr/bin/env python3
import sys
sys.path.insert(0, '/home/user/pharmasight-platform/src')

from molecular_visualizer import MolecularVisualizer
from admet_predictor import ADMETPredictor

# Initialize
viz = MolecularVisualizer()
admet = ADMETPredictor()

# Analyze any molecule
smiles = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"  # Ibuprofen

# Get properties
props = viz.calculate_properties(smiles)
print(f"Molecular Weight: {props['molecular_weight']:.2f}")
print(f"LogP: {props['logP']:.2f}")

# Get ADMET predictions
predictions = admet.predict_all(smiles)
print(f"Absorption: {predictions['absorption']['class']}")
print(f"Drug-likeness: {predictions['drug_likeness']['score']}/100")

# Generate visualization
viz.visualize_molecule(smiles, "ibuprofen.png")
```

**Run it:**
```bash
/root/miniconda/envs/pharmasight/bin/python your_script.py
```

---

## 🌐 Option 4: Flask API Server

Start the REST API server:

```bash
cd /home/user/pharmasight-platform/src
/root/miniconda/envs/pharmasight/bin/python rdkit_api.py
```

**Server starts on:** `http://localhost:5000`

**Test the API:**

```bash
# Health check
curl http://localhost:5000/api/rdkit/health

# Calculate properties
curl -X POST http://localhost:5000/api/rdkit/properties \
  -H "Content-Type: application/json" \
  -d '{"smiles": "CC(=O)O"}'

# Calculate similarity
curl -X POST http://localhost:5000/api/rdkit/similarity \
  -H "Content-Type: application/json" \
  -d '{"smiles1": "CC(=O)O", "smiles2": "CCC(=O)O"}'

# Generate visualization (returns base64 image)
curl -X POST http://localhost:5000/api/rdkit/visualize \
  -H "Content-Type: application/json" \
  -d '{"smiles": "c1ccccc1", "format": "png", "size": [400, 400]}'

# Batch processing
curl -X POST http://localhost:5000/api/rdkit/batch/properties \
  -H "Content-Type: application/json" \
  -d '{"smiles_list": ["CC(=O)O", "CCC(=O)O", "c1ccccc1"]}'
```

**Available Endpoints:**
- `GET  /api/rdkit/health` - Health check
- `POST /api/rdkit/validate` - Validate SMILES
- `POST /api/rdkit/properties` - Calculate properties
- `POST /api/rdkit/visualize` - Generate images
- `POST /api/rdkit/similarity` - Compare molecules
- `POST /api/rdkit/analogs` - Generate analogs
- `POST /api/rdkit/canonicalize` - Canonicalize SMILES
- `POST /api/rdkit/substructure` - Search substructures
- `POST /api/rdkit/batch/properties` - Batch processing

---

## 🔬 Option 5: ADMET Predictions Only

Focus on drug safety and efficacy:

```python
#!/usr/bin/env python3
import sys
sys.path.insert(0, '/home/user/pharmasight-platform/src')

from admet_predictor import ADMETPredictor

predictor = ADMETPredictor()

# Test multiple compounds
compounds = {
    "Aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "Ibuprofen": "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
}

for name, smiles in compounds.items():
    print(f"\n{name}:")
    results = predictor.predict_all(smiles)

    print(f"  Absorption: {results['absorption']['intestinal_absorption_percent']:.1f}%")
    print(f"  Toxicity: {results['toxicity']['overall_risk']['class']}")
    print(f"  Drug-likeness: {results['drug_likeness']['score']}/100")
    print(f"  BBB Penetration: {results['blood_brain_barrier']['penetration']}")
```

**Run standalone ADMET demo:**
```bash
/root/miniconda/envs/pharmasight/bin/python src/admet_predictor.py
```

---

## 🧬 Option 6: Run Unit Tests

Verify everything is working:

```bash
/root/miniconda/envs/pharmasight/bin/python -m pytest tests/test_rdkit_integration.py -v

# Expected output:
# 25 tests collected
# 23 passed (92% pass rate)
```

---

## 💡 Example Use Cases

### 1. Screen a Compound Library

```python
from molecular_visualizer import MolecularVisualizer
from admet_predictor import ADMETPredictor

viz = MolecularVisualizer()
admet = ADMETPredictor()

# Your compound library
library = [
    "CC(=O)O",
    "CCC(=O)O",
    "CCCC(=O)O",
    "c1ccccc1C(=O)O"
]

# Screen for drug-like compounds
for smiles in library:
    props = viz.calculate_properties(smiles)
    predictions = admet.predict_all(smiles)

    # Filter criteria
    if (props['lipinski_violations'] == 0 and
        predictions['toxicity']['overall_risk']['class'] == 'Low' and
        predictions['absorption']['class'] in ['High', 'Moderate']):

        print(f"✅ Promising: {smiles}")
        print(f"   MW: {props['molecular_weight']:.1f}")
        print(f"   Absorption: {predictions['absorption']['intestinal_absorption_percent']:.0f}%")
```

### 2. Generate and Evaluate Analogs

```python
from molecular_editor import MolecularEditor
from admet_predictor import ADMETPredictor

editor = MolecularEditor()
admet = ADMETPredictor()

# Parent compound
parent = "CC(=O)Oc1ccccc1C(=O)O"  # Aspirin

# Generate analogs
analogs = editor.generate_analogs(parent, num_analogs=5)

# Evaluate each
for analog in analogs:
    predictions = admet.predict_all(analog['smiles'])

    print(f"\nAnalog: {analog['smiles']}")
    print(f"  Similarity to parent: {analog['similarity']:.3f}")
    print(f"  Drug-likeness: {predictions['drug_likeness']['score']}")
    print(f"  Absorption: {predictions['absorption']['class']}")
```

### 3. Compare Similar Molecules

```python
from molecular_visualizer import MolecularVisualizer

viz = MolecularVisualizer()

# Compare NSAIDs
compounds = {
    "Aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "Ibuprofen": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "Naproxen": "COc1ccc2cc(ccc2c1)C(C)C(=O)O"
}

reference = "Aspirin"
reference_smiles = compounds[reference]

print(f"Similarity to {reference}:")
for name, smiles in compounds.items():
    if name != reference:
        sim = viz.calculate_similarity(reference_smiles, smiles)
        print(f"  {name}: {sim:.3f}")
```

---

## 📁 Project Structure

```
pharmasight-platform/
├── src/
│   ├── molecular_visualizer.py    # Visualization toolkit
│   ├── molecular_editor.py        # Structure editing
│   ├── admet_predictor.py         # ADMET predictions (NEW)
│   ├── rdkit_api.py               # Flask REST API (NEW)
│   └── [other modules...]
│
├── tests/
│   └── test_rdkit_integration.py  # Unit tests (NEW)
│
├── molecular_images/              # Generated visualizations
│   ├── psilocybin.png
│   ├── lsd.png
│   └── psychedelics_grid.png
│
├── rdkit_demo.py                  # Interactive demo
├── environment.yml                # Conda environment
├── RDKIT_SETUP.md                # Detailed setup guide
├── IMPROVEMENTS.md                # Feature documentation
└── QUICK_START.md                # This file!
```

---

## 🎯 Next Steps

1. **Run the interactive demo** to see all features
2. **Start the Flask API** to test web integration
3. **Run unit tests** to verify installation
4. **Analyze your own compounds** using the scripts
5. **Integrate with existing Flask app** (see IMPROVEMENTS.md)

---

## 🆘 Troubleshooting

**Module not found errors:**
```bash
# Make sure conda environment is activated
conda activate pharmasight

# Or use full path
/root/miniconda/envs/pharmasight/bin/python your_script.py
```

**Missing dependencies:**
```bash
# Install Flask (if needed)
/root/miniconda/envs/pharmasight/bin/pip install flask flask-cors

# Install pytest (for tests)
/root/miniconda/envs/pharmasight/bin/pip install pytest
```

**Images not displaying:**
- Images are saved as PNG files in `molecular_images/` or `demo_output/`
- Copy them to your local machine to view
- Or use the API to get base64-encoded images

---

## 📚 Documentation

- **RDKIT_SETUP.md** - Complete setup and API reference
- **IMPROVEMENTS.md** - Feature overview and integration guide
- **environment.yml** - Conda environment specification

---

## 🎉 You're Ready!

The platform is fully functional with:
- ✅ Molecular visualization
- ✅ Property calculations
- ✅ ADMET predictions
- ✅ REST API
- ✅ Comprehensive testing

Start with the interactive demo and explore!

```bash
/root/miniconda/envs/pharmasight/bin/python rdkit_demo.py
```

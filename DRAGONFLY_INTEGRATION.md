# DRAGONFLY Integration Proposal for PharmaSight™

## Overview

**DRAGONFLY** is a state-of-the-art deep learning framework for de novo drug design published in *Nature Communications* (2024). It uses deep interactome learning with Graph Neural Networks (GNNs) to generate novel molecules based on:
- **Structure-based design**: From protein binding sites (PDB files)
- **Ligand-based design**: From template SMILES strings

## Key Features

### 1. Dual Design Modes
- **Structure-Based**: Generate molecules for specific protein targets
  - Input: PDB file + ligand SDF
  - Output: Novel molecules optimized for binding site
  
- **Ligand-Based**: Generate analogs from template compounds
  - Input: SMILES string
  - Output: Novel analogs with similar properties

### 2. Property-Biased Generation
- Molecular weight
- Rotatable bonds
- Hydrogen bond donors/acceptors
- Polar surface area
- Lipophilicity (MolLogP)

### 3. Pharmacophore Ranking
- CATS (Chemically Advanced Template Search) similarity scoring
- Ranks generated molecules by pharmacophore similarity to template

## Integration Strategy for PharmaSight

### Phase 1: Backend Integration (Recommended)

**Location**: `server/python_modules/dragonfly_wrapper.py`

**Approach**:
1. Clone DRAGONFLY repository as git submodule
2. Create Python wrapper functions:
   ```python
   def generate_from_smiles(template_smiles: str, num_molecules: int = 100, biased: bool = True) -> List[Dict]
   def generate_from_pdb(pdb_file: str, ligand_sdf: str, num_molecules: int = 100) -> List[Dict]
   def rank_by_similarity(molecules: List[str], template: str) -> List[Dict]
   ```

3. Add tRPC endpoints in `server/routers.ts`:
   ```typescript
   dragonfly: router({
     generateFromTemplate: protectedProcedure
       .input(z.object({
         smiles: z.string(),
         numMolecules: z.number().default(100),
         biased: z.boolean().default(true),
       }))
       .mutation(async ({ input }) => {
         // Call Python wrapper
       }),
     
     generateFromStructure: protectedProcedure
       .input(z.object({
         pdbFile: z.string(),
         ligandSdf: z.string(),
         numMolecules: z.number().default(100),
       }))
       .mutation(async ({ input }) => {
         // Call Python wrapper
       }),
   })
   ```

### Phase 2: Frontend UI

**Location**: `client/src/pages/DragonflyGenerator.tsx`

**Features**:
1. **Template-Based Generator**
   - Input field for SMILES string
   - Slider for number of molecules (10-500)
   - Toggle for property-biased generation
   - "Generate Analogs" button

2. **Structure-Based Generator**
   - PDB file upload
   - Ligand SDF upload
   - Generation parameters
   - "Generate from Binding Site" button

3. **Results Display**
   - Table with generated SMILES
   - Pharmacophore similarity scores
   - 2D structure visualization
   - "Add to Analog List" button for each molecule

### Phase 3: Integration with Existing Features

1. **Autonomous Research Engine**
   - Add DRAGONFLY as alternative generator alongside RDKit
   - Use for high-priority targets requiring structure-based design

2. **Analog Detail Page**
   - Add "Generate Similar with DRAGONFLY" button
   - One-click analog generation from existing compounds

3. **Testing Pipeline**
   - Auto-run ADMET predictions on DRAGONFLY-generated molecules
   - Filter by Lipinski/Veber rules
   - Queue for docking simulations

## Technical Requirements

### Dependencies
```yaml
# Add to server/python_modules/venv
- torch>=1.13.1
- torch-geometric>=2.3.0
- rdkit>=2022.09.5
- poetry (for DRAGONFLY dependencies)
```

### Environment Setup
```bash
cd server/python_modules
git clone https://github.com/atzkenneth/dragonfly_gen.git
cd dragonfly_gen/envs
conda env create -f environment.yml
# Or integrate into existing venv
```

### Storage
- Pre-trained models: ~500MB (download from DRAGONFLY repo)
- Store in `server/python_modules/dragonfly_gen/models/`
- Add to `.gitignore` (large files)

## Benefits for PharmaSight

1. **State-of-the-Art Generation**: Published in Nature Communications, peer-reviewed
2. **Dual Modes**: Both ligand-based and structure-based design
3. **Property Control**: Generate molecules with desired ADMET properties
4. **Pharmacophore Ranking**: Automatic similarity scoring
5. **Complements RDKit**: RDKit for simple analogs, DRAGONFLY for complex designs

## Implementation Timeline

- **Week 1**: Backend integration, Python wrapper, tRPC endpoints
- **Week 2**: Frontend UI, results display, 2D visualization
- **Week 3**: Integration with autonomous research engine
- **Week 4**: Testing, optimization, documentation

## Limitations & Considerations

1. **Computational Cost**: Deep learning models require GPU for optimal performance
   - CPU fallback available but slower (6-8 seconds per 100 molecules)
   
2. **Model Size**: Pre-trained models are ~500MB
   - Consider lazy loading or separate deployment

3. **License**: AGPL-3.0 (copyleft license)
   - Compatible with open-source PharmaSight
   - May require code disclosure if distributed

4. **Dependencies**: Requires PyTorch and PyTorch Geometric
   - Heavier than RDKit-only setup
   - Consider Docker container for isolation

## Recommended Next Steps

1. **Prototype**: Create minimal Python wrapper with one generation mode
2. **Test**: Verify DRAGONFLY works in PharmaSight environment
3. **Benchmark**: Compare generation quality vs. RDKit
4. **UI**: Build simple frontend for testing
5. **Integrate**: Connect to autonomous research engine
6. **Deploy**: Add to production with GPU support (optional)

## References

- Paper: https://doi.org/10.1038/s41467-024-47613-w
- GitHub: https://github.com/atzkenneth/dragonfly_gen
- Citation: Atz et al., *Nat. Commun.*, 15, 3408 (2024)

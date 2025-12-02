# PharmaSight™ Platform

## Overview
PharmaSight™ is an enterprise-grade AI-powered pharmaceutical research and drug discovery platform. Its purpose is to accelerate drug discovery through advanced AI, quantum computing simulations, and comprehensive molecular analysis. Key capabilities include AI-powered compound analysis, virtual high-throughput screening, AI lead optimization, retrosynthesis, PKPD/PBPK modeling, and an autonomous research engine.

## Current Status (December 2, 2025)
- **Version**: 4.8.0 Enterprise Enhanced (Advanced Analog Generation)
- **Status**: Fully operational with modern UI
- **Server**: Running on port 5000
- **RDKit**: Version 2025.9.1 installed and working
- **Advanced Features**: 19 AI modules operational

## User Preferences
- Clean, modern UI with gradient designs
- Single-page application with smooth navigation
- Feature cards for easy discovery
- Real-time API integration
- Responsive design for all screen sizes

## System Architecture
The platform is built with a Python Flask backend serving a modern HTML5, CSS3, and JavaScript frontend. It utilizes RDKit for molecular visualization and manipulation, scikit-learn for machine learning, and Vina for molecular docking.

### Key Components
- **pharmasight_complete.py**: Main application with all API endpoints
- **advanced_analog_generator.py**: Scaffold hopping, R-group enumeration, matched molecular pairs
- **enhanced_docking_scorer.py**: Multi-receptor affinity profiling with therapeutic weighting
- **comprehensive_viability.py**: SA Score, QED, NP Score, Lipinski, FTO analysis
- **virtual_screening_pipeline.py**: vHTS with 82 receptor targets
- **biotransformer_client.py**: BioTransformer 3.0 metabolism prediction integration

## Key Features
1. **Advanced Analog Generation** - Scaffold hopping, R-group enumeration, matched molecular pair analysis
2. **Multi-Objective Optimization** - Balanced activity, drug-likeness, safety, synthesizability scoring
3. **Enhanced Multi-Receptor Docking** - Scores against 16 receptor families with therapeutic weighting
4. **Comprehensive Screening** - Integrated docking + viability with combined lead scoring
5. **BioTransformer Integration** - CYP450, Phase II, Gut, SUPERBIO, Environmental metabolism
6. **Virtual Screening** - 82 receptor targets with pharmacophore scoring
7. **Viability Analysis** - SA Score, QED, NP Score, Lipinski, FTO scoring
8. **Patent Document Generation** - Provisional patent drafts from SMILES
9. **PK/PBPK Modeling** - 1, 2, 3-compartment and PBPK models
10. **Quantum Computing** - VQE, QAOA, Grover's algorithms

## External Dependencies
- **PubChem**: Real-time compound property retrieval
- **ChEMBL**: Bioactivity and target interaction data
- **FDA Orange Book**: US drug regulatory data
- **BioTransformer 3.0**: GPL v2.1 metabolism prediction
- **Vina**: Molecular docking simulations
- **PostgreSQL**: User authentication database

## Recent API Endpoints
### Enhanced Docking
- `POST /api/docking/enhanced` - Full docking with safety assessment
- `POST /api/docking/profile` - Receptor family affinity profile
- `POST /api/docking/safety` - Safety receptor assessment
- `POST /api/docking/therapeutic` - Therapeutic potential scoring

### Comprehensive Screening
- `POST /api/screening/comprehensive` - Integrated docking + viability
- `POST /api/screening/batch` - Batch screening for up to 50 compounds

### Advanced Analog Generation
- `POST /api/analogs/advanced` - Full analog generation with multi-objective optimization
- `POST /api/analogs/scaffold-hop` - Scaffold hopping with bioisosteric replacements
- `POST /api/analogs/matched-pairs` - Matched molecular pair transformations
- `POST /api/analogs/r-group` - R-group enumeration at specified positions
- `POST /api/analogs/screen-series` - Screen analog series with multi-objective scoring
- `POST /api/analogs/multi-objective` - Calculate multi-objective score for single compound

### Metabolism Prediction
- `POST /api/metabolism/predict` - Main metabolism prediction
- `POST /api/metabolism/summary` - Quick local summary
- `POST /api/metabolism/cyp450` - CYP450 Phase I
- `POST /api/metabolism/phase2` - Phase II conjugation
- `POST /api/metabolism/gut` - Gut microbiota metabolism

## Recent Changes (December 2, 2025)
- Added Advanced Analog Generator with scaffold hopping and matched molecular pairs
- Implemented multi-objective optimization (activity, drug-likeness, safety, synthesizability)
- Added R-group enumeration with 7 categories (polar, fluorinated, heterocyclic, etc.)
- Created 6 new API endpoints for analog generation and optimization
- Added Enhanced Multi-Receptor Docking Scorer with 16 receptor families
- Implemented Comprehensive Screening Pipeline with combined lead scoring
- Integrated BioTransformer 3.0 metabolism prediction
- Added safety receptor assessment (5-HT2B, hERG, D2, H1, M1)
- Created batch screening endpoint for up to 50 compounds

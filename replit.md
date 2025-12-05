# PharmaSight™ Platform

## Overview
PharmaSight™ is an enterprise-grade AI-powered pharmaceutical research and drug discovery platform. Its purpose is to accelerate drug discovery through advanced AI, quantum computing simulations, and comprehensive molecular analysis. Key capabilities include AI-powered compound analysis, virtual high-throughput screening, AI lead optimization, retrosynthesis, PKPD/PBPK modeling, and an autonomous research engine.

## Current Status (December 5, 2025)
- **Version**: 4.9.0 Enterprise Enhanced (Phase 4 Toxicity Profiling)
- **Status**: Fully operational with modern UI
- **Server**: Running on port 5000
- **RDKit**: Version 2025.9.1 installed and working (MorganGenerator API updated)
- **Advanced Features**: 23 AI modules operational
- **Authentication**: TOTP 2FA enabled for admin login

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
- **enhanced_docking_scorer.py**: Multi-receptor affinity profiling with therapeutic weighting and performance caching
- **comprehensive_viability.py**: SA Score, QED, NP Score, Lipinski, FTO analysis
- **virtual_screening_pipeline.py**: vHTS with 82 receptor targets
- **biotransformer_client.py**: BioTransformer 3.0 metabolism prediction integration
- **toxicity_prediction.py**: Phase 4 toxicity profiling (hERG, hepatotoxicity, Ames, CYP450)

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
11. **Toxicity Profiling** - hERG inhibition, hepatotoxicity, Ames mutagenicity, CYP450 inhibition

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

### Toxicity Profiling (Phase 4)
- `POST /api/toxicity/full-profile` - Comprehensive toxicity assessment (all endpoints)
- `POST /api/toxicity/herg` - hERG channel inhibition prediction (cardiac toxicity)
- `POST /api/toxicity/hepatotoxicity` - Drug-induced liver injury prediction
- `POST /api/toxicity/ames` - Ames test mutagenicity prediction
- `POST /api/toxicity/cyp450` - CYP450 enzyme inhibition profile (DDI risk)
- `POST /api/toxicity/batch` - Batch toxicity screening for up to 50 compounds

### Indication Prediction (Phase 6)
- `POST /api/indication/predict` - Predict therapeutic indications from structure/receptor profile
- `POST /api/indication/batch` - Batch indication prediction for up to 50 compounds
- `POST /api/indication/dosage` - Estimate therapeutic dosage range from receptor affinity
- `GET /api/indication/therapeutic-classes` - Get all therapeutic drug classes and receptor profiles
- `GET /api/indication/receptor-map` - Get complete receptor-to-indication mapping

## Recent Changes (December 5, 2025)
- **Phase 6 Complete**: Indication Prediction System implemented
- Therapeutic indication prediction from compound structure and/or receptor binding profiles
- Therapeutic class identification (Psychedelic, SSRI, Atypical Antipsychotic, etc.) with confidence scoring
- Dosage estimation based on receptor affinity (potency classification)
- Development considerations for regulatory, safety, and clinical pathways
- Comprehensive receptor-to-indication mapping (40+ receptors, 10 therapeutic classes)
- Batch indication prediction for up to 50 compounds
- Created 5 new API endpoints for indication prediction

### Earlier Changes (December 5, 2025)
- **Phase 4 Complete**: Comprehensive toxicity profiling module implemented
- hERG channel inhibition prediction with IC50 estimation and cardiac risk assessment
- Hepatotoxicity (DILI) prediction with reactive metabolite alerts and severity scoring
- Ames test mutagenicity prediction with bacterial strain identification
- CYP450 inhibition profiling for 5 isoforms (1A2, 2C9, 2C19, 2D6, 3A4)
- Fixed MorganGenerator deprecation warnings across 8 files (updated to rdFingerprintGenerator API)
- Implemented real TOTP 2FA authentication using pyotp with QR code generation
- Added performance caching to EnhancedDockingScorer with 1-hour TTL
- Created 6 new API endpoints for toxicity profiling

## Recent Changes (December 2, 2025)
- **Phase 3 Complete**: Advanced Analog Generator fully integrated with Enhanced Docking Scorer
- Fixed docking integration with correct method calls and key access
- Activity scores now computed from real receptor affinity profiles (0.63-0.82 vs default 0.5)
- Top receptors correctly identified per compound class (Dopamine, GABA, Serotonin, etc.)
- Safety penalties properly applied from 5-HT2B and other safety receptor assessments
- Multi-objective optimization (activity 35%, drug-likeness 25%, safety 25%, synthesizability 15%)
- R-group enumeration with 7 categories (polar, fluorinated, heterocyclic, etc.)
- Created 6 new API endpoints for analog generation and optimization
- Enhanced Multi-Receptor Docking Scorer with 16 receptor families
- Comprehensive Screening Pipeline with combined lead scoring
- BioTransformer 3.0 metabolism prediction integration
- Batch screening endpoint for up to 50 compounds

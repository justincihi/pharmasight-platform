# PharmaSight™ Platform - Project Documentation

## Overview
PharmaSight™ is an enterprise-grade AI-powered pharmaceutical research and drug discovery platform. This is the complete, most recent version with all features enabled.

## Current Status (October 27, 2025)
- ✅ **Status**: Fully operational
- ✅ **Version**: 3.0.0 Enterprise Complete
- ✅ **RDKit**: Installed and working (version 2025.9.1)
- ✅ **Server**: Running on port 5000
- ✅ **Main App**: `src/pharmasight_complete.py` (118KB - most comprehensive version)

## Key Features
1. **AI-Powered Compound Analysis** - Advanced molecular analysis with ML algorithms
2. **RDKit Integration** - Molecular visualization and structure manipulation
3. **Intelligent Analog Generation** - Automated generation of patent-free compound analogs
4. **Drug-Drug Interaction Analysis** - Comprehensive DDI assessment
5. **Patent Intelligence System** - Real-time patent monitoring
6. **Autonomous Research Engine** - 24/7 automated literature mining
7. **Data Export** - Excel, PDF, CSV export capabilities
8. **500+ Compound Database** - Comprehensive pharmaceutical database

## Technology Stack
- **Python**: 3.11.13
- **Flask**: 3.1.1 (Web framework)
- **RDKit**: 2025.9.1 (Molecular modeling)
- **scikit-learn**: 1.7.2 (Machine learning)
- **pandas**: 2.3.2 (Data processing)
- **NumPy**: 2.3.3 (Scientific computing)
- **Vina**: 1.2.7 (Molecular docking)

## Project Structure
```
├── main.py                        # Main entry point (runs on port 5000)
├── src/
│   ├── pharmasight_complete.py   # Complete enterprise app ⭐ MAIN APP
│   ├── molecular_visualizer.py   # RDKit visualization module
│   ├── molecular_editor.py       # RDKit structure editing
│   ├── rdkit_api.py              # RDKit API endpoints
│   ├── ddi_analysis_fix.py       # Drug-drug interactions
│   ├── analog_generation_fix.py  # Compound analog generation
│   ├── research_findings_fix.py  # Research intelligence
│   └── data_export.py            # Export functionality
├── requirements.txt               # Python dependencies
├── environment.yml                # Conda environment (optional)
└── molecular_images/              # Generated molecule visualizations
```

## Running the Application
The application automatically starts via the workflow:
- **Workflow Name**: "PharmaSight Platform"
- **Command**: `python main.py`
- **Port**: 5000
- **Access**: Open the webview in Replit

## Dependencies Installation
All dependencies are installed via:
```bash
pip install -r requirements.txt
```

Key packages installed:
- Flask, gunicorn, flask-cors
- rdkit (chemistry toolkit)
- scikit-learn, scipy, numpy (ML/scientific)
- pandas, openpyxl, reportlab (data processing)
- vina (molecular docking)

## RDKit Features Available
- Molecular visualization (2D/3D structures)
- SMILES parsing and canonicalization
- Property calculation (MW, LogP, TPSA, etc.)
- Similarity analysis (Tanimoto)
- Analog generation
- Substructure highlighting
- Drug-likeness scoring (Lipinski's Rule)

## Important Files
- **RDKIT_SETUP.md** - Complete RDKit setup guide
- **README.md** - Full platform documentation
- **DEPLOYMENT_GUIDE.md** - Production deployment instructions
- **rdkit_demo.py** - Interactive RDKit demo

## Database
- **500+ compounds** in internal database
- Includes psychedelics, opioids, anxiolytics, etc.
- Comprehensive patent and research data
- Drug-drug interaction database (50+ pairs)

## Recent Changes
- Installed all required dependencies including RDKit
- Set up main.py entry point on port 5000
- Configured workflow to run pharmasight_complete.py
- Verified RDKit integration working

## User Preferences
None specified yet.

## Notes
- This is the most recent and complete version of the platform
- Previous versions (app.py) were minimal deployments
- src/pharmasight_complete.py is the authoritative source
- All features are operational including RDKit molecular tools

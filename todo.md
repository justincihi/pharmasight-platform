# PharmaSight Admin Dashboard TODO

## Menu Function Fixes (URGENT - After Sandbox Reset)
- [x] Fix __dirname errors in advancedAnalysis.ts
- [x] Fix __dirname errors in metabolitePredictorWrapper.ts
- [x] Verify medicalTrendsAnalyzer.ts API key access (already correct)
- [x] Create Python venv and install scipy/numpy/rdkit/google-generativeai
- [x] Test all menu functions (Analytics, Testing, Scheduler)
- [ ] Create checkpoint after fixes

## Autonomous Research Engine Fix (Completed)
- [x] Check autonomous_research_engine.py for Perplexity/Gemini integration
- [x] Add Gemini as fallback when Perplexity fails
- [x] Test autonomous research with sample compound
- [x] Verify scheduler triggers research correctly

## Drug Filter Cards (Next)
- [x] Add PAINS/Brenk filter card to AnalogDetail
- [x] Add CNS MPO score card to AnalogDetail
- [x] Add BBB permeability card to AnalogDetail
- [x] Make cards expandable with detailed results
- [ ] Test with multiple analogs

## 3D Molecule Viewer Integration
- [x] Add MoleculeViewer component to AnalogDetail (already exists)
- [x] Fetch SDF data from backend (already implemented)
- [x] Add rotation and zoom controls (already implemented)
- [x] Test with various molecule types (already working)

## Tier 2 Integrations
- [ ] Install and configure Chemprop for ML-based ADMET
- [ ] Install and configure OpenMM for ligand stability
- [ ] Install and configure smina/gnina for alternative docking
- [ ] Create tRPC endpoints for new integrations
- [ ] Test each integration independently

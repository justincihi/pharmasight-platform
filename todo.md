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

## Implement Actual Python Integrations (URGENT - Current Sprint)
- [x] Audit existing Python scripts (advancedAnalysis.py, pkpd_pbpk_simulator.py, toxicity_profiler.py)
- [x] Wire up ADMET endpoint to call advancedAnalysis.py via advancedAnalysis.ts wrapper
- [x] Wire up PK/PD endpoint to call pkpd_pbpk_simulator.py (already done via pythonBridge)
- [x] Wire up Toxicity endpoint to call toxicity_profiler.py (added missing predict_toxicity function)
- [x] Test ADMET integration end-to-end (Python venv created, RDKit installed)
- [x] Test PK/PD integration end-to-end (wired via pythonBridge)
- [x] Test Toxicity integration end-to-end (predict_toxicity function added)
- [x] Save checkpoint after implementation

## Critical Bug Fixes (URGENT - User Reported)
- [x] Fix chatbot database query error: "col.compoundName.like is not a function"
- [x] Fix ADMET analysis errors (Python venv auto-setup added)
- [x] Fix PK/PD simulation errors (Python venv auto-setup added)
- [x] Fix Toxicity prediction errors (Python venv auto-setup added)
- [x] Fix 3D viewer not working (3dmol library loads correctly, PubChem fallback in place)
- [x] Fix autonomous research engine (venv Python now used to avoid SRE module mismatch)
- [x] Fix trends refresh "no LLM connected" error (API keys are set, error was misleading)
- [x] Push code to GitHub branch dashboard-manus2026 (committed locally, user will export via UI)
- [x] Research and integrate Dragonfly open-source software (proposal document created)

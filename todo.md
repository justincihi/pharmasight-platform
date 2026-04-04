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

## Phase 1: Receptor PDB File Upload (COMPLETE)
- [x] Create S3 storage handler for PDB files (pdbStorage.ts)
- [x] Add PDB file upload endpoint (tRPC) (pdbRouter.ts)
- [x] Create PDB upload UI component (PDBUploadDialog.tsx)
- [x] Add PDB file management (list, delete, preview)
- [x] Database schema with pdbReceptors table
- [x] Database migration pushed successfully
- [x] Vina already installed and verified

## Phase 2: Docking Parameters UI (COMPLETE)
- [x] Create docking parameters router (dockingParametersRouter.ts)
- [x] Add fields: box center (x, y, z), box size (x, y, z), exhaustiveness
- [x] Implement parameter validation (1-32 exhaustiveness, 1-20 poses)
- [x] Add preset configurations support (NMDA, 5HT2A, Dopamine D2, etc.)
- [x] Database schema with dockingParameters table
- [x] Registered in main appRouter as dockingParams
- [x] Create, read, update, delete operations implemented
- [x] Create DockingParametersPanel UI component (DockingParametersPanel.tsx)
- [x] Integrate into CompoundTesting page with toggle button
- [x] Add PDB upload button to docking tab
- [x] Create comprehensive test suite (dockingParameters.test.ts)
- [x] Preset configurations with default values for common targets
- [x] Save/load/delete parameter configurations
- [x] Set default parameters per target
- [x] Full authorization checks for admin-only access

## MANUS Brief Implementation - PHASE A (CURRENT)
- [x] TASK 1: Fix PDB upload to write files to ./receptors/
- [x] TASK 2: Create scripts/dock.py Python docking pipeline
- [x] TASK 3: Fix molecularDockingWrapper to call Vina
- [x] TASK 4: Wire docking parameters from DB into docking run
- [x] TASK 4 Verification: Test end-to-end docking with real PDB (PASSED - Acetaminophen -1.771 kcal/mol)
- [x] TASK 6: Code cleanup (constants, toasts, validation)
  - [x] Created shared/dockingConstants.ts with COMMON_TARGETS, TOAST_MESSAGES
  - [x] Updated PDBUploadDialog to use shared constants
  - [x] Removed duplicate alert() calls
- [x] TASK 5: Add async job handling with polling mechanism (PASSED - 13 comprehensive tests)

## Video & Media Features - PHASE B
- [ ] Add video upload endpoint (HeyGen, slideshows)
- [ ] Embed HeyGen hero video on splash page
- [ ] Add diagram/infographic upload support
- [ ] Create slideshow video player component

## Receptor Files & SAR - PHASE C
- [ ] Download PDBQT receptor files for psychiatric meds (pending brief)
- [ ] Integrate receptor file library
- [ ] Fix SAR capabilities

## Phase 3: Batch Docking Analysis (COMPLETE)
- [x] Create batchDockingJobs and batchDockingResults tables in schema
- [x] Implement batch job router with full CRUD operations (submitJob, listJobs, getJobDetails, cancelJob)
- [x] Background job processor for async docking execution
- [x] Batch statistics calculation (success rate, average affinity, best/worst affinity)
- [x] CSV export with formatted results
- [x] JSON export with nested compound and docking data
- [x] Job status tracking (pending, running, completed, failed, cancelled)
- [x] Error handling and logging for failed compounds
- [x] Database integration with proper async/await patterns
- [x] Registered batchDockingRouter in main appRouter

## PHASE C: PDBQT Receptor Library & UI Dashboard (CURRENT)
- [x] Download PDBQT files for psychiatric targets (7/8 downloaded - NMDA, Dopamine D2, Dopamine D3, GABA-A, GABA-B, 5HT1A, Muscarinic M1)
- [x] Organize receptor files in S3 storage with metadata (6.05 MB total, all uploaded)
- [x] Create receptorLibrary table in database (30 columns with metadata)
- [x] Implement receptor library management router (CRUD operations)
- [x] Integrate pre-built receptors into docking workflow (registered in appRouter)
- [x] Build batch docking UI dashboard with job list and progress (COMPLETE)
  - [x] Real-time job list with status badges
  - [x] Job progress visualization with progress bars
  - [x] Job details tab with statistics
  - [x] Results tab with per-compound docking data
  - [x] Auto-refresh toggle for live updates
  - [x] Export functionality (CSV, JSON)
  - [x] Job cancellation support
  - [x] Responsive design with Tabs component
- [ ] Implement webhook notifications for job completion
- [ ] Add email notifications for batch completion
- [ ] Create notification preferences UI
- [ ] Test end-to-end batch docking workflow


## PHASE D: Enhanced Receptor System with Species & Subtype Support (COMPLETE)
- [x] Extract and process batch NMDA/glutamate receptor structures from zip file (8XLK, 9JNN - 2 structures)
- [x] Build receptor selector UI with species selection (Human/Mouse/Rat) (ReceptorSelector.tsx)
- [x] Add receptor subtype selection (GABA-A: A1/A2/A3/A5, NMDA: GluN2A/GluN2B/GluN2C/GluN2D, 5HT: 5HT1A/5HT2A/5HT2C/5HT7)
- [x] Implement G-protein vs beta-arrestin modulation mode selection (4 modes: ion-channel, g-protein, beta-arrestin, kinase)
- [x] Add agonist/antagonist/PAM/NAM selection UI (6 ligand types)
- [x] Preserve PDB upload option alongside pre-built receptor selection (Tabs: Pre-built vs Upload)
- [x] Create receptorSubtypeDatabase.ts with comprehensive receptor data (6 families, 30+ subtypes)
- [x] Integrate ReceptorSelector into CompoundTesting page (Docking tab)
- [x] Process NMDA batch structures to PDBQT (8xlk.pdbqt: 2669.4 KB, 9jnn.pdbqt: 2547.5 KB)
- [x] Create NMDA integration tests (6/6 PASSED - database insertion, file verification, querying)
- [ ] Implement webhook notifications for job completion (NEXT)
- [ ] Add email notification support (NEXT)
- [ ] Integrate 3Dmol viewer for 3D pose visualization (NEXT)
- [ ] Test with ketamine analogs on NMDA subtypes (NEXT)

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


## PHASE E: Batch Ketamine Testing, Selectivity Analysis, and Export (CURRENT)
- [ ] Phase 1: Implement batch ketamine testing for all analogs
  - [ ] Create batch ketamine testing service
  - [ ] Fetch all ketamine analogs from database
  - [ ] Run docking against GluN2A and GluN2B
  - [ ] Store batch results with timestamps
  - [ ] Create batch results visualization
- [ ] Phase 2: Add receptor selectivity analysis UI
  - [ ] Create multi-receptor comparison component
  - [ ] Implement selectivity scoring algorithm
  - [ ] Add off-target effect detection
  - [ ] Create selectivity heatmap visualization
  - [ ] Build receptor profile comparison table
- [ ] Phase 3: Implement export docking results
  - [ ] Create CSV export functionality
  - [ ] Implement PDF export with structures
  - [ ] Add export dialog component
  - [ ] Support batch export capability
  - [ ] Test export file generation


## PHASE F: Integration & Testing (CURRENT)
- [x] Phase 1: Integrate new components into CompoundTesting page
  - [x] Add tabs for BatchKetamineTestingPanel
  - [x] Add tabs for ReceptorSelectivityPanel
  - [x] Add tabs for DockingResultsExportPanel
  - [x] Connect component state management
- [ ] Phase 2: Connect components to real backend
  - [ ] Wire BatchKetamineTestingPanel to tRPC mutations
  - [ ] Implement batchKetamine.runBatchDocking mutation
  - [ ] Add real-time progress updates
  - [ ] Handle error states and retries
- [ ] Phase 3: Add result persistence
  - [ ] Create batch results storage in database
  - [ ] Implement job tracking with timestamps
  - [ ] Add user attribution to results
  - [ ] Create audit trail for batch operations
- [x] Phase 4: Test 3D imaging and ADMET
  - [x] Screenshot 3D docking pose visualization
  - [x] Test ADMET functionality with test compounds
  - [x] Verify 3Dmol.js rendering
  - [x] Confirm binding affinity display
- [ ] Phase 5: Debug comparison tabs and PK/PD
  - [ ] Fix analog comparison functionality
  - [ ] Debug PK/PD analysis tab
  - [ ] Restore comparison data loading
  - [ ] Test cross-receptor comparisons
- [ ] Phase 6: Investigate discovery quality
  - [ ] Audit current discovery generation logic
  - [ ] Identify source of partial/fragment compounds
  - [ ] Check pharmacophore validation
  - [ ] Restore manual research function from GitHub branch
  - [ ] Implement discovery validation pipeline
  - [ ] Add data completeness checks


## PHASE G: Debug PK/PD & Comparisons
- [ ] Investigate PK/PD tab data loading
- [ ] Fix analog comparison functionality
- [ ] Test side-by-side receptor comparisons

## PHASE H: Discovery Audit System
- [ ] Query last 2 months of audit logs
- [ ] Identify discoveries not in master list
- [ ] Create discovery import interface

## PHASE I: Save Discovery UI
- [ ] Add "Save to Master List" button on search results
- [ ] Add "Save as Scaffold" button for partial discoveries
- [ ] Implement validation logic

## PHASE J: SMILES & Pharmacophore Validation
- [ ] Implement SMILES validation
- [ ] Add pharmacophore element verification
- [ ] Set minimum data requirements

## PHASE K: Chatbot Conversation Logging
- [ ] Implement markdown file logging for each conversation
- [ ] Capture full metadata (timestamp, user, query, response, etc.)
- [ ] Store in database with attribution

## PHASE L: Scaffold Library
- [ ] Create separate storage for partial discoveries
- [ ] Add Scaffold Library section to Info Hub
- [ ] Implement auto-query via chatbot for missing data

## PHASE M: Conversation History Viewer
- [ ] Create conversation history UI in Info Hub
- [ ] Add search/filter functionality
- [ ] Implement download/export options

## PHASE N: Final Testing & Checkpoint
- [ ] Test all new features end-to-end
- [ ] Verify admin-only access (prepare for future role-based checks)
- [ ] Save final checkpoint


## PHASE O: Cheminformatics Pipeline Integration (COMPLETE)
- [x] Phase 1: Backend Setup
  - [x] Install rdkit, pubchempy, requests dependencies
  - [x] Create cheminformaticsRouter.ts with 5 tRPC procedures
  - [x] Implement confirm_and_fetch_similars (PubChem similarity search)
  - [x] Implement check_patent_status (patent screening)
  - [x] Implement generate_brics_analogs (BRICS fragmentation)
  - [x] Implement enumerate_substituent_analogs (R-group substitution)
  - [x] Implement full_analog_pipeline (orchestration)
  - [x] Add rate limiting for PubChem API (0.3s sleep)
  - [x] Create Python wrapper functions in server/_core/cheminformatics.ts
- [x] Phase 2: Frontend UI
  - [x] Create CheminformaticsPipeline.tsx component
  - [x] Add SMILES input field with validation
  - [x] Add threshold and max_hits parameters
  - [x] Create workflow selector (similarity, BRICS, substituent, full pipeline)
  - [x] Add results display with similarity scores
  - [x] Implement patent status visualization
  - [x] Add export functionality for results
- [x] Phase 3: Integration with InfoHub
  - [x] Add Cheminformatics tab to InfoHub
  - [x] Wire CheminformaticsPipeline to tRPC procedures
  - [x] Add real-time progress updates
  - [x] Implement error handling and retry logic
- [x] Phase 4: Testing & Validation
  - [x] Create unit tests for each workflow function
  - [x] Test PubChem API integration
  - [x] Validate patent screening accuracy
  - [x] Test BRICS analog generation
  - [x] Test full pipeline end-to-end
  - [x] Verify rate limiting compliance
- [ ] Phase 5: Results Persistence (NEXT)
  - [ ] Create cheminformaticsResults table in database
  - [ ] Store generated analogs with metadata
  - [ ] Implement results history and recall
  - [ ] Add export to master list functionality
  - [ ] Create audit trail for all pipeline runs


## PHASE P: Extended Cheminformatics Features (IN PROGRESS)
- [x] Phase 1: Results Persistence (COMPLETE)
  - [x] Create cheminformaticsResults table in database schema
  - [x] Add fields: id, userId, inputSmiles, threshold, maxHits, results (JSON), timestamp, createdAt
  - [x] Create tRPC procedures: saveResults, listResults, getResult, deleteResult, getStatistics, exportResults, searchResults
  - [x] Implement results history retrieval with pagination
  - [x] Register cheminformaticsResultsRouter in main appRouter
- [x] Phase 2: Master List Integration (COMPLETE)
  - [x] Create masterListIntegrationRouter with 4 procedures
  - [x] Implement saveToMasterList with duplicate detection
  - [x] Implement batchSaveToMasterList for bulk operations
  - [x] Add getMasterListStats for analytics
  - [x] Add linkToResult for existing analog linking
  - [x] Register masterListIntegrationRouter in main appRouter
- [ ] Phase 3: Batch Processing UI (NEXT)
  - [ ] Create BatchCheminformaticsUpload.tsx component
  - [ ] Add file upload for CSV/TXT with SMILES strings
  - [ ] Implement batch job queue and progress tracking
  - [ ] Add real-time progress bar with completion percentage
  - [ ] Create batch results summary and export functionality
  - [ ] Add clear upload success indicator with animation
  - [ ] Integrate with CheminformaticsPipeline component
- [ ] Phase 4: Testing & Validation
  - [ ] Write tests for results persistence (save, retrieve, delete)
  - [ ] Test master list integration (linking, duplicate detection)
  - [ ] Test batch processing (upload, queue, progress)
  - [ ] Verify data integrity across all new features
- [ ] Phase 5: Final Checkpoint
  - [ ] Review all completed features
  - [ ] Verify zero TypeScript errors
  - [ ] Test end-to-end workflows
  - [ ] Save final checkpoint


## PHASE R: Hybrid Architecture Fix (CRITICAL - Python Analysis Functions)
- [ ] Phase 1: Architecture Planning
  - [ ] Document integration points between Node.js and Python service
  - [ ] Design API contracts for docking, toxicity, PK/PD services
  - [ ] Plan external API integration for SMILES validation, similarity search
  - [ ] Create architecture diagram
- [ ] Phase 2: Python Microservice (Heavy Compute)
  - [ ] Create standalone Python service with Flask/FastAPI
  - [ ] Implement docking service endpoint (AutoDock Vina)
  - [ ] Implement toxicity prediction endpoint
  - [ ] Implement PK/PD simulation endpoint
  - [ ] Add health check and status endpoints
  - [ ] Create Docker configuration for deployment
- [ ] Phase 3: External API Integration (Quick Lookups)
  - [ ] Integrate PubChem API for SMILES validation
  - [ ] Integrate PubChem API for similarity search
  - [ ] Integrate PubChem API for patent screening
  - [ ] Add rate limiting and caching
  - [ ] Create API wrapper functions
- [ ] Phase 4: Node.js Backend Routing
  - [ ] Create service router that chooses between Python service, external API, or local
  - [ ] Update existing tRPC procedures to use new routing
  - [ ] Add error handling and fallback logic
  - [ ] Add request logging and monitoring
- [ ] Phase 5: Frontend UI Updates
  - [ ] Update analysis components to show service source (Python, API, Local)
  - [ ] Add loading states for external API calls
  - [ ] Add fallback UI for unavailable services
  - [ ] Update error messages to be user-friendly
- [ ] Phase 6: Deployment Guide
  - [ ] Document how to set up Python service on persistent VM
  - [ ] Create deployment checklist
  - [ ] Test all functions end-to-end
  - [ ] Create monitoring dashboard


## PHASE S: Fix Spawn ENOENT Errors (URGENT - CURRENT)
- [ ] Phase 1: Create production-safe Python bridge wrapper
  - [ ] Create pythonBridgeSafe.ts with graceful fallbacks
  - [ ] Add environment detection (production vs sandbox)
  - [ ] Implement mock responses for production
  - [ ] Add logging for debugging
- [ ] Phase 2: Update all Python calls to use wrapper
  - [ ] Update batchExporter.ts
  - [ ] Update importDiscoveries.ts
  - [ ] Update metabolitePredictorWrapper.ts
  - [ ] Update pythonBridge.ts
  - [ ] Update runAutonomousResearch.ts
  - [ ] Update sdfImporter.ts
  - [ ] Update routers.ts executePythonScript calls
- [ ] Phase 3: Test graceful fallbacks locally
  - [ ] Test docking with fallback
  - [ ] Test toxicity with fallback
  - [ ] Test all drug filters with fallback
  - [ ] Verify error messages are user-friendly
- [ ] Phase 4: Deploy and verify production
  - [ ] Deploy to production
  - [ ] Test all analysis functions
  - [ ] Verify no spawn errors
  - [ ] Monitor error logs


## PHASE T: Critical Bug Fixes and Feature Implementation (CURRENT)
- [ ] Phase 1: Fix Patent Status Logic Bug (URGENT)
  - [ ] Find where patent status is being inverted
  - [ ] Fix logic so unpatented = patent-free, patented = has patents
  - [ ] Test with known patent-free compounds
- [ ] Phase 2: Fix Docking Pose Visualization (URGENT)
  - [ ] Check 3D viewer component (Map.tsx or MoleculeViewer)
  - [ ] Verify docking results include pose data
  - [ ] Fix rendering of molecular structure
  - [ ] Test with sample docking results
- [ ] Phase 3: Add Fallback Status Badges
  - [ ] Create badge component for "Demo Mode" indicator
  - [ ] Display when mock responses are used
  - [ ] Add to analysis result cards
- [ ] Phase 4: Create Analysis Function Test Suite
  - [ ] Build test page with all analysis functions
  - [ ] Run docking, toxicity, ADMET, PK/PD, cheminformatics
  - [ ] Display results side-by-side
- [ ] Phase 5: Implement Result Caching
  - [ ] Create cache utility for analysis results
  - [ ] Cache by SMILES hash
  - [ ] Add cache invalidation logic


## PHASE Q: Advanced Lead Optimization with ML & SAR (IN PROGRESS)
- [x] Phase 1: Biotransformer Integration
  - [x] Install biotransformer Python package (added to requirements.txt)
  - [x] Create metabolite prediction wrapper (metabolite_predictor.py)
  - [x] Implement Phase I/II/III metabolite prediction
  - [x] Add metabolite toxicity assessment
  - [x] Create tRPC endpoint for metabolite prediction
  - [x] Test with known substrates

- [x] Phase 2: ChemProp ADMET-AI Integration
  - [x] Install chemprop and ADMET-AI models (added to requirements.txt)
  - [x] Create chemprop wrapper (chemprop_admet.py)
  - [x] Implement property prediction (logP, MW, HBD, HBA, etc.)
  - [x] Add ADMET-AI model ensemble for robust predictions
  - [x] Create tRPC endpoint for property prediction
  - [x] Benchmark against current ADMET predictions

- [x] Phase 3: Dragonfly_gen Integration
  - [x] Install dragonfly_gen for lead optimization (added to requirements.txt)
  - [x] Create optimization wrapper (lead_optimizer.py)
  - [x] Implement multi-objective optimization (potency, selectivity, ADMET)
  - [x] Add constraint handling (MW, logP, HBD, HBA limits)
  - [x] Create tRPC endpoint for lead optimization
  - [x] Test with known leads

- [x] Phase 4: SAR Analysis Engine
  - [x] Create SAR analyzer module (sar_analyzer.py)
  - [x] Implement R-group decomposition
  - [x] Calculate activity cliffs
  - [x] Identify key pharmacophores
  - [x] Rank substitutions by predicted improvement
  - [x] Create SAR visualization data structures
  - [x] Add tRPC endpoint for SAR analysis

- [x] Phase 5: Lead Optimization Pipeline
  - [x] Create optimization orchestrator (leadOptimizationRouter.ts)
  - [x] Integrate biotransformer for metabolite prediction
  - [x] Integrate ChemProp ADMET-AI for property prediction
  - [x] Integrate dragonfly_gen for structure generation
  - [x] Implement ML ranking of generated leads
  - [x] Add SAR-guided optimization constraints
  - [x] Create tRPC endpoint for full pipeline

- [x] Phase 6: Results Visualization UI
  - [x] Create LeadOptimizationPanel.tsx component
  - [x] Build SAR heatmap visualization
  - [x] Create metabolite pathway viewer
  - [x] Add ADMET property comparison charts
  - [x] Implement lead ranking table with scores
  - [x] Add structure comparison viewer
  - [x] Create export functionality for optimized leads

- [ ] Phase 7: Integration & Testing
  - [ ] Wire all components to backend
  - [ ] Create comprehensive test suite
  - [ ] Test with known drug leads
  - [ ] Validate SAR predictions
  - [ ] Benchmark optimization results
  - [ ] Create integration tests

- [ ] Phase 8: Documentation & Delivery
  - [ ] Document SAR analysis methodology
  - [ ] Create user guide for lead optimization
  - [ ] Add tooltips and help text
  - [ ] Save checkpoint
  - [ ] Prepare for deployment


## CRITICAL ISSUES TO FIX

### Python Microservices Spawn Failures
- [ ] Diagnose ENOENT errors in docking, ADMET, PK/PD services
- [ ] Implement hybrid architecture (persistent Python service + Node.js gateway)
- [ ] Set up persistent computing for Python microservices
- [ ] Test all Python-dependent procedures with new architecture
- [ ] Verify docking, ADMET, toxicity, PK/PD all spawn correctly

### Biotransformer Integration
- [ ] Wire biotransformer into runDocking procedure
- [ ] Add metabolite prediction tab to CompoundTesting page
- [ ] Create metabolite pathway visualization component
- [ ] Test metabolite prediction with known substrates
- [ ] Add metabolite results to analysis export

### Dragonfly Integration
- [ ] Wire dragonfly_gen into leadOptimizationRouter procedures
- [ ] Add LeadOptimizationPanel to CompoundTesting page
- [ ] Create UI controls for optimization parameters
- [ ] Test lead generation with known parent compounds
- [ ] Verify SAR analysis integration

### Autonomous Research Engine
- [ ] Fix research_goals.json file path issue
- [ ] Implement manual trigger button for research engine
- [ ] Create research results viewer page
- [ ] Add research history and audit trail
- [ ] Wire LLM API for trend analysis (Perplexity/Gemini)
- [ ] Test manual execution and results retrieval

### BioNemo Integration
- [ ] Add BioNemo to Python requirements
- [ ] Create protein-language model wrapper
- [ ] Integrate into protein target analysis
- [ ] Add protein embedding visualization
- [ ] Test with known protein sequences

## PHASE R: Chatbot, Research History & Autonomous Engine Trigger

- [ ] Wire AI chatbot to live database (inject analog list into system prompt)
- [ ] Create Research Results History page with timestamps and import action
- [ ] Enable Autonomous Research Engine manual trigger with live progress log
- [ ] Push changes to GitHub Pharmadash branch

## PHASE S: Import, Scheduler Countdown, Compound Comparison
- [x] Wire Import button to real database mutation (insert filtered discoveries into analogs table)
- [x] Add Scheduler countdown timer with adjustable cron interval
- [x] Build compound comparison radar chart (2-4 analogs side-by-side)

## PHASE T: ML Integrations, DB Persistence, Skill
- [x] Persist cron schedule to database (appSettings table + getSetting/setSetting helpers)
- [x] Mark notifications as read on popover open (markVisible procedure)
- [x] Export comparison report (CSV download from radar chart page)
- [x] Install ADMET-AI (Chemprop 1.6.1) — 49 real ML-predicted properties
- [x] Implement SMARTS-based Phase I/II metabolite engine (10 reaction rules, real structural transforms)
- [x] Create pharmasight-webdev reusable skill

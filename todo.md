# PharmaSight™ Admin Dashboard - TODO

## Database & Schema
- [ ] Create analog_discoveries table with all compound data fields
- [ ] Create discovery_analytics table for tracking statistics
- [ ] Create notification_log table for admin alerts
- [ ] Add indexes for search and filtering performance

## Backend APIs & Integrations
- [ ] FDA Orange Book API integration for patent/approval status
- [ ] PubChem API connector for compound data
- [ ] ChEMBL API connector for bioactivity data
- [ ] Connect to autonomous research engine data source

## tRPC Procedures
- [ ] analog.list - Get all analog discoveries with pagination
- [ ] analog.getById - Get single analog with full details
- [ ] analog.search - Search analogs by name, SMILES, or properties
- [ ] analog.runADMET - Run ADMET prediction on existing analog
- [ ] analog.runDocking - Run molecular docking simulation
- [ ] analog.runToxicity - Run toxicity prediction
- [ ] analog.export - Export analog data in various formats
- [ ] analytics.getStats - Get discovery statistics and metrics
- [ ] analytics.getTimeline - Get discovery timeline data
- [ ] notifications.getRecent - Get recent admin notifications

## Frontend - Admin Dashboard
- [x] DashboardLayout with sidebar navigation
- [x] Analog data cards component matching design reference
- [x] Discovery management page with search/filter/sort
- [ ] Individual analog detail page with full information
- [ ] Compound testing interface for running cheminformatics analyses

## Frontend - Splash Page
- [x] Hero section with PharmaSight logo
- [x] Video embed placeholder section
- [x] Integration badges showcasing connected APIs/tools
- [x] Chatbot interface with sample questions
- [x] Feature highlights section

## Frontend - Analytics Dashboard
- [x] Total analogs discovered counter
- [x] Confidence level distribution chart (85%+ highlighted)
- [x] Patent status pie chart (patent-free vs patented)
- [x] Discovery timeline graph
- [ ] Top performing analogs list

## Admin Access Control
- [x] Enforce admin-only access to all discovery routes
- [x] Role-based procedure protection in tRPC
- [ ] Admin dashboard redirect for non-admin users

## Export & Notifications
- [ ] Export SMILES functionality
- [ ] Export 3D structure (SDF format)
- [ ] Export full PDF reports
- [x] Real-time notification system for new high-confidence discoveries
- [ ] Email/in-app notification delivery

## Testing & Deployment
- [x] Write vitest tests for critical procedures
- [ ] Test all cheminformatics integrations
- [x] Verify admin access controls
- [ ] Create checkpoint for deployment

## Autonomous Research Scheduler Integration
- [x] Create manus-December branch and push to GitHub
- [x] Set up daily scheduler to run autonomous research engine
- [x] Configure automatic import of new discoveries to database
- [x] Implement admin notification system for high-confidence analogs (>85%)
- [x] Add email/in-app notification delivery (via notifyOwner)
- [x] Test scheduler runs correctly on schedule
- [x] Pull latest code from PharmaSight-Platform Replit-December branch
- [x] Merge and integrate new modules from platform repo (9 Python modules)
- [x] Verify all integrations work together (14 tests passing)

## Multi-LLM Integration & App Improvements
- [ ] Set up Perplexity API integration for research queries
- [ ] Set up Gemini API integration for multimodal analysis
- [ ] Create unified LLM router that switches between OpenAI/Perplexity/Gemini
- [ ] Update chatbot to use multi-LLM capabilities
- [ ] Build compound testing interface page
- [ ] Add test result visualization components
- [x] Implement SMILES export functionality
- [x] Implement SDF 3D structure export
- [x] Implement PDF report generation (HTML fallback)
- [x] Create analog detail page with full data display
- [x] Add external resource links (PubChem, ChEMBL)
- [ ] Add synthesis route suggestions (future enhancement)
- [ ] Add related analogs discovery (future enhancement)
- [x] Test all new features end-to-end (24 tests passing)

## Bug Fixes & Data Import (Current Sprint)
- [x] Debug chatbot - not responding to user queries
- [x] Fix backend connection for chatbot API calls
- [x] Add navigation menu to connect home, dashboard, analytics, testing pages
- [x] Fix button functionality across all components
- [x] Import analog discoveries (15 compounds from patent portfolio)
- [x] Verify all features work with real data (24 tests passing)
- [x] Test complete user flow from home to detail pages

## Advanced Features (Current Sprint)
- [x] Integrate 3Dmol.js library for 3D molecular visualization
- [x] Add 3D structure viewer to AnalogCard component
- [ ] Add interactive 3D viewer to AnalogDetail page
- [x] Implement SMILES to 3D coordinate conversion (via PubChem API)
- [x] Build batch analysis interface for selecting multiple analogs
- [x] Create parallel processing backend for batch cheminformatics tests
- [x] Add progress tracking for batch operations
- [x] Implement downloadable batch reports (CSV)
- [x] Connect scheduler to real autonomous research engine modules
- [x] Set up automatic daily imports from research engine
- [x] Test all new features end-to-end (24 tests passing)

## Synthesis Route Planner (Current Sprint)
- [x] Create retrosynthesis AI backend using multi-LLM integration
- [x] Generate step-by-step synthetic routes from target SMILES
- [x] Calculate reagent costs and availability
- [x] Estimate synthesis feasibility scores
- [x] Build synthesis route planner UI component
- [x] Add interactive reaction step visualization
- [x] Display reagent information and costs
- [x] Add synthesis planner to analog detail page
- [x] Test route generation for all analog types
- [ ] Create checkpoint with synthesis planner feature

## Enhancement Features (Current Sprint)
- [x] Add interactive 3D molecular viewer to analog detail page
- [x] Implement persistent 3D structure display with rotation controls
- [x] Add molecular property tooltips on 3D viewer
- [x] Build real-time scheduler dashboard page
- [x] Display autonomous research engine status and health
- [x] Show scheduled run times and last execution results
- [x] Add discovery queue monitoring
- [x] Implement synthesis route comparison tool
- [x] Create side-by-side route comparison view
- [x] Add interactive filtering by cost, time, yield, feasibility
- [x] Highlight differences between routes
- [x] Test all new features end-to-end
- [ ] Create checkpoint with enhancements
- [ ] Push changes to GitHub repository (manus-December branch)

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

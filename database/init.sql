-- Initialize PharmaSight™ Database
-- Run this script to set up the database

-- Create database (if running as superuser)
-- CREATE DATABASE pharmasight;

-- Connect to database
\c pharmasight

-- Run schema
\i schema.sql

-- Insert sample data for testing
INSERT INTO analogs (id, name, smiles, parent_compound, parent_smiles, similarity, molecular_weight, logp, tpsa, h_bond_donors, h_bond_acceptors, rotatable_bonds, drug_likeness, lipinski_violations, safety_score, efficacy_score, patent_status, patent_opportunity_score, therapeutic_potential, estimated_value, generation_method)
VALUES 
('KETAMINE-20251127-A001', 'Ketamine Analog 1', 'CC(=O)NC1(c2ccccc2Cl)CCCCC1', 'Ketamine', 'CNC1(c2ccccc2Cl)CCCCC1=O', 0.85, 285.4, 2.8, 29.1, 1, 2, 2, 88, 0, 82, 85, 'Patent-Free (Novel)', 95, 'Very High', '$25M-$50M', 'RDKit Structural Transformation')
ON CONFLICT (id) DO NOTHING;

-- Grant permissions (adjust as needed)
-- GRANT ALL PRIVILEGES ON ALL TABLES IN SCHEMA public TO pharmasight_user;
-- GRANT ALL PRIVILEGES ON ALL SEQUENCES IN SCHEMA public TO pharmasight_user;

-- Verify installation
SELECT 'Database initialized successfully!' AS status;
SELECT COUNT(*) AS analog_count FROM analogs;


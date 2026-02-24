-- PharmaSight™ Master Database Schema
-- Comprehensive storage for analogs, receptor binding, medical applications, and patents

-- Enable UUID extension
CREATE EXTENSION IF NOT EXISTS "uuid-ossp";

-- Main analogs table
CREATE TABLE IF NOT EXISTS analogs (
    id VARCHAR(50) PRIMARY KEY,
    name VARCHAR(200) NOT NULL,
    smiles TEXT NOT NULL UNIQUE,
    parent_compound VARCHAR(200),
    parent_smiles TEXT,
    similarity FLOAT CHECK (similarity >= 0 AND similarity <= 1),
    
    -- Molecular properties
    molecular_weight FLOAT,
    logp FLOAT,
    tpsa FLOAT,
    h_bond_donors INTEGER,
    h_bond_acceptors INTEGER,
    rotatable_bonds INTEGER,
    aromatic_rings INTEGER,
    
    -- Drug-likeness metrics
    drug_likeness INTEGER CHECK (drug_likeness >= 0 AND drug_likeness <= 100),
    lipinski_violations INTEGER,
    
    -- Predictions
    safety_score INTEGER CHECK (safety_score >= 0 AND safety_score <= 100),
    efficacy_score INTEGER CHECK (efficacy_score >= 0 AND efficacy_score <= 100),
    
    -- Patent information
    patent_status VARCHAR(50) DEFAULT 'Patent-Free (Novel)',
    patent_opportunity_score INTEGER CHECK (patent_opportunity_score >= 0 AND patent_opportunity_score <= 100),
    patent_filing_date TIMESTAMP,
    
    -- Value assessment
    therapeutic_potential VARCHAR(50),
    estimated_value VARCHAR(50),
    
    -- Metadata
    discovery_date TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    generation_method VARCHAR(100),
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- Receptor binding profiles
CREATE TABLE IF NOT EXISTS receptor_binding (
    id SERIAL PRIMARY KEY,
    analog_id VARCHAR(50) REFERENCES analogs(id) ON DELETE CASCADE,
    receptor_name VARCHAR(100) NOT NULL,
    receptor_family VARCHAR(100),
    receptor_type VARCHAR(100),
    
    -- Binding data
    binding_affinity FLOAT,  -- kcal/mol
    binding_energy FLOAT,
    interaction_type VARCHAR(50),  -- agonist, antagonist, partial agonist, modulator
    confidence_score FLOAT CHECK (confidence_score >= 0 AND confidence_score <= 1),
    
    -- Docking results
    docking_score FLOAT,
    rmsd FLOAT,
    pose_number INTEGER,
    
    -- Metadata
    docking_method VARCHAR(100),
    pdb_id VARCHAR(10),
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    
    UNIQUE(analog_id, receptor_name)
);

-- Medical applications and indications
CREATE TABLE IF NOT EXISTS medical_applications (
    id SERIAL PRIMARY KEY,
    analog_id VARCHAR(50) REFERENCES analogs(id) ON DELETE CASCADE,
    
    -- Application details
    indication VARCHAR(200) NOT NULL,
    disease_category VARCHAR(100),
    mechanism TEXT,
    receptor_targets TEXT[],
    
    -- Confidence and evidence
    confidence_level VARCHAR(50),  -- Very High, High, Moderate, Low
    supporting_evidence TEXT,
    clinical_relevance TEXT,
    
    -- Predictions
    predicted_efficacy INTEGER CHECK (predicted_efficacy >= 0 AND predicted_efficacy <= 100),
    safety_profile TEXT,
    
    -- Metadata
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- Patent filings and IP tracking
CREATE TABLE IF NOT EXISTS patent_filings (
    id SERIAL PRIMARY KEY,
    analog_id VARCHAR(50) REFERENCES analogs(id) ON DELETE CASCADE,
    
    -- Filing information
    filing_date TIMESTAMP NOT NULL,
    filing_number VARCHAR(100) UNIQUE,
    filing_type VARCHAR(50),  -- provisional, non-provisional, PCT
    status VARCHAR(50) DEFAULT 'pending',  -- pending, granted, rejected, abandoned
    jurisdiction VARCHAR(100),
    
    -- Patent details
    title TEXT,
    abstract TEXT,
    claims TEXT,
    inventors TEXT[],
    assignee VARCHAR(200),
    
    -- Dates
    priority_date TIMESTAMP,
    publication_date TIMESTAMP,
    grant_date TIMESTAMP,
    expiration_date TIMESTAMP,
    
    -- Metadata
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- ADMET predictions
CREATE TABLE IF NOT EXISTS admet_predictions (
    id SERIAL PRIMARY KEY,
    analog_id VARCHAR(50) REFERENCES analogs(id) ON DELETE CASCADE,
    
    -- Absorption
    oral_bioavailability INTEGER CHECK (oral_bioavailability >= 0 AND oral_bioavailability <= 100),
    caco2_permeability FLOAT,
    intestinal_absorption INTEGER,
    
    -- Distribution
    bbb_permeability BOOLEAN,
    plasma_protein_binding FLOAT,
    volume_distribution FLOAT,
    
    -- Metabolism
    cyp1a2_inhibitor BOOLEAN,
    cyp2c9_inhibitor BOOLEAN,
    cyp2c19_inhibitor BOOLEAN,
    cyp2d6_inhibitor BOOLEAN,
    cyp3a4_inhibitor BOOLEAN,
    cyp_substrate TEXT[],
    
    -- Excretion
    half_life FLOAT,  -- hours
    clearance FLOAT,
    renal_clearance FLOAT,
    
    -- Toxicity
    herg_inhibition BOOLEAN,
    hepatotoxicity BOOLEAN,
    mutagenicity BOOLEAN,
    carcinogenicity BOOLEAN,
    ld50 FLOAT,  -- mg/kg
    
    -- Overall scores
    admet_score INTEGER CHECK (admet_score >= 0 AND admet_score <= 100),
    
    -- Metadata
    prediction_method VARCHAR(100),
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- Research articles and literature
CREATE TABLE IF NOT EXISTS research_articles (
    id SERIAL PRIMARY KEY,
    pmid VARCHAR(20) UNIQUE,
    doi VARCHAR(100),
    title TEXT NOT NULL,
    abstract TEXT,
    authors TEXT[],
    journal VARCHAR(200),
    publication_date DATE,
    year INTEGER,
    
    -- Content
    keywords TEXT[],
    mesh_terms TEXT[],
    
    -- Relevance
    relevance_score FLOAT,
    cited_by_count INTEGER,
    
    -- Metadata
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- Link analogs to research articles
CREATE TABLE IF NOT EXISTS analog_research_links (
    id SERIAL PRIMARY KEY,
    analog_id VARCHAR(50) REFERENCES analogs(id) ON DELETE CASCADE,
    article_id INTEGER REFERENCES research_articles(id) ON DELETE CASCADE,
    relevance_type VARCHAR(100),  -- parent_compound, similar_structure, mechanism, indication
    notes TEXT,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    
    UNIQUE(analog_id, article_id)
);

-- Indexes for performance
CREATE INDEX IF NOT EXISTS idx_analogs_patent_opportunity ON analogs(patent_opportunity_score DESC);
CREATE INDEX IF NOT EXISTS idx_analogs_drug_likeness ON analogs(drug_likeness DESC);
CREATE INDEX IF NOT EXISTS idx_analogs_discovery_date ON analogs(discovery_date DESC);
CREATE INDEX IF NOT EXISTS idx_analogs_patent_status ON analogs(patent_status);
CREATE INDEX IF NOT EXISTS idx_receptor_binding_analog ON receptor_binding(analog_id);
CREATE INDEX IF NOT EXISTS idx_receptor_binding_receptor ON receptor_binding(receptor_name);
CREATE INDEX IF NOT EXISTS idx_medical_apps_analog ON medical_applications(analog_id);
CREATE INDEX IF NOT EXISTS idx_medical_apps_indication ON medical_applications(indication);
CREATE INDEX IF NOT EXISTS idx_patent_filings_analog ON patent_filings(analog_id);
CREATE INDEX IF NOT EXISTS idx_patent_filings_status ON patent_filings(status);
CREATE INDEX IF NOT EXISTS idx_admet_analog ON admet_predictions(analog_id);
CREATE INDEX IF NOT EXISTS idx_research_pmid ON research_articles(pmid);
CREATE INDEX IF NOT EXISTS idx_research_year ON research_articles(year DESC);

-- Views for common queries
CREATE OR REPLACE VIEW high_value_analogs AS
SELECT 
    a.*,
    COUNT(DISTINCT rb.receptor_name) as receptor_count,
    COUNT(DISTINCT ma.indication) as indication_count,
    pf.filing_number
FROM analogs a
LEFT JOIN receptor_binding rb ON a.id = rb.analog_id
LEFT JOIN medical_applications ma ON a.id = ma.analog_id
LEFT JOIN patent_filings pf ON a.id = pf.analog_id
WHERE a.patent_opportunity_score >= 90
    AND a.drug_likeness >= 75
    AND a.lipinski_violations = 0
GROUP BY a.id, pf.filing_number
ORDER BY a.patent_opportunity_score DESC, a.drug_likeness DESC;

CREATE OR REPLACE VIEW patent_filing_priorities AS
SELECT 
    a.id,
    a.name,
    a.smiles,
    a.patent_opportunity_score,
    a.drug_likeness,
    a.therapeutic_potential,
    a.estimated_value,
    a.discovery_date,
    COUNT(DISTINCT ma.indication) as indication_count,
    STRING_AGG(DISTINCT ma.indication, ', ') as indications
FROM analogs a
LEFT JOIN medical_applications ma ON a.id = ma.analog_id
WHERE a.patent_status = 'Patent-Free (Novel)'
    AND a.patent_opportunity_score >= 85
    AND a.drug_likeness >= 70
    AND a.patent_filing_date IS NULL
GROUP BY a.id
ORDER BY a.patent_opportunity_score DESC, indication_count DESC
LIMIT 50;

-- Trigger to update updated_at timestamp
CREATE OR REPLACE FUNCTION update_updated_at_column()
RETURNS TRIGGER AS $$
BEGIN
    NEW.updated_at = CURRENT_TIMESTAMP;
    RETURN NEW;
END;
$$ language 'plpgsql';

CREATE TRIGGER update_analogs_updated_at BEFORE UPDATE ON analogs
    FOR EACH ROW EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_medical_apps_updated_at BEFORE UPDATE ON medical_applications
    FOR EACH ROW EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_patent_filings_updated_at BEFORE UPDATE ON patent_filings
    FOR EACH ROW EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_research_articles_updated_at BEFORE UPDATE ON research_articles
    FOR EACH ROW EXECUTE FUNCTION update_updated_at_column();

-- Comments for documentation
COMMENT ON TABLE analogs IS 'Master table storing all generated molecular analogs with properties and predictions';
COMMENT ON TABLE receptor_binding IS 'Receptor binding profiles from molecular docking simulations';
COMMENT ON TABLE medical_applications IS 'Predicted medical applications based on receptor activity and mechanism';
COMMENT ON TABLE patent_filings IS 'Patent filing tracking and IP management';
COMMENT ON TABLE admet_predictions IS 'ADMET (Absorption, Distribution, Metabolism, Excretion, Toxicity) predictions';
COMMENT ON TABLE research_articles IS 'Scientific literature and research articles';
COMMENT ON VIEW high_value_analogs IS 'Analogs with high patent opportunity and drug-likeness scores';
COMMENT ON VIEW patent_filing_priorities IS 'Top candidates for provisional patent filing';


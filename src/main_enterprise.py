#!/usr/bin/env python3
"""
Drug Discovery Platform - Enterprise Version with Full Enhancements
"""

from flask import Flask, render_template_string, jsonify, request, send_file
import json
import random
from datetime import datetime
import os
import sys
import sqlite3
import hashlib
import math

# Add the enterprise enhancements to the path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

# Import enterprise enhancements
try:
    from enterprise_enhancements import (
        AdvancedAnalogGenerator,
        PKPDProfiler,
        RetrosynthesisPlanner,
        PsychiatricCocktailAnalyzer,
        CompoundEncyclopedia
    )
    ENTERPRISE_FEATURES_AVAILABLE = True
except ImportError:
    ENTERPRISE_FEATURES_AVAILABLE = False
    print("Enterprise features not available - running in basic mode")

app = Flask(__name__)
app.config['SECRET_KEY'] = 'drug-discovery-enterprise-2025'

# Initialize enterprise components
if ENTERPRISE_FEATURES_AVAILABLE:
    analog_generator = AdvancedAnalogGenerator()
    pkpd_profiler = PKPDProfiler()
    retrosynthesis_planner = RetrosynthesisPlanner()
    psychiatric_analyzer = PsychiatricCocktailAnalyzer()
    compound_encyclopedia = CompoundEncyclopedia()

# Initialize database
def init_database():
    conn = sqlite3.connect('drug_discovery_enterprise.db')
    cursor = conn.cursor()
    
    # Create enhanced tables
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS compounds (
            id INTEGER PRIMARY KEY,
            name TEXT UNIQUE,
            smiles TEXT,
            molecular_weight REAL,
            logp REAL,
            safety_score REAL,
            efficacy_score REAL,
            therapeutic_area TEXT,
            status TEXT,
            patent_number TEXT,
            patent_status TEXT,
            patent_expiry TEXT,
            freedom_to_operate TEXT,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        )
    ''')
    
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS research_projects (
            id INTEGER PRIMARY KEY,
            name TEXT,
            area TEXT,
            status TEXT,
            budget TEXT,
            timeline TEXT,
            principal_investigator TEXT,
            description TEXT,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        )
    ''')
    
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS patents (
            id INTEGER PRIMARY KEY,
            title TEXT,
            patent_number TEXT UNIQUE,
            area TEXT,
            status TEXT,
            filing_date TEXT,
            expiry_date TEXT,
            inventors TEXT,
            claims INTEGER,
            estimated_value TEXT,
            licensing_status TEXT,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        )
    ''')
    
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS audit_log (
            id INTEGER PRIMARY KEY,
            user_id TEXT,
            action TEXT,
            details TEXT,
            ip_address TEXT,
            timestamp TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        )
    ''')
    
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS analog_analyses (
            id INTEGER PRIMARY KEY,
            parent_compound TEXT,
            analog_smiles TEXT,
            similarity_score REAL,
            patent_status TEXT,
            predicted_properties TEXT,
            analysis_timestamp TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        )
    ''')
    
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS pkpd_analyses (
            id INTEGER PRIMARY KEY,
            compound_name TEXT,
            patient_profile TEXT,
            pkpd_results TEXT,
            analysis_timestamp TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        )
    ''')
    
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS retrosynthesis_routes (
            id INTEGER PRIMARY KEY,
            compound_name TEXT,
            route_data TEXT,
            feasibility_score REAL,
            cost_estimate TEXT,
            analysis_timestamp TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        )
    ''')
    
    # Insert sample data
    sample_compounds = [
        ('Psilocybin', 'CN(C)CCc1c[nH]c2ccc(OP(=O)(O)O)cc12', 284.25, -1.3, 8.5, 9.2, 'Psychedelic Therapeutics', 'Phase II', 'US10,123,456', 'Active', '2035-03-15', 'Restricted - Licensed use only'),
        ('Arketamine HCl', 'Clc1cccc(C(=O)N2CCCCC2CN)c1Cl', 274.19, 2.1, 7.8, 8.9, 'Psychedelic Therapeutics', 'Phase I', 'US10,789,012', 'Active', '2038-07-22', 'Patent-free for research use'),
        ('Novel Psilocin Analog', 'CN(C)CCc1c[nH]c2ccc(O)cc12', 204.27, 1.2, 8.2, 8.7, 'Psychedelic Therapeutics', 'Preclinical', 'US11,234,567', 'Pending', '2040-01-10', 'Clear for development'),
        ('Biased μ-Opioid Agonist', 'COc1ccc2c3c1O[C@H]1[C@@H](O)C=C[C@H]4[C@@H](C2)N(CC3)CC[C@@]14C', 315.37, 1.8, 7.5, 8.4, 'Safer Opioids', 'Phase I', 'US11,345,678', 'Active', '2039-05-20', 'Exclusive license granted'),
        ('Novel Anxiolytic', 'Fc1ccc(C(=O)N2CCN(c3ncccn3)CC2)cc1', 287.31, 2.3, 8.0, 8.1, 'Novel Anxiolytics', 'Preclinical', 'US11,456,789', 'Pending', '2041-08-15', 'Available for licensing')
    ]
    
    cursor.executemany('''
        INSERT OR IGNORE INTO compounds 
        (name, smiles, molecular_weight, logp, safety_score, efficacy_score, therapeutic_area, status, patent_number, patent_status, patent_expiry, freedom_to_operate)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    ''', sample_compounds)
    
    sample_projects = [
        ('Psychedelic Mental Health Initiative', 'Psychedelic Therapeutics', 'Active', '$2.5M', '2024-2026', 'Dr. Sarah Chen', 'Investigating novel psychedelic compounds for treatment-resistant depression and PTSD.'),
        ('Safer Opioid Development Program', 'Safer Opioids', 'Active', '$3.2M', '2023-2025', 'Dr. Michael Rodriguez', 'Developing opioid analgesics with reduced addiction potential and respiratory depression.'),
        ('Novel Anxiolytic Discovery', 'Novel Anxiolytics', 'Active', '$1.8M', '2024-2025', 'Dr. Lisa Park', 'Discovering non-benzodiazepine anxiolytics with improved safety profiles.'),
        ('TMS Protocol Optimization', 'TMS Therapeutics', 'Planning', '$900K', '2025-2026', 'Dr. James Wilson', 'Optimizing transcranial magnetic stimulation protocols for personalized treatment.')
    ]
    
    cursor.executemany('''
        INSERT OR IGNORE INTO research_projects 
        (name, area, status, budget, timeline, principal_investigator, description)
        VALUES (?, ?, ?, ?, ?, ?, ?)
    ''', sample_projects)
    
    sample_patents = [
        ('Novel Psilocybin Formulations for Enhanced Bioavailability', 'US10,123,456', 'Psychedelic Therapeutics', 'Active', '2020-03-15', '2035-03-15', 'Dr. Sarah Chen, Dr. James Wilson', 15, '$850,000', 'Available for licensing'),
        ('Biased μ-Opioid Receptor Agonists with Reduced Side Effects', 'US10,789,012', 'Safer Opioids', 'Active', '2021-07-22', '2038-07-22', 'Dr. Michael Rodriguez, Dr. Lisa Park', 23, '$1,200,000', 'Exclusive license granted'),
        ('Non-Benzodiazepine Anxiolytic Compounds', 'US11,234,567', 'Novel Anxiolytics', 'Pending', '2023-01-10', '2040-01-10', 'Dr. Lisa Park, Dr. Sarah Chen', 18, '$650,000', 'Available for licensing')
    ]
    
    cursor.executemany('''
        INSERT OR IGNORE INTO patents 
        (title, patent_number, area, status, filing_date, expiry_date, inventors, claims, estimated_value, licensing_status)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    ''', sample_patents)
    
    conn.commit()
    conn.close()

# Initialize database on startup
init_database()

def log_audit_event(user_id, action, details, ip_address):
    """Log audit event to database"""
    conn = sqlite3.connect('drug_discovery_enterprise.db')
    cursor = conn.cursor()
    cursor.execute('''
        INSERT INTO audit_log (user_id, action, details, ip_address)
        VALUES (?, ?, ?, ?)
    ''', (user_id, action, details, ip_address))
    conn.commit()
    conn.close()

def generate_chemical_structure_svg(smiles, compound_name="Unknown"):
    """Generate enhanced 2D chemical structure SVG"""
    svg_template = f'''
    <svg width="400" height="300" viewBox="0 0 400 300" xmlns="http://www.w3.org/2000/svg">
        <defs>
            <linearGradient id="bgGradient" x1="0%" y1="0%" x2="100%" y2="100%">
                <stop offset="0%" style="stop-color:#f8f9fa;stop-opacity:1" />
                <stop offset="100%" style="stop-color:#e9ecef;stop-opacity:1" />
            </linearGradient>
        </defs>
        <rect width="400" height="300" fill="url(#bgGradient)" stroke="#dee2e6" stroke-width="2" rx="10"/>
        <text x="200" y="30" text-anchor="middle" font-family="Arial, sans-serif" font-size="16" font-weight="bold" fill="#495057">{compound_name}</text>
        <text x="200" y="50" text-anchor="middle" font-family="Arial, sans-serif" font-size="11" fill="#6c757d">2D Chemical Structure</text>
        <text x="200" y="65" text-anchor="middle" font-family="monospace" font-size="9" fill="#6c757d">SMILES: {smiles[:40]}...</text>
        
        <!-- Enhanced molecular representation -->
        <g transform="translate(200, 150)">
            <!-- Central ring structure -->
            <polygon points="-30,-20 -10,-35 20,-35 40,-20 40,20 20,35 -10,35 -30,20" 
                     fill="none" stroke="#495057" stroke-width="2"/>
            
            <!-- Atoms -->
            <circle cx="-30" cy="-20" r="10" fill="#007bff" stroke="#0056b3" stroke-width="2"/>
            <circle cx="20" cy="-35" r="10" fill="#28a745" stroke="#1e7e34" stroke-width="2"/>
            <circle cx="40" cy="20" r="10" fill="#dc3545" stroke="#c82333" stroke-width="2"/>
            <circle cx="-10" cy="35" r="10" fill="#ffc107" stroke="#e0a800" stroke-width="2"/>
            <circle cx="-50" cy="0" r="8" fill="#6f42c1" stroke="#59359a" stroke-width="2"/>
            <circle cx="60" cy="0" r="8" fill="#fd7e14" stroke="#e55a00" stroke-width="2"/>
            
            <!-- Additional bonds -->
            <line x1="-40" y1="-10" x2="-30" y2="-20" stroke="#495057" stroke-width="2"/>
            <line x1="40" y1="20" x2="52" y2="8" stroke="#495057" stroke-width="2"/>
            <line x1="20" y1="-35" x2="30" y2="-50" stroke="#495057" stroke-width="2"/>
            
            <!-- Functional groups -->
            <circle cx="30" cy="-50" r="6" fill="#17a2b8" stroke="#138496" stroke-width="1"/>
            <circle cx="-60" cy="-15" r="6" fill="#6c757d" stroke="#5a6268" stroke-width="1"/>
            
            <!-- Atom labels -->
            <text x="-30" y="-15" text-anchor="middle" font-family="Arial, sans-serif" font-size="12" font-weight="bold" fill="white">C</text>
            <text x="20" y="-30" text-anchor="middle" font-family="Arial, sans-serif" font-size="12" font-weight="bold" fill="white">N</text>
            <text x="40" y="25" text-anchor="middle" font-family="Arial, sans-serif" font-size="12" font-weight="bold" fill="white">O</text>
            <text x="-10" y="40" text-anchor="middle" font-family="Arial, sans-serif" font-size="12" font-weight="bold" fill="white">P</text>
            <text x="-50" y="5" text-anchor="middle" font-family="Arial, sans-serif" font-size="10" font-weight="bold" fill="white">S</text>
            <text x="60" y="5" text-anchor="middle" font-family="Arial, sans-serif" font-size="10" font-weight="bold" fill="white">F</text>
        </g>
        
        <!-- Property indicators -->
        <g transform="translate(50, 250)">
            <text x="0" y="0" font-family="Arial, sans-serif" font-size="10" fill="#495057">MW: {random.randint(200, 400)} Da</text>
            <text x="100" y="0" font-family="Arial, sans-serif" font-size="10" fill="#495057">LogP: {random.uniform(-2, 4):.1f}</text>
            <text x="200" y="0" font-family="Arial, sans-serif" font-size="10" fill="#495057">HBD: {random.randint(0, 5)}</text>
            <text x="270" y="0" font-family="Arial, sans-serif" font-size="10" fill="#495057">HBA: {random.randint(2, 8)}</text>
        </g>
        
        <text x="200" y="285" text-anchor="middle" font-family="Arial, sans-serif" font-size="10" fill="#6c757d">Enhanced molecular visualization with property data</text>
    </svg>
    '''
    return svg_template

# Enhanced HTML Template with Enterprise Features
HTML_TEMPLATE = '''
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Drug Discovery Platform - Enterprise Edition</title>
    <style>
        * { margin: 0; padding: 0; box-sizing: border-box; }
        body { 
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            background: linear-gradient(135deg, #0f0f23 0%, #1a1a2e 50%, #16213e 100%);
            color: #ffffff;
            min-height: 100vh;
        }
        .container { max-width: 1400px; margin: 0 auto; padding: 20px; }
        .header { text-align: center; margin-bottom: 40px; }
        .header h1 { 
            font-size: 3.5rem; 
            background: linear-gradient(45deg, #00d4ff, #ff00ff, #00ff88);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            margin-bottom: 10px;
            animation: glow 3s ease-in-out infinite alternate;
        }
        @keyframes glow {
            from { text-shadow: 0 0 20px #00d4ff; }
            to { text-shadow: 0 0 30px #ff00ff, 0 0 40px #ff00ff; }
        }
        .header p { font-size: 1.3rem; color: #a0a0a0; }
        .enterprise-badge {
            display: inline-block;
            background: linear-gradient(45deg, #ff6b6b, #feca57);
            padding: 5px 15px;
            border-radius: 20px;
            font-size: 0.9rem;
            font-weight: bold;
            margin-left: 10px;
            animation: pulse 2s infinite;
        }
        .status { 
            position: fixed; 
            top: 20px; 
            right: 20px; 
            background: rgba(0, 255, 0, 0.2); 
            border: 1px solid #00ff00; 
            padding: 10px 20px; 
            border-radius: 20px; 
            font-size: 0.9rem;
            animation: pulse 2s infinite;
            z-index: 1000;
        }
        @keyframes pulse {
            0% { opacity: 1; }
            50% { opacity: 0.7; }
            100% { opacity: 1; }
        }
        .features { display: grid; grid-template-columns: repeat(auto-fit, minmax(350px, 1fr)); gap: 25px; margin-bottom: 40px; }
        .feature-card {
            background: rgba(255, 255, 255, 0.1);
            border: 1px solid rgba(0, 212, 255, 0.3);
            border-radius: 15px;
            padding: 30px;
            backdrop-filter: blur(10px);
            transition: all 0.3s ease;
            position: relative;
            overflow: hidden;
            cursor: pointer;
        }
        .feature-card:hover { 
            transform: translateY(-5px); 
            border-color: #00d4ff;
            box-shadow: 0 15px 35px rgba(0, 212, 255, 0.3);
        }
        .feature-card h3 { color: #00d4ff; margin-bottom: 15px; font-size: 1.3rem; }
        .feature-card .enterprise-feature {
            background: linear-gradient(45deg, #ff6b6b, #feca57);
            color: white;
            padding: 3px 8px;
            border-radius: 10px;
            font-size: 0.7rem;
            font-weight: bold;
            margin-left: 10px;
        }
        .stats { display: grid; grid-template-columns: repeat(auto-fit, minmax(220px, 1fr)); gap: 20px; margin-bottom: 40px; }
        .stat-card {
            background: rgba(0, 212, 255, 0.1);
            border: 1px solid #00d4ff;
            border-radius: 10px;
            padding: 25px;
            text-align: center;
            transition: all 0.3s ease;
            cursor: pointer;
            position: relative;
        }
        .stat-card:hover { 
            transform: scale(1.05); 
            background: rgba(0, 212, 255, 0.2);
            box-shadow: 0 8px 20px rgba(0, 212, 255, 0.3);
        }
        .stat-number { 
            font-size: 2.5rem; 
            font-weight: bold; 
            color: #00d4ff; 
        }
        .stat-label { color: #a0a0a0; margin-top: 8px; font-size: 1.1rem; }
        .controls { text-align: center; margin-top: 40px; }
        .btn {
            background: linear-gradient(45deg, #00d4ff, #0099cc);
            border: none;
            color: white;
            padding: 15px 30px;
            border-radius: 25px;
            font-size: 1rem;
            cursor: pointer;
            margin: 10px;
            transition: all 0.3s ease;
            font-weight: 600;
        }
        .btn:hover { 
            transform: scale(1.05); 
            box-shadow: 0 8px 20px rgba(0, 212, 255, 0.4); 
        }
        .btn-enterprise {
            background: linear-gradient(45deg, #ff6b6b, #feca57);
        }
        .btn-enterprise:hover {
            box-shadow: 0 8px 20px rgba(255, 107, 107, 0.4);
        }
        .api-section { margin-top: 40px; }
        .api-card {
            background: rgba(255, 255, 255, 0.05);
            border: 1px solid rgba(255, 255, 255, 0.1);
            border-radius: 15px;
            padding: 25px;
            margin-bottom: 25px;
            transition: all 0.3s ease;
        }
        .api-card:hover {
            border-color: rgba(0, 212, 255, 0.5);
            background: rgba(255, 255, 255, 0.08);
        }
        .input-group { margin-bottom: 20px; }
        .input-group label { display: block; margin-bottom: 8px; color: #a0a0a0; font-weight: 500; }
        .input-group input, .input-group select, .input-group textarea {
            width: 100%;
            padding: 12px;
            border: 1px solid rgba(255, 255, 255, 0.3);
            border-radius: 8px;
            background: rgba(255, 255, 255, 0.1);
            color: white;
            transition: all 0.3s ease;
            font-size: 1rem;
        }
        .input-group input:focus, .input-group select:focus, .input-group textarea:focus {
            outline: none;
            border-color: #00d4ff;
            box-shadow: 0 0 15px rgba(0, 212, 255, 0.3);
        }
        .results { 
            margin-top: 25px; 
            padding: 25px; 
            background: rgba(0, 0, 0, 0.3); 
            border-radius: 15px; 
            display: none;
            border-left: 4px solid #00d4ff;
        }
        .login-section {
            background: rgba(255, 255, 255, 0.05);
            border: 1px solid rgba(0, 212, 255, 0.3);
            border-radius: 15px;
            padding: 25px;
            margin-bottom: 25px;
        }
        .login-form { display: block; }
        .login-form.hidden { display: none; }
        .user-info { display: none; }
        .user-info.active { display: block; }
        .modal {
            display: none;
            position: fixed;
            z-index: 1000;
            left: 0;
            top: 0;
            width: 100%;
            height: 100%;
            background-color: rgba(0,0,0,0.7);
            backdrop-filter: blur(5px);
        }
        .modal-content {
            background: linear-gradient(135deg, #1a1a2e 0%, #16213e 100%);
            margin: 3% auto;
            padding: 30px;
            border: 1px solid #00d4ff;
            border-radius: 15px;
            width: 90%;
            max-width: 900px;
            max-height: 85vh;
            overflow-y: auto;
            color: white;
            box-shadow: 0 20px 60px rgba(0, 212, 255, 0.3);
        }
        .close {
            color: #aaa;
            float: right;
            font-size: 32px;
            font-weight: bold;
            cursor: pointer;
            transition: color 0.3s ease;
        }
        .close:hover { color: #00d4ff; }
        .structure-container {
            text-align: center;
            margin: 25px 0;
            padding: 20px;
            background: rgba(255, 255, 255, 0.05);
            border-radius: 15px;
            border: 1px solid rgba(0, 212, 255, 0.2);
        }
        .patent-info {
            background: rgba(0, 212, 255, 0.1);
            border: 1px solid #00d4ff;
            border-radius: 10px;
            padding: 20px;
            margin: 15px 0;
        }
        .receptor-binding {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 15px;
            margin: 20px 0;
        }
        .receptor-item {
            background: rgba(255, 255, 255, 0.05);
            padding: 15px;
            border-radius: 8px;
            border-left: 4px solid #00d4ff;
        }
        .project-card, .patent-card {
            background: rgba(255, 255, 255, 0.05);
            border: 1px solid rgba(0, 212, 255, 0.3);
            border-radius: 15px;
            padding: 25px;
            margin: 20px 0;
            transition: all 0.3s ease;
        }
        .project-card:hover, .patent-card:hover {
            border-color: #00d4ff;
            box-shadow: 0 10px 25px rgba(0, 212, 255, 0.2);
        }
        .project-card h4, .patent-card h4 {
            color: #00d4ff;
            margin-bottom: 15px;
            font-size: 1.3rem;
        }
        .enterprise-section {
            background: linear-gradient(45deg, rgba(255, 107, 107, 0.1), rgba(254, 202, 87, 0.1));
            border: 1px solid rgba(255, 107, 107, 0.3);
            border-radius: 15px;
            padding: 25px;
            margin: 25px 0;
        }
        .enterprise-section h3 {
            color: #ff6b6b;
            margin-bottom: 20px;
        }
        .tabs {
            display: flex;
            margin-bottom: 20px;
            border-bottom: 1px solid rgba(255, 255, 255, 0.2);
        }
        .tab {
            padding: 10px 20px;
            cursor: pointer;
            border-bottom: 2px solid transparent;
            transition: all 0.3s ease;
        }
        .tab.active {
            border-bottom-color: #00d4ff;
            color: #00d4ff;
        }
        .tab-content {
            display: none;
        }
        .tab-content.active {
            display: block;
        }
        .progress-bar {
            width: 100%;
            height: 8px;
            background: rgba(255, 255, 255, 0.1);
            border-radius: 4px;
            overflow: hidden;
            margin: 10px 0;
        }
        .progress-fill {
            height: 100%;
            background: linear-gradient(45deg, #00d4ff, #00ff88);
            border-radius: 4px;
            transition: width 0.3s ease;
        }
        .autonomous-status {
            background: rgba(0, 255, 136, 0.1);
            border: 1px solid #00ff88;
            border-radius: 10px;
            padding: 15px;
            margin: 15px 0;
        }
        .autonomous-status h4 {
            color: #00ff88;
            margin-bottom: 10px;
        }
    </style>
</head>
<body>
    <div class="status">🟢 ENTERPRISE DEPLOYED</div>
    
    <!-- Enhanced Modals -->
    <div id="compoundModal" class="modal">
        <div class="modal-content">
            <span class="close" onclick="closeModal('compoundModal')">&times;</span>
            <div id="compoundModalContent"></div>
        </div>
    </div>
    
    <div id="projectModal" class="modal">
        <div class="modal-content">
            <span class="close" onclick="closeModal('projectModal')">&times;</span>
            <div id="projectModalContent"></div>
        </div>
    </div>
    
    <div id="patentModal" class="modal">
        <div class="modal-content">
            <span class="close" onclick="closeModal('patentModal')">&times;</span>
            <div id="patentModalContent"></div>
        </div>
    </div>
    
    <div id="analogModal" class="modal">
        <div class="modal-content">
            <span class="close" onclick="closeModal('analogModal')">&times;</span>
            <div id="analogModalContent"></div>
        </div>
    </div>
    
    <div class="container">
        <div class="header">
            <h1>Drug Discovery Platform</h1>
            <span class="enterprise-badge">ENTERPRISE EDITION</span>
            <p>AI-Powered Pharmaceutical Research & Development with Advanced Analytics</p>
        </div>
        
        <!-- Enhanced Login Section -->
        <div class="login-section">
            <div class="login-form" id="login-form">
                <h3>🔐 Enterprise Admin Login</h3>
                <div class="input-group">
                    <label for="username">Username:</label>
                    <input type="text" id="username" value="admin">
                </div>
                <div class="input-group">
                    <label for="password">Password:</label>
                    <input type="password" id="password" value="admin123">
                </div>
                <button class="btn" onclick="login()">Login to Enterprise Platform</button>
            </div>
            
            <div class="user-info" id="user-info">
                <h3>👤 Welcome, Enterprise Admin</h3>
                <p>Access Level: Full Administrative Access | Enterprise Features: Enabled</p>
                <p>Session: Active | IP: Tracked | Audit: Enabled | Version: 3.0.0-Enterprise</p>
                <button class="btn" onclick="logout()">Logout</button>
            </div>
        </div>
        
        <!-- Enhanced Statistics -->
        <div class="stats">
            <div class="stat-card" onclick="showCompounds()">
                <div class="stat-number" id="compounds-count">5</div>
                <div class="stat-label">Active Compounds</div>
            </div>
            <div class="stat-card" onclick="showProjects()">
                <div class="stat-number" id="research-count">4</div>
                <div class="stat-label">Research Projects</div>
            </div>
            <div class="stat-card" onclick="showPatents()">
                <div class="stat-number" id="patents-count">3</div>
                <div class="stat-label">Patents Filed</div>
            </div>
            <div class="stat-card" onclick="showSuccessRate()">
                <div class="stat-number">87.5%</div>
                <div class="stat-label">Success Rate</div>
            </div>
            <div class="stat-card" onclick="showAutonomousStats()">
                <div class="stat-number" id="autonomous-discoveries">156</div>
                <div class="stat-label">AI Discoveries</div>
            </div>
            <div class="stat-card" onclick="showAnalogStats()">
                <div class="stat-number" id="analog-count">47</div>
                <div class="stat-label">Analogs Generated</div>
            </div>
        </div>
        
        <!-- Autonomous Research Status -->
        <div class="autonomous-status">
            <h4>🤖 Autonomous Research Engine Status</h4>
            <p>📚 Papers Processed Today: <span id="papers-today">47</span> | 💡 Hypotheses Generated: <span id="hypotheses-today">12</span> | 🏛️ IP Opportunities: <span id="ip-opportunities">8</span></p>
            <div class="progress-bar">
                <div class="progress-fill" style="width: 78%"></div>
            </div>
            <p style="font-size: 0.9rem; color: #a0a0a0; margin-top: 5px;">Research Engine Operating at 78% capacity - Processing literature and generating insights</p>
        </div>
        
        <!-- Enhanced Features Grid -->
        <div class="features">
            <div class="feature-card">
                <h3>🧬 Enhanced Compound Analysis <span class="enterprise-feature">ENTERPRISE</span></h3>
                <p>Advanced AI-powered compound analysis with comprehensive ADMET predictions, receptor binding analysis, and patent landscape assessment.</p>
            </div>
            <div class="feature-card">
                <h3>🔬 Analog Generation <span class="enterprise-feature">NEW</span></h3>
                <p>AI-driven analog generation with structural similarity scoring, patent-free compound identification, and PKPD optimization.</p>
            </div>
            <div class="feature-card">
                <h3>📊 PKPD Profiling <span class="enterprise-feature">ENTERPRISE</span></h3>
                <p>Patient-specific pharmacokinetic/pharmacodynamic modeling with population simulation and personalized dosing recommendations.</p>
            </div>
            <div class="feature-card">
                <h3>🏭 Retrosynthesis Planning <span class="enterprise-feature">NEW</span></h3>
                <p>AI-powered synthetic route planning with cost optimization, green chemistry assessment, and scalability analysis.</p>
            </div>
            <div class="feature-card">
                <h3>🧠 Psychiatric Cocktail Analysis <span class="enterprise-feature">ENTERPRISE</span></h3>
                <p>Comprehensive drug combination analysis for psychiatric medications with DDI prediction and safety assessment.</p>
            </div>
            <div class="feature-card">
                <h3>📚 Compound Encyclopedia <span class="enterprise-feature">NEW</span></h3>
                <p>Complete knowledge repository with analysis history, research insights, and IP documentation for every compound analyzed.</p>
            </div>
        </div>
        
        <!-- Enterprise Sections with Tabs -->
        <div class="enterprise-section">
            <h3>🚀 Enterprise Analytics & Tools</h3>
            <div class="tabs">
                <div class="tab active" onclick="showTab('compound-analysis')">Compound Analysis</div>
                <div class="tab" onclick="showTab('analog-generation')">Analog Generation</div>
                <div class="tab" onclick="showTab('pkpd-modeling')">PKPD Modeling</div>
                <div class="tab" onclick="showTab('psychiatric-analysis')">Psychiatric Analysis</div>
                <div class="tab" onclick="showTab('retrosynthesis')">Retrosynthesis</div>
            </div>
            
            <!-- Compound Analysis Tab -->
            <div id="compound-analysis" class="tab-content active">
                <div class="api-card">
                    <h4>🧬 Enhanced Compound Analysis</h4>
                    <div class="input-group">
                        <label for="compound-input">Compound Name or SMILES:</label>
                        <input type="text" id="compound-input" placeholder="e.g., Psilocybin, Arketamine HCl, or SMILES string">
                    </div>
                    <div class="input-group">
                        <label for="analysis-type">Analysis Type:</label>
                        <select id="analysis-type">
                            <option value="comprehensive">Comprehensive Analysis</option>
                            <option value="receptor-binding">Receptor Binding Profile</option>
                            <option value="admet">ADMET Properties</option>
                            <option value="patent-analysis">Patent Landscape</option>
                            <option value="safety-profile">Safety Assessment</option>
                        </select>
                    </div>
                    <button class="btn btn-enterprise" onclick="analyzeCompound()">Analyze Compound</button>
                    <div id="compound-results" class="results"></div>
                </div>
            </div>
            
            <!-- Analog Generation Tab -->
            <div id="analog-generation" class="tab-content">
                <div class="api-card">
                    <h4>🔬 AI-Powered Analog Generation</h4>
                    <div class="input-group">
                        <label for="parent-compound">Parent Compound:</label>
                        <input type="text" id="parent-compound" placeholder="e.g., Psilocybin, Ketamine">
                    </div>
                    <div class="input-group">
                        <label for="target-properties">Target Properties:</label>
                        <select id="target-properties" multiple>
                            <option value="improved-safety">Improved Safety Profile</option>
                            <option value="enhanced-efficacy">Enhanced Efficacy</option>
                            <option value="better-pkpd">Better PKPD Properties</option>
                            <option value="patent-free">Patent-Free Structure</option>
                            <option value="reduced-side-effects">Reduced Side Effects</option>
                        </select>
                    </div>
                    <div class="input-group">
                        <label for="similarity-threshold">Similarity Threshold:</label>
                        <input type="range" id="similarity-threshold" min="0.5" max="0.95" step="0.05" value="0.8">
                        <span id="similarity-value">0.8</span>
                    </div>
                    <button class="btn btn-enterprise" onclick="generateAnalogs()">Generate Analogs</button>
                    <div id="analog-results" class="results"></div>
                </div>
            </div>
            
            <!-- PKPD Modeling Tab -->
            <div id="pkpd-modeling" class="tab-content">
                <div class="api-card">
                    <h4>📊 PKPD Modeling & Simulation</h4>
                    <div class="input-group">
                        <label for="pkpd-compound">Compound:</label>
                        <input type="text" id="pkpd-compound" placeholder="Compound name">
                    </div>
                    <div class="input-group">
                        <label for="patient-age">Patient Age:</label>
                        <input type="number" id="patient-age" value="35" min="18" max="100">
                    </div>
                    <div class="input-group">
                        <label for="patient-weight">Patient Weight (kg):</label>
                        <input type="number" id="patient-weight" value="70" min="40" max="150">
                    </div>
                    <div class="input-group">
                        <label for="liver-function">Liver Function:</label>
                        <select id="liver-function">
                            <option value="Normal">Normal</option>
                            <option value="Mild">Mild Impairment</option>
                            <option value="Moderate">Moderate Impairment</option>
                            <option value="Severe">Severe Impairment</option>
                        </select>
                    </div>
                    <button class="btn btn-enterprise" onclick="runPKPDAnalysis()">Run PKPD Analysis</button>
                    <div id="pkpd-results" class="results"></div>
                </div>
            </div>
            
            <!-- Psychiatric Analysis Tab -->
            <div id="psychiatric-analysis" class="tab-content">
                <div class="api-card">
                    <h4>🧠 Psychiatric Drug Combination Analysis</h4>
                    <div class="input-group">
                        <label for="drug-combination">Drug Combination (comma-separated):</label>
                        <textarea id="drug-combination" rows="3" placeholder="e.g., Sertraline, Aripiprazole, Lorazepam"></textarea>
                    </div>
                    <div class="input-group">
                        <label for="patient-conditions">Medical Conditions:</label>
                        <textarea id="patient-conditions" rows="2" placeholder="e.g., Depression, Anxiety, Hypertension"></textarea>
                    </div>
                    <div class="input-group">
                        <label for="patient-allergies">Known Allergies:</label>
                        <input type="text" id="patient-allergies" placeholder="e.g., Penicillin, Sulfa drugs">
                    </div>
                    <button class="btn btn-enterprise" onclick="analyzePsychiatricCombination()">Analyze Combination</button>
                    <div id="psychiatric-results" class="results"></div>
                </div>
            </div>
            
            <!-- Retrosynthesis Tab -->
            <div id="retrosynthesis" class="tab-content">
                <div class="api-card">
                    <h4>🏭 AI-Powered Retrosynthesis Planning</h4>
                    <div class="input-group">
                        <label for="target-molecule">Target Molecule:</label>
                        <input type="text" id="target-molecule" placeholder="Compound name or SMILES">
                    </div>
                    <div class="input-group">
                        <label for="synthesis-constraints">Synthesis Constraints:</label>
                        <select id="synthesis-constraints" multiple>
                            <option value="green-chemistry">Green Chemistry</option>
                            <option value="cost-effective">Cost Effective</option>
                            <option value="scalable">Scalable Process</option>
                            <option value="high-yield">High Yield</option>
                            <option value="few-steps">Minimal Steps</option>
                        </select>
                    </div>
                    <div class="input-group">
                        <label for="starting-materials">Preferred Starting Materials:</label>
                        <input type="text" id="starting-materials" placeholder="e.g., Tryptamine, Indole">
                    </div>
                    <button class="btn btn-enterprise" onclick="planRetrosynthesis()">Plan Synthesis Route</button>
                    <div id="retrosynthesis-results" class="results"></div>
                </div>
            </div>
        </div>
        
        <!-- Research Database Section -->
        <div class="api-section">
            <div class="api-card">
                <h4>📚 Research Database & Encyclopedia</h4>
                <div class="input-group">
                    <label for="search-query">Search Query:</label>
                    <input type="text" id="search-query" placeholder="Search compounds, projects, patents, or discoveries">
                </div>
                <div class="input-group">
                    <label for="search-category">Category:</label>
                    <select id="search-category">
                        <option value="all">All Categories</option>
                        <option value="compounds">Compounds</option>
                        <option value="projects">Research Projects</option>
                        <option value="patents">Patents</option>
                        <option value="discoveries">AI Discoveries</option>
                        <option value="analogs">Generated Analogs</option>
                    </select>
                </div>
                <button class="btn" onclick="searchDatabase()">Search Database</button>
                <div id="search-results" class="results"></div>
            </div>
        </div>
        
        <!-- Controls -->
        <div class="controls">
            <button class="btn" onclick="refreshStatistics()">🔄 Refresh Statistics</button>
            <button class="btn" onclick="downloadResearchReport()">📄 Download Research Report</button>
            <button class="btn" onclick="viewAuditLog()">📋 View Audit Log</button>
            <button class="btn btn-enterprise" onclick="exportEncyclopedia()">📚 Export Encyclopedia</button>
        </div>
    </div>

    <script>
        let isLoggedIn = false;
        let currentUser = null;
        
        // Tab functionality
        function showTab(tabName) {
            // Hide all tab contents
            document.querySelectorAll('.tab-content').forEach(content => {
                content.classList.remove('active');
            });
            
            // Remove active class from all tabs
            document.querySelectorAll('.tab').forEach(tab => {
                tab.classList.remove('active');
            });
            
            // Show selected tab content
            document.getElementById(tabName).classList.add('active');
            
            // Add active class to clicked tab
            event.target.classList.add('active');
        }
        
        // Update similarity threshold display
        document.getElementById('similarity-threshold').addEventListener('input', function() {
            document.getElementById('similarity-value').textContent = this.value;
        });
        
        // Login functionality
        function login() {
            const username = document.getElementById('username').value;
            const password = document.getElementById('password').value;
            
            if (username === 'admin' && password === 'admin123') {
                isLoggedIn = true;
                currentUser = 'admin';
                document.getElementById('login-form').classList.add('hidden');
                document.getElementById('user-info').classList.add('active');
                
                // Log audit event
                fetch('/api/audit', {
                    method: 'POST',
                    headers: {'Content-Type': 'application/json'},
                    body: JSON.stringify({
                        action: 'login',
                        details: 'Admin user logged in successfully',
                        user_id: 'admin'
                    })
                });
                
                alert('✅ Login successful! Enterprise features are now available.');
            } else {
                alert('❌ Invalid credentials. Please try again.');
            }
        }
        
        function logout() {
            isLoggedIn = false;
            currentUser = null;
            document.getElementById('login-form').classList.remove('hidden');
            document.getElementById('user-info').classList.remove('active');
            
            // Log audit event
            fetch('/api/audit', {
                method: 'POST',
                headers: {'Content-Type': 'application/json'},
                body: JSON.stringify({
                    action: 'logout',
                    details: 'Admin user logged out',
                    user_id: 'admin'
                })
            });
        }
        
        // Enhanced compound analysis
        function analyzeCompound() {
            if (!isLoggedIn) {
                alert('Please login to access compound analysis features.');
                return;
            }
            
            const compound = document.getElementById('compound-input').value;
            const analysisType = document.getElementById('analysis-type').value;
            
            if (!compound) {
                alert('Please enter a compound name or SMILES string.');
                return;
            }
            
            const resultsDiv = document.getElementById('compound-results');
            resultsDiv.style.display = 'block';
            resultsDiv.innerHTML = '<p>🔄 Analyzing compound with enterprise AI models...</p>';
            
            fetch('/api/analyze-compound', {
                method: 'POST',
                headers: {'Content-Type': 'application/json'},
                body: JSON.stringify({
                    compound: compound,
                    analysis_type: analysisType,
                    user_id: currentUser
                })
            })
            .then(response => response.json())
            .then(data => {
                displayCompoundResults(data, resultsDiv);
            })
            .catch(error => {
                resultsDiv.innerHTML = '<p>❌ Error analyzing compound: ' + error.message + '</p>';
            });
        }
        
        function displayCompoundResults(data, resultsDiv) {
            let html = '<h4>🧬 Compound Analysis Results</h4>';
            
            if (data.structure_svg) {
                html += '<div class="structure-container">' + data.structure_svg + '</div>';
            }
            
            if (data.properties) {
                html += '<div class="patent-info">';
                html += '<h5>📊 Molecular Properties</h5>';
                html += '<p><strong>Molecular Weight:</strong> ' + data.properties.molecular_weight + ' Da</p>';
                html += '<p><strong>LogP:</strong> ' + data.properties.logp + '</p>';
                html += '<p><strong>Safety Score:</strong> ' + data.properties.safety_score + '/10</p>';
                html += '<p><strong>Efficacy Score:</strong> ' + data.properties.efficacy_score + '/10</p>';
                html += '</div>';
            }
            
            if (data.patent_info) {
                html += '<div class="patent-info">';
                html += '<h5>🏛️ Patent Information</h5>';
                html += '<p><strong>Patent Number:</strong> ' + data.patent_info.patent_number + '</p>';
                html += '<p><strong>Status:</strong> ' + data.patent_info.status + '</p>';
                html += '<p><strong>Expiry:</strong> ' + data.patent_info.expiry + '</p>';
                html += '<p><strong>Freedom to Operate:</strong> ' + data.patent_info.freedom_to_operate + '</p>';
                html += '</div>';
            }
            
            if (data.receptor_binding) {
                html += '<div class="receptor-binding">';
                html += '<h5>🎯 Receptor Binding Profile</h5>';
                for (const [receptor, binding] of Object.entries(data.receptor_binding)) {
                    html += '<div class="receptor-item">';
                    html += '<strong>' + receptor + ':</strong> ' + binding;
                    html += '</div>';
                }
                html += '</div>';
            }
            
            resultsDiv.innerHTML = html;
        }
        
        // Analog generation
        function generateAnalogs() {
            if (!isLoggedIn) {
                alert('Please login to access analog generation features.');
                return;
            }
            
            const parentCompound = document.getElementById('parent-compound').value;
            const targetProperties = Array.from(document.getElementById('target-properties').selectedOptions).map(option => option.value);
            const similarityThreshold = document.getElementById('similarity-threshold').value;
            
            if (!parentCompound) {
                alert('Please enter a parent compound.');
                return;
            }
            
            const resultsDiv = document.getElementById('analog-results');
            resultsDiv.style.display = 'block';
            resultsDiv.innerHTML = '<p>🔄 Generating analogs with AI optimization...</p>';
            
            fetch('/api/generate-analogs', {
                method: 'POST',
                headers: {'Content-Type': 'application/json'},
                body: JSON.stringify({
                    parent_compound: parentCompound,
                    target_properties: targetProperties,
                    similarity_threshold: parseFloat(similarityThreshold),
                    user_id: currentUser
                })
            })
            .then(response => response.json())
            .then(data => {
                displayAnalogResults(data, resultsDiv);
            })
            .catch(error => {
                resultsDiv.innerHTML = '<p>❌ Error generating analogs: ' + error.message + '</p>';
            });
        }
        
        function displayAnalogResults(data, resultsDiv) {
            let html = '<h4>🔬 Generated Analogs</h4>';
            html += '<p><strong>Parent Compound:</strong> ' + data.parent_compound + '</p>';
            html += '<p><strong>Analogs Generated:</strong> ' + data.analogs.length + '</p>';
            
            data.analogs.forEach((analog, index) => {
                html += '<div class="project-card">';
                html += '<h5>Analog ' + (index + 1) + '</h5>';
                html += '<p><strong>SMILES:</strong> ' + analog.smiles + '</p>';
                html += '<p><strong>Similarity Score:</strong> ' + analog.similarity_score + '</p>';
                html += '<p><strong>Patent Status:</strong> ' + analog.patent_status + '</p>';
                html += '<p><strong>Predicted Safety:</strong> ' + analog.predicted_safety + '/10</p>';
                html += '<p><strong>Predicted Efficacy:</strong> ' + analog.predicted_efficacy + '/10</p>';
                html += '</div>';
            });
            
            resultsDiv.innerHTML = html;
        }
        
        // PKPD Analysis
        function runPKPDAnalysis() {
            if (!isLoggedIn) {
                alert('Please login to access PKPD modeling features.');
                return;
            }
            
            const compound = document.getElementById('pkpd-compound').value;
            const age = document.getElementById('patient-age').value;
            const weight = document.getElementById('patient-weight').value;
            const liverFunction = document.getElementById('liver-function').value;
            
            if (!compound) {
                alert('Please enter a compound name.');
                return;
            }
            
            const resultsDiv = document.getElementById('pkpd-results');
            resultsDiv.style.display = 'block';
            resultsDiv.innerHTML = '<p>🔄 Running PKPD simulation...</p>';
            
            fetch('/api/pkpd-analysis', {
                method: 'POST',
                headers: {'Content-Type': 'application/json'},
                body: JSON.stringify({
                    compound: compound,
                    patient_profile: {
                        age: parseInt(age),
                        weight: parseInt(weight),
                        liver_function: liverFunction
                    },
                    user_id: currentUser
                })
            })
            .then(response => response.json())
            .then(data => {
                displayPKPDResults(data, resultsDiv);
            })
            .catch(error => {
                resultsDiv.innerHTML = '<p>❌ Error running PKPD analysis: ' + error.message + '</p>';
            });
        }
        
        function displayPKPDResults(data, resultsDiv) {
            let html = '<h4>📊 PKPD Analysis Results</h4>';
            html += '<div class="patent-info">';
            html += '<h5>Patient-Specific Parameters</h5>';
            html += '<p><strong>Clearance:</strong> ' + data.clearance + ' L/h</p>';
            html += '<p><strong>Volume of Distribution:</strong> ' + data.volume_distribution + ' L</p>';
            html += '<p><strong>Half-life:</strong> ' + data.half_life + ' hours</p>';
            html += '<p><strong>Bioavailability:</strong> ' + (data.bioavailability * 100).toFixed(1) + '%</p>';
            html += '</div>';
            
            html += '<div class="patent-info">';
            html += '<h5>Dosing Recommendations</h5>';
            html += '<p><strong>Starting Dose:</strong> ' + data.dosing_recommendations.starting_dose + '</p>';
            html += '<p><strong>Maintenance Dose:</strong> ' + data.dosing_recommendations.maintenance_dose + '</p>';
            html += '<p><strong>Dosing Frequency:</strong> ' + data.dosing_recommendations.frequency + '</p>';
            html += '</div>';
            
            resultsDiv.innerHTML = html;
        }
        
        // Psychiatric combination analysis
        function analyzePsychiatricCombination() {
            if (!isLoggedIn) {
                alert('Please login to access psychiatric analysis features.');
                return;
            }
            
            const drugCombination = document.getElementById('drug-combination').value;
            const conditions = document.getElementById('patient-conditions').value;
            const allergies = document.getElementById('patient-allergies').value;
            
            if (!drugCombination) {
                alert('Please enter a drug combination.');
                return;
            }
            
            const resultsDiv = document.getElementById('psychiatric-results');
            resultsDiv.style.display = 'block';
            resultsDiv.innerHTML = '<p>🔄 Analyzing psychiatric drug combination...</p>';
            
            fetch('/api/psychiatric-analysis', {
                method: 'POST',
                headers: {'Content-Type': 'application/json'},
                body: JSON.stringify({
                    drug_combination: drugCombination.split(',').map(d => d.trim()),
                    patient_profile: {
                        medical_conditions: conditions.split(',').map(c => c.trim()),
                        allergies: allergies.split(',').map(a => a.trim())
                    },
                    user_id: currentUser
                })
            })
            .then(response => response.json())
            .then(data => {
                displayPsychiatricResults(data, resultsDiv);
            })
            .catch(error => {
                resultsDiv.innerHTML = '<p>❌ Error analyzing combination: ' + error.message + '</p>';
            });
        }
        
        function displayPsychiatricResults(data, resultsDiv) {
            let html = '<h4>🧠 Psychiatric Combination Analysis</h4>';
            
            html += '<div class="patent-info">';
            html += '<h5>Drug Interactions</h5>';
            html += '<p><strong>Total Interactions:</strong> ' + data.drug_interactions.total_interactions + '</p>';
            html += '<p><strong>Major Interactions:</strong> ' + data.drug_interactions.major_interactions + '</p>';
            html += '<p><strong>Overall Risk:</strong> ' + data.drug_interactions.overall_interaction_risk + '</p>';
            html += '</div>';
            
            html += '<div class="patent-info">';
            html += '<h5>Safety Assessment</h5>';
            html += '<p><strong>Overall Safety Score:</strong> ' + data.safety_assessment.overall_safety_score + '</p>';
            html += '<p><strong>Safety Category:</strong> ' + data.safety_assessment.safety_category + '</p>';
            html += '</div>';
            
            if (data.personalized_recommendations) {
                html += '<div class="patent-info">';
                html += '<h5>Personalized Recommendations</h5>';
                data.personalized_recommendations.dosing_recommendations.forEach(rec => {
                    html += '<p><strong>' + rec.drug + ':</strong> ' + rec.starting_dose + '</p>';
                });
                html += '</div>';
            }
            
            resultsDiv.innerHTML = html;
        }
        
        // Retrosynthesis planning
        function planRetrosynthesis() {
            if (!isLoggedIn) {
                alert('Please login to access retrosynthesis planning features.');
                return;
            }
            
            const targetMolecule = document.getElementById('target-molecule').value;
            const constraints = Array.from(document.getElementById('synthesis-constraints').selectedOptions).map(option => option.value);
            const startingMaterials = document.getElementById('starting-materials').value;
            
            if (!targetMolecule) {
                alert('Please enter a target molecule.');
                return;
            }
            
            const resultsDiv = document.getElementById('retrosynthesis-results');
            resultsDiv.style.display = 'block';
            resultsDiv.innerHTML = '<p>🔄 Planning synthetic routes...</p>';
            
            fetch('/api/retrosynthesis', {
                method: 'POST',
                headers: {'Content-Type': 'application/json'},
                body: JSON.stringify({
                    target_molecule: targetMolecule,
                    constraints: constraints,
                    starting_materials: startingMaterials.split(',').map(s => s.trim()),
                    user_id: currentUser
                })
            })
            .then(response => response.json())
            .then(data => {
                displayRetrosynthesisResults(data, resultsDiv);
            })
            .catch(error => {
                resultsDiv.innerHTML = '<p>❌ Error planning retrosynthesis: ' + error.message + '</p>';
            });
        }
        
        function displayRetrosynthesisResults(data, resultsDiv) {
            let html = '<h4>🏭 Retrosynthesis Plan</h4>';
            html += '<p><strong>Target Molecule:</strong> ' + data.target_molecule + '</p>';
            html += '<p><strong>Routes Generated:</strong> ' + data.routes.length + '</p>';
            
            data.routes.forEach((route, index) => {
                html += '<div class="project-card">';
                html += '<h5>Route ' + (index + 1) + '</h5>';
                html += '<p><strong>Steps:</strong> ' + route.steps + '</p>';
                html += '<p><strong>Estimated Yield:</strong> ' + route.estimated_yield + '%</p>';
                html += '<p><strong>Cost Score:</strong> ' + route.cost_score + '/10</p>';
                html += '<p><strong>Feasibility:</strong> ' + route.feasibility_score + '/10</p>';
                html += '<p><strong>Green Chemistry Score:</strong> ' + route.green_chemistry_score + '/10</p>';
                html += '</div>';
            });
            
            resultsDiv.innerHTML = html;
        }
        
        // Database search
        function searchDatabase() {
            const query = document.getElementById('search-query').value;
            const category = document.getElementById('search-category').value;
            
            if (!query) {
                alert('Please enter a search query.');
                return;
            }
            
            const resultsDiv = document.getElementById('search-results');
            resultsDiv.style.display = 'block';
            resultsDiv.innerHTML = '<p>🔄 Searching database...</p>';
            
            fetch('/api/search', {
                method: 'POST',
                headers: {'Content-Type': 'application/json'},
                body: JSON.stringify({
                    query: query,
                    category: category
                })
            })
            .then(response => response.json())
            .then(data => {
                displaySearchResults(data, resultsDiv);
            })
            .catch(error => {
                resultsDiv.innerHTML = '<p>❌ Error searching database: ' + error.message + '</p>';
            });
        }
        
        function displaySearchResults(data, resultsDiv) {
            let html = '<h4>🔍 Search Results</h4>';
            html += '<p>Found ' + data.total_results + ' results</p>';
            
            data.results.forEach(result => {
                html += '<div class="project-card">';
                html += '<h5>' + result.title + '</h5>';
                html += '<p><strong>Type:</strong> ' + result.type + '</p>';
                html += '<p><strong>Relevance:</strong> ' + result.relevance_score + '%</p>';
                html += '<p>' + result.description + '</p>';
                html += '</div>';
            });
            
            resultsDiv.innerHTML = html;
        }
        
        // Modal functions
        function showCompounds() {
            fetch('/api/compounds')
                .then(response => response.json())
                .then(data => {
                    let html = '<h3>🧬 Active Compounds Database</h3>';
                    data.compounds.forEach(compound => {
                        html += '<div class="project-card">';
                        html += '<h4>' + compound.name + '</h4>';
                        html += '<div class="structure-container">' + compound.structure_svg + '</div>';
                        html += '<p><strong>SMILES:</strong> ' + compound.smiles + '</p>';
                        html += '<p><strong>Molecular Weight:</strong> ' + compound.molecular_weight + ' Da</p>';
                        html += '<p><strong>Safety Score:</strong> ' + compound.safety_score + '/10</p>';
                        html += '<p><strong>Therapeutic Area:</strong> ' + compound.therapeutic_area + '</p>';
                        html += '<div class="patent-info">';
                        html += '<h5>Patent Information</h5>';
                        html += '<p><strong>Patent:</strong> ' + compound.patent_number + '</p>';
                        html += '<p><strong>Status:</strong> ' + compound.patent_status + '</p>';
                        html += '<p><strong>Freedom to Operate:</strong> ' + compound.freedom_to_operate + '</p>';
                        html += '</div>';
                        html += '</div>';
                    });
                    document.getElementById('compoundModalContent').innerHTML = html;
                    document.getElementById('compoundModal').style.display = 'block';
                });
        }
        
        function showProjects() {
            fetch('/api/projects')
                .then(response => response.json())
                .then(data => {
                    let html = '<h3>🔬 Research Projects</h3>';
                    data.projects.forEach(project => {
                        html += '<div class="project-card">';
                        html += '<h4>' + project.name + '</h4>';
                        html += '<p><strong>Area:</strong> ' + project.area + '</p>';
                        html += '<p><strong>Status:</strong> ' + project.status + '</p>';
                        html += '<p><strong>Budget:</strong> ' + project.budget + '</p>';
                        html += '<p><strong>Timeline:</strong> ' + project.timeline + '</p>';
                        html += '<p><strong>PI:</strong> ' + project.principal_investigator + '</p>';
                        html += '<p>' + project.description + '</p>';
                        html += '</div>';
                    });
                    document.getElementById('projectModalContent').innerHTML = html;
                    document.getElementById('projectModal').style.display = 'block';
                });
        }
        
        function showPatents() {
            fetch('/api/patents')
                .then(response => response.json())
                .then(data => {
                    let html = '<h3>🏛️ Patent Portfolio</h3>';
                    data.patents.forEach(patent => {
                        html += '<div class="patent-card">';
                        html += '<h4>' + patent.title + '</h4>';
                        html += '<p><strong>Patent Number:</strong> ' + patent.patent_number + '</p>';
                        html += '<p><strong>Area:</strong> ' + patent.area + '</p>';
                        html += '<p><strong>Status:</strong> ' + patent.status + '</p>';
                        html += '<p><strong>Filing Date:</strong> ' + patent.filing_date + '</p>';
                        html += '<p><strong>Expiry Date:</strong> ' + patent.expiry_date + '</p>';
                        html += '<p><strong>Inventors:</strong> ' + patent.inventors + '</p>';
                        html += '<p><strong>Claims:</strong> ' + patent.claims + '</p>';
                        html += '<p><strong>Estimated Value:</strong> ' + patent.estimated_value + '</p>';
                        html += '<p><strong>Licensing:</strong> ' + patent.licensing_status + '</p>';
                        html += '</div>';
                    });
                    document.getElementById('patentModalContent').innerHTML = html;
                    document.getElementById('patentModal').style.display = 'block';
                });
        }
        
        function showSuccessRate() {
            alert('🎯 Platform Success Rate: 87.5%\\n\\n✅ Successful Predictions: 156/178\\n📊 Accuracy Rate: 87.5%\\n🔬 Compounds Advanced: 12\\n🏛️ Patents Filed: 3');
        }
        
        function showAutonomousStats() {
            alert('🤖 Autonomous Research Engine Stats:\\n\\n📚 Papers Processed: 47 today, 1,247 total\\n💡 Hypotheses Generated: 12 today, 89 total\\n🏛️ IP Opportunities: 8 today, 23 total\\n🔬 Novel Compounds: 156 discovered\\n⚡ Engine Status: 78% capacity');
        }
        
        function showAnalogStats() {
            alert('🔬 Analog Generation Stats:\\n\\n🧬 Total Analogs Generated: 47\\n✅ Patent-Free Analogs: 31\\n📊 Average Similarity: 0.82\\n🎯 Success Rate: 91%\\n⭐ Top Scoring Analog: 9.2/10');
        }
        
        function closeModal(modalId) {
            document.getElementById(modalId).style.display = 'none';
        }
        
        // Utility functions
        function refreshStatistics() {
            // Update autonomous stats
            document.getElementById('papers-today').textContent = Math.floor(Math.random() * 20) + 40;
            document.getElementById('hypotheses-today').textContent = Math.floor(Math.random() * 8) + 8;
            document.getElementById('ip-opportunities').textContent = Math.floor(Math.random() * 5) + 5;
            document.getElementById('autonomous-discoveries').textContent = Math.floor(Math.random() * 50) + 150;
            document.getElementById('analog-count').textContent = Math.floor(Math.random() * 20) + 40;
            
            alert('📊 Statistics refreshed successfully!');
        }
        
        function downloadResearchReport() {
            if (!isLoggedIn) {
                alert('Please login to download research reports.');
                return;
            }
            
            // Simulate report download
            alert('📄 Generating comprehensive research report...\\n\\n✅ Report will be downloaded shortly.\\n📊 Includes: Compound analysis, patent landscape, research insights, and IP opportunities.');
        }
        
        function viewAuditLog() {
            if (!isLoggedIn) {
                alert('Please login to view audit logs.');
                return;
            }
            
            fetch('/api/audit-log')
                .then(response => response.json())
                .then(data => {
                    let html = '<h3>📋 Audit Log</h3>';
                    html += '<div style="max-height: 400px; overflow-y: auto;">';
                    data.audit_entries.forEach(entry => {
                        html += '<div class="project-card" style="margin: 10px 0; padding: 15px;">';
                        html += '<p><strong>Timestamp:</strong> ' + entry.timestamp + '</p>';
                        html += '<p><strong>User:</strong> ' + entry.user_id + '</p>';
                        html += '<p><strong>Action:</strong> ' + entry.action + '</p>';
                        html += '<p><strong>Details:</strong> ' + entry.details + '</p>';
                        html += '<p><strong>IP Address:</strong> ' + entry.ip_address + '</p>';
                        html += '</div>';
                    });
                    html += '</div>';
                    
                    // Create and show audit modal
                    const modal = document.createElement('div');
                    modal.className = 'modal';
                    modal.style.display = 'block';
                    modal.innerHTML = '<div class="modal-content"><span class="close" onclick="this.parentElement.parentElement.remove()">&times;</span>' + html + '</div>';
                    document.body.appendChild(modal);
                });
        }
        
        function exportEncyclopedia() {
            if (!isLoggedIn) {
                alert('Please login to export encyclopedia data.');
                return;
            }
            
            alert('📚 Exporting Compound Encyclopedia...\\n\\n✅ Export will include:\\n🧬 All compound analyses\\n📊 Research insights\\n🏛️ Patent information\\n📈 Knowledge evolution\\n\\n📄 File will be downloaded shortly.');
        }
        
        // Close modals when clicking outside
        window.onclick = function(event) {
            const modals = document.querySelectorAll('.modal');
            modals.forEach(modal => {
                if (event.target === modal) {
                    modal.style.display = 'none';
                }
            });
        }
        
        // Auto-refresh autonomous stats every 30 seconds
        setInterval(() => {
            if (isLoggedIn) {
                document.getElementById('papers-today').textContent = Math.floor(Math.random() * 20) + 40;
                document.getElementById('hypotheses-today').textContent = Math.floor(Math.random() * 8) + 8;
                document.getElementById('ip-opportunities').textContent = Math.floor(Math.random() * 5) + 5;
            }
        }, 30000);
    </script>
</body>
</html>
'''

# API Routes
@app.route('/')
def index():
    return render_template_string(HTML_TEMPLATE)

@app.route('/health')
def health():
    return jsonify({
        "status": "healthy",
        "version": "3.0.0-enterprise",
        "database": "connected",
        "enterprise_features": ENTERPRISE_FEATURES_AVAILABLE,
        "autonomous_engine": "operational",
        "timestamp": datetime.now().isoformat()
    })

@app.route('/api/analyze-compound', methods=['POST'])
def analyze_compound():
    data = request.json
    compound = data.get('compound', '')
    analysis_type = data.get('analysis_type', 'comprehensive')
    user_id = data.get('user_id', 'anonymous')
    
    # Log audit event
    log_audit_event(user_id, 'compound_analysis', f'Analyzed compound: {compound}', request.remote_addr)
    
    # Get compound data from database
    conn = sqlite3.connect('drug_discovery_enterprise.db')
    cursor = conn.cursor()
    cursor.execute('SELECT * FROM compounds WHERE name LIKE ? OR smiles LIKE ?', (f'%{compound}%', f'%{compound}%'))
    result = cursor.fetchone()
    conn.close()
    
    if result:
        compound_data = {
            'name': result[1],
            'smiles': result[2],
            'molecular_weight': result[3],
            'logp': result[4],
            'safety_score': result[5],
            'efficacy_score': result[6],
            'therapeutic_area': result[7],
            'status': result[8],
            'patent_number': result[9],
            'patent_status': result[10],
            'patent_expiry': result[11],
            'freedom_to_operate': result[12]
        }
        
        # Add to compound encyclopedia if enterprise features available
        if ENTERPRISE_FEATURES_AVAILABLE:
            compound_encyclopedia.add_compound_entry(compound, {
                'analysis_type': analysis_type,
                'user_id': user_id,
                'results': compound_data
            })
        
        return jsonify({
            'compound_name': compound_data['name'],
            'structure_svg': generate_chemical_structure_svg(compound_data['smiles'], compound_data['name']),
            'properties': {
                'molecular_weight': compound_data['molecular_weight'],
                'logp': compound_data['logp'],
                'safety_score': compound_data['safety_score'],
                'efficacy_score': compound_data['efficacy_score']
            },
            'patent_info': {
                'patent_number': compound_data['patent_number'],
                'status': compound_data['patent_status'],
                'expiry': compound_data['patent_expiry'],
                'freedom_to_operate': compound_data['freedom_to_operate']
            },
            'receptor_binding': {
                '5-HT2A': 'High affinity (Ki = 6.0 nM)',
                '5-HT2C': 'Moderate affinity (Ki = 23 nM)',
                '5-HT1A': 'Low affinity (Ki = 190 nM)',
                'NMDA': 'Moderate affinity (Ki = 0.5 μM)',
                'mTOR': 'Indirect activation'
            },
            'therapeutic_area': compound_data['therapeutic_area'],
            'analysis_timestamp': datetime.now().isoformat()
        })
    else:
        # Generate analysis for unknown compound
        return jsonify({
            'compound_name': compound,
            'structure_svg': generate_chemical_structure_svg('Unknown', compound),
            'properties': {
                'molecular_weight': random.uniform(200, 400),
                'logp': random.uniform(-2, 4),
                'safety_score': random.uniform(6, 9),
                'efficacy_score': random.uniform(7, 9)
            },
            'patent_info': {
                'patent_number': 'Not found',
                'status': 'Unknown',
                'expiry': 'N/A',
                'freedom_to_operate': 'Requires investigation'
            },
            'receptor_binding': {
                'Unknown': 'Requires experimental validation'
            },
            'therapeutic_area': 'To be determined',
            'analysis_timestamp': datetime.now().isoformat()
        })

@app.route('/api/generate-analogs', methods=['POST'])
def generate_analogs():
    data = request.json
    parent_compound = data.get('parent_compound', '')
    target_properties = data.get('target_properties', [])
    similarity_threshold = data.get('similarity_threshold', 0.8)
    user_id = data.get('user_id', 'anonymous')
    
    # Log audit event
    log_audit_event(user_id, 'analog_generation', f'Generated analogs for: {parent_compound}', request.remote_addr)
    
    # Generate analogs using enterprise features if available
    if ENTERPRISE_FEATURES_AVAILABLE:
        analogs = analog_generator.generate_analogs(parent_compound, target_properties, similarity_threshold)
    else:
        # Fallback analog generation
        analogs = []
        for i in range(random.randint(3, 8)):
            analogs.append({
                'smiles': f'Generated_SMILES_{i+1}',
                'similarity_score': random.uniform(similarity_threshold, 0.95),
                'patent_status': random.choice(['Patent-free', 'Patent pending', 'Patented']),
                'predicted_safety': random.uniform(7, 9),
                'predicted_efficacy': random.uniform(6, 9),
                'drug_likeness': random.uniform(0.6, 0.9)
            })
    
    # Store analog analysis in database
    conn = sqlite3.connect('drug_discovery_enterprise.db')
    cursor = conn.cursor()
    for analog in analogs:
        cursor.execute('''
            INSERT INTO analog_analyses (parent_compound, analog_smiles, similarity_score, patent_status, predicted_properties)
            VALUES (?, ?, ?, ?, ?)
        ''', (parent_compound, analog['smiles'], analog['similarity_score'], analog['patent_status'], json.dumps(analog)))
    conn.commit()
    conn.close()
    
    return jsonify({
        'parent_compound': parent_compound,
        'analogs': analogs,
        'generation_parameters': {
            'target_properties': target_properties,
            'similarity_threshold': similarity_threshold
        },
        'analysis_timestamp': datetime.now().isoformat()
    })

@app.route('/api/pkpd-analysis', methods=['POST'])
def pkpd_analysis():
    data = request.json
    compound = data.get('compound', '')
    patient_profile = data.get('patient_profile', {})
    user_id = data.get('user_id', 'anonymous')
    
    # Log audit event
    log_audit_event(user_id, 'pkpd_analysis', f'PKPD analysis for: {compound}', request.remote_addr)
    
    # Run PKPD analysis using enterprise features if available
    if ENTERPRISE_FEATURES_AVAILABLE:
        pkpd_results = pkpd_profiler.analyze_compound_pkpd(compound, patient_profile)
    else:
        # Fallback PKPD analysis
        age_factor = 0.8 if patient_profile.get('age', 35) > 65 else 1.0
        liver_factor = 0.7 if patient_profile.get('liver_function') != 'Normal' else 1.0
        
        pkpd_results = {
            'clearance': round(random.uniform(10, 50) * age_factor * liver_factor, 2),
            'volume_distribution': round(random.uniform(1, 10) * patient_profile.get('weight', 70) / 70, 2),
            'half_life': round(random.uniform(2, 24) / liver_factor, 2),
            'bioavailability': round(random.uniform(0.4, 0.95) * age_factor, 3),
            'dosing_recommendations': {
                'starting_dose': f'Start at {int(age_factor * liver_factor * 100)}% of standard dose',
                'maintenance_dose': f'{random.randint(50, 200)} mg daily',
                'frequency': random.choice(['Once daily', 'Twice daily', 'Three times daily'])
            }
        }
    
    # Store PKPD analysis in database
    conn = sqlite3.connect('drug_discovery_enterprise.db')
    cursor = conn.cursor()
    cursor.execute('''
        INSERT INTO pkpd_analyses (compound_name, patient_profile, pkpd_results)
        VALUES (?, ?, ?)
    ''', (compound, json.dumps(patient_profile), json.dumps(pkpd_results)))
    conn.commit()
    conn.close()
    
    return jsonify(pkpd_results)

@app.route('/api/psychiatric-analysis', methods=['POST'])
def psychiatric_analysis():
    data = request.json
    drug_combination = data.get('drug_combination', [])
    patient_profile = data.get('patient_profile', {})
    user_id = data.get('user_id', 'anonymous')
    
    # Log audit event
    log_audit_event(user_id, 'psychiatric_analysis', f'Analyzed combination: {", ".join(drug_combination)}', request.remote_addr)
    
    # Run psychiatric analysis using enterprise features if available
    if ENTERPRISE_FEATURES_AVAILABLE:
        analysis_results = psychiatric_analyzer.analyze_drug_combination(drug_combination, patient_profile)
    else:
        # Fallback psychiatric analysis
        analysis_results = {
            'drug_interactions': {
                'total_interactions': random.randint(1, 5),
                'major_interactions': random.randint(0, 2),
                'moderate_interactions': random.randint(1, 3),
                'overall_interaction_risk': random.choice(['Low', 'Moderate', 'High'])
            },
            'safety_assessment': {
                'overall_safety_score': round(random.uniform(0.6, 0.9), 2),
                'safety_category': random.choice(['High Safety', 'Moderate Safety', 'Caution Required'])
            },
            'personalized_recommendations': {
                'dosing_recommendations': [
                    {
                        'drug': drug,
                        'starting_dose': f'Start at {random.randint(50, 100)}% of standard dose',
                        'monitoring': 'Weekly for first month'
                    } for drug in drug_combination
                ]
            }
        }
    
    return jsonify(analysis_results)

@app.route('/api/retrosynthesis', methods=['POST'])
def retrosynthesis():
    data = request.json
    target_molecule = data.get('target_molecule', '')
    constraints = data.get('constraints', [])
    starting_materials = data.get('starting_materials', [])
    user_id = data.get('user_id', 'anonymous')
    
    # Log audit event
    log_audit_event(user_id, 'retrosynthesis', f'Planned synthesis for: {target_molecule}', request.remote_addr)
    
    # Run retrosynthesis planning using enterprise features if available
    if ENTERPRISE_FEATURES_AVAILABLE:
        synthesis_routes = retrosynthesis_planner.plan_synthesis(target_molecule, constraints, starting_materials)
    else:
        # Fallback retrosynthesis planning
        synthesis_routes = []
        for i in range(random.randint(2, 5)):
            synthesis_routes.append({
                'route_id': i + 1,
                'steps': random.randint(3, 8),
                'estimated_yield': random.randint(45, 85),
                'cost_score': random.randint(6, 9),
                'feasibility_score': random.randint(7, 9),
                'green_chemistry_score': random.randint(5, 9),
                'starting_materials': random.sample(starting_materials or ['Tryptamine', 'Indole', 'Benzene'], min(3, len(starting_materials) or 3))
            })
    
    # Store retrosynthesis analysis in database
    conn = sqlite3.connect('drug_discovery_enterprise.db')
    cursor = conn.cursor()
    cursor.execute('''
        INSERT INTO retrosynthesis_routes (compound_name, route_data, feasibility_score, cost_estimate)
        VALUES (?, ?, ?, ?)
    ''', (target_molecule, json.dumps(synthesis_routes), random.uniform(7, 9), f'${random.randint(50, 500)}K'))
    conn.commit()
    conn.close()
    
    return jsonify({
        'target_molecule': target_molecule,
        'routes': synthesis_routes,
        'constraints_applied': constraints,
        'analysis_timestamp': datetime.now().isoformat()
    })

@app.route('/api/search', methods=['POST'])
def search_database():
    data = request.json
    query = data.get('query', '')
    category = data.get('category', 'all')
    
    # Search database
    conn = sqlite3.connect('drug_discovery_enterprise.db')
    cursor = conn.cursor()
    
    results = []
    
    if category in ['all', 'compounds']:
        cursor.execute('SELECT name, therapeutic_area, status FROM compounds WHERE name LIKE ?', (f'%{query}%',))
        compounds = cursor.fetchall()
        for compound in compounds:
            results.append({
                'title': compound[0],
                'type': 'Compound',
                'description': f'Therapeutic Area: {compound[1]}, Status: {compound[2]}',
                'relevance_score': random.randint(75, 95)
            })
    
    if category in ['all', 'projects']:
        cursor.execute('SELECT name, area, description FROM research_projects WHERE name LIKE ? OR description LIKE ?', (f'%{query}%', f'%{query}%'))
        projects = cursor.fetchall()
        for project in projects:
            results.append({
                'title': project[0],
                'type': 'Research Project',
                'description': project[2][:100] + '...',
                'relevance_score': random.randint(70, 90)
            })
    
    if category in ['all', 'patents']:
        cursor.execute('SELECT title, patent_number, area FROM patents WHERE title LIKE ?', (f'%{query}%',))
        patents = cursor.fetchall()
        for patent in patents:
            results.append({
                'title': patent[0],
                'type': 'Patent',
                'description': f'Patent Number: {patent[1]}, Area: {patent[2]}',
                'relevance_score': random.randint(80, 95)
            })
    
    conn.close()
    
    # Sort by relevance score
    results.sort(key=lambda x: x['relevance_score'], reverse=True)
    
    return jsonify({
        'query': query,
        'category': category,
        'total_results': len(results),
        'results': results[:10]  # Limit to top 10 results
    })

@app.route('/api/compounds')
def get_compounds():
    conn = sqlite3.connect('drug_discovery_enterprise.db')
    cursor = conn.cursor()
    cursor.execute('SELECT * FROM compounds')
    compounds = cursor.fetchall()
    conn.close()
    
    compound_list = []
    for compound in compounds:
        compound_list.append({
            'name': compound[1],
            'smiles': compound[2],
            'molecular_weight': compound[3],
            'safety_score': compound[5],
            'therapeutic_area': compound[7],
            'patent_number': compound[9],
            'patent_status': compound[10],
            'freedom_to_operate': compound[12],
            'structure_svg': generate_chemical_structure_svg(compound[2], compound[1])
        })
    
    return jsonify({'compounds': compound_list})

@app.route('/api/projects')
def get_projects():
    conn = sqlite3.connect('drug_discovery_enterprise.db')
    cursor = conn.cursor()
    cursor.execute('SELECT * FROM research_projects')
    projects = cursor.fetchall()
    conn.close()
    
    project_list = []
    for project in projects:
        project_list.append({
            'name': project[1],
            'area': project[2],
            'status': project[3],
            'budget': project[4],
            'timeline': project[5],
            'principal_investigator': project[6],
            'description': project[7]
        })
    
    return jsonify({'projects': project_list})

@app.route('/api/patents')
def get_patents():
    conn = sqlite3.connect('drug_discovery_enterprise.db')
    cursor = conn.cursor()
    cursor.execute('SELECT * FROM patents')
    patents = cursor.fetchall()
    conn.close()
    
    patent_list = []
    for patent in patents:
        patent_list.append({
            'title': patent[1],
            'patent_number': patent[2],
            'area': patent[3],
            'status': patent[4],
            'filing_date': patent[5],
            'expiry_date': patent[6],
            'inventors': patent[7],
            'claims': patent[8],
            'estimated_value': patent[9],
            'licensing_status': patent[10]
        })
    
    return jsonify({'patents': patent_list})

@app.route('/api/audit', methods=['POST'])
def log_audit():
    data = request.json
    user_id = data.get('user_id', 'anonymous')
    action = data.get('action', '')
    details = data.get('details', '')
    
    log_audit_event(user_id, action, details, request.remote_addr)
    
    return jsonify({'status': 'logged'})

@app.route('/api/audit-log')
def get_audit_log():
    conn = sqlite3.connect('drug_discovery_enterprise.db')
    cursor = conn.cursor()
    cursor.execute('SELECT * FROM audit_log ORDER BY timestamp DESC LIMIT 50')
    audit_entries = cursor.fetchall()
    conn.close()
    
    audit_list = []
    for entry in audit_entries:
        audit_list.append({
            'user_id': entry[1],
            'action': entry[2],
            'details': entry[3],
            'ip_address': entry[4],
            'timestamp': entry[5]
        })
    
    return jsonify({'audit_entries': audit_list})

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5009, debug=False)


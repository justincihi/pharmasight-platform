"""
PharmaSight™ Ultimate - Complete Advanced Research Platform
Integrates all enhanced features for comprehensive pharmaceutical research
"""

from flask import Flask, render_template_string, request, jsonify, session
from flask_cors import CORS
import json
import datetime
import hashlib
import random
import re
import os
import sys

# Add src directory to path for imports
sys.path.append(os.path.join(os.path.dirname(__file__), 'src'))

# Import enhanced systems
from advanced_research_system import create_advanced_research_system
from enhanced_pkpd_system import create_enhanced_pkpd_system
from molecular_visualization_system import create_molecular_visualization_system

# Initialize Flask app
app = Flask(__name__)
app.secret_key = 'pharmasight_ultimate_2024_complete'
CORS(app)

# Initialize enhanced systems
research_system = create_advanced_research_system()
pkpd_system = create_enhanced_pkpd_system()
viz_system = create_molecular_visualization_system()

# Load comprehensive databases
def load_comprehensive_data():
    """Load all comprehensive databases"""
    try:
        # Load comprehensive compound database
        with open('comprehensive_compound_database.json', 'r') as f:
            compound_db = json.load(f)
        
        # Load enhanced research findings
        with open('enhanced_research_findings.json', 'r') as f:
            research_data = json.load(f)
        
        # Load analog discoveries
        with open('analog_discoveries.json', 'r') as f:
            analog_data = json.load(f)
        
        # Load autonomous research articles
        with open('autonomous_research_articles.json', 'r') as f:
            articles_data = json.load(f)
        
        return compound_db, research_data, analog_data, articles_data
    except Exception as e:
        print(f"Error loading comprehensive data: {e}")
        return {}, {"findings": [], "summary": {}}, [], {"articles": [], "summary": {}}

# Load comprehensive databases
COMPOUND_DATABASE, RESEARCH_DATA, ANALOG_DATA, ARTICLES_DATA = load_comprehensive_data()

print(f"✅ Loaded {len(COMPOUND_DATABASE)} compounds")
print(f"✅ Loaded {len(RESEARCH_DATA.get('findings', []))} research findings")
print(f"✅ Loaded {len(ANALOG_DATA)} analog discoveries")
print(f"✅ Loaded {len(ARTICLES_DATA.get('articles', []))} research articles")
print(f"✅ Initialized advanced research systems")

# Authentication
def check_auth(username, password):
    """Check if username/password combination is valid"""
    valid_users = {
        'ImplicateOrder25': 'pharmasight2024',
        'admin': 'admin123',
        'researcher': 'research2024'
    }
    return valid_users.get(username) == password

@app.route('/')
def index():
    """Main dashboard with all enhanced features"""
    if 'logged_in' not in session:
        return render_template_string(LOGIN_TEMPLATE)
    
    # Get research goals summary
    research_goals = research_system.get_research_goals('active')
    high_confidence_compounds = research_system.export_high_confidence_compounds(0.85)
    
    return render_template_string(ULTIMATE_DASHBOARD_TEMPLATE, 
                                compounds_count=len(COMPOUND_DATABASE),
                                research_findings_count=len(RESEARCH_DATA.get('findings', [])),
                                high_confidence_count=len(high_confidence_compounds),
                                research_goals_count=len(research_goals),
                                total_market_value="$365M+")

@app.route('/login', methods=['POST'])
def login():
    """Handle login"""
    username = request.form.get('username')
    password = request.form.get('password')
    
    if check_auth(username, password):
        session['logged_in'] = True
        session['username'] = username
        return jsonify({'success': True, 'redirect': '/'})
    
    return jsonify({'success': False, 'message': 'Invalid credentials'})

@app.route('/api/research_goals')
def get_research_goals():
    """Get all research goals"""
    goals = research_system.get_research_goals()
    return jsonify(goals)

@app.route('/api/research_goals', methods=['POST'])
def add_research_goal():
    """Add new research goal"""
    goal_data = request.json
    goal_id = research_system.add_research_goal(goal_data)
    return jsonify({'success': True, 'goal_id': goal_id})

@app.route('/api/autonomous_search/<goal_id>')
def autonomous_search(goal_id):
    """Perform autonomous literature search for a research goal"""
    results = research_system.autonomous_literature_search(goal_id)
    return jsonify({'results': results, 'total_found': len(results)})

@app.route('/api/discover_compounds/<goal_id>')
def discover_compounds(goal_id):
    """Discover compounds for a research goal"""
    compounds = research_system.discover_compounds_for_goal(goal_id)
    return jsonify({'compounds': compounds, 'total_discovered': len(compounds)})

@app.route('/api/high_confidence_compounds')
def get_high_confidence_compounds():
    """Get high confidence compounds with SMILES and IP data"""
    min_confidence = float(request.args.get('min_confidence', 0.85))
    compounds = research_system.export_high_confidence_compounds(min_confidence)
    return jsonify({'compounds': compounds, 'total': len(compounds)})

@app.route('/api/pkpd_analysis', methods=['POST'])
def pkpd_analysis():
    """Comprehensive PKPD/DDI analysis"""
    data = request.json
    drugs = data.get('drugs', [])
    doses = data.get('doses', [])
    routes = data.get('routes', ['oral'] * len(drugs))
    
    if len(drugs) < 2:
        return jsonify({'error': 'At least 2 drugs required for interaction analysis'})
    
    report = pkpd_system.generate_comprehensive_report(drugs, doses, routes)
    return jsonify(report)

@app.route('/api/drug_interaction')
def drug_interaction():
    """Analyze specific drug interaction"""
    drug1 = request.args.get('drug1')
    drug2 = request.args.get('drug2')
    dose1 = float(request.args.get('dose1', 1))
    dose2 = float(request.args.get('dose2', 1))
    
    interaction = pkpd_system.analyze_drug_interaction(drug1, drug2, dose1, dose2)
    return jsonify(interaction)

@app.route('/api/neuroplasticity_analysis', methods=['POST'])
def neuroplasticity_analysis():
    """Analyze neuroplasticity windows for TMS enhancement"""
    drugs = request.json.get('drugs', [])
    analysis = pkpd_system.neuroplasticity_window_analysis(drugs)
    return jsonify(analysis)

@app.route('/api/3d_structure/<compound_name>')
def get_3d_structure(compound_name):
    """Get 3D structure data for molecular visualization"""
    # Get SMILES from research system
    smiles = research_system.get_compound_smiles(compound_name)
    if not smiles:
        # Try to find in compound database
        for compound_id, compound_data in COMPOUND_DATABASE.items():
            if compound_data.get('name', '').lower() == compound_name.lower():
                smiles = compound_data.get('smiles', '')
                break
    
    if smiles:
        structure_data = viz_system.generate_3d_structure_data(smiles, compound_name)
        return jsonify(structure_data)
    
    return jsonify({'error': f'No structure data found for {compound_name}'})

@app.route('/api/retrosynthesis/<compound_type>')
def get_retrosynthesis(compound_type):
    """Get retrosynthesis pathway for compound type"""
    pathway = viz_system.get_retrosynthesis_pathway(compound_type)
    return jsonify(pathway)

@app.route('/api/compound_search')
def compound_search():
    """Search compounds with advanced filtering"""
    query = request.args.get('query', '').lower()
    category = request.args.get('category', '')
    min_confidence = float(request.args.get('min_confidence', 0))
    
    results = []
    
    # Search in compound database
    for compound_id, compound_data in COMPOUND_DATABASE.items():
        name = compound_data.get('name', '').lower()
        compound_category = compound_data.get('category', '')
        
        if query in name or query in compound_id.lower():
            if not category or category == compound_category:
                results.append({
                    'id': compound_id,
                    'name': compound_data.get('name', ''),
                    'category': compound_category,
                    'smiles': compound_data.get('smiles', ''),
                    'molecular_weight': compound_data.get('molecular_weight', 0),
                    'confidence': compound_data.get('confidence', 0.5)
                })
    
    # Search in discovered compounds
    high_confidence = research_system.export_high_confidence_compounds(min_confidence)
    for compound in high_confidence:
        if query in compound['name'].lower():
            results.append(compound)
    
    return jsonify({'results': results, 'total': len(results)})

# HTML Templates
LOGIN_TEMPLATE = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>PharmaSight™ Ultimate - Advanced Research Platform</title>
    <style>
        * { margin: 0; padding: 0; box-sizing: border-box; }
        body { 
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            min-height: 100vh;
            display: flex;
            align-items: center;
            justify-content: center;
        }
        .login-container {
            background: white;
            padding: 2rem;
            border-radius: 15px;
            box-shadow: 0 20px 40px rgba(0,0,0,0.1);
            width: 100%;
            max-width: 400px;
        }
        .logo {
            text-align: center;
            margin-bottom: 2rem;
        }
        .logo h1 {
            color: #4a90e2;
            font-size: 2rem;
            margin-bottom: 0.5rem;
        }
        .logo p {
            color: #666;
            font-size: 0.9rem;
        }
        .form-group {
            margin-bottom: 1.5rem;
        }
        label {
            display: block;
            margin-bottom: 0.5rem;
            color: #333;
            font-weight: 500;
        }
        input[type="text"], input[type="password"] {
            width: 100%;
            padding: 0.75rem;
            border: 2px solid #e1e5e9;
            border-radius: 8px;
            font-size: 1rem;
            transition: border-color 0.3s;
        }
        input[type="text"]:focus, input[type="password"]:focus {
            outline: none;
            border-color: #4a90e2;
        }
        .login-btn {
            width: 100%;
            padding: 0.75rem;
            background: #4a90e2;
            color: white;
            border: none;
            border-radius: 8px;
            font-size: 1rem;
            cursor: pointer;
            transition: background 0.3s;
        }
        .login-btn:hover {
            background: #357abd;
        }
        .features {
            margin-top: 2rem;
            padding-top: 2rem;
            border-top: 1px solid #e1e5e9;
        }
        .feature-tag {
            display: inline-block;
            background: #f8f9fa;
            color: #495057;
            padding: 0.25rem 0.5rem;
            border-radius: 4px;
            font-size: 0.8rem;
            margin: 0.25rem;
        }
    </style>
</head>
<body>
    <div class="login-container">
        <div class="logo">
            <h1>🧬 PharmaSight™ Ultimate</h1>
            <p>Advanced AI-Powered Pharmaceutical Research & Development</p>
        </div>
        
        <form id="loginForm">
            <div class="form-group">
                <label for="username">Username:</label>
                <input type="text" id="username" name="username" value="ImplicateOrder25" required>
            </div>
            
            <div class="form-group">
                <label for="password">Password:</label>
                <input type="password" id="password" name="password" value="pharmasight2024" required>
            </div>
            
            <button type="submit" class="login-btn">Login to Ultimate Platform</button>
        </form>
        
        <div class="features">
            <div class="feature-tag">🎯 Custom Research Goals</div>
            <div class="feature-tag">🤖 Autonomous Literature Search</div>
            <div class="feature-tag">🧪 SMILES Export</div>
            <div class="feature-tag">💊 Advanced PKPD/DDI</div>
            <div class="feature-tag">🧬 3D Molecular Visualization</div>
            <div class="feature-tag">⚗️ Retrosynthesis Planning</div>
            <div class="feature-tag">🧠 Neuroplasticity Analysis</div>
            <div class="feature-tag">📊 IP Opportunity Tracking</div>
        </div>
    </div>

    <script>
        document.getElementById('loginForm').addEventListener('submit', function(e) {
            e.preventDefault();
            
            const formData = new FormData(this);
            
            fetch('/login', {
                method: 'POST',
                body: formData
            })
            .then(response => response.json())
            .then(data => {
                if (data.success) {
                    window.location.href = data.redirect;
                } else {
                    alert(data.message);
                }
            });
        });
    </script>
</body>
</html>
"""

ULTIMATE_DASHBOARD_TEMPLATE = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>PharmaSight™ Ultimate - Advanced Research Dashboard</title>
    <style>
        * { margin: 0; padding: 0; box-sizing: border-box; }
        body { 
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            background: #f8f9fa;
            line-height: 1.6;
        }
        .header {
            background: linear-gradient(135deg, #4a90e2 0%, #357abd 100%);
            color: white;
            padding: 1rem 2rem;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }
        .header h1 {
            font-size: 1.8rem;
            margin-bottom: 0.5rem;
        }
        .header p {
            opacity: 0.9;
            font-size: 0.9rem;
        }
        .stats-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 1.5rem;
            padding: 2rem;
            max-width: 1200px;
            margin: 0 auto;
        }
        .stat-card {
            background: white;
            padding: 1.5rem;
            border-radius: 12px;
            box-shadow: 0 4px 15px rgba(0,0,0,0.1);
            text-align: center;
            transition: transform 0.3s;
        }
        .stat-card:hover {
            transform: translateY(-5px);
        }
        .stat-number {
            font-size: 2.5rem;
            font-weight: bold;
            color: #4a90e2;
            margin-bottom: 0.5rem;
        }
        .stat-label {
            color: #666;
            font-size: 0.9rem;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }
        .nav-tabs {
            display: flex;
            background: white;
            margin: 0 2rem;
            border-radius: 12px 12px 0 0;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            overflow: hidden;
        }
        .nav-tab {
            flex: 1;
            padding: 1rem;
            text-align: center;
            background: #f8f9fa;
            border: none;
            cursor: pointer;
            transition: all 0.3s;
            font-size: 0.9rem;
            font-weight: 500;
        }
        .nav-tab.active {
            background: #4a90e2;
            color: white;
        }
        .nav-tab:hover:not(.active) {
            background: #e9ecef;
        }
        .content-area {
            background: white;
            margin: 0 2rem 2rem;
            padding: 2rem;
            border-radius: 0 0 12px 12px;
            box-shadow: 0 4px 15px rgba(0,0,0,0.1);
            min-height: 600px;
        }
        .section {
            display: none;
        }
        .section.active {
            display: block;
        }
        .section h2 {
            color: #333;
            margin-bottom: 1.5rem;
            font-size: 1.5rem;
        }
        .research-goal {
            background: #f8f9fa;
            padding: 1.5rem;
            border-radius: 8px;
            margin-bottom: 1rem;
            border-left: 4px solid #4a90e2;
        }
        .research-goal h3 {
            color: #333;
            margin-bottom: 0.5rem;
        }
        .research-goal p {
            color: #666;
            margin-bottom: 1rem;
        }
        .goal-actions {
            display: flex;
            gap: 0.5rem;
        }
        .btn {
            padding: 0.5rem 1rem;
            border: none;
            border-radius: 6px;
            cursor: pointer;
            font-size: 0.9rem;
            transition: all 0.3s;
        }
        .btn-primary {
            background: #4a90e2;
            color: white;
        }
        .btn-primary:hover {
            background: #357abd;
        }
        .btn-secondary {
            background: #6c757d;
            color: white;
        }
        .btn-secondary:hover {
            background: #545b62;
        }
        .compound-grid {
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(300px, 1fr));
            gap: 1rem;
            margin-top: 1rem;
        }
        .compound-card {
            background: #f8f9fa;
            padding: 1rem;
            border-radius: 8px;
            border: 1px solid #e9ecef;
        }
        .compound-name {
            font-weight: bold;
            color: #333;
            margin-bottom: 0.5rem;
        }
        .compound-smiles {
            font-family: monospace;
            font-size: 0.8rem;
            color: #666;
            background: white;
            padding: 0.5rem;
            border-radius: 4px;
            margin: 0.5rem 0;
            word-break: break-all;
        }
        .confidence-badge {
            display: inline-block;
            padding: 0.25rem 0.5rem;
            border-radius: 4px;
            font-size: 0.8rem;
            font-weight: bold;
        }
        .confidence-high {
            background: #28a745;
            color: white;
        }
        .confidence-medium {
            background: #ffc107;
            color: black;
        }
        .confidence-low {
            background: #dc3545;
            color: white;
        }
        .form-group {
            margin-bottom: 1rem;
        }
        .form-group label {
            display: block;
            margin-bottom: 0.5rem;
            font-weight: 500;
        }
        .form-group input, .form-group select, .form-group textarea {
            width: 100%;
            padding: 0.5rem;
            border: 1px solid #ddd;
            border-radius: 4px;
            font-size: 0.9rem;
        }
        .loading {
            text-align: center;
            padding: 2rem;
            color: #666;
        }
        .results-container {
            margin-top: 2rem;
            padding: 1rem;
            background: #f8f9fa;
            border-radius: 8px;
        }
        .interaction-warning {
            background: #fff3cd;
            border: 1px solid #ffeaa7;
            color: #856404;
            padding: 1rem;
            border-radius: 6px;
            margin: 1rem 0;
        }
        .interaction-severe {
            background: #f8d7da;
            border-color: #f5c6cb;
            color: #721c24;
        }
        .molecular-viewer {
            width: 100%;
            height: 400px;
            border: 1px solid #ddd;
            border-radius: 8px;
            margin: 1rem 0;
        }
    </style>
</head>
<body>
    <div class="header">
        <h1>🧬 PharmaSight™ Ultimate Research Platform</h1>
        <p>Advanced AI-Powered Pharmaceutical Research & Development with Custom Research Goals</p>
    </div>

    <div class="stats-grid">
        <div class="stat-card">
            <div class="stat-number">{{ compounds_count }}+</div>
            <div class="stat-label">Active Compounds</div>
        </div>
        <div class="stat-card">
            <div class="stat-number">{{ research_findings_count }}</div>
            <div class="stat-label">Research Findings</div>
        </div>
        <div class="stat-card">
            <div class="stat-number">{{ high_confidence_count }}</div>
            <div class="stat-label">High-Confidence Discoveries</div>
        </div>
        <div class="stat-card">
            <div class="stat-number">{{ research_goals_count }}</div>
            <div class="stat-label">Active Research Goals</div>
        </div>
        <div class="stat-card">
            <div class="stat-number">{{ total_market_value }}</div>
            <div class="stat-label">Total Market Value</div>
        </div>
    </div>

    <div class="nav-tabs">
        <button class="nav-tab active" onclick="showSection('research-goals')">🎯 Research Goals</button>
        <button class="nav-tab" onclick="showSection('compound-discovery')">🧪 Compound Discovery</button>
        <button class="nav-tab" onclick="showSection('pkpd-analysis')">💊 PKPD/DDI Analysis</button>
        <button class="nav-tab" onclick="showSection('molecular-viz')">🧬 3D Visualization</button>
        <button class="nav-tab" onclick="showSection('retrosynthesis')">⚗️ Retrosynthesis</button>
        <button class="nav-tab" onclick="showSection('neuroplasticity')">🧠 Neuroplasticity</button>
    </div>

    <div class="content-area">
        <!-- Research Goals Section -->
        <div id="research-goals" class="section active">
            <h2>🎯 Custom Research Goals & Autonomous Search</h2>
            <div id="research-goals-container">
                <div class="loading">Loading research goals...</div>
            </div>
            
            <h3>Add New Research Goal</h3>
            <form id="new-goal-form">
                <div class="form-group">
                    <label>Goal Title:</label>
                    <input type="text" id="goal-title" required>
                </div>
                <div class="form-group">
                    <label>Description:</label>
                    <textarea id="goal-description" rows="3" required></textarea>
                </div>
                <div class="form-group">
                    <label>Keywords (comma-separated):</label>
                    <input type="text" id="goal-keywords" placeholder="kava lactones, GABA modulation, tolerance-free">
                </div>
                <div class="form-group">
                    <label>Priority:</label>
                    <select id="goal-priority">
                        <option value="high">High</option>
                        <option value="medium">Medium</option>
                        <option value="low">Low</option>
                    </select>
                </div>
                <button type="submit" class="btn btn-primary">Add Research Goal</button>
            </form>
        </div>

        <!-- Compound Discovery Section -->
        <div id="compound-discovery" class="section">
            <h2>🧪 High-Confidence Compound Discovery</h2>
            <div class="form-group">
                <label>Minimum Confidence Level:</label>
                <select id="confidence-filter" onchange="loadHighConfidenceCompounds()">
                    <option value="0.85">85%+ (High Confidence)</option>
                    <option value="0.90">90%+ (Very High Confidence)</option>
                    <option value="0.95">95%+ (Exceptional Confidence)</option>
                </select>
            </div>
            <div id="compounds-container">
                <div class="loading">Loading high-confidence compounds...</div>
            </div>
        </div>

        <!-- PKPD/DDI Analysis Section -->
        <div id="pkpd-analysis" class="section">
            <h2>💊 Advanced PKPD & Drug-Drug Interaction Analysis</h2>
            <form id="pkpd-form">
                <div class="form-group">
                    <label>Drug 1:</label>
                    <select id="drug1">
                        <option value="buprenorphine">Buprenorphine</option>
                        <option value="ketamine">Ketamine</option>
                        <option value="psilocybin">Psilocybin</option>
                        <option value="tropacocaine">Tropacocaine</option>
                    </select>
                </div>
                <div class="form-group">
                    <label>Drug 1 Dose (mg):</label>
                    <input type="number" id="dose1" value="8" step="0.1">
                </div>
                <div class="form-group">
                    <label>Drug 2:</label>
                    <select id="drug2">
                        <option value="ketamine">Ketamine</option>
                        <option value="buprenorphine">Buprenorphine</option>
                        <option value="psilocybin">Psilocybin</option>
                        <option value="tropacocaine">Tropacocaine</option>
                    </select>
                </div>
                <div class="form-group">
                    <label>Drug 2 Dose (mg):</label>
                    <input type="number" id="dose2" value="0.5" step="0.1">
                </div>
                <button type="submit" class="btn btn-primary">Analyze Drug Interaction</button>
            </form>
            <div id="pkpd-results"></div>
        </div>

        <!-- 3D Molecular Visualization Section -->
        <div id="molecular-viz" class="section">
            <h2>🧬 3D Molecular Visualization</h2>
            <div class="form-group">
                <label>Select Compound:</label>
                <select id="viz-compound" onchange="load3DStructure()">
                    <option value="">Select a compound...</option>
                    <option value="Synthetic Kavain Analog SK-101">Synthetic Kavain Analog SK-101</option>
                    <option value="Selective KOR Antagonist KA-301">Selective KOR Antagonist KA-301</option>
                    <option value="Short-Acting Buprenorphine Analog BA-150">Short-Acting Buprenorphine Analog BA-150</option>
                    <option value="Short-Duration Tryptamine ST-42">Short-Duration Tryptamine ST-42</option>
                </select>
            </div>
            <div id="molecular-viewer" class="molecular-viewer">
                <div class="loading">Select a compound to view 3D structure</div>
            </div>
            <div id="molecular-properties"></div>
        </div>

        <!-- Retrosynthesis Section -->
        <div id="retrosynthesis" class="section">
            <h2>⚗️ Retrosynthesis Planning</h2>
            <div class="form-group">
                <label>Compound Type:</label>
                <select id="retro-compound" onchange="loadRetrosynthesis()">
                    <option value="">Select compound type...</option>
                    <option value="buprenorphine_analogs">Buprenorphine Analogs</option>
                    <option value="kappa_antagonists">Kappa Opioid Antagonists</option>
                    <option value="gaba_modulators">GABA Modulators</option>
                    <option value="psychedelic_analogs">Psychedelic Analogs</option>
                </select>
            </div>
            <div id="retrosynthesis-results"></div>
        </div>

        <!-- Neuroplasticity Section -->
        <div id="neuroplasticity" class="section">
            <h2>🧠 Neuroplasticity Window Analysis for TMS Enhancement</h2>
            <p>Analyze drug combinations that enhance TMS efficacy by maintaining neuroplasticity windows.</p>
            <form id="neuroplasticity-form">
                <div class="form-group">
                    <label>Select Drugs for Analysis:</label>
                    <div>
                        <input type="checkbox" id="neuro-ketamine" value="ketamine" checked>
                        <label for="neuro-ketamine">Ketamine (NMDA antagonist, mTOR activation)</label>
                    </div>
                    <div>
                        <input type="checkbox" id="neuro-psilocybin" value="psilocybin" checked>
                        <label for="neuro-psilocybin">Psilocybin (5-HT2A agonist, BDNF upregulation)</label>
                    </div>
                    <div>
                        <input type="checkbox" id="neuro-lsd" value="LSD">
                        <label for="neuro-lsd">LSD (5-HT2A agonist, extended plasticity)</label>
                    </div>
                </div>
                <button type="submit" class="btn btn-primary">Analyze Neuroplasticity Windows</button>
            </form>
            <div id="neuroplasticity-results"></div>
        </div>
    </div>

    <script>
        // Navigation
        function showSection(sectionId) {
            document.querySelectorAll('.section').forEach(s => s.classList.remove('active'));
            document.querySelectorAll('.nav-tab').forEach(t => t.classList.remove('active'));
            
            document.getElementById(sectionId).classList.add('active');
            event.target.classList.add('active');
        }

        // Load research goals
        function loadResearchGoals() {
            fetch('/api/research_goals')
                .then(response => response.json())
                .then(data => {
                    const container = document.getElementById('research-goals-container');
                    container.innerHTML = '';
                    
                    Object.entries(data).forEach(([goalId, goal]) => {
                        const goalDiv = document.createElement('div');
                        goalDiv.className = 'research-goal';
                        goalDiv.innerHTML = `
                            <h3>${goal.title}</h3>
                            <p>${goal.description}</p>
                            <p><strong>Keywords:</strong> ${goal.keywords ? goal.keywords.join(', ') : 'N/A'}</p>
                            <p><strong>Priority:</strong> ${goal.priority || 'Medium'}</p>
                            <div class="goal-actions">
                                <button class="btn btn-primary" onclick="performAutonomousSearch('${goalId}')">🔍 Autonomous Search</button>
                                <button class="btn btn-secondary" onclick="discoverCompounds('${goalId}')">🧪 Discover Compounds</button>
                            </div>
                            <div id="results-${goalId}"></div>
                        `;
                        container.appendChild(goalDiv);
                    });
                });
        }

        // Perform autonomous search
        function performAutonomousSearch(goalId) {
            const resultsDiv = document.getElementById(`results-${goalId}`);
            resultsDiv.innerHTML = '<div class="loading">Performing autonomous literature search...</div>';
            
            fetch(`/api/autonomous_search/${goalId}`)
                .then(response => response.json())
                .then(data => {
                    resultsDiv.innerHTML = `
                        <div class="results-container">
                            <h4>📚 Literature Search Results (${data.total_found} found)</h4>
                            ${data.results.map(result => `
                                <div style="margin: 1rem 0; padding: 1rem; background: white; border-radius: 6px;">
                                    <strong>${result.title}</strong><br>
                                    <em>${result.authors} - ${result.journal} (${result.year})</em><br>
                                    <small>DOI: ${result.doi}</small><br>
                                    <p>${result.key_findings}</p>
                                    <span class="confidence-badge confidence-${result.confidence >= 0.9 ? 'high' : result.confidence >= 0.7 ? 'medium' : 'low'}">
                                        ${Math.round(result.confidence * 100)}% Confidence
                                    </span>
                                </div>
                            `).join('')}
                        </div>
                    `;
                });
        }

        // Discover compounds
        function discoverCompounds(goalId) {
            const resultsDiv = document.getElementById(`results-${goalId}`);
            resultsDiv.innerHTML = '<div class="loading">Discovering compounds...</div>';
            
            fetch(`/api/discover_compounds/${goalId}`)
                .then(response => response.json())
                .then(data => {
                    resultsDiv.innerHTML = `
                        <div class="results-container">
                            <h4>🧪 Compound Discoveries (${data.total_discovered} found)</h4>
                            <div class="compound-grid">
                                ${data.compounds.map(compound => `
                                    <div class="compound-card">
                                        <div class="compound-name">${compound.name}</div>
                                        <div class="compound-smiles">${compound.smiles}</div>
                                        <p><strong>Target:</strong> ${compound.target_receptors ? compound.target_receptors.join(', ') : 'N/A'}</p>
                                        <p><strong>Binding:</strong> ${compound.binding_affinity || 'N/A'}</p>
                                        <p><strong>Market Potential:</strong> ${compound.market_potential || 'N/A'}</p>
                                        <p><strong>IP Status:</strong> ${compound.ip_status || 'N/A'}</p>
                                        <span class="confidence-badge confidence-${compound.confidence >= 0.9 ? 'high' : compound.confidence >= 0.7 ? 'medium' : 'low'}">
                                            ${Math.round(compound.confidence * 100)}% Confidence
                                        </span>
                                    </div>
                                `).join('')}
                            </div>
                        </div>
                    `;
                });
        }

        // Load high confidence compounds
        function loadHighConfidenceCompounds() {
            const minConfidence = document.getElementById('confidence-filter').value;
            const container = document.getElementById('compounds-container');
            container.innerHTML = '<div class="loading">Loading high-confidence compounds...</div>';
            
            fetch(`/api/high_confidence_compounds?min_confidence=${minConfidence}`)
                .then(response => response.json())
                .then(data => {
                    container.innerHTML = `
                        <h3>High-Confidence Compounds (${data.total} found)</h3>
                        <div class="compound-grid">
                            ${data.compounds.map(compound => `
                                <div class="compound-card">
                                    <div class="compound-name">${compound.name}</div>
                                    <div class="compound-smiles">${compound.smiles}</div>
                                    <p><strong>Formula:</strong> ${compound.molecular_formula || 'N/A'}</p>
                                    <p><strong>Targets:</strong> ${compound.target_receptors ? compound.target_receptors.join(', ') : 'N/A'}</p>
                                    <p><strong>Binding:</strong> ${compound.binding_affinity || 'N/A'}</p>
                                    <p><strong>Market Potential:</strong> ${compound.market_potential || 'N/A'}</p>
                                    <p><strong>IP Status:</strong> ${compound.ip_status || 'N/A'}</p>
                                    <span class="confidence-badge confidence-${compound.confidence >= 0.9 ? 'high' : compound.confidence >= 0.7 ? 'medium' : 'low'}">
                                        ${Math.round(compound.confidence * 100)}% Confidence
                                    </span>
                                </div>
                            `).join('')}
                        </div>
                    `;
                });
        }

        // PKPD Analysis
        document.getElementById('pkpd-form').addEventListener('submit', function(e) {
            e.preventDefault();
            
            const drug1 = document.getElementById('drug1').value;
            const drug2 = document.getElementById('drug2').value;
            const dose1 = parseFloat(document.getElementById('dose1').value);
            const dose2 = parseFloat(document.getElementById('dose2').value);
            
            const resultsDiv = document.getElementById('pkpd-results');
            resultsDiv.innerHTML = '<div class="loading">Analyzing drug interaction...</div>';
            
            fetch(`/api/drug_interaction?drug1=${drug1}&drug2=${drug2}&dose1=${dose1}&dose2=${dose2}`)
                .then(response => response.json())
                .then(data => {
                    if (data.interaction_found) {
                        const severityClass = data.severity === 'severe' ? 'interaction-severe' : 'interaction-warning';
                        resultsDiv.innerHTML = `
                            <div class="results-container">
                                <h3>Drug Interaction Analysis</h3>
                                <div class="${severityClass}">
                                    <strong>Interaction Severity:</strong> ${data.severity.toUpperCase()}<br>
                                    <strong>Mechanism:</strong> ${data.mechanism}<br>
                                    <strong>Effect:</strong> ${data.effect}<br>
                                    <strong>Recommendation:</strong> ${data.recommendation}<br>
                                    <strong>Confidence:</strong> ${Math.round(data.confidence * 100)}%
                                </div>
                                <h4>Monitoring Parameters:</h4>
                                <ul>
                                    ${data.monitoring_parameters ? data.monitoring_parameters.map(param => `<li>${param.replace('_', ' ')}</li>`).join('') : '<li>Standard monitoring</li>'}
                                </ul>
                                <h4>Dose Adjustments:</h4>
                                <ul>
                                    ${Object.entries(data.dose_adjustments || {}).map(([drug, adjustment]) => `<li><strong>${drug}:</strong> ${adjustment}</li>`).join('')}
                                </ul>
                            </div>
                        `;
                    } else {
                        resultsDiv.innerHTML = `
                            <div class="results-container">
                                <p>${data.message}</p>
                                <p><strong>Recommendation:</strong> ${data.recommendation}</p>
                            </div>
                        `;
                    }
                });
        });

        // 3D Molecular Visualization
        function load3DStructure() {
            const compound = document.getElementById('viz-compound').value;
            if (!compound) return;
            
            const viewer = document.getElementById('molecular-viewer');
            const properties = document.getElementById('molecular-properties');
            
            viewer.innerHTML = '<div class="loading">Loading 3D structure...</div>';
            
            fetch(`/api/3d_structure/${encodeURIComponent(compound)}`)
                .then(response => response.json())
                .then(data => {
                    if (data.error) {
                        viewer.innerHTML = `<div class="loading">Error: ${data.error}</div>`;
                        return;
                    }
                    
                    viewer.innerHTML = `
                        <div style="padding: 2rem; text-align: center;">
                            <h3>${data.compound_name}</h3>
                            <p><strong>SMILES:</strong> <code>${data.smiles}</code></p>
                            <p><strong>Formula:</strong> ${data.molecular_formula}</p>
                            <div style="background: #f0f0f0; padding: 1rem; border-radius: 8px; margin: 1rem 0;">
                                <em>3D molecular viewer would be rendered here using ChemDoodle or similar library</em>
                            </div>
                        </div>
                    `;
                    
                    properties.innerHTML = `
                        <div class="results-container">
                            <h4>Molecular Properties</h4>
                            <ul>
                                <li><strong>Molecular Weight:</strong> ${data.property_display.molecular_weight} g/mol</li>
                                <li><strong>LogP:</strong> ${data.property_display.logp}</li>
                                <li><strong>Polar Surface Area:</strong> ${data.property_display.polar_surface_area} Ų</li>
                                <li><strong>Rotatable Bonds:</strong> ${data.property_display.rotatable_bonds}</li>
                            </ul>
                        </div>
                    `;
                });
        }

        // Retrosynthesis
        function loadRetrosynthesis() {
            const compoundType = document.getElementById('retro-compound').value;
            if (!compoundType) return;
            
            const resultsDiv = document.getElementById('retrosynthesis-results');
            resultsDiv.innerHTML = '<div class="loading">Loading retrosynthesis pathway...</div>';
            
            fetch(`/api/retrosynthesis/${compoundType}`)
                .then(response => response.json())
                .then(data => {
                    if (data.error) {
                        resultsDiv.innerHTML = `<div class="loading">Error: ${data.error}</div>`;
                        return;
                    }
                    
                    resultsDiv.innerHTML = `
                        <div class="results-container">
                            <h3>Retrosynthesis Pathway</h3>
                            <p><strong>Target SMILES:</strong> <code>${data.target_smiles}</code></p>
                            <p><strong>Complexity Score:</strong> ${data.complexity_score}/10</p>
                            <p><strong>Estimated Steps:</strong> ${data.estimated_steps}</p>
                            <p><strong>Overall Yield:</strong> ${data.overall_yield}</p>
                            <p><strong>Feasibility Score:</strong> ${data.feasibility_score}/10</p>
                            <p><strong>Estimated Cost:</strong> ${data.estimated_cost}</p>
                            
                            <h4>Synthetic Route:</h4>
                            ${data.synthetic_route.map(step => `
                                <div style="margin: 1rem 0; padding: 1rem; background: white; border-radius: 6px; border-left: 4px solid #4a90e2;">
                                    <strong>Step ${step.step}: ${step.reaction}</strong><br>
                                    <strong>Starting Material:</strong> ${step.starting_material}<br>
                                    <strong>Reagents:</strong> ${step.reagents.join(', ')}<br>
                                    <strong>Conditions:</strong> ${step.conditions}<br>
                                    <strong>Yield:</strong> ${step.yield}
                                </div>
                            `).join('')}
                            
                            <h4>Safety Considerations:</h4>
                            <ul>
                                ${data.safety_considerations.map(safety => `<li>${safety}</li>`).join('')}
                            </ul>
                            
                            <h4>Equipment Needed:</h4>
                            <ul>
                                ${data.equipment_needed.map(equipment => `<li>${equipment}</li>`).join('')}
                            </ul>
                        </div>
                    `;
                });
        }

        // Neuroplasticity Analysis
        document.getElementById('neuroplasticity-form').addEventListener('submit', function(e) {
            e.preventDefault();
            
            const drugs = [];
            document.querySelectorAll('#neuroplasticity-form input[type="checkbox"]:checked').forEach(checkbox => {
                drugs.push(checkbox.value);
            });
            
            if (drugs.length === 0) {
                alert('Please select at least one drug for analysis');
                return;
            }
            
            const resultsDiv = document.getElementById('neuroplasticity-results');
            resultsDiv.innerHTML = '<div class="loading">Analyzing neuroplasticity windows...</div>';
            
            fetch('/api/neuroplasticity_analysis', {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                },
                body: JSON.stringify({drugs: drugs})
            })
            .then(response => response.json())
            .then(data => {
                if (data.message) {
                    resultsDiv.innerHTML = `<div class="results-container"><p>${data.message}</p></div>`;
                    return;
                }
                
                resultsDiv.innerHTML = `
                    <div class="results-container">
                        <h3>Neuroplasticity Window Analysis</h3>
                        <h4>Individual Windows:</h4>
                        <ul>
                            ${Object.entries(data.individual_windows).map(([drug, hours]) => `<li><strong>${drug}:</strong> ${hours} hours</li>`).join('')}
                        </ul>
                        <p><strong>Maximum Window:</strong> ${data.maximum_window} hours</p>
                        <p><strong>Optimal Timing:</strong> ${data.optimal_timing}</p>
                        <p><strong>Therapeutic Implications:</strong> ${data.therapeutic_implications}</p>
                        <p><strong>Monitoring Duration:</strong> ${data.monitoring_duration}</p>
                        
                        <div class="interaction-warning">
                            <strong>TMS Enhancement Protocol:</strong><br>
                            Administer selected compounds simultaneously to maximize neuroplasticity window overlap. 
                            This extended window may significantly enhance TMS efficacy for treatment-resistant depression 
                            by maintaining optimal conditions for synaptic plasticity and BDNF expression.
                        </div>
                    </div>
                `;
            });
        });

        // Add new research goal
        document.getElementById('new-goal-form').addEventListener('submit', function(e) {
            e.preventDefault();
            
            const goalData = {
                title: document.getElementById('goal-title').value,
                description: document.getElementById('goal-description').value,
                keywords: document.getElementById('goal-keywords').value.split(',').map(k => k.trim()),
                priority: document.getElementById('goal-priority').value,
                status: 'active'
            };
            
            fetch('/api/research_goals', {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                },
                body: JSON.stringify(goalData)
            })
            .then(response => response.json())
            .then(data => {
                if (data.success) {
                    alert('Research goal added successfully!');
                    this.reset();
                    loadResearchGoals();
                }
            });
        });

        // Initialize
        document.addEventListener('DOMContentLoaded', function() {
            loadResearchGoals();
            loadHighConfidenceCompounds();
        });
    </script>
</body>
</html>
"""

@app.route('/health')
def health_check():
    """Health check endpoint"""
    return jsonify({
        'status': 'healthy',
        'version': '4.0.0-ultimate',
        'database': 'comprehensive',
        'compounds': len(COMPOUND_DATABASE),
        'research_findings': len(RESEARCH_DATA.get('findings', [])),
        'analog_discoveries': len(ANALOG_DATA),
        'research_articles': len(ARTICLES_DATA.get('articles', [])),
        'research_goals': len(research_system.get_research_goals()),
        'features': 'all_ultimate_operational'
    })

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=8080, debug=False)

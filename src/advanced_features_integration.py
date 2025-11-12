#!/usr/bin/env python3
"""
Advanced Features Integration Module
Integrates all new drug discovery features into the main PharmaSight platform
"""

from flask import Blueprint, request, jsonify, render_template_string
import json
import traceback

# Import all new modules
from virtual_screening_pipeline import VirtualScreeningPipeline
from ai_lead_optimization import AILeadOptimizer
from sar_explorer import SARExplorer
from pharmacophore_modeling import PharmacophoreModeler
from off_target_prediction import OffTargetPredictor
from quantum_computing_module import create_quantum_routes

# Create Flask blueprint for advanced features
advanced_bp = Blueprint('advanced_features', __name__)

# Initialize modules
vhts_pipeline = VirtualScreeningPipeline()
lead_optimizer = AILeadOptimizer()
sar_explorer = SARExplorer()
pharmacophore_modeler = PharmacophoreModeler()
off_target_predictor = OffTargetPredictor()

# Virtual Screening Routes
@advanced_bp.route('/api/vhts/screen_compound', methods=['POST'])
def vhts_screen_compound():
    """Screen a single compound against all 82 receptors"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')
        compound_name = data.get('name')
        
        if not smiles:
            return jsonify({'error': 'SMILES string required'}), 400
        
        result = vhts_pipeline.screen_compound(smiles, compound_name)
        return jsonify(result)
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@advanced_bp.route('/api/vhts/batch_screen', methods=['POST'])
def vhts_batch_screen():
    """Screen multiple compounds in parallel"""
    try:
        data = request.get_json()
        compounds = data.get('compounds', [])
        
        if not compounds:
            return jsonify({'error': 'Compounds list required'}), 400
        
        # Convert to required format
        compound_list = [(c['smiles'], c.get('name', f'CPD_{i}')) 
                        for i, c in enumerate(compounds)]
        
        results = vhts_pipeline.batch_screen(compound_list, parallel=True)
        report = vhts_pipeline.generate_screening_report(results)
        
        return jsonify({
            'results': results[:10],  # First 10 results
            'report': report
        })
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@advanced_bp.route('/api/vhts/find_selective', methods=['POST'])
def vhts_find_selective():
    """Find compounds selective for a specific receptor"""
    try:
        data = request.get_json()
        target = data.get('target_receptor')
        avoid = data.get('avoid_receptors', [])
        
        if not target:
            return jsonify({'error': 'Target receptor required'}), 400
        
        results = vhts_pipeline.find_selective_ligands(target, avoid)
        return jsonify({'selective_compounds': results})
    except Exception as e:
        return jsonify({'error': str(e)}), 500

# AI Lead Optimization Routes
@advanced_bp.route('/api/lead_opt/optimize', methods=['POST'])
def lead_optimize():
    """Generate optimization suggestions for a lead compound"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')
        target_profile = data.get('target_profile', {})
        
        if not smiles:
            return jsonify({'error': 'SMILES string required'}), 400
        
        suggestions = lead_optimizer.optimize_lead(smiles, target_profile)
        return jsonify(suggestions)
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@advanced_bp.route('/api/lead_opt/generate_analogs', methods=['POST'])
def generate_analogs():
    """Generate optimized analog series"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')
        num_analogs = data.get('num_analogs', 10)
        
        if not smiles:
            return jsonify({'error': 'SMILES string required'}), 400
        
        analogs = lead_optimizer.generate_analog_series(smiles, num_analogs)
        return jsonify({'analogs': analogs})
    except Exception as e:
        return jsonify({'error': str(e)}), 500

# SAR Explorer Routes
@advanced_bp.route('/api/sar/analyze', methods=['POST'])
def sar_analyze():
    """Analyze SAR for a compound series"""
    try:
        data = request.get_json()
        compounds = data.get('compounds', [])
        
        if not compounds:
            return jsonify({'error': 'Compound series required'}), 400
        
        analysis = sar_explorer.analyze_sar(compounds)
        return jsonify(analysis)
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@advanced_bp.route('/api/sar/predict_activity', methods=['POST'])
def sar_predict():
    """Predict activity based on SAR model"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')
        
        if not smiles:
            return jsonify({'error': 'SMILES string required'}), 400
        
        prediction = sar_explorer.predict_activity(smiles)
        return jsonify(prediction)
    except Exception as e:
        return jsonify({'error': str(e)}), 500

# Pharmacophore Modeling Routes
@advanced_bp.route('/api/pharmacophore/extract', methods=['POST'])
def extract_pharmacophore():
    """Extract pharmacophore from active compounds"""
    try:
        data = request.get_json()
        compounds = data.get('active_compounds', [])
        
        if not compounds:
            return jsonify({'error': 'Active compounds required'}), 400
        
        pharmacophore = pharmacophore_modeler.extract_pharmacophore(compounds)
        return jsonify(pharmacophore)
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@advanced_bp.route('/api/pharmacophore/screen', methods=['POST'])
def screen_pharmacophore():
    """Screen compound against pharmacophore"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')
        pharmacophore = data.get('pharmacophore')
        
        if not smiles or not pharmacophore:
            return jsonify({'error': 'SMILES and pharmacophore required'}), 400
        
        result = pharmacophore_modeler.screen_compound(smiles, pharmacophore)
        return jsonify(result)
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@advanced_bp.route('/api/pharmacophore/get_model/<receptor>')
def get_pharmacophore_model(receptor):
    """Get predefined pharmacophore for a receptor"""
    try:
        model = pharmacophore_modeler.generate_pharmacophore_query(receptor)
        return jsonify(model)
    except Exception as e:
        return jsonify({'error': str(e)}), 500

# Off-Target Prediction Routes
@advanced_bp.route('/api/off_target/predict', methods=['POST'])
def predict_off_targets():
    """Predict off-target interactions and side effects"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')
        primary_target = data.get('primary_target')
        
        if not smiles:
            return jsonify({'error': 'SMILES string required'}), 400
        
        prediction = off_target_predictor.predict_off_targets(smiles, primary_target)
        report = off_target_predictor.generate_safety_report(prediction)
        
        return jsonify({
            'prediction': prediction,
            'safety_report': report
        })
    except Exception as e:
        return jsonify({'error': str(e)}), 500

# Combined Analysis Route
@advanced_bp.route('/api/advanced/full_analysis', methods=['POST'])
def full_compound_analysis():
    """Perform comprehensive analysis using all modules"""
    try:
        data = request.get_json()
        smiles = data.get('smiles')
        compound_name = data.get('name', 'Compound')
        
        if not smiles:
            return jsonify({'error': 'SMILES string required'}), 400
        
        # Run all analyses
        results = {
            'compound': compound_name,
            'smiles': smiles,
            'analyses': {}
        }
        
        # Virtual screening
        try:
            results['analyses']['receptor_screening'] = vhts_pipeline.screen_compound(
                smiles, compound_name
            )
        except:
            results['analyses']['receptor_screening'] = {'error': 'Failed to screen'}
        
        # Lead optimization
        try:
            results['analyses']['optimization_suggestions'] = lead_optimizer.optimize_lead(smiles)
        except:
            results['analyses']['optimization_suggestions'] = {'error': 'Failed to optimize'}
        
        # Off-target prediction
        try:
            results['analyses']['safety_assessment'] = off_target_predictor.predict_off_targets(smiles)
        except:
            results['analyses']['safety_assessment'] = {'error': 'Failed to predict'}
        
        return jsonify(results)
    except Exception as e:
        return jsonify({'error': str(e)}), 500

# Advanced Features UI Route
@advanced_bp.route('/advanced_features')
def advanced_features_ui():
    """Render advanced features interface"""
    return render_template_string(ADVANCED_FEATURES_HTML)

# HTML Template for Advanced Features
ADVANCED_FEATURES_HTML = '''
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>PharmaSight™ - Advanced Drug Discovery Features</title>
    <style>
        * { margin: 0; padding: 0; box-sizing: border-box; }
        
        body {
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            min-height: 100vh;
            color: #333;
        }
        
        .header {
            background: rgba(255, 255, 255, 0.95);
            padding: 20px;
            box-shadow: 0 2px 10px rgba(0, 0, 0, 0.1);
        }
        
        .header h1 {
            color: #764ba2;
            font-size: 32px;
            display: flex;
            align-items: center;
            gap: 15px;
        }
        
        .container {
            max-width: 1400px;
            margin: 30px auto;
            padding: 0 20px;
        }
        
        .feature-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(400px, 1fr));
            gap: 25px;
            margin-top: 30px;
        }
        
        .feature-card {
            background: white;
            border-radius: 15px;
            padding: 25px;
            box-shadow: 0 4px 15px rgba(0, 0, 0, 0.1);
            transition: transform 0.3s, box-shadow 0.3s;
        }
        
        .feature-card:hover {
            transform: translateY(-5px);
            box-shadow: 0 8px 25px rgba(0, 0, 0, 0.2);
        }
        
        .feature-card h3 {
            color: #667eea;
            font-size: 24px;
            margin-bottom: 15px;
            display: flex;
            align-items: center;
            gap: 10px;
        }
        
        .feature-card p {
            color: #666;
            line-height: 1.6;
            margin-bottom: 20px;
        }
        
        .feature-stats {
            display: grid;
            grid-template-columns: repeat(2, 1fr);
            gap: 10px;
            margin-top: 15px;
        }
        
        .stat {
            background: #f0f4ff;
            padding: 10px;
            border-radius: 8px;
            text-align: center;
        }
        
        .stat-value {
            font-size: 24px;
            font-weight: bold;
            color: #667eea;
        }
        
        .stat-label {
            font-size: 12px;
            color: #666;
            text-transform: uppercase;
        }
        
        .input-group {
            margin-bottom: 20px;
        }
        
        .input-group label {
            display: block;
            margin-bottom: 5px;
            color: #555;
            font-weight: 500;
        }
        
        .input-group input, .input-group textarea {
            width: 100%;
            padding: 10px;
            border: 2px solid #e0e0e0;
            border-radius: 8px;
            font-size: 14px;
            transition: border-color 0.3s;
        }
        
        .input-group input:focus, .input-group textarea:focus {
            outline: none;
            border-color: #667eea;
        }
        
        .btn {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            border: none;
            padding: 12px 30px;
            border-radius: 8px;
            font-size: 16px;
            cursor: pointer;
            transition: opacity 0.3s;
        }
        
        .btn:hover {
            opacity: 0.9;
        }
        
        .results {
            margin-top: 20px;
            padding: 15px;
            background: #f8f9fa;
            border-radius: 8px;
            max-height: 400px;
            overflow-y: auto;
        }
        
        .loading {
            display: none;
            text-align: center;
            padding: 20px;
        }
        
        .loading.active {
            display: block;
        }
        
        .spinner {
            border: 3px solid #f3f3f3;
            border-top: 3px solid #667eea;
            border-radius: 50%;
            width: 40px;
            height: 40px;
            animation: spin 1s linear infinite;
            margin: 0 auto;
        }
        
        @keyframes spin {
            0% { transform: rotate(0deg); }
            100% { transform: rotate(360deg); }
        }
        
        .tabs {
            display: flex;
            gap: 10px;
            margin-bottom: 20px;
            border-bottom: 2px solid #e0e0e0;
        }
        
        .tab {
            padding: 10px 20px;
            background: none;
            border: none;
            color: #666;
            cursor: pointer;
            font-size: 16px;
            transition: color 0.3s;
        }
        
        .tab.active {
            color: #667eea;
            border-bottom: 3px solid #667eea;
        }
        
        .tab-content {
            display: none;
        }
        
        .tab-content.active {
            display: block;
        }
    </style>
</head>
<body>
    <div class="header">
        <h1>🚀 PharmaSight™ Advanced Drug Discovery Suite</h1>
        <p style="color: #666; margin-top: 10px;">
            Powered by AI • 82 Receptor Targets • Real-time Analysis
        </p>
    </div>

    <div class="container">
        <!-- Statistics Overview -->
        <div class="feature-grid">
            <div class="feature-card">
                <h3>📊 Platform Statistics</h3>
                <div class="feature-stats">
                    <div class="stat">
                        <div class="stat-value">82</div>
                        <div class="stat-label">Receptor Targets</div>
                    </div>
                    <div class="stat">
                        <div class="stat-value">5</div>
                        <div class="stat-label">AI Modules</div>
                    </div>
                    <div class="stat">
                        <div class="stat-value">500+</div>
                        <div class="stat-label">Compounds</div>
                    </div>
                    <div class="stat">
                        <div class="stat-value">Real-time</div>
                        <div class="stat-label">Analysis</div>
                    </div>
                </div>
            </div>
        </div>

        <!-- Feature Cards -->
        <div class="feature-grid">
            <!-- Virtual Screening -->
            <div class="feature-card">
                <h3>🔬 Virtual High-Throughput Screening</h3>
                <p>Screen compounds against all 82 receptor targets simultaneously. Identify multi-target drugs and repurposing opportunities.</p>
                
                <div class="input-group">
                    <label>SMILES String</label>
                    <input type="text" id="vhts-smiles" placeholder="Enter SMILES (e.g., CC(C)NCC(O)COc1ccccc1)" value="CN1CCC[C@H]1c2cccnc2">
                </div>
                
                <button class="btn" onclick="runVirtualScreening()">Run Screening</button>
                
                <div id="vhts-loading" class="loading">
                    <div class="spinner"></div>
                    <p>Screening against 82 receptors...</p>
                </div>
                
                <div id="vhts-results" class="results" style="display: none;"></div>
            </div>

            <!-- Lead Optimization -->
            <div class="feature-card">
                <h3>🧬 AI-Powered Lead Optimization</h3>
                <p>Get ML-based suggestions to improve binding affinity, selectivity, and drug-likeness of your lead compounds.</p>
                
                <div class="input-group">
                    <label>Lead Compound SMILES</label>
                    <input type="text" id="lead-smiles" placeholder="Enter SMILES" value="Cc1ccc(cc1Nc2nccc(n2)c3cccnc3)NC(=O)c4ccc(cc4)CN5CCN(CC5)C">
                </div>
                
                <button class="btn" onclick="optimizeLead()">Optimize Lead</button>
                
                <div id="lead-loading" class="loading">
                    <div class="spinner"></div>
                    <p>Generating optimization strategies...</p>
                </div>
                
                <div id="lead-results" class="results" style="display: none;"></div>
            </div>

            <!-- SAR Explorer -->
            <div class="feature-card">
                <h3>📈 Structure-Activity Relationship Explorer</h3>
                <p>Visualize how molecular changes affect receptor binding and identify activity cliffs in your compound series.</p>
                
                <div class="input-group">
                    <label>Compound Series (JSON)</label>
                    <textarea id="sar-compounds" rows="4" placeholder='[{"smiles": "...", "activity": 0.8}, ...]'>[
  {"smiles": "CC(C)NCC(O)COc1ccccc1", "activity": 0.75, "id": "CPD1"},
  {"smiles": "CC(C)NCC(O)COc1ccccc1F", "activity": 0.85, "id": "CPD2"},
  {"smiles": "CC(C)NCC(O)COc1ccccc1Cl", "activity": 0.65, "id": "CPD3"}
]</textarea>
                </div>
                
                <button class="btn" onclick="analyzeSAR()">Analyze SAR</button>
                
                <div id="sar-loading" class="loading">
                    <div class="spinner"></div>
                    <p>Analyzing structure-activity relationships...</p>
                </div>
                
                <div id="sar-results" class="results" style="display: none;"></div>
            </div>

            <!-- Pharmacophore Modeling -->
            <div class="feature-card">
                <h3>🎯 Pharmacophore Modeling</h3>
                <p>Extract 3D pharmacophore models from active compounds and screen new molecules for key binding features.</p>
                
                <div class="tabs">
                    <button class="tab active" onclick="switchTab('pharm', 'extract')">Extract Model</button>
                    <button class="tab" onclick="switchTab('pharm', 'screen')">Screen Compound</button>
                </div>
                
                <div id="pharm-extract" class="tab-content active">
                    <div class="input-group">
                        <label>Active Compounds (JSON)</label>
                        <textarea id="pharm-actives" rows="3" placeholder='[{"smiles": "..."}, ...]'>[
  {"smiles": "CN1CCC[C@H]1c2cccnc2", "activity": 0.9},
  {"smiles": "CC(C)NCC(O)COc1ccccc1", "activity": 0.8}
]</textarea>
                    </div>
                    <button class="btn" onclick="extractPharmacophore()">Extract Pharmacophore</button>
                </div>
                
                <div id="pharm-screen" class="tab-content">
                    <div class="input-group">
                        <label>Receptor Target</label>
                        <select id="pharm-receptor">
                            <option value="5-HT2A">5-HT2A Serotonin</option>
                            <option value="D2">D2 Dopamine</option>
                            <option value="MOR">μ-Opioid (MOR)</option>
                            <option value="GABA-A">GABA-A</option>
                            <option value="Beta2">β2-Adrenergic</option>
                        </select>
                    </div>
                    <button class="btn" onclick="getPharmacophoreModel()">Get Model</button>
                </div>
                
                <div id="pharm-loading" class="loading">
                    <div class="spinner"></div>
                    <p>Processing pharmacophore...</p>
                </div>
                
                <div id="pharm-results" class="results" style="display: none;"></div>
            </div>

            <!-- Off-Target Prediction -->
            <div class="feature-card">
                <h3>⚠️ Off-Target & Safety Prediction</h3>
                <p>Predict potential side effects and safety liabilities by analyzing off-target receptor interactions.</p>
                
                <div class="input-group">
                    <label>Compound SMILES</label>
                    <input type="text" id="safety-smiles" placeholder="Enter SMILES" value="O=C(C)Oc1ccccc1C(=O)O">
                </div>
                
                <div class="input-group">
                    <label>Primary Target (Optional)</label>
                    <input type="text" id="safety-target" placeholder="e.g., COX-2" value="COX-2">
                </div>
                
                <button class="btn" onclick="predictSafety()">Predict Safety</button>
                
                <div id="safety-loading" class="loading">
                    <div class="spinner"></div>
                    <p>Analyzing safety profile...</p>
                </div>
                
                <div id="safety-results" class="results" style="display: none;"></div>
            </div>

            <!-- Full Analysis -->
            <div class="feature-card">
                <h3>🔍 Comprehensive Analysis</h3>
                <p>Run all analyses simultaneously: receptor screening, lead optimization, and safety assessment in one click.</p>
                
                <div class="input-group">
                    <label>Compound SMILES</label>
                    <input type="text" id="full-smiles" placeholder="Enter SMILES" value="CN1C=NC2=C1C(=O)N(C(=O)N2C)C">
                </div>
                
                <div class="input-group">
                    <label>Compound Name</label>
                    <input type="text" id="full-name" placeholder="Enter name" value="Caffeine">
                </div>
                
                <button class="btn" onclick="runFullAnalysis()">Run Full Analysis</button>
                
                <div id="full-loading" class="loading">
                    <div class="spinner"></div>
                    <p>Running comprehensive analysis...</p>
                </div>
                
                <div id="full-results" class="results" style="display: none;"></div>
            </div>
        </div>
    </div>

    <script>
        // Virtual Screening
        async function runVirtualScreening() {
            const smiles = document.getElementById('vhts-smiles').value;
            const loading = document.getElementById('vhts-loading');
            const results = document.getElementById('vhts-results');
            
            loading.classList.add('active');
            results.style.display = 'none';
            
            try {
                const response = await fetch('/api/vhts/screen_compound', {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({ smiles: smiles, name: 'Test Compound' })
                });
                
                const data = await response.json();
                
                results.innerHTML = `
                    <h4>Screening Results</h4>
                    <p><strong>Compound:</strong> ${data.compound_id}</p>
                    <p><strong>Hits:</strong> ${data.receptor_hits ? data.receptor_hits.length : 0} receptors</p>
                    <p><strong>Selectivity:</strong> ${data.selectivity_profile?.status || 'N/A'}</p>
                    <p><strong>Polypharmacology Score:</strong> ${data.polypharmacology_score || 0}</p>
                    ${data.safety_flags?.length ? `<p style="color: red;"><strong>Safety Flags:</strong> ${data.safety_flags.join(', ')}</p>` : ''}
                    <h5>Top Receptor Hits:</h5>
                    <ul>${(data.receptor_hits || []).slice(0, 5).map(h => 
                        `<li>${h.receptor} - ${h.binding_score} (${h.predicted_ki})</li>`
                    ).join('')}</ul>
                `;
                
                results.style.display = 'block';
            } catch (error) {
                results.innerHTML = `<p style="color: red;">Error: ${error.message}</p>`;
                results.style.display = 'block';
            } finally {
                loading.classList.remove('active');
            }
        }
        
        // Lead Optimization
        async function optimizeLead() {
            const smiles = document.getElementById('lead-smiles').value;
            const loading = document.getElementById('lead-loading');
            const results = document.getElementById('lead-results');
            
            loading.classList.add('active');
            results.style.display = 'none';
            
            try {
                const response = await fetch('/api/lead_opt/optimize', {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({ smiles: smiles })
                });
                
                const data = await response.json();
                
                results.innerHTML = `
                    <h4>Lead Optimization Suggestions</h4>
                    <p><strong>Current Properties:</strong></p>
                    <ul>
                        <li>MW: ${data.current_properties?.mw}</li>
                        <li>LogP: ${data.current_properties?.logp}</li>
                        <li>QED: ${data.current_properties?.qed}</li>
                    </ul>
                    <h5>Top Modifications:</h5>
                    <ol>${(data.modified_structures || []).slice(0, 3).map(m => 
                        `<li>${m.modification} (Score: ${m.improvement_score})</li>`
                    ).join('')}</ol>
                `;
                
                results.style.display = 'block';
            } catch (error) {
                results.innerHTML = `<p style="color: red;">Error: ${error.message}</p>`;
                results.style.display = 'block';
            } finally {
                loading.classList.remove('active');
            }
        }
        
        // SAR Analysis
        async function analyzeSAR() {
            const compounds = JSON.parse(document.getElementById('sar-compounds').value);
            const loading = document.getElementById('sar-loading');
            const results = document.getElementById('sar-results');
            
            loading.classList.add('active');
            results.style.display = 'none';
            
            try {
                const response = await fetch('/api/sar/analyze', {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({ compounds: compounds })
                });
                
                const data = await response.json();
                
                results.innerHTML = `
                    <h4>SAR Analysis</h4>
                    <p><strong>Series Size:</strong> ${data.series_size} compounds</p>
                    <p><strong>Activity Cliffs:</strong> ${data.activity_cliffs?.length || 0} identified</p>
                    <h5>Key Features:</h5>
                    <ul>${(data.key_features || []).slice(0, 5).map(f => 
                        `<li>${f.fragment} - ${f.importance} (contribution: ${f.contribution})</li>`
                    ).join('')}</ul>
                `;
                
                results.style.display = 'block';
            } catch (error) {
                results.innerHTML = `<p style="color: red;">Error: ${error.message}</p>`;
                results.style.display = 'block';
            } finally {
                loading.classList.remove('active');
            }
        }
        
        // Pharmacophore
        async function extractPharmacophore() {
            const compounds = JSON.parse(document.getElementById('pharm-actives').value);
            const loading = document.getElementById('pharm-loading');
            const results = document.getElementById('pharm-results');
            
            loading.classList.add('active');
            results.style.display = 'none';
            
            try {
                const response = await fetch('/api/pharmacophore/extract', {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({ active_compounds: compounds })
                });
                
                const data = await response.json();
                
                results.innerHTML = `
                    <h4>Pharmacophore Model</h4>
                    <p><strong>Confidence:</strong> ${data.confidence}%</p>
                    <p><strong>Common Features:</strong></p>
                    <ul>${(data.common_features || []).map(f => 
                        `<li>${f.type} ${f.required ? '(Required)' : '(Optional)'}</li>`
                    ).join('')}</ul>
                `;
                
                results.style.display = 'block';
            } catch (error) {
                results.innerHTML = `<p style="color: red;">Error: ${error.message}</p>`;
                results.style.display = 'block';
            } finally {
                loading.classList.remove('active');
            }
        }
        
        async function getPharmacophoreModel() {
            const receptor = document.getElementById('pharm-receptor').value;
            const loading = document.getElementById('pharm-loading');
            const results = document.getElementById('pharm-results');
            
            loading.classList.add('active');
            results.style.display = 'none';
            
            try {
                const response = await fetch(`/api/pharmacophore/get_model/${receptor}`);
                const data = await response.json();
                
                results.innerHTML = `
                    <h4>${receptor} Pharmacophore Model</h4>
                    <p><strong>Features:</strong> ${data.num_features}</p>
                    <p><strong>Constraints:</strong> ${data.num_constraints}</p>
                    <p>${data.description}</p>
                `;
                
                results.style.display = 'block';
            } catch (error) {
                results.innerHTML = `<p style="color: red;">Error: ${error.message}</p>`;
                results.style.display = 'block';
            } finally {
                loading.classList.remove('active');
            }
        }
        
        // Safety Prediction
        async function predictSafety() {
            const smiles = document.getElementById('safety-smiles').value;
            const target = document.getElementById('safety-target').value;
            const loading = document.getElementById('safety-loading');
            const results = document.getElementById('safety-results');
            
            loading.classList.add('active');
            results.style.display = 'none';
            
            try {
                const response = await fetch('/api/off_target/predict', {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({ smiles: smiles, primary_target: target })
                });
                
                const data = await response.json();
                const prediction = data.prediction;
                
                results.innerHTML = `
                    <h4>Safety Assessment</h4>
                    <p><strong>Risk Score:</strong> ${prediction.risk_score}/100</p>
                    <p><strong>Classification:</strong> <span style="color: ${prediction.risk_score > 50 ? 'red' : 'green'};">${prediction.safety_classification}</span></p>
                    <p><strong>DDI Risk:</strong> ${prediction.drug_drug_interaction_risk}</p>
                    ${prediction.safety_alerts?.length ? 
                        `<h5>Safety Alerts:</h5><ul>${prediction.safety_alerts.map(a => 
                            `<li style="color: ${a.level === 'CRITICAL' ? 'red' : 'orange'};">${a.message}</li>`
                        ).join('')}</ul>` : ''}
                    ${prediction.predicted_side_effects?.length ?
                        `<h5>Predicted Side Effects:</h5><ul>${prediction.predicted_side_effects.slice(0, 5).map(e => 
                            `<li>${e}</li>`
                        ).join('')}</ul>` : ''}
                `;
                
                results.style.display = 'block';
            } catch (error) {
                results.innerHTML = `<p style="color: red;">Error: ${error.message}</p>`;
                results.style.display = 'block';
            } finally {
                loading.classList.remove('active');
            }
        }
        
        // Full Analysis
        async function runFullAnalysis() {
            const smiles = document.getElementById('full-smiles').value;
            const name = document.getElementById('full-name').value;
            const loading = document.getElementById('full-loading');
            const results = document.getElementById('full-results');
            
            loading.classList.add('active');
            results.style.display = 'none';
            
            try {
                const response = await fetch('/api/advanced/full_analysis', {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({ smiles: smiles, name: name })
                });
                
                const data = await response.json();
                
                let html = `<h4>Comprehensive Analysis: ${data.compound}</h4>`;
                
                // Receptor Screening
                if (data.analyses.receptor_screening) {
                    const rs = data.analyses.receptor_screening;
                    html += `
                        <h5>Receptor Screening:</h5>
                        <p>Hits: ${rs.receptor_hits?.length || 0} | 
                           Selectivity: ${rs.selectivity_profile?.status || 'N/A'} | 
                           Polypharmacology: ${rs.polypharmacology_score || 0}</p>
                    `;
                }
                
                // Safety
                if (data.analyses.safety_assessment) {
                    const sa = data.analyses.safety_assessment;
                    html += `
                        <h5>Safety Profile:</h5>
                        <p>Risk Score: ${sa.risk_score}/100 | 
                           Classification: ${sa.safety_classification} | 
                           Off-targets: ${sa.off_target_hits?.length || 0}</p>
                    `;
                }
                
                // Optimization
                if (data.analyses.optimization_suggestions) {
                    const os = data.analyses.optimization_suggestions;
                    html += `
                        <h5>Optimization Potential:</h5>
                        <p>Strategies: ${os.optimization_strategies?.length || 0} | 
                           Suggested modifications: ${os.modified_structures?.length || 0}</p>
                    `;
                }
                
                results.innerHTML = html;
                results.style.display = 'block';
            } catch (error) {
                results.innerHTML = `<p style="color: red;">Error: ${error.message}</p>`;
                results.style.display = 'block';
            } finally {
                loading.classList.remove('active');
            }
        }
        
        // Tab switching
        function switchTab(feature, tab) {
            const tabs = document.querySelectorAll(`#${feature}-${tab}`).parentElement.querySelectorAll('.tab');
            const contents = document.querySelectorAll(`#${feature}-${tab}`).parentElement.parentElement.querySelectorAll('.tab-content');
            
            tabs.forEach(t => t.classList.remove('active'));
            contents.forEach(c => c.classList.remove('active'));
            
            event.target.classList.add('active');
            document.getElementById(`${feature}-${tab}`).classList.add('active');
        }
    </script>
</body>
</html>
'''

def register_advanced_features(app):
    """Register advanced features blueprint with the main app"""
    app.register_blueprint(advanced_bp)
    
    # Register quantum computing routes
    quantum_bp = create_quantum_routes(app)
    app.register_blueprint(quantum_bp)
    
    print("✅ Advanced Drug Discovery Features Integrated")
    print("   - Virtual High-Throughput Screening (82 receptors)")
    print("   - AI-Powered Lead Optimization")
    print("   - Structure-Activity Relationship Explorer")
    print("   - Pharmacophore Modeling")
    print("   - Off-Target Prediction System")
    print("   - Quantum Computing Module")
    return app
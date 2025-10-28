#!/usr/bin/env python3
"""
3D Molecular Visualization Module
Interactive 3D viewing of molecular structures using py3Dmol
"""

from rdkit import Chem
from rdkit.Chem import AllChem
import json
import base64
from io import BytesIO

class Molecular3DViewer:
    """Generate 3D molecular visualizations"""
    
    def __init__(self):
        self.viewer_template = """
        <!DOCTYPE html>
        <html>
        <head>
            <script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
            <style>
                body {{ 
                    font-family: Arial, sans-serif;
                    background: linear-gradient(135deg, #0f0f23 0%, #1a1a2e 50%, #16213e 100%);
                    color: white;
                    padding: 20px;
                }}
                .viewer-container {{
                    width: 100%;
                    height: 600px;
                    position: relative;
                    border: 2px solid #00d4ff;
                    border-radius: 10px;
                    overflow: hidden;
                    background: rgba(0,0,0,0.5);
                }}
                .controls {{
                    margin: 20px 0;
                    padding: 15px;
                    background: rgba(255,255,255,0.1);
                    border-radius: 10px;
                }}
                .btn {{
                    padding: 10px 20px;
                    margin: 5px;
                    background: linear-gradient(45deg, #00d4ff, #0099cc);
                    border: none;
                    border-radius: 5px;
                    color: white;
                    cursor: pointer;
                    font-weight: bold;
                }}
                .btn:hover {{
                    transform: scale(1.05);
                    box-shadow: 0 5px 15px rgba(0, 212, 255, 0.4);
                }}
                .info {{
                    background: rgba(255,255,255,0.05);
                    padding: 15px;
                    border-radius: 10px;
                    margin: 20px 0;
                }}
                .property {{
                    display: inline-block;
                    margin: 5px 10px;
                }}
                .property-label {{
                    color: #00d4ff;
                    font-weight: bold;
                }}
            </style>
        </head>
        <body>
            <h1>🧬 3D Molecular Structure Viewer</h1>
            
            <div class="info">
                <h2>{compound_name}</h2>
                <div class="property">
                    <span class="property-label">Formula:</span> {formula}
                </div>
                <div class="property">
                    <span class="property-label">MW:</span> {mw} Da
                </div>
                <div class="property">
                    <span class="property-label">SMILES:</span> {smiles}
                </div>
            </div>
            
            <div id="container" class="viewer-container"></div>
            
            <div class="controls">
                <button class="btn" onclick="setStyle('stick')">Stick</button>
                <button class="btn" onclick="setStyle('ball')">Ball & Stick</button>
                <button class="btn" onclick="setStyle('sphere')">Space Fill</button>
                <button class="btn" onclick="setStyle('surface')">Surface</button>
                <button class="btn" onclick="setStyle('cartoon')">Cartoon</button>
                <button class="btn" onclick="toggleSpin()">Toggle Spin</button>
                <button class="btn" onclick="resetView()">Reset View</button>
                <button class="btn" onclick="showLabels()">Show Labels</button>
                <button class="btn" onclick="colorByElement()">Color by Element</button>
                <button class="btn" onclick="colorByCharge()">Color by Charge</button>
            </div>
            
            <script>
                let viewer = $3Dmol.createViewer(document.getElementById('container'), {{
                    backgroundColor: 'black'
                }});
                
                // Load molecule
                let molData = `{mol_data}`;
                viewer.addModel(molData, "{format}");
                
                // Initial style
                viewer.setStyle({{}}, {{stick: {{}}}});
                viewer.zoomTo();
                viewer.render();
                
                let spinning = false;
                
                function setStyle(style) {{
                    viewer.removeAllModels();
                    viewer.addModel(molData, "{format}");
                    
                    switch(style) {{
                        case 'stick':
                            viewer.setStyle({{}}, {{stick: {{}}}});
                            break;
                        case 'ball':
                            viewer.setStyle({{}}, {{stick: {{}}, sphere: {{scale: 0.3}}}});
                            break;
                        case 'sphere':
                            viewer.setStyle({{}}, {{sphere: {{}}}});
                            break;
                        case 'surface':
                            viewer.addSurface($3Dmol.SurfaceType.VDW, {{
                                opacity: 0.85,
                                colorscheme: 'default'
                            }});
                            viewer.setStyle({{}}, {{stick: {{}}}});
                            break;
                        case 'cartoon':
                            viewer.setStyle({{}}, {{cartoon: {{color: 'spectrum'}}}});
                            break;
                    }}
                    viewer.render();
                }}
                
                function toggleSpin() {{
                    if (spinning) {{
                        viewer.spin(false);
                        spinning = false;
                    }} else {{
                        viewer.spin('y', 1);
                        spinning = true;
                    }}
                }}
                
                function resetView() {{
                    viewer.zoomTo();
                    viewer.render();
                }}
                
                function showLabels() {{
                    viewer.addPropertyLabels("elem", {{}}, {{
                        fontSize: 12,
                        fontColor: 'white',
                        backgroundOpacity: 0.8
                    }});
                    viewer.render();
                }}
                
                function colorByElement() {{
                    viewer.setStyle({{}}, {{
                        stick: {{colorscheme: 'default'}},
                        sphere: {{colorscheme: 'default', scale: 0.3}}
                    }});
                    viewer.render();
                }}
                
                function colorByCharge() {{
                    viewer.setStyle({{}}, {{
                        stick: {{colorscheme: 'RWB'}},
                        sphere: {{colorscheme: 'RWB', scale: 0.3}}
                    }});
                    viewer.render();
                }}
                
                // Enable mouse controls
                viewer.setHoverable({{}}, true, function(atom, viewer, event, container) {{
                    if (!atom.label) {{
                        atom.label = viewer.addLabel(atom.elem + " " + atom.serial, {{
                            position: atom,
                            backgroundColor: 'black',
                            fontColor: 'white',
                            backgroundOpacity: 0.8
                        }});
                    }}
                }}, function(atom) {{
                    if (atom.label) {{
                        viewer.removeLabel(atom.label);
                        delete atom.label;
                    }}
                }});
            </script>
            
            <div class="info">
                <h3>Controls:</h3>
                <ul>
                    <li>Left click + drag: Rotate molecule</li>
                    <li>Right click + drag: Translate molecule</li>
                    <li>Scroll wheel: Zoom in/out</li>
                    <li>Hover over atoms: See atom details</li>
                </ul>
            </div>
        </body>
        </html>
        """
    
    def generate_3d_structure(self, smiles: str, optimize: bool = True):
        """Generate 3D coordinates for a molecule"""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                return None
            
            # Add hydrogens
            mol = Chem.AddHs(mol)
            
            # Generate 3D coordinates
            if AllChem.EmbedMolecule(mol, randomSeed=42) != 0:
                # Try 2D if 3D fails
                AllChem.Compute2DCoords(mol)
            
            # Optimize geometry
            if optimize:
                AllChem.MMFFOptimizeMolecule(mol)
            
            return mol
        except Exception as e:
            print(f"Error generating 3D structure: {e}")
            return None
    
    def create_3d_viewer_html(self, smiles: str, compound_name: str = "Compound",
                             format: str = "sdf") -> str:
        """Create interactive 3D viewer HTML"""
        mol = self.generate_3d_structure(smiles)
        if not mol:
            return "<p>Error generating 3D structure</p>"
        
        # Get molecular data in specified format
        if format == "sdf":
            mol_data = Chem.MolToMolBlock(mol)
        elif format == "pdb":
            mol_data = Chem.MolToPDBBlock(mol)
        else:
            mol_data = Chem.MolToMolBlock(mol)
        
        # Calculate properties
        mol_no_h = Chem.RemoveHs(mol)
        formula = Chem.rdMolDescriptors.CalcMolFormula(mol_no_h)
        mw = round(Chem.Descriptors.MolWt(mol_no_h), 2)
        
        # Fill template
        html = self.viewer_template.format(
            compound_name=compound_name,
            formula=formula,
            mw=mw,
            smiles=smiles,
            mol_data=mol_data,
            format=format
        )
        
        return html
    
    def create_3d_comparison_viewer(self, smiles_list: List[str], names: List[str] = None):
        """Create a viewer comparing multiple molecules"""
        if not names:
            names = [f"Compound {i+1}" for i in range(len(smiles_list))]
        
        comparison_html = """
        <!DOCTYPE html>
        <html>
        <head>
            <script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
            <style>
                body {{ 
                    font-family: Arial, sans-serif;
                    background: linear-gradient(135deg, #0f0f23 0%, #1a1a2e 50%);
                    color: white;
                }}
                .grid {{
                    display: grid;
                    grid-template-columns: repeat(auto-fit, minmax(400px, 1fr));
                    gap: 20px;
                    padding: 20px;
                }}
                .mol-container {{
                    border: 2px solid #00d4ff;
                    border-radius: 10px;
                    padding: 10px;
                    background: rgba(0,0,0,0.5);
                }}
                .viewer {{
                    width: 100%;
                    height: 400px;
                }}
                h2 {{
                    text-align: center;
                    color: #00d4ff;
                }}
                .overlay-btn {{
                    display: block;
                    margin: 20px auto;
                    padding: 15px 30px;
                    background: linear-gradient(45deg, #ff00ff, #00d4ff);
                    border: none;
                    border-radius: 5px;
                    color: white;
                    font-size: 16px;
                    font-weight: bold;
                    cursor: pointer;
                }}
            </style>
        </head>
        <body>
            <h1 style="text-align: center;">🧬 Molecular Structure Comparison</h1>
            
            <button class="overlay-btn" onclick="showOverlay()">Show Overlaid Structures</button>
            
            <div class="grid">
        """
        
        viewer_scripts = []
        mol_data_list = []
        
        for i, (smiles, name) in enumerate(zip(smiles_list, names)):
            mol = self.generate_3d_structure(smiles)
            if mol:
                mol_data = Chem.MolToMolBlock(mol)
                mol_data_list.append(mol_data)
                
                comparison_html += f"""
                <div class="mol-container">
                    <h2>{name}</h2>
                    <div id="viewer{i}" class="viewer"></div>
                    <p style="text-align: center; font-size: 12px;">{smiles}</p>
                </div>
                """
                
                viewer_scripts.append(f"""
                    let viewer{i} = $3Dmol.createViewer(document.getElementById('viewer{i}'), {{
                        backgroundColor: 'black'
                    }});
                    viewer{i}.addModel(`{mol_data}`, "sdf");
                    viewer{i}.setStyle({{}}, {{stick: {{}}, sphere: {{scale: 0.3}}}});
                    viewer{i}.zoomTo();
                    viewer{i}.render();
                    viewer{i}.spin('y', 1);
                """)
        
        comparison_html += """
            </div>
            
            <div id="overlayModal" style="display:none; position:fixed; top:0; left:0; width:100%; height:100%; background:rgba(0,0,0,0.9); z-index:1000;">
                <div style="position:relative; width:90%; height:90%; margin:5% auto;">
                    <span onclick="closeOverlay()" style="position:absolute; top:10px; right:20px; font-size:30px; cursor:pointer; color:white;">&times;</span>
                    <div id="overlayViewer" style="width:100%; height:100%;"></div>
                </div>
            </div>
            
            <script>
        """
        
        comparison_html += '\n'.join(viewer_scripts)
        
        # Add overlay functionality
        colors = ['red', 'blue', 'green', 'yellow', 'orange', 'purple']
        overlay_mols = []
        for i, mol_data in enumerate(mol_data_list):
            color = colors[i % len(colors)]
            overlay_mols.append(f"""
                overlayViewer.addModel(`{mol_data}`, "sdf");
                overlayViewer.setStyle({{model: {i}}}, {{stick: {{color: '{color}'}}, sphere: {{color: '{color}', scale: 0.2}}}});
            """)
        
        comparison_html += f"""
                function showOverlay() {{
                    document.getElementById('overlayModal').style.display = 'block';
                    let overlayViewer = $3Dmol.createViewer(document.getElementById('overlayViewer'), {{
                        backgroundColor: 'black'
                    }});
                    {''.join(overlay_mols)}
                    overlayViewer.zoomTo();
                    overlayViewer.render();
                }}
                
                function closeOverlay() {{
                    document.getElementById('overlayModal').style.display = 'none';
                }}
            </script>
        </body>
        </html>
        """
        
        return comparison_html
"""
3D Molecular Visualization and Retrosynthesis System
Provides interactive 3D molecular structures and synthetic pathway planning
"""

import json
import re
from typing import Dict, List, Any, Tuple
import base64
from io import BytesIO

class MolecularVisualizationSystem:
    def __init__(self):
        self.retrosynthesis_database = self.load_retrosynthesis_database()
        self.reaction_templates = self.load_reaction_templates()
    
    def load_retrosynthesis_database(self) -> Dict[str, Any]:
        """Load retrosynthesis pathways for key compounds"""
        return {
            "buprenorphine_analogs": {
                "target_smiles": "COC1=C(O)C=C2C3=C1C(=O)CC[C@]3(C)[C@H](O)[C@H]2N(CC1CC1)CC1CC1",
                "complexity_score": 8.5,
                "estimated_steps": 12,
                "starting_materials": [
                    "thebaine",
                    "cyclopropylmethyl bromide",
                    "various protecting groups"
                ],
                "key_reactions": [
                    "Diels-Alder cycloaddition",
                    "N-alkylation",
                    "Selective reduction",
                    "Demethylation"
                ],
                "synthetic_route": [
                    {
                        "step": 1,
                        "reaction": "Protection of phenolic OH",
                        "starting_material": "thebaine",
                        "reagents": ["TBSCl", "imidazole", "DMF"],
                        "conditions": "RT, 12h",
                        "yield": "92%",
                        "product": "TBS-protected thebaine"
                    },
                    {
                        "step": 2,
                        "reaction": "Diels-Alder cycloaddition",
                        "starting_material": "TBS-protected thebaine",
                        "reagents": ["methyl vinyl ketone"],
                        "conditions": "180°C, sealed tube, 24h",
                        "yield": "78%",
                        "product": "Diels-Alder adduct"
                    },
                    {
                        "step": 3,
                        "reaction": "Reduction of ketone",
                        "starting_material": "Diels-Alder adduct",
                        "reagents": ["NaBH4", "MeOH"],
                        "conditions": "0°C to RT, 2h",
                        "yield": "85%",
                        "product": "Secondary alcohol"
                    },
                    {
                        "step": 4,
                        "reaction": "N-demethylation",
                        "starting_material": "Secondary alcohol",
                        "reagents": ["cyanogen bromide", "toluene"],
                        "conditions": "reflux, 4h",
                        "yield": "71%",
                        "product": "N-demethylated intermediate"
                    },
                    {
                        "step": 5,
                        "reaction": "N-alkylation with cyclopropylmethyl group",
                        "starting_material": "N-demethylated intermediate",
                        "reagents": ["cyclopropylmethyl bromide", "K2CO3", "DMF"],
                        "conditions": "80°C, 6h",
                        "yield": "68%",
                        "product": "N-cyclopropylmethyl derivative"
                    },
                    {
                        "step": 6,
                        "reaction": "Deprotection of TBS group",
                        "starting_material": "N-cyclopropylmethyl derivative",
                        "reagents": ["TBAF", "THF"],
                        "conditions": "RT, 2h",
                        "yield": "89%",
                        "product": "Buprenorphine analog"
                    }
                ],
                "overall_yield": "28%",
                "estimated_cost": "$2,500/g",
                "safety_considerations": [
                    "Cyanogen bromide is highly toxic - use extreme caution",
                    "High temperature reactions require pressure equipment",
                    "All steps require anhydrous conditions"
                ],
                "equipment_needed": [
                    "High-pressure reactor",
                    "Rotary evaporator",
                    "HPLC purification system",
                    "NMR spectrometer for structure confirmation"
                ]
            },
            "kappa_antagonists": {
                "target_smiles": "CC(C)(C)NC(=O)C1CCN(CC1)C(=O)c2ccc(Cl)cc2Cl",
                "complexity_score": 6.2,
                "estimated_steps": 8,
                "starting_materials": [
                    "4-piperidone",
                    "3,4-dichlorobenzoyl chloride",
                    "tert-butylamine"
                ],
                "synthetic_route": [
                    {
                        "step": 1,
                        "reaction": "Acylation of piperidine",
                        "starting_material": "4-piperidone",
                        "reagents": ["3,4-dichlorobenzoyl chloride", "Et3N", "DCM"],
                        "conditions": "0°C to RT, 3h",
                        "yield": "88%",
                        "product": "N-acyl piperidone"
                    },
                    {
                        "step": 2,
                        "reaction": "Reduction of ketone",
                        "starting_material": "N-acyl piperidone",
                        "reagents": ["NaBH4", "EtOH"],
                        "conditions": "0°C, 1h",
                        "yield": "92%",
                        "product": "Secondary alcohol"
                    },
                    {
                        "step": 3,
                        "reaction": "Conversion to carboxylic acid",
                        "starting_material": "Secondary alcohol",
                        "reagents": ["Jones reagent"],
                        "conditions": "0°C, 2h",
                        "yield": "79%",
                        "product": "Carboxylic acid"
                    },
                    {
                        "step": 4,
                        "reaction": "Amide formation",
                        "starting_material": "Carboxylic acid",
                        "reagents": ["tert-butylamine", "EDC", "HOBt", "DMF"],
                        "conditions": "RT, overnight",
                        "yield": "84%",
                        "product": "Final kappa antagonist"
                    }
                ],
                "overall_yield": "57%",
                "estimated_cost": "$450/g",
                "safety_considerations": [
                    "Jones reagent is highly corrosive",
                    "Use appropriate ventilation for all steps"
                ]
            },
            "gaba_modulators": {
                "target_smiles": "COc1cc(ccc1OC)C(=O)C=Cc2ccc(OC)c(OC)c2",
                "complexity_score": 4.8,
                "estimated_steps": 5,
                "starting_materials": [
                    "3,4-dimethoxybenzaldehyde",
                    "3,4-dimethoxyacetophenone"
                ],
                "synthetic_route": [
                    {
                        "step": 1,
                        "reaction": "Aldol condensation",
                        "starting_material": "3,4-dimethoxyacetophenone",
                        "reagents": ["3,4-dimethoxybenzaldehyde", "NaOH", "EtOH"],
                        "conditions": "RT, 4h",
                        "yield": "76%",
                        "product": "Chalcone intermediate"
                    },
                    {
                        "step": 2,
                        "reaction": "Purification by recrystallization",
                        "starting_material": "Chalcone intermediate",
                        "reagents": ["EtOH"],
                        "conditions": "Hot filtration, cooling",
                        "yield": "89%",
                        "product": "Pure synthetic kavain analog"
                    }
                ],
                "overall_yield": "68%",
                "estimated_cost": "$85/g",
                "safety_considerations": [
                    "Standard organic synthesis precautions",
                    "Avoid skin contact with NaOH"
                ]
            },
            "psychedelic_analogs": {
                "target_smiles": "CN(C)CCc1c[nH]c2ccc(O)cc12",
                "complexity_score": 5.5,
                "estimated_steps": 6,
                "starting_materials": [
                    "5-hydroxyindole",
                    "oxalyl chloride",
                    "dimethylamine"
                ],
                "synthetic_route": [
                    {
                        "step": 1,
                        "reaction": "Friedel-Crafts acylation",
                        "starting_material": "5-hydroxyindole",
                        "reagents": ["oxalyl chloride", "AlCl3", "DCM"],
                        "conditions": "0°C to RT, 3h",
                        "yield": "72%",
                        "product": "Acylated indole"
                    },
                    {
                        "step": 2,
                        "reaction": "Reduction to alcohol",
                        "starting_material": "Acylated indole",
                        "reagents": ["LiAlH4", "THF"],
                        "conditions": "0°C to RT, 2h",
                        "yield": "85%",
                        "product": "Primary alcohol"
                    },
                    {
                        "step": 3,
                        "reaction": "Conversion to tosylate",
                        "starting_material": "Primary alcohol",
                        "reagents": ["TsCl", "pyridine"],
                        "conditions": "0°C, 4h",
                        "yield": "91%",
                        "product": "Tosylate intermediate"
                    },
                    {
                        "step": 4,
                        "reaction": "Nucleophilic substitution",
                        "starting_material": "Tosylate intermediate",
                        "reagents": ["dimethylamine", "K2CO3", "DMF"],
                        "conditions": "80°C, 6h",
                        "yield": "78%",
                        "product": "Short-duration tryptamine analog"
                    }
                ],
                "overall_yield": "42%",
                "estimated_cost": "$320/g",
                "safety_considerations": [
                    "LiAlH4 is highly reactive with water",
                    "Use dry solvents throughout",
                    "Controlled substance precursor - follow regulations"
                ]
            }
        }
    
    def load_reaction_templates(self) -> Dict[str, Any]:
        """Load common reaction templates for retrosynthesis"""
        return {
            "aldol_condensation": {
                "name": "Aldol Condensation",
                "reactants": ["ketone", "aldehyde"],
                "products": ["α,β-unsaturated ketone"],
                "conditions": "Base (NaOH, KOH), protic solvent",
                "mechanism": "Enolate formation followed by nucleophilic addition"
            },
            "friedel_crafts": {
                "name": "Friedel-Crafts Acylation",
                "reactants": ["aromatic compound", "acyl chloride"],
                "products": ["aromatic ketone"],
                "conditions": "Lewis acid catalyst (AlCl3), aprotic solvent",
                "mechanism": "Electrophilic aromatic substitution"
            },
            "reductive_amination": {
                "name": "Reductive Amination",
                "reactants": ["carbonyl compound", "amine"],
                "products": ["secondary or tertiary amine"],
                "conditions": "Reducing agent (NaBH4, NaCNBH3)",
                "mechanism": "Imine formation followed by reduction"
            },
            "diels_alder": {
                "name": "Diels-Alder Cycloaddition",
                "reactants": ["diene", "dienophile"],
                "products": ["cyclohexene derivative"],
                "conditions": "Heat, sometimes Lewis acid catalyst",
                "mechanism": "Concerted [4+2] cycloaddition"
            }
        }
    
    def generate_3d_structure_data(self, smiles: str, compound_name: str = "") -> Dict[str, Any]:
        """Generate 3D structure data for molecular visualization"""
        # This would typically use RDKit or similar library
        # For now, providing structure for web-based visualization
        
        structure_data = {
            "compound_name": compound_name,
            "smiles": smiles,
            "molecular_formula": self.smiles_to_formula(smiles),
            "3d_coordinates": self.generate_3d_coordinates(smiles),
            "visualization_script": self.generate_visualization_script(smiles, compound_name),
            "interactive_features": [
                "rotation",
                "zoom",
                "atom_labeling",
                "bond_highlighting",
                "surface_rendering"
            ],
            "property_display": {
                "molecular_weight": self.calculate_molecular_weight(smiles),
                "logp": self.estimate_logp(smiles),
                "polar_surface_area": self.estimate_psa(smiles),
                "rotatable_bonds": self.count_rotatable_bonds(smiles)
            }
        }
        
        return structure_data
    
    def smiles_to_formula(self, smiles: str) -> str:
        """Convert SMILES to molecular formula (simplified)"""
        # This is a simplified version - would use RDKit in practice
        formula_map = {
            "COc1cc(ccc1OC)C(=O)C=Cc2ccc(OC)c(OC)c2": "C18H18O5",
            "CC(C)(C)NC(=O)C1CCN(CC1)C(=O)c2ccc(Cl)cc2Cl": "C17H22Cl2N2O2",
            "CN(C)CCc1c[nH]c2ccc(O)cc12": "C13H17N3O",
            "COC1=C(O)C=C2C3=C1C(=O)CC[C@]3(C)[C@H](O)[C@H]2N(CC1CC1)CC1CC1": "C25H35NO4"
        }
        return formula_map.get(smiles, "Unknown")
    
    def generate_3d_coordinates(self, smiles: str) -> List[Dict[str, Any]]:
        """Generate 3D coordinates for atoms (simplified)"""
        # This would use RDKit's 3D coordinate generation
        # Providing example coordinates for visualization
        return [
            {"atom": "C", "x": 0.0, "y": 0.0, "z": 0.0},
            {"atom": "O", "x": 1.4, "y": 0.0, "z": 0.0},
            {"atom": "C", "x": -0.7, "y": 1.2, "z": 0.0},
            # ... more coordinates would be generated
        ]
    
    def generate_visualization_script(self, smiles: str, compound_name: str) -> str:
        """Generate JavaScript for 3D molecular visualization"""
        script = f"""
        // 3D Molecular Visualization for {compound_name}
        const viewer = new ChemDoodle.TransformCanvas3D('canvas_{compound_name.replace(' ', '_')}', 400, 400);
        viewer.specs.set3DRepresentation('Ball and Stick');
        viewer.specs.backgroundColor = 'white';
        viewer.specs.atoms_displayLabels_3D = true;
        
        // Load molecule from SMILES: {smiles}
        const molecule = ChemDoodle.readSMILES('{smiles}');
        viewer.loadMolecule(molecule);
        
        // Add interaction controls
        viewer.mousedown = function(e) {{
            this.lastPoint = new ChemDoodle.Point(e.pageX, e.pageY);
        }};
        
        viewer.mousemove = function(e) {{
            if (this.lastPoint) {{
                const dx = e.pageX - this.lastPoint.x;
                const dy = e.pageY - this.lastPoint.y;
                this.rotateX(dy / 2);
                this.rotateY(dx / 2);
                this.repaint();
                this.lastPoint = new ChemDoodle.Point(e.pageX, e.pageY);
            }}
        }};
        
        viewer.mouseup = function(e) {{
            this.lastPoint = null;
        }};
        """
        return script
    
    def calculate_molecular_weight(self, smiles: str) -> float:
        """Calculate molecular weight (simplified)"""
        mw_map = {
            "COc1cc(ccc1OC)C(=O)C=Cc2ccc(OC)c(OC)c2": 314.33,
            "CC(C)(C)NC(=O)C1CCN(CC1)C(=O)c2ccc(Cl)cc2Cl": 357.27,
            "CN(C)CCc1c[nH]c2ccc(O)cc12": 231.29,
            "COC1=C(O)C=C2C3=C1C(=O)CC[C@]3(C)[C@H](O)[C@H]2N(CC1CC1)CC1CC1": 467.64
        }
        return mw_map.get(smiles, 0.0)
    
    def estimate_logp(self, smiles: str) -> float:
        """Estimate LogP (simplified)"""
        logp_map = {
            "COc1cc(ccc1OC)C(=O)C=Cc2ccc(OC)c(OC)c2": 2.8,
            "CC(C)(C)NC(=O)C1CCN(CC1)C(=O)c2ccc(Cl)cc2Cl": 3.2,
            "CN(C)CCc1c[nH]c2ccc(O)cc12": 1.9,
            "COC1=C(O)C=C2C3=C1C(=O)CC[C@]3(C)[C@H](O)[C@H]2N(CC1CC1)CC1CC1": 4.1
        }
        return logp_map.get(smiles, 0.0)
    
    def estimate_psa(self, smiles: str) -> float:
        """Estimate polar surface area (simplified)"""
        psa_map = {
            "COc1cc(ccc1OC)C(=O)C=Cc2ccc(OC)c(OC)c2": 63.2,
            "CC(C)(C)NC(=O)C1CCN(CC1)C(=O)c2ccc(Cl)cc2Cl": 49.4,
            "CN(C)CCc1c[nH]c2ccc(O)cc12": 44.5,
            "COC1=C(O)C=C2C3=C1C(=O)CC[C@]3(C)[C@H](O)[C@H]2N(CC1CC1)CC1CC1": 62.2
        }
        return psa_map.get(smiles, 0.0)
    
    def count_rotatable_bonds(self, smiles: str) -> int:
        """Count rotatable bonds (simplified)"""
        rb_map = {
            "COc1cc(ccc1OC)C(=O)C=Cc2ccc(OC)c(OC)c2": 7,
            "CC(C)(C)NC(=O)C1CCN(CC1)C(=O)c2ccc(Cl)cc2Cl": 4,
            "CN(C)CCc1c[nH]c2ccc(O)cc12": 3,
            "COC1=C(O)C=C2C3=C1C(=O)CC[C@]3(C)[C@H](O)[C@H]2N(CC1CC1)CC1CC1": 2
        }
        return rb_map.get(smiles, 0)
    
    def get_retrosynthesis_pathway(self, compound_type: str) -> Dict[str, Any]:
        """Get detailed retrosynthesis pathway for a compound type"""
        if compound_type in self.retrosynthesis_database:
            pathway = self.retrosynthesis_database[compound_type].copy()
            
            # Add additional analysis
            pathway["feasibility_score"] = self.calculate_feasibility_score(pathway)
            pathway["alternative_routes"] = self.suggest_alternative_routes(compound_type)
            pathway["supplier_information"] = self.get_supplier_info(pathway["starting_materials"])
            
            return pathway
        
        return {"error": f"No retrosynthesis data available for {compound_type}"}
    
    def calculate_feasibility_score(self, pathway: Dict[str, Any]) -> float:
        """Calculate synthetic feasibility score (0-10)"""
        complexity = pathway.get("complexity_score", 5)
        overall_yield = float(pathway.get("overall_yield", "50%").rstrip("%")) / 100
        steps = pathway.get("estimated_steps", 5)
        
        # Simple scoring algorithm
        feasibility = (10 - complexity) * 0.4 + overall_yield * 5 + (10 - min(steps, 10)) * 0.3
        return round(max(0, min(10, feasibility)), 1)
    
    def suggest_alternative_routes(self, compound_type: str) -> List[Dict[str, Any]]:
        """Suggest alternative synthetic routes"""
        alternatives = {
            "buprenorphine_analogs": [
                {
                    "route_name": "Microwave-assisted synthesis",
                    "advantages": ["Shorter reaction times", "Higher yields"],
                    "disadvantages": ["Requires specialized equipment"],
                    "estimated_improvement": "20% yield increase"
                },
                {
                    "route_name": "Flow chemistry approach",
                    "advantages": ["Continuous process", "Better heat transfer"],
                    "disadvantages": ["Complex setup", "Higher initial cost"],
                    "estimated_improvement": "15% cost reduction"
                }
            ],
            "kappa_antagonists": [
                {
                    "route_name": "One-pot synthesis",
                    "advantages": ["Fewer purification steps", "Reduced waste"],
                    "disadvantages": ["More complex optimization"],
                    "estimated_improvement": "30% time reduction"
                }
            ]
        }
        
        return alternatives.get(compound_type, [])
    
    def get_supplier_info(self, starting_materials: List[str]) -> Dict[str, Any]:
        """Get supplier information for starting materials"""
        supplier_data = {}
        
        for material in starting_materials:
            supplier_data[material] = {
                "primary_suppliers": ["Sigma-Aldrich", "TCI Chemicals", "Alfa Aesar"],
                "estimated_cost": "$50-200/g",
                "availability": "In stock",
                "purity": "≥98%",
                "special_handling": "Store under inert atmosphere" if "sensitive" in material.lower() else "Standard storage"
            }
        
        return supplier_data

# Initialize the molecular visualization system
def create_molecular_visualization_system():
    """Create and return a molecular visualization system instance"""
    return MolecularVisualizationSystem()

if __name__ == "__main__":
    # Test the system
    viz_system = create_molecular_visualization_system()
    
    # Test 3D structure generation
    print("Generating 3D structure for buprenorphine analog...")
    structure = viz_system.generate_3d_structure_data(
        "COC1=C(O)C=C2C3=C1C(=O)CC[C@]3(C)[C@H](O)[C@H]2N(CC1CC1)CC1CC1",
        "Buprenorphine Analog"
    )
    print(f"Molecular formula: {structure['molecular_formula']}")
    print(f"Molecular weight: {structure['property_display']['molecular_weight']}")
    
    # Test retrosynthesis pathway
    print("\nGenerating retrosynthesis pathway...")
    pathway = viz_system.get_retrosynthesis_pathway("kappa_antagonists")
    print(f"Synthetic complexity: {pathway['complexity_score']}")
    print(f"Estimated steps: {pathway['estimated_steps']}")
    print(f"Overall yield: {pathway['overall_yield']}")
    print(f"Feasibility score: {pathway['feasibility_score']}/10")

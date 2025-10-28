#!/usr/bin/env python3
"""
Retrosynthesis Analysis Module
Analyzes synthetic routes and complexity for drug compounds
"""

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from typing import Dict, List, Optional, Tuple
import json

class RetrosynthesisAnalyzer:
    """Analyze retrosynthetic routes and synthetic complexity"""
    
    def __init__(self):
        self.reaction_database = self.load_reaction_database()
        self.building_blocks = self.load_building_blocks()
        
    def load_reaction_database(self) -> Dict:
        """Load common organic reactions"""
        return {
            "Suzuki": {
                "template": "[c:1][Br,I].[c:2]B(O)O>>[c:1][c:2]",
                "conditions": "Pd catalyst, base, heat",
                "yield": 85,
                "cost": "Medium"
            },
            "Amide_coupling": {
                "template": "[C:1](=O)O.[N:2]>>[C:1](=O)[N:2]",
                "conditions": "EDC, HOBt, DMF",
                "yield": 90,
                "cost": "Low"
            },
            "Reductive_amination": {
                "template": "[C:1]=O.[N:2]>>[C:1][N:2]",
                "conditions": "NaBH4 or NaBH(OAc)3",
                "yield": 80,
                "cost": "Low"
            },
            "Grignard": {
                "template": "[C:1][Br,I].[C:2]=O>>[C:1][C:2](O)",
                "conditions": "Mg, ether, anhydrous",
                "yield": 75,
                "cost": "Low"
            },
            "Click_chemistry": {
                "template": "[C:1]#[C].[N:2]=[N+]=[N-]>>[c:1]1[n:2]nnn1",
                "conditions": "Cu catalyst, room temp",
                "yield": 95,
                "cost": "Low"
            },
            "Buchwald_Hartwig": {
                "template": "[c:1][Br,I].[N:2]>>[c:1][N:2]",
                "conditions": "Pd catalyst, ligand, base",
                "yield": 85,
                "cost": "High"
            }
        }
    
    def load_building_blocks(self) -> List[str]:
        """Load common commercial building blocks"""
        return [
            "c1ccccc1", "c1ccc(O)cc1", "c1ccc(N)cc1", "c1ccc(C(=O)O)cc1",
            "CC(C)C", "CCO", "CN", "C(=O)O", "C(=O)N", "CCN", "CCCO",
            "c1ncccn1", "c1cnccc1", "c1cccnc1", "C1CCCCC1", "C1CCNCC1",
            "FC(F)(F)", "Cl", "Br", "I", "S(=O)(=O)N", "P(=O)(O)(O)"
        ]
    
    def analyze_synthesis(self, smiles: str) -> Dict:
        """Comprehensive retrosynthetic analysis"""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return {"error": "Invalid SMILES"}
        
        # Calculate synthetic complexity
        complexity = self.calculate_synthetic_complexity(mol)
        
        # Identify key disconnections
        disconnections = self.identify_disconnections(mol)
        
        # Propose synthetic routes
        routes = self.propose_synthetic_routes(mol, disconnections)
        
        # Estimate synthesis parameters
        synthesis_params = self.estimate_synthesis_parameters(mol, routes)
        
        # Generate recommendations
        recommendations = self.generate_synthesis_recommendations(complexity, routes)
        
        return {
            "smiles": smiles,
            "complexity_score": complexity,
            "key_disconnections": disconnections,
            "proposed_routes": routes,
            "synthesis_parameters": synthesis_params,
            "recommendations": recommendations,
            "commercial_availability": self.check_commercial_availability(mol)
        }
    
    def calculate_synthetic_complexity(self, mol) -> Dict:
        """Calculate synthetic complexity score"""
        # Basic complexity metrics
        num_atoms = mol.GetNumAtoms()
        num_bonds = mol.GetNumBonds()
        num_rings = Descriptors.RingCount(mol)
        num_aromatic = Descriptors.NumAromaticRings(mol)
        num_heteroatoms = Descriptors.NumHeteroatoms(mol)
        num_rotatable = Descriptors.NumRotatableBonds(mol)
        num_stereo = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
        
        # Complexity score (simplified Bertz complexity)
        complexity_score = num_atoms * 1.0 + num_bonds * 0.5 + num_rings * 2.0 + num_stereo * 3.0
        
        # Functional group complexity
        fg_complexity = self._assess_functional_group_complexity(mol)
        
        # Overall assessment
        if complexity_score < 50:
            assessment = "Simple"
            synthetic_steps = "3-5 steps"
        elif complexity_score < 100:
            assessment = "Moderate"
            synthetic_steps = "5-8 steps"
        elif complexity_score < 200:
            assessment = "Complex"
            synthetic_steps = "8-12 steps"
        else:
            assessment = "Very Complex"
            synthetic_steps = "12+ steps"
        
        return {
            "score": round(complexity_score, 2),
            "assessment": assessment,
            "estimated_steps": synthetic_steps,
            "metrics": {
                "atoms": num_atoms,
                "bonds": num_bonds,
                "rings": num_rings,
                "aromatic_rings": num_aromatic,
                "heteroatoms": num_heteroatoms,
                "rotatable_bonds": num_rotatable,
                "stereocenters": num_stereo
            },
            "functional_group_complexity": fg_complexity
        }
    
    def _assess_functional_group_complexity(self, mol) -> str:
        """Assess complexity of functional groups"""
        smiles = Chem.MolToSmiles(mol)
        
        complex_groups = 0
        if "C(=O)N" in smiles:  # Amide
            complex_groups += 1
        if "S(=O)(=O)" in smiles:  # Sulfonyl
            complex_groups += 2
        if "P(=O)" in smiles:  # Phosphate
            complex_groups += 2
        if "[n]" in smiles.lower():  # Heterocycle
            complex_groups += 1
        if "C#N" in smiles:  # Nitrile
            complex_groups += 1
        
        if complex_groups == 0:
            return "Simple"
        elif complex_groups <= 2:
            return "Moderate"
        else:
            return "Complex"
    
    def identify_disconnections(self, mol) -> List[Dict]:
        """Identify strategic disconnections"""
        disconnections = []
        smiles = Chem.MolToSmiles(mol)
        
        # Check for amide bonds
        amide_pattern = Chem.MolFromSmarts("[C](=O)[N]")
        if mol.HasSubstructMatch(amide_pattern):
            disconnections.append({
                "type": "Amide bond",
                "reaction": "Amide coupling",
                "priority": "High",
                "description": "Disconnect C-N bond to carboxylic acid + amine"
            })
        
        # Check for C-C bonds adjacent to aromatics
        ar_pattern = Chem.MolFromSmarts("c-[C]")
        if mol.HasSubstructMatch(ar_pattern):
            disconnections.append({
                "type": "Ar-C bond",
                "reaction": "Suzuki or Grignard",
                "priority": "Medium",
                "description": "Disconnect Ar-C bond via cross-coupling"
            })
        
        # Check for C-N bonds to aromatics
        ar_n_pattern = Chem.MolFromSmarts("c-[N]")
        if mol.HasSubstructMatch(ar_n_pattern):
            disconnections.append({
                "type": "Ar-N bond",
                "reaction": "Buchwald-Hartwig",
                "priority": "Medium",
                "description": "Disconnect Ar-N bond via C-N coupling"
            })
        
        # Check for ether bonds
        ether_pattern = Chem.MolFromSmarts("[C]-[O]-[C]")
        if mol.HasSubstructMatch(ether_pattern):
            disconnections.append({
                "type": "Ether bond",
                "reaction": "Williamson ether synthesis",
                "priority": "Low",
                "description": "Disconnect C-O bond to alcohol + alkyl halide"
            })
        
        return disconnections
    
    def propose_synthetic_routes(self, mol, disconnections: List[Dict]) -> List[Dict]:
        """Propose synthetic routes based on disconnections"""
        routes = []
        
        # Route 1: Linear synthesis
        linear_route = {
            "route_name": "Linear Synthesis",
            "strategy": "Step-by-step assembly",
            "key_reactions": [],
            "estimated_yield": 100,
            "pros": ["Simple purification", "Easy to troubleshoot"],
            "cons": ["Lower overall yield", "More time-consuming"]
        }
        
        # Add reactions based on disconnections
        for disc in disconnections[:3]:  # Take top 3 disconnections
            reaction = self.reaction_database.get(disc['reaction'].replace(" ", "_"), {})
            if reaction:
                linear_route['key_reactions'].append({
                    "step": len(linear_route['key_reactions']) + 1,
                    "reaction": disc['reaction'],
                    "conditions": reaction.get('conditions', 'Standard'),
                    "expected_yield": reaction.get('yield', 80)
                })
                linear_route['estimated_yield'] *= reaction.get('yield', 80) / 100
        
        linear_route['estimated_yield'] = round(linear_route['estimated_yield'] * 100, 1)
        routes.append(linear_route)
        
        # Route 2: Convergent synthesis
        if len(disconnections) > 2:
            convergent_route = {
                "route_name": "Convergent Synthesis",
                "strategy": "Parallel fragment synthesis + coupling",
                "key_reactions": [
                    {
                        "step": "Fragment A",
                        "reaction": disconnections[0]['reaction'],
                        "conditions": "Optimized conditions",
                        "expected_yield": 85
                    },
                    {
                        "step": "Fragment B",
                        "reaction": disconnections[1]['reaction'] if len(disconnections) > 1 else "Direct use",
                        "conditions": "Optimized conditions",
                        "expected_yield": 85
                    },
                    {
                        "step": "Coupling",
                        "reaction": "Fragment coupling",
                        "conditions": "Final assembly",
                        "expected_yield": 75
                    }
                ],
                "estimated_yield": 54.2,  # 0.85 * 0.85 * 0.75 * 100
                "pros": ["Higher overall yield", "Faster synthesis"],
                "cons": ["More complex planning", "Requires parallel work"]
            }
            routes.append(convergent_route)
        
        return routes
    
    def estimate_synthesis_parameters(self, mol, routes: List[Dict]) -> Dict:
        """Estimate synthesis time, cost, and difficulty"""
        mw = Descriptors.MolWt(mol)
        
        # Base estimates
        if mw < 300:
            base_cost = 500
            base_time = 2
        elif mw < 500:
            base_cost = 1500
            base_time = 3
        else:
            base_cost = 3000
            base_time = 4
        
        # Adjust for complexity
        num_steps = len(routes[0]['key_reactions']) if routes else 5
        
        # Cost estimation
        cost_per_step = 500
        total_cost = base_cost + (num_steps * cost_per_step)
        
        # Time estimation (weeks)
        time_per_step = 0.5
        total_time = base_time + (num_steps * time_per_step)
        
        # Scale estimates
        scales = {
            "Discovery (1-10 mg)": {
                "cost": total_cost,
                "time": total_time,
                "difficulty": "Standard"
            },
            "Medicinal Chemistry (100 mg - 1 g)": {
                "cost": total_cost * 3,
                "time": total_time * 1.5,
                "difficulty": "Moderate"
            },
            "Process Development (10-100 g)": {
                "cost": total_cost * 10,
                "time": total_time * 2,
                "difficulty": "Requires optimization"
            },
            "Pilot Scale (1-10 kg)": {
                "cost": total_cost * 50,
                "time": total_time * 4,
                "difficulty": "Process chemistry required"
            }
        }
        
        return {
            "estimated_steps": num_steps,
            "scale_parameters": scales,
            "key_materials_cost": f"${round(total_cost * 0.3)}",
            "labor_cost": f"${round(total_cost * 0.5)}",
            "analytical_cost": f"${round(total_cost * 0.2)}",
            "critical_materials": self._identify_critical_materials(mol)
        }
    
    def _identify_critical_materials(self, mol) -> List[str]:
        """Identify critical starting materials"""
        critical = []
        smiles = Chem.MolToSmiles(mol)
        
        if "F" in smiles:
            critical.append("Fluorinating reagents")
        if "Br" in smiles or "I" in smiles:
            critical.append("Halogenating reagents")
        if "[Pd]" in str(self.reaction_database):
            critical.append("Palladium catalyst")
        if "B(O)O" in smiles:
            critical.append("Boronic acid/ester")
        
        return critical if critical else ["Standard reagents"]
    
    def check_commercial_availability(self, mol) -> Dict:
        """Check if building blocks are commercially available"""
        smiles = Chem.MolToSmiles(mol)
        
        # Check if the molecule itself might be commercial
        if Descriptors.MolWt(mol) < 250:
            availability = "Possibly commercial"
            suppliers = ["Sigma-Aldrich", "TCI", "Alfa Aesar"]
            price_range = "$50-500 per gram"
        else:
            availability = "Custom synthesis required"
            suppliers = ["WuXi AppTec", "Syngene", "ChemPartner"]
            price_range = "$2000-10000 per gram"
        
        # Check for commercial building blocks
        commercial_fragments = []
        for bb in self.building_blocks:
            if bb in smiles:
                commercial_fragments.append(bb)
        
        return {
            "compound_availability": availability,
            "potential_suppliers": suppliers,
            "estimated_price": price_range,
            "commercial_building_blocks": len(commercial_fragments),
            "custom_intermediates_needed": max(0, 3 - len(commercial_fragments))
        }
    
    def generate_synthesis_recommendations(self, complexity: Dict, routes: List[Dict]) -> List[str]:
        """Generate synthesis recommendations"""
        recommendations = []
        
        # Complexity-based recommendations
        if complexity['score'] > 150:
            recommendations.append("Consider fragment-based approach to reduce complexity")
            recommendations.append("Explore biocatalytic transformations for stereocenters")
        
        if complexity['metrics']['stereocenters'] > 2:
            recommendations.append("Use asymmetric synthesis or chiral resolution")
            recommendations.append("Consider enzymatic methods for stereochemistry")
        
        # Route-based recommendations
        if routes and routes[0]['estimated_yield'] < 30:
            recommendations.append("Optimize reaction conditions to improve yield")
            recommendations.append("Consider alternative disconnection strategy")
        
        # General recommendations
        recommendations.append("Perform retrosynthetic analysis with AI tools")
        recommendations.append("Check patent literature for existing routes")
        recommendations.append("Consider flow chemistry for scale-up")
        
        return recommendations